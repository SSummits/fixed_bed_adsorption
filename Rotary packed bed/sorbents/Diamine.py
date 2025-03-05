from pyomo.environ import Var, Param, Block, exp, log, units, Reals, NonNegativeReals
from torch import P

def add_diamine_parameters(RPB):
    RPB.DA = Block()
    DA = RPB.DA

    DA.eb = Param(initialize=(0.68), doc="bed voidage")
    DA.ep = Param(initialize=(0.68), doc="particle porosity")
    DA.dp = Param(
        initialize=(0.000525), units=units.m, doc="particle diameter [m]"
    )
    DA.Cp_sol = Param(
        initialize=(1.457),
        units=units.kJ / units.kg / units.K,
        doc="solid heat capacity [kJ/kg/K]",
    )
    DA.rho_sol = Param(
        initialize=(1000),
        units=units.kg / units.m**3,
        mutable=True,
        doc="solid particle densitry [kg/m^3]",
    )

    @DA.Expression(doc="particle radius [m]")
    def rp(b):
        return b.dp / 2

    @DA.Expression(
        doc="specific particle area for mass transfer, bed voidage"
        "included [m^2/m^3 bed]"
    )
    def a_s(b):
        return 6 / b.dp * (1 - b.eb)
    
    # Isotherm Parameters
    DA.b_chem_0 = Param(
        initialize=28.56, units=1/units.bar, doc="isotherm parameter"
    )
    DA.Q_st_chem = Param(
        initialize=72.56, units=units.kJ/units.mol, doc="isotherm parameter"
    )
    DA.n_chem_0 = Param(
        initialize=0.21, doc="isotherm parameter"
    )
    DA.b_phys_0 = Param(
        initialize=0.62, units=1/units.bar, doc="isotherm parameter"
    )
    DA.Q_st_phys = Param(
        initialize=43.83, units=units.kJ/units.mol, doc="isotherm parameter"
    )
    DA.n_phys = Param(
        initialize=1.46, doc="isotherm parameter"
    )
    DA.N_phys = Param(
        initialize=3.52, units=units.mmol/units.g, doc="isotherm parameter"
    )
    DA.N_chem = Param(
        initialize=3.82, units=units.mmol/units.g, doc="isotherm parameter"
    )
    DA.K_a = Param(
        initialize=-0.92, doc="isotherm parameter"
    )
    DA.K_b = Param(
        initialize=324.86, units=units.K, doc="isotherm parameter"
    )
    DA.E_n = Param(
        initialize=11.29, units=units.kJ/units.mol, doc="isotherm parameter"
    )
    DA.K_c = Param(
        initialize=-71.14, doc="isotherm parameter"
    )
    DA.K_d = Param(
        initialize=2.84e4, units=units.K, doc="isotherm parameter"
    )
    DA.k_chem_0 = Param(
        initialize=0.0136, units=1/units.s, doc="mass transfer parameter"
    )
    DA.E_chem = Param(
        initialize=23.21, units=units.kJ/units.mol, doc="mass transfer parameter"
    )
    DA.k_phys_0 = Param(
        initialize=0.0823, units=1/units.s, doc="mass transfer parameter"
    )
    DA.E_phys = Param(
        initialize=7.18, units=units.kJ/units.mol, doc="mass transfer parameter"
    )

    # Mass transfer parameters
    DA.C1 = Param(
        initialize=(4.11e-12),
        units=units.m**2 / units.K**0.5 / units.s,
        doc="lumped MT parameter [m^2/K^0.5/s]",
    )

def add_diamine_isotherm(blk, initial_guesses):
    RPB = blk.parent_block()
    DA_param = RPB.DA
    blk.DA = Block()
    DA = blk.DA

    if initial_guesses == "adsorption":
        qCO2_in_init = 1
        Ts_in_init = 100 + 273
    elif initial_guesses == "desorption":
        qCO2_in_init = 2.5
        Ts_in_init = 110 + 273
    else:
        qCO2_in_init = 1
        Ts_in_init = 100 + 273

    blk.DA.qCO2 = Var(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        ['chem', 'phys'],
        initialize=qCO2_in_init,
        domain=NonNegativeReals,
        bounds=(0, 5),
        doc="CO2 loading [mol/kg]",
        units=units.mol / units.kg,
    )

    def b_chem(T):
        return DA_param.b_chem_0 * exp(
            (DA_param.Q_st_chem / RPB.R / DA_param.T0)*(DA_param.T0/T - 1)
        )
    def b_phys(T):
        return DA_param.b_phys_0 * exp(
            (DA_param.Q_st_phys / RPB.R / DA_param.T0)*(DA_param.T0/T - 1)
        )
    
    def q_chem_inf(T):
        return DA_param.N_chem * (
            exp(DA_param.K_a + DA_param.K_b / T) / 
            (1 + exp(DA_param.K_a + DA_param.K_b / T))
        )
    def q_phys_inf(T):
        return DA_param.N_phys * (
            exp(DA_param.K_c + DA_param.K_d / T) / 
            (1 + exp(DA_param.K_c + DA_param.K_d / T))
        )
    
    def n_chem(T):
        return DA_param.n_chem_0 * exp(
            (DA_param.E_n / RPB.R / DA_param.T0)*(DA_param.T0/T - 1)
        )
    
    def q_star_chem(T,P):
        return q_chem_inf(T) * ((b_chem(T) * P)**(1/n_chem(T)) / (1 + (b_chem(T) * P)**(1/n_chem(T))))
    def q_star_phys(T,P):
        return q_phys_inf(T) * ((b_phys(T) * P)**(1/DA_param.n_phys) / (1 + (b_phys(T) * P)**(1/DA_param.n_phys))) 
    
    @DA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="isotherm loading expression [mol/kg]",
    )
    def qCO2_eq(b, t, z, o):
        return q_star_chem(blk.Ts[t, z, o], blk.P_surf[t, z, o]) + q_star_phys(blk.Ts[t, z, o], blk.P_surf[t, z, o])
    

    a1_FL = 0.02
    a2_FL = 0.98
    sig_FL = 0.01

    def FL(z):
        def FL_1(z):
            return exp((z - a1_FL) / sig_FL) / (1 + exp((z - a1_FL) / sig_FL))

        def FL_2(z):
            return exp((z - a2_FL) / sig_FL) / (1 + exp((z - a2_FL) / sig_FL))

        return FL_1(z) - FL_2(z)
    
    def k_chem(T):
        return DA_param.k_chem_0 * exp(
            -DA_param.E_chem / RPB.R / DA_param.T0 * (DA_param.T0 / T - 1)
        )
    def k_phys(T):
        return DA_param.k_phys_0 * exp(
            -DA_param.E_phys / RPB.R / DA_param.T0 * (DA_param.T0 / T - 1)
        )
    
    @blk.DA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="effective diffusion in solids [m^2/s]",
    )
    def Deff(b, t, z, o):
        return DA_param.C1 * blk.Ts[t, z, o] ** 0.5

    @blk.DA.Expression(
        RPB.flowsheet().time, blk.z, blk.o, doc="internal MT coeff. [1/s]"
    )
    def k_I(b, t, z, o):
        return (
            blk.R_MT_coeff * (15 * DA_param.ep * b.Deff[t, z, o] / DA_param.rp**2)
            + (1 - blk.R_MT_coeff) * 0.001 / units.s
        )
    
    blk.DA.k_0 = Var(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        ['chem', 'phys'],
        initialize=1,
        units=1/units.s,
        doc="Mass transfer coefficient for chemical and physical adsorption"
    )
    @blk.DA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        ['chem', 'phys'],
        doc="Mass transfer coefficient equation for chemical and physical adsorption"
    )
    def k_0_eqn(b, t, z, o, i):
        if i == 'chem':
            k = k_chem(b.Ts[t, z, o])
        elif i == 'phys':
            k = k_phys(b.Ts[t, z, o])
        
        return k == b.k_0[t, z, o, i] * (k / b.k_I[t, z, o] + 1)

    blk.DA.Rs_CO2 = Var(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        initialize=0,
        domain=Reals,
        units=units.mol / units.s / units.m**3,
        doc="solids mass transfer rate [mol/s/m^3 bed]",
    )

    @blk.DA.Constraint(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="solids mass transfer rate [mol/s/m^3 bed]",
    )
    def Rs_CO2_eq(b, t, z, o):
        flux_lim = FL(z)

        if 0 < z < 1 and 0 < o < 1:
            Ts = blk.Ts[t, z, o]
            P = blk.P_surf[t, z, o]
            return (
                b.Rs_CO2[t, z, o]
                == flux_lim
                * (k_chem(b.Ts[t, z, o]) * (q_star_chem(blk.Ts[t, z, o], blk.P_surf[t, z, o]) - b.qCO2[t, z, o, 'chem'])
                   + k_phys(b.Ts[t, z, o]) * (q_star_phys(blk.Ts[t, z, o], blk.P_surf[t, z, o]) - b.qCO2[t, z, o, 'phys'])
                )
                * (1 - RPB.eb[z])
                * DA_param.rho_sol * RPB.sorbent_weight[z, 'DA']
            )
        else:
            return b.Rs_CO2[t, z, o] == 0 * units.mol / units.s / units.m**3
    
def add_diamine_adsortion_heat(DA):
    RPB = DA.parent_block().parent_block()
    blk = DA.parent_block()
    DA_param = RPB.DA
    @DA.Expression(
        RPB.flowsheet().time, blk.z, blk.o, doc="heat of adsorption [kJ/mol]"
    )
    def delH_CO2(b, t, z, o):
        return -65 * units.kJ / units.mol
        # return -(
        #     TA_param.delH_1
        #     - (TA_param.delH_1 - TA_param.delH_2)
        #     * exp(TA_param.delH_a1 * (b.qCO2_eq[t, z, o] - TA_param.delH_b1))
        #     / (1 + exp(TA_param.delH_a1 * (b.qCO2_eq[t, z, o] - TA_param.delH_b1)))
        #     - (TA_param.delH_2 - TA_param.delH_3)
        #     * exp(TA_param.delH_a2 * (b.qCO2_eq[t, z, o] - TA_param.delH_b2))
        #     / (1 + exp(TA_param.delH_a2 * (b.qCO2_eq[t, z, o] - TA_param.delH_b2)))
        # )