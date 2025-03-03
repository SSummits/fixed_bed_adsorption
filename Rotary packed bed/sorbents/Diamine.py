from pyomo.environ import Param, Block, exp, log, units

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

    # Mass transfer parameters
    DA.C1 = Param(
        initialize=(4.11e-12),
        units=units.m**2 / units.K**0.5 / units.s,
        doc="lumped MT parameter [m^2/K^0.5/s]",
    )

def add_diamine_isotherm(blk):
    RPB = blk.parent_block()
    DA_param = RPB.DA
    blk.DA = Block()
    DA = blk.DA

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