from pyomo.environ import Param, Block, exp, log, units

def add_tetraamine_parameters(RPB):
    RPB.TA = Block()
    TA = RPB.TA

    TA.eb = Param(initialize=(0.73), doc="bed voidage")
    TA.ep = Param(initialize=(0.68), doc="particle porosity")
    TA.dp = Param(
        initialize=(0.000525), units=units.m, doc="particle diameter [m]"
    )
    TA.Cp_sol = Param(
        initialize=(1.457),
        units=units.kJ / units.kg / units.K,
        doc="solid heat capacity [kJ/kg/K]",
    )
    TA.rho_sol = Param(
        initialize=(1144),
        units=units.kg / units.m**3,
        mutable=True,
        doc="solid particle densitry [kg/m^3]",
    )

    @TA.Expression(doc="particle radius [m]")
    def rp(b):
        return b.dp / 2

    @TA.Expression(
        doc="specific particle area for mass transfer, bed voidage"
        "included [m^2/m^3 bed]"
    )
    def a_s(b):
        return 6 / b.dp * (1 - b.eb)

    # Isotherm Parameters
    TA.q_inf_1 = Param(
        initialize=2.87e-02, units=units.mol / units.kg, doc="isotherm parameter"
    )
    TA.q_inf_2 = Param(
        initialize=1.95, units=units.mol / units.kg, doc="isotherm parameter"
    )
    TA.q_inf_3 = Param(
        initialize=3.45, units=units.mol / units.kg, doc="isotherm parameter"
    )

    TA.d_inf_1 = Param(
        initialize=1670.31, units=units.bar**-1, doc="isotherm parameter"
    )
    TA.d_inf_2 = Param(
        initialize=789.01, units=units.bar**-1, doc="isotherm parameter"
    )
    TA.d_inf_3 = Param(
        initialize=10990.67, units=units.bar**-1, doc="isotherm parameter"
    )
    TA.d_inf_4 = Param(
        initialize=0.28,
        units=units.mol / units.kg / units.bar,
        doc="isotherm parameter",
    )

    TA.E_1 = Param(
        initialize=-76.15, units=units.kJ / units.mol, doc="isotherm parameter"
    )
    TA.E_2 = Param(
        initialize=-77.44, units=units.kJ / units.mol, doc="isotherm parameter"
    )
    TA.E_3 = Param(
        initialize=-194.48, units=units.kJ / units.mol, doc="isotherm parameter"
    )
    TA.E_4 = Param(
        initialize=-6.76, units=units.kJ / units.mol, doc="isotherm parameter"
    )

    TA.X_11 = Param(initialize=4.20e-2, doc="isotherm parameter")
    TA.X_21 = Param(initialize=2.97, units=units.K, doc="isotherm parameter")
    TA.X_12 = Param(initialize=7.74e-2, doc="isotherm parameter")
    TA.X_22 = Param(initialize=1.66, units=units.K, doc="isotherm parameter")

    TA.P_step_01 = Param(
        initialize=1.85e-03, units=units.bar, doc="isotherm parameter"
    )
    TA.P_step_02 = Param(
        initialize=1.78e-02, units=units.bar, doc="isotherm parameter"
    )

    @TA.Expression()
    def ln_P0_1(b):
        return log(b.P_step_01 / units.bar)

    @TA.Expression()
    def ln_P0_2(b):
        return log(b.P_step_02 / units.bar)

    TA.H_step_1 = Param(
        initialize=-99.64, units=units.kJ / units.mol, doc="isotherm parameter"
    )
    TA.H_step_2 = Param(
        initialize=-78.19, units=units.kJ / units.mol, doc="isotherm parameter"
    )

    TA.gamma_1 = Param(initialize=894.67, doc="isotherm parameter")
    TA.gamma_2 = Param(initialize=95.22, doc="isotherm parameter")

    TA.T0 = Param(initialize=363.15, units=units.K, doc="isotherm parameter")

    # Mass transfer parameters
    TA.C1 = Param(
        initialize=(2.562434e-12),
        units=units.m**2 / units.K**0.5 / units.s,
        doc="lumped MT parameter [m^2/K^0.5/s]",
    )

    # heat of adsorption parameters
    TA.delH_a1 = Param(
        initialize=21.68,
        units=units.kg / units.mol,
        doc="heat of adsorption parameter",
    )
    TA.delH_a2 = Param(
        initialize=29.10,
        units=units.kg / units.mol,
        doc="heat of adsorption parameter",
    )
    TA.delH_b1 = Param(
        initialize=1.59,
        units=units.mol / units.kg,
        doc="heat of adsorption parameter",
    )
    TA.delH_b2 = Param(
        initialize=3.39,
        units=units.mol / units.kg,
        doc="heat of adsorption parameter",
    )
    TA.delH_1 = Param(
        initialize=98.76,
        units=units.kJ / units.mol,
        doc="heat of adsorption parameter",
    )
    TA.delH_2 = Param(
        initialize=77.11,
        units=units.kJ / units.mol,
        doc="heat of adsorption parameter",
    )
    TA.delH_3 = Param(
        initialize=21.25,
        units=units.kJ / units.mol,
        doc="heat of adsorption parameter",
    )

def add_tetraamine_isotherm(blk):
    RPB = blk.parent_block()
    TA_param = RPB.TA
    blk.TA = Block()
    TA = blk.TA

    def d_1(T):
        return TA_param.d_inf_1 * exp(
            -TA_param.E_1 / (RPB.R * TA_param.T0) * (TA_param.T0 / T - 1)
        )

    def d_2(T):
        return TA_param.d_inf_2 * exp(
            -TA_param.E_2 / (RPB.R * TA_param.T0) * (TA_param.T0 / T - 1)
        )

    def d_3(T):
        return TA_param.d_inf_3 * exp(
            -TA_param.E_3 / (RPB.R * TA_param.T0) * (TA_param.T0 / T - 1)
        )

    def d_4(T):
        return TA_param.d_inf_4 * exp(
            -TA_param.E_4 / (RPB.R * TA_param.T0) * (TA_param.T0 / T - 1)
        )

    def sigma_1(T):
        return TA_param.X_11 * exp(TA_param.X_21 * (1 / TA_param.T0 - 1 / T))

    def sigma_2(T):
        return TA_param.X_12 * exp(TA_param.X_22 * (1 / TA_param.T0 - 1 / T))

    def ln_pstep1(T):
        return TA_param.ln_P0_1 + (-TA_param.H_step_1 / RPB.R * (1 / TA_param.T0 - 1 / T))

    def ln_pstep2(T):
        return TA_param.ln_P0_2 + (-TA_param.H_step_2 / RPB.R * (1 / TA_param.T0 - 1 / T))

    def q_star_1(P, T):
        return TA_param.q_inf_1 * d_1(T) * P / (1 + d_1(T) * P)

    def q_star_2(P, T):
        return TA_param.q_inf_2 * d_2(T) * P / (1 + d_2(T) * P)

    def q_star_3(P, T):
        return TA_param.q_inf_3 * d_3(T) * P / (1 + d_3(T) * P) + d_4(T) * P

    @TA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="Partial pressure of CO2 at particle surface [bar] (ideal gas law)",
    )
    def P_surf(b, t, z, o):
        # smooth max operator: max(0, x) = 0.5*(x + (x^2 + eps)^0.5)
        eps = 1e-8
        Cs_r_smooth_max = 0.5 * (
            blk.Cs_r[t, z, o]
            + (blk.Cs_r[t, z, o] ** 2 + eps * (units.mol / units.m**3) ** 2) ** 0.5
        )
        return Cs_r_smooth_max * RPB.Rg * blk.Ts[t, z, o]
        # return smooth_max(0,m.Cs_r[z, o]) * RPB.Rg * m.Ts[z, o] #idaes smooth_max doesn't carry units through

    @TA.Expression(RPB.flowsheet().time, blk.z, blk.o, doc="log(Psurf)")
    def ln_Psurf(b, t, z, o):
        return log(b.P_surf[t, z, o] / units.bar)  # must make dimensionless

    @TA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="weighting function term1: (ln_Psurf-ln_Pstep)/sigma",
    )
    def iso_w_term1(b, t, z, o):
        return (b.ln_Psurf[t, z, o] - ln_pstep1(blk.Ts[t, z, o])) / sigma_1(
            blk.Ts[t, z, o]
        )

    @TA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="weighting function term2: (ln_Psurf-ln_Pstep)/sigma",
    )
    def iso_w_term2(b, t, z, o):
        return (b.ln_Psurf[t, z, o] - ln_pstep2(blk.Ts[t, z, o])) / sigma_2(
            blk.Ts[t, z, o]
        )

    @TA.Expression(
        RPB.flowsheet().time, blk.z, blk.o, doc="log of weighting function 1"
    )
    def ln_w1(b, t, z, o):
        # return gamma_1*log(exp(m.iso_w_term1[z,o])/(1+exp(m.iso_w_term1[z,o])))
        # return gamma_1*(log(exp(m.iso_w_term1[z,o])) - log(1+exp(m.iso_w_term1[z,o])))
        return TA_param.gamma_1 * (
            b.iso_w_term1[t, z, o] - log(1 + exp(b.iso_w_term1[t, z, o]))
        )

    @TA.Expression(
        RPB.flowsheet().time, blk.z, blk.o, doc="log of weighting function 2"
    )
    def ln_w2(b, t, z, o):
        # return gamma_2*log(exp(m.iso_w_term2[z,o])/(1+exp(m.iso_w_term2[z,o])))
        # return gamma_2*(log(exp(m.iso_w_term2[z,o])) - log(1+exp(m.iso_w_term2[z,o])))
        return TA_param.gamma_2 * (
            b.iso_w_term2[t, z, o] - log(1 + exp(b.iso_w_term2[t, z, o]))
        )

    @TA.Expression(RPB.flowsheet().time, blk.z, blk.o, doc="weighting function 1")
    def iso_w1(b, t, z, o):
        return exp(b.ln_w1[t, z, o])

    @TA.Expression(RPB.flowsheet().time, blk.z, blk.o, doc="weighting function 2")
    def iso_w2(b, t, z, o):
        return exp(b.ln_w2[t, z, o])

    @TA.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="isotherm loading expression [mol/kg]",
    )
    def qCO2_eq(b, t, z, o):
        return (
            (1 - b.iso_w1[t, z, o]) * q_star_1(b.P_surf[t, z, o], blk.Ts[t, z, o])
            + (b.iso_w1[t, z, o] - b.iso_w2[t, z, o])
            * q_star_2(b.P_surf[t, z, o], blk.Ts[t, z, o])
            + b.iso_w2[t, z, o] * q_star_3(b.P_surf[t, z, o], blk.Ts[t, z, o])
        )

def add_tetraamine_adsortion_heat(TA):
    RPB = TA.parent_block().parent_block()
    blk = TA.parent_block()
    TA_param = RPB.TA
    @TA.Expression(
        RPB.flowsheet().time, blk.z, blk.o, doc="heat of adsorption [kJ/mol]"
    )
    def delH_CO2(b, t, z, o):
        return -(
            TA_param.delH_1
            - (TA_param.delH_1 - TA_param.delH_2)
            * exp(TA_param.delH_a1 * (b.qCO2_eq[t, z, o] - TA_param.delH_b1))
            / (1 + exp(TA_param.delH_a1 * (b.qCO2_eq[t, z, o] - TA_param.delH_b1)))
            - (TA_param.delH_2 - TA_param.delH_3)
            * exp(TA_param.delH_a2 * (b.qCO2_eq[t, z, o] - TA_param.delH_b2))
            / (1 + exp(TA_param.delH_a2 * (b.qCO2_eq[t, z, o] - TA_param.delH_b2)))
        )