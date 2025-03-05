from pyomo.environ import Var, Constraint, Reals, units
from .Tetraamine import add_tetraamine_parameters, add_tetraamine_isotherm, add_tetraamine_adsortion_heat
from .Diamine import add_diamine_parameters, add_diamine_isotherm, add_diamine_adsortion_heat
import idaes.core.util.scaling as iscale

def mix_sorbent_params(RPB, sorbent_list):
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        sorbent_set.append('TA')
        add_tetraamine_parameters(RPB)
    if 'Diamine' in sorbent_list:
        sorbent_set.append('DA')
        add_diamine_parameters(RPB)

    RPB.sorbent_weight = Var(RPB.z, sorbent_set, initialize=0.5, doc='sorbent bed composition')

    @RPB.Expression(RPB.z,
                    doc = "bed voidage")
    def eb(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).eb for s in sorbent_set)
    
    @RPB.Expression(RPB.z,
                    doc = "particle porosity")
    def ep(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).ep for s in sorbent_set)
    
    @RPB.Expression(RPB.z,
                    doc = "particle diameter [m]")
    def dp(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).dp for s in sorbent_set)
    
    @RPB.Expression(RPB.z,
                    doc = "solid particle density [kg/m^3]")
    def rho_sol(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).rho_sol for s in sorbent_set)
    
    @RPB.Expression(RPB.z,
                    doc = "solid heat capacity [kJ/kg/K]")
    def Cp_sol(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).Cp_sol for s in sorbent_set)
    
    @RPB.Expression(RPB.z,
                    doc = "specific particle area for mass transfer, "
                    "bed voidage included [m^2/m^3 bed]")
    def a_s(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).a_s for s in sorbent_set)
    
    @RPB.Expression(RPB.z,
                    doc = "particle radius [m]")
    def rp(b, z):
        return b.dp[z] / 2
    
    @RPB.Expression(RPB.z,
                    doc = "lumped MT parameter [m^2/K^0.5/s]")
    def C1(b, z):
        return sum(b.sorbent_weight[z, s] * getattr(b, s).C1 for s in sorbent_set)
    
def mix_sorbent_isotherm(blk, sorbent_list, initial_guesses):
    RPB = blk.parent_block()
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        add_tetraamine_isotherm(blk, initial_guesses)
        sorbent_set.append('TA')
    if 'Diamine' in sorbent_list:
        add_diamine_isotherm(blk, initial_guesses)
        sorbent_set.append('DA')

    @blk.Constraint(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="Total mass transfer rate from all sorbents"
    )
    def Rs_CO2_mix(b, t, z, o):
        return b.Rs_CO2[t, z, o] == sum(
            getattr(b, s).Rs_CO2[t, z, o] for s in sorbent_set
        )
    
    @blk.Constraint(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="Total loading from all sorbents"
    )
    def qCO2_mix(b, t, z, o):
        lhs = b.qCO2[t, z, o] * RPB.rho_sol[z]
        rhs = 0
        if 'TA' in sorbent_set:
            rhs += b.TA.qCO2[t, z, o] * RPB.TA.rho_sol * RPB.sorbent_weight[z, 'TA']
        if 'DA' in sorbent_set:
            rhs += sum(b.DA.qCO2[t, z, o, i] for i in ['chem', 'phys']) * RPB.DA.rho_sol * RPB.sorbent_weight[z, 'DA']
        return lhs == rhs

def mix_sorbent_adsorption_heat(blk, sorbent_list):
    RPB = blk.parent_block()
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        add_tetraamine_adsortion_heat(blk.TA)
        sorbent_set.append('TA')
    if 'Diamine' in sorbent_list:
        add_diamine_adsortion_heat(blk.DA)
        sorbent_set.append('DA')

    blk.Q_delH = Var(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        initialize=0,
        domain=Reals,
        units=units.kJ / units.s / units.m**3,
        doc="adsorption/desorption heat rate [kJ/s/m^3 bed]",
    )

    @blk.Constraint(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="adsorption/desorption heat rate [kJ/s/m^3 bed]",
    )
    def Q_delH_eq(b, t, z, o):
        Q = 0
        if 'TA' in sorbent_set:
            Q += b.TA.delH_CO2[t, z, o] * b.TA.Rs_CO2[t, z, o]
        if 'DA' in sorbent_set:
            Q += b.DA.delH_CO2[t, z, o] * b.DA.Rs_CO2[t, z, o]
        return b.Q_delH[t, z, o] == blk.R_delH * Q

    @blk.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="mixed sorbent adsorption heat [kJ/kg/K]"
    )
    def delH_CO2(b, t, z, o):
        return sum(
            RPB.sorbent_weight[z, s] * getattr(b, s).delH_CO2[t, z, o]
            for s in sorbent_set
        )
    
def scale_sorbent_mix(blk, sorbent_list):
    RPB = blk.parent_block()
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        sorbent_set.append('TA')
    if 'Diamine' in sorbent_list:
        sorbent_set.append('DA')
    
    for s in sorbent_set:
        sorbent = getattr(blk, s)

        for t in RPB.flowsheet().time:
            for z in blk.z:
                
                for o in blk.o:
                    if s == 'DA':
                        iscale.set_scaling_factor(sorbent.qCO2[t, z, o, 'chem'], 10)
                        iscale.set_scaling_factor(sorbent.qCO2[t, z, o, 'phys'], 10)
                    else:
                        iscale.set_scaling_factor(sorbent.qCO2[t, z, o], 10)
                    
                    if 0 < z < 1 and 0 < o < 1:
                        if s == 'DA':
                            iscale.set_scaling_factor(sorbent.pde_solidMB[t, z, o, 'chem'], 1e-3)
                            iscale.set_scaling_factor(sorbent.pde_solidMB[t, z, o, 'phys'], 1e-3)
                        else:
                            iscale.set_scaling_factor(sorbent.pde_solidMB[t, z, o], 1e-3)

                        iscale.set_scaling_factor(sorbent.Rs_CO2[t, z, o], 0.5)

def connect_mixed_sorbent_sections(RPB, sorbent_list):
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        sorbent_set.append('TA')
    if 'Diamine' in sorbent_list:
        sorbent_set.append('DA')

    if 'TA' in sorbent_set:
        for t in RPB.flowsheet().time:
            for z in [0, 1]:
                RPB.des.TA.qCO2[t, z, 0].fix()
                RPB.ads.TA.qCO2[t, z, 0].fix()

        @RPB.Constraint(
            RPB.flowsheet().time,
            RPB.z,
            doc="Tetraamine rich loading constraint"
        )
        def TA_rich_loading_constraint(b, t, z):
            if 0 < z < 1:
                return b.des.TA.qCO2[t, z, 0] == b.ads.TA.qCO2[t, z, 1]
            else:
                return Constraint.Skip
        @RPB.Constraint(
            RPB.flowsheet().time,
            RPB.z,
            doc="Tetraamine lean loading constraint"
        )
        def TA_lean_loading_constraint(b, t, z):
            if 0 < z < 1:
                return b.ads.TA.qCO2[t, z, 0] == b.des.TA.qCO2[t, z, 1]
            else:
                return Constraint.Skip
            
    if 'DA' in sorbent_set:
        for t in RPB.flowsheet().time:
            for z in [0, 1]:
                for i in ['chem', 'phys']:
                    RPB.des.DA.qCO2[t, z, 0, i].fix()
                    RPB.ads.DA.qCO2[t, z, 0, i].fix()

        @RPB.Constraint(
            RPB.flowsheet().time,
            RPB.z,
            ['chem', 'phys'],
            doc="Diamine rich loading constraint"
        )
        def DA_rich_loading_constraint(b, t, z, i):
            if 0 < z < 1:
                return b.des.DA.qCO2[t, z, 0, i] == b.ads.DA.qCO2[t, z, 1, i]
            else:
                return Constraint.Skip
        @RPB.Constraint(
            RPB.flowsheet().time,
            RPB.z,
            ['chem', 'phys'],
            doc="Diamine lean loading constraint"
        )
        def DA_lean_loading_constraint(b, t, z, i):
            if 0 < z < 1:
                return b.ads.DA.qCO2[t, z, 0, i] == b.des.DA.qCO2[t, z, 1, i]
            else:
                return Constraint.Skip