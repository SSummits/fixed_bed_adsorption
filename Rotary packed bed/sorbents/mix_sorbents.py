from pydoc import doc
from matplotlib.pylab import f
from pyomo.environ import Var, Constraint
from .Tetraamine import add_tetraamine_parameters, add_tetraamine_isotherm, add_tetraamine_adsortion_heat

def mix_sorbent_params(RPB, sorbent_list):
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        sorbent_set.append('TA')
        add_tetraamine_parameters(RPB)
    if 'Diamine' in sorbent_list:
        sorbent_set.append('DA')

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
    
def mix_sorbent_isotherm(blk, sorbent_list):
    RPB = blk.parent_block()
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        add_tetraamine_isotherm(blk)
        sorbent_set.append('TA')
    if 'Diamine' in sorbent_list:
        sorbent_set.append('DA')
        pass

    @blk.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="mixed sorbent isothermal loading [mol/kg]"
    )
    def qCO2_eq(b, t, z, o):
        return sum(
            RPB.sorbent_weight[z, s] * getattr(b, s).qCO2_eq[t, z, o]
            for s in sorbent_set
        )

def mix_sorbent_adsorption_heat(blk, sorbent_list):
    RPB = blk.parent_block()
    sorbent_set = []
    if 'Tetraamine' in sorbent_list:
        add_tetraamine_adsortion_heat(blk.TA)
        sorbent_set = []
    if 'Diamine' in sorbent_list:
        sorbent_set.append('DA')
        pass

    @blk.Expression(
        RPB.flowsheet().time,
        blk.z,
        blk.o,
        doc="mixed sorbent adsorption heat [kJ/kg/K]"
    )
    def delH_CO2(b, t, z, o):
        return sum(
            RPB.sorbent_weight[z, s] * getattr(b, s).Cp_s_eq[t, z, o]
            for s in sorbent_set
        )
    