from math import pi
from re import L

from cycler import V
from pyomo.environ import Block, Var, Param, units as pyunits

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
    )

def RPB_Polishing_Costing(fs):

    fs.RPB_cost = Block()

    @fs.RPB_cost.Integral(
        fs.RPB.z,
        wrt=fs.RPB.z
    )
    def tetraamine_volume(b, z):
        D = fs.RPB.D
        L = fs.RPB.L
        V = pi * D**2 * L / 4
        return V * (1 - fs.RPB.eb[z]) * fs.RPB.sorbent_weight[z, 'TA']
    
    fs.RPB_cost.tetraamine_cost = Var(initialize = 200,
                                   units = pyunits.USD_2018/pyunits.ft**3)
    fs.RPB_cost.tetraamine_cost.fix()

    fs.RPB_cost.tetraamine_bare_erected_cost = Var(initialize=1,
                                        units = pyunits.MUSD_2018)
    @fs.RPB_cost.Constraint()
    def tetraamine_bare_erected_cost_eqn(b):
        return b.bare_erected_cost == pyunits.convert(
            b.tetraamine_volume * b.tetraamine_cost,
            pyunits.MUSD_2018
        )
    
    @fs.RPB_cost.Integral(
        fs.RPB.z,
        wrt=fs.RPB.z
    )
    def diamine_volume(b, z):
        D = fs.RPB.D
        L = fs.RPB.L
        V = pi * D**2 * L / 4
        return V * (1 - fs.RPB.eb[z]) * fs.RPB.sorbent_weight[z, 'DA']
    
    fs.RPB_cost.diamine_cost = Var(initialize = 200,
                                   units = pyunits.USD_2018/pyunits.ft**3)
    fs.RPB_cost.diamine_cost.fix()

    fs.RPB_cost.diamine_bare_erected_cost = Var(initialize=1,
                                        units = pyunits.MUSD_2018)
    @fs.RPB_cost.Constraint()
    def diamine_bare_erected_cost_eqn(b):
        return b.bare_erected_cost == pyunits.convert(
            b.diamine_volume * b.diamine_cost,
            pyunits.MUSD_2018
        )
    
    fs.RPB_cost.eng_fee = Param(initialize = 0)
    fs.RPB_cost.process_conting = Param(initialize = 0)
    fs.RPB_cost.project_conting = Param(initialize = 0)

    fs.RPB_cost.total_plant_cost = Var(initialize = 1,
                                       units = pyunits.MUSD_2018,
                                       bounds = (0, 10000))
    @fs.RPB_cost.Constraint()
    def total_plant_cost_eqn(b):
        return b.total_plant_cost == (b.tetraamine_bare_erected_cost+b.diamaine_bare_erected_cost) * (
                1 + b.eng_fee + b.process_conting
            ) * (1 + b.project_conting)
    
    fs.tonneCO2cap = Var(initialize = 1,
                         units = pyunits.tonne/pyunits.yr,
                         bounds = (0,None),
                         doc = "Yearly capture rate in metric tonnes")
    
    @fs.Constraint()
    def tonneCO2cap_eqn(b):
        co2_stream = b.RPB.des_gas_outlet
        co2_flow_mol = co2_stream.flow_mol[0] * co2_stream.mole_frac_comp[0,"CO2"]
        MW = 0.044 * pyunits.kg / pyunits.mol
        return b.tonneCO2cap == pyunits.convert(co2_flow_mol * MW,
                                                pyunits.tonne/pyunits.yr)
    
    resources = ['des_steam', 'ads_cw']

    fs.des_steam_rate = Var(fs.time,
                            initialize = 1,
                            units = pyunits.kW,
                            bounds = (0,None))
    fs.ads_cw_rate = Var(fs.time,
                         initialize = 1,
                         units = pyunits.kW,
                         bounds = (0,None))
    rates = [fs.des_steam_rate, fs.ads_cw_rate]
    
    @fs.Constraint(fs.time)
    def des_steam_rate_eqn(b, t):
        return b.des_steam_rate[t] == pyunits.convert(b.RPB.total_thermal_energy[t],
                                              pyunits.kW)
    @fs.Constraint(fs.time)
    def ads_cw_rate_eqn(b, t):
        return b.ads_cw_rate[t] == pyunits.convert(b.RPB.ads.Q_ghx_tot_kW[t],
                                           pyunits.kW)
    
    prices = {
        'des_steam': 3.33 * pyunits.USD_2018 / pyunits.GJ,
        'ads_cw': 0.354 * pyunits.USD_2018 / pyunits.GJ,
    }

    fs.costing = QGESSCosting()
    fs.costing.build_process_costs(
        total_plant_cost=fs.RPB_cost.total_plant_cost,
        labor_rate=0, #38.5,
        labor_burden=0,#30,
        operators_per_shift=2,
        tech=6,
        fixed_OM=True,
        variable_OM=True,
        resources=resources,
        rates=rates,
        prices=prices,
        # transport_cost=fs.costing.transport_cost,
        tonne_CO2_capture=fs.tonneCO2cap,
        CE_index_year="2018",
    )
    @fs.Constraint()
    def total_TPC_eqn(b):
        return b.costing.total_TPC == pyunits.convert(b.RPB_cost.total_plant_cost,
                                                      pyunits.MUSD_2018)
    fs.costing.total_TPC.unfix()

    capacity_factor = fs.costing.capacity_factor
    capital_cost_orig = fs.costing.annualized_cost

    capital_cost = capital_cost_orig
    fixed_cost = fs.costing.total_fixed_OM_cost
    variable_cost = fs.costing.total_variable_OM_cost[0] * capacity_factor * pyunits.a

    capturecost = capital_cost + fixed_cost + variable_cost
    tonnescaptured = fs.costing.tonne_CO2_capture

    @fs.costing.Expression()
    def LCOC(b):
        return pyunits.convert(capturecost, pyunits.USD_2018) / tonnescaptured

    