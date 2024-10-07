import json

# import pytest
from numpy import number
import pandas as pd
from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units,
    value,
    Var,
    Expression,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigValue

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom


from idaes.models.properties import iapws95


from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)
import pyomo.environ as pyo

def build_RPB_costing(fs):
    with open('costing/dac_eb_costing_data.json', "r") as f:
        costing_data = json.load(f)

    cost_params = costing_data
    fs.costing = QGESSCosting()
    CE_index_year = "2018"
    CE_index_units = getattr(units, "MUSD_" + CE_index_year)

    fs.number_of_units = Var(initialize=1, bounds=(0, None), doc="Number of units")
    fs.number_of_units.fix(1)


    # Parameters
    RPB_steam_power = fs.RPB.total_thermal_energy[0]
    fs.steam_flow_mass = Var(fs.time, initialize=0, units=units.lb / units.h)
    @fs.Constraint(fs.time)
    def steam_flow_mass_eqn(b, t):
        steam_heat_of_vaporization = 2257.92 * units.kJ / units.kg
        return b.steam_flow_mass[t] == units.convert(b.RPB.total_thermal_energy[t] / steam_heat_of_vaporization,
                                                     units.lb / units.h)

    # fs.raw_water_system = UnitModelBlock()
    # fs.raw_water_system.costing = UnitModelCostingBlock(
    #     flowsheet_costing_block=fs.costing,
    #     costing_method=QGESSCostingData.get_PP_costing,
    #     costing_method_arguments={
    #         "cost_accounts": ["3.2", "3.4", "9.5", "14.6"],
    #         "scaled_param": fs.raw_water_withdrawal[0] * fs.number_of_units,
    #         "tech": 8,
    #         "ccs": "B",
    #         "additional_costing_params": cost_params,
    #     },
    # )

    # Electric Boiler 15.9
    fs.electric_boiler = UnitModelBlock()
    fs.electric_boiler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.9"],
            "scaled_param": fs.steam_flow_mass[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": cost_params,
        },
    )

    fs.vessels = UnitModelBlock()
    fs.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.1"],
            "scaled_param": units.convert(fs.RPB.ads.vol_tot, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": cost_params,
        },
    )

    # Building plant cost
    distributed_systems = [
        fs.vessels.costing,
        # fs.duct_dampers.costing,
        # fs.feed_fans.costing,
        # fs.desorption_gas_handling.costing,
        # fs.controls_equipment.costing,
    ]

    centralized_TPCs = []

    for b in fs.costing._registered_unit_costing:
        if b not in distributed_systems:
            for key in b.total_plant_cost.keys():
                centralized_TPCs.append(b.total_plant_cost[key])

    fs.costing.total_TPC = Var(
        initialize=100,
        bounds=(0, 1e4),
        units=CE_index_units,
        doc="total TPC in $MM",
    )

    @fs.costing.Constraint()
    def total_TPC_eqn(b):
        return b.total_TPC == (
            sum(centralized_TPCs) +
            sum(sum(block.total_plant_cost[key] for key in block.total_plant_cost.keys()) for block in distributed_systems) * fs.number_of_units
        )


    resources = [
        "sorbent"
    ]

    capacity_factor = 0.85

    @fs.costing.Expression(fs.time)
    def sorbent_rate(b, t):
        RPB = b.parent_block().RPB
        sorbent_lifespan = 2 * units.year
        return (
            units.convert(RPB.ads.vol_solids_tot, units.ft**3)
            / sorbent_lifespan
            * b.parent_block().number_of_units
        )
    
    rates = [
        fs.costing.sorbent_rate
    ]

    prices = {
        "sorbent": 200 * units.USD_2018 / units.ft**3,

    }

    fs.costing.tonne_CO2_capture = Var(
        initialize=1,
        units=units.tonne / units.yr,
        bounds=(0, None),
        doc="Yearly capture rate in metric tonnes"
    )
    
    @fs.costing.Constraint()
    def tonne_CO2_capture_eqn(b):
        co2_stream = b.parent_block().RPB.des_gas_outlet
        co2_flow_mol = co2_stream.flow_mol[0] * co2_stream.mole_frac_comp[0,"CO2"]
        MW = 0.044 * units.kg / units.mol
        return b.tonne_CO2_capture == units.convert(co2_flow_mol * MW,
                                                units.tonne/units.yr)

    # Initialize costing
    fs.costing.build_process_costs(
        net_power = None,
        total_plant_cost = True,
        labor_rate = 38.50,
        labor_burden=30,
        operators_per_shift=0,
        tech=6,
        fixed_OM=True,
        # arguments related owners costs
        variable_OM=True,
        capacity_factor=capacity_factor,
        land_cost=None,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
        waste=None,
        chemicals=None,
        tonne_CO2_capture=fs.costing.tonne_CO2_capture,
    )


    fs.costing.costing_initialization()

    @fs.costing.Expression()
    def LCOC(b):
        return units.convert(b.costing.cost_of_capture, units.USD / units.tonne)