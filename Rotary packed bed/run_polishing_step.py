# Import RPB model along with other utility functions

from attr import mutable
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Feed, Product
from RPB_model import RotaryPackedBed
import RPB_util
# from RPB_costing import RPB_Polishing_Costing

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TransformationFactory,
    Reference,
    units as pyunits,
    Param,
)

import idaes.core.util as iutil
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state

from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)

from pyomo.network import Arc

from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import numpy as np
import pandas as pd
import solver_methods

from costing.rpb_costing import build_RPB_costing

# create Flowsheet block
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic = False)


# create gas phase properties block
flue_species={"H2O", "CO2", "N2"}
prop_config = get_prop(flue_species, ["Vap"], eos=EosType.IDEAL)
prop_config["state_bounds"]["pressure"] = (0.99*1e5,1.02*1e5,2.5*1e5, pyunits.Pa)
prop_config["state_bounds"]["temperature"] = (25+273.15,90+273.15,180+273.15, pyunits.K)

m.fs.gas_props = GenericParameterBlock(
    **prop_config,
    doc = "Flue gas properties",
)

m.fs.gas_props.set_default_scaling("temperature", 1e-2)
m.fs.gas_props.set_default_scaling("pressure", 1e-4)


# create feed and product blocks
m.fs.flue_gas_in = Feed(property_package = m.fs.gas_props)
m.fs.flue_gas_out = Product(property_package = m.fs.gas_props)
m.fs.steam_sweep_feed = Feed(property_package = m.fs.gas_props)
m.fs.regeneration_prod = Product(property_package = m.fs.gas_props)


# limited discretization, much faster
m.fs.RPB = RotaryPackedBed(
    property_package = m.fs.gas_props,
    z_init_points = (0.01,0.99),
    o_init_points = (0.01,0.99),
    mixed_sorbent_list = ['Tetraamine', 'Diamine'],
)

# increased number of discretization points, lower mass balance error
# z_init_points=tuple(np.geomspace(0.01, 0.5, 9)[:-1]) + tuple((1 - np.geomspace(0.01, 0.5, 9))[::-1])
# o_init_points=tuple(np.geomspace(0.005, 0.1, 8)) + tuple(np.linspace(0.1, 0.995, 10)[1:])
# z_nfe=20
# o_nfe=20
# m.fs.RPB = RotaryPackedBed(
#     property_package = m.fs.gas_props,
#     z_init_points=z_init_points,
#     o_init_points=o_init_points,
#     z_nfe=z_nfe,
#     o_nfe=o_nfe,
# )
# z_init_points=tuple(np.geomspace(0.01, 0.5, 7)[:-1]) + tuple((1 - np.geomspace(0.01, 0.5, 7))[::-1])
# o_init_points=tuple(np.geomspace(0.005, 0.1, 6)) + tuple(np.linspace(0.1, 0.995, 8)[1:])
# m.fs.RPB = RotaryPackedBed(
#     property_package = m.fs.gas_props,
#     z_init_points=z_init_points,
#     o_init_points=o_init_points,
# )

# add stream connections
m.fs.s_flue_gas = Arc(source=m.fs.flue_gas_in.outlet, destination=m.fs.RPB.ads_gas_inlet)
m.fs.s_cleaned_flue_gas = Arc(source=m.fs.RPB.ads_gas_outlet, destination=m.fs.flue_gas_out.inlet)
m.fs.s_steam_feed = Arc(source=m.fs.steam_sweep_feed.outlet, destination=m.fs.RPB.des_gas_inlet)
m.fs.s_regeneration_prod = Arc(source=m.fs.RPB.des_gas_outlet, destination=m.fs.regeneration_prod.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)


# fix state variables in feed and product blocks
# ads side
m.fs.flue_gas_in.pressure.fix(1.5*1e5)
m.fs.flue_gas_in.temperature.fix(90+273.15)
m.fs.flue_gas_out.pressure.fix(1.03*1e5)
m.fs.flue_gas_in.mole_frac_comp[0,"CO2"].fix(0.004)
m.fs.flue_gas_in.mole_frac_comp[0,"H2O"].fix(0.07)
m.fs.flue_gas_in.mole_frac_comp[0,"N2"].fix(1-0.004-0.07)

#des side
m.fs.steam_sweep_feed.pressure.fix(1.03*1e5)
m.fs.steam_sweep_feed.temperature.fix(120+273.15)
m.fs.regeneration_prod.pressure.fix(1.01325*1e5)
m.fs.steam_sweep_feed.mole_frac_comp[0,"CO2"].fix(1e-5)
m.fs.steam_sweep_feed.mole_frac_comp[0,"N2"].fix(1e-3)
m.fs.steam_sweep_feed.mole_frac_comp[0,"H2O"].fix(1-1e-5-1e-3)

# fix design variables of the RPB
m.fs.RPB.ads.Tx.fix(348.04)
m.fs.RPB.des.Tx.fix(433)
m.fs.RPB.w_rpm.fix(1)
m.fs.RPB.L.fix(6)
m.fs.RPB.ads.theta.fix(0.60)
m.fs.RPB.des.theta = 1-0.60

# Fix sorbent composition
for z in m.fs.RPB.z:
    if z < 0.5:
        m.fs.RPB.sorbent_weight[z, 'TA'].fix(1e-8)
        m.fs.RPB.sorbent_weight[z, 'DA'].fix(1)
    else:
        m.fs.RPB.sorbent_weight[z, 'TA'].fix(1)
        m.fs.RPB.sorbent_weight[z, 'DA'].fix(1e-8)


# initialize feed and product blocks
m.fs.flue_gas_in.initialize()
m.fs.flue_gas_out.initialize()
m.fs.steam_sweep_feed.initialize()
m.fs.regeneration_prod.initialize()

# propagate feed and product blocks (for initial RPB guesses)
propagate_state(arc = m.fs.s_flue_gas, direction="forward")
propagate_state(arc = m.fs.s_steam_feed, direction="forward")
propagate_state(arc = m.fs.s_cleaned_flue_gas, direction="backward")
propagate_state(arc = m.fs.s_regeneration_prod, direction="backward")

# Initialize RPB
optarg = {
    # "halt_on_ampl_error": "yes",
    "max_iter": 5000,
    "bound_push": 1e-22,
    # "mu_init": 1e-3,
    "nlp_scaling_method": "user-scaling",
    'halt_on_ampl_error': 'yes',
}
init_points = [1e-5,1e-3,1e-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
# init_points = [1]
# init_points = [1e-5, 0.25, 0.5, 0.75, 1]
RPB_util.set_bounds(m.fs)

# m.fs.RPB.initialize(outlvl=idaeslog.DEBUG, optarg=optarg, initialization_points=init_points)

iutil.from_json(m.fs, fname='layered_polish_90_init5.json.gz')
# iutil.from_json(m.fs, fname='HighNfe_TA_NGCC_init.json.gz')
# iutil.from_json(m, fname='json_files/archive_polish/90_PCC_80_RPB.json.gz')

design_variables = [
    m.fs.flue_gas_in.pressure,
    m.fs.steam_sweep_feed.pressure,
    m.fs.RPB.ads.Tx,
    m.fs.RPB.des.Tx,
    m.fs.RPB.w_rpm,
    m.fs.RPB.ads.theta,
    m.fs.RPB.L,
    # m.fs.RPB.D,
]
m.fs.RPB.ads.CO2_capture.unfix()
m.fs.RPB.ads.flow_mol_inlet.unfix()
for v in design_variables:
    v.fix()

# full solve with IPOPT
Solver = get_solver("ipopt", optarg)
Solver.solve(m, tee=True).write()

build_RPB_costing(m.fs)

Solver.solve(m, tee=True)

m.fs.RPB.ads.flow_mol_inlet.fix()



m.fs.alpha_obj = Param(initialize=0.1, mutable=True)
@m.fs.Objective()
def min_energy(b):
    # return b.alpha_obj * b.RPB.energy_requirement[0] - (1-b.alpha_obj)*b.RPB.productivity[0]
    return b.costing.LCOC

for v in design_variables:
    v.unfix()

RPB_util.set_bounds(m.fs)
# iutil.from_json(m, fname='RPB_init_try2.json.gz')


res_df = RPB_util.make_results_table(m.fs)
cap_init = m.fs.RPB.ads.CO2_capture[0]()
for cap in np.linspace(cap_init, 0.99, 10):
    m.fs.RPB.ads.CO2_capture.fix(cap)
    # Solver.solve(m, tee=True)
    # opt_res = solver_methods.NEOS_solver(m.fs)
    res_df = res_df.join(RPB_util.make_results_table(m.fs), rsuffix=f'_{cap}')
    # print(m.fs.RPB.report_custom())
    # print(f'================\n{m.fs.RPB.ads.F_in[0]()}\n{m.fs.RPB.des.Tx[0]()}\n================')


