# Import RPB model along with other utility functions

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Feed, Product
from RPB_model import RotaryPackedBed
from RPB_costing import RPB_Polishing_Costing

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
import solver_methods

# create Flowsheet block
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic = False)


# create gas phase properties block
flue_species={"H2O", "CO2", "N2"}
prop_config = get_prop(flue_species, ["Vap"], eos=EosType.IDEAL)
prop_config["state_bounds"]["pressure"] = (0.99*1e5,1.02*1e5,1.5*1e5, pyunits.Pa)
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
# m.fs.RPB = RotaryPackedBed(
#     property_package = m.fs.gas_props,
#     z_init_points = (0.01,0.99),
#     o_init_points = (0.01,0.99),
# )

# increased number of discretization points, lower mass balance error
z_init_points=tuple(np.geomspace(0.01, 0.5, 9)[:-1]) + tuple((1 - np.geomspace(0.01, 0.5, 9))[::-1])
o_init_points=tuple(np.geomspace(0.005, 0.1, 8)) + tuple(np.linspace(0.1, 0.995, 10)[1:])
z_nfe=20
o_nfe=20
m.fs.RPB = RotaryPackedBed(
    property_package = m.fs.gas_props,
    z_init_points=z_init_points,
    o_init_points=o_init_points,
    z_nfe=z_nfe,
    o_nfe=o_nfe,
)

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
m.fs.flue_gas_out.pressure.fix(1.01325*1e5)
m.fs.flue_gas_in.mole_frac_comp[0,"CO2"].fix(0.0022)
m.fs.flue_gas_in.mole_frac_comp[0,"H2O"].fix(0.09)
m.fs.flue_gas_in.mole_frac_comp[0,"N2"].fix(1-0.0022-0.09)

#des side
m.fs.steam_sweep_feed.pressure.fix(1.015*1e5)
m.fs.steam_sweep_feed.temperature.fix(120+273.15)
m.fs.regeneration_prod.pressure.fix(1.01325*1e5)
m.fs.steam_sweep_feed.mole_frac_comp[0,"CO2"].fix(1e-5)
m.fs.steam_sweep_feed.mole_frac_comp[0,"N2"].fix(1e-3)
m.fs.steam_sweep_feed.mole_frac_comp[0,"H2O"].fix(1-1e-5-1e-3)

# fix design variables of the RPB
m.fs.RPB.ads.Tx.fix(298)
m.fs.RPB.des.Tx.fix(433)
m.fs.RPB.w_rpm.fix(0.002)
m.fs.RPB.L.fix(3.9)

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
    "max_iter": 1000,
    "bound_push": 1e-22,
    # "mu_init": 1e-3,
    "nlp_scaling_method": "user-scaling",
}
init_points = [1e-5,1e-3,1e-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
# init_points = [1e-5, 0.25, 0.5, 0.75, 1]
# m.fs.RPB.initialize(outlvl=idaeslog.DEBUG, optarg=optarg, initialization_points=init_points)

# iutil.from_json(m, fname='RPB_polishing_init.json.gz')
iutil.from_json(m, fname='95PCC_high_elems.json.gz')

# full solve with IPOPT
Solver = get_solver("ipopt", optarg)
# Solver.solve(m, tee=True).write()

RPB_Polishing_Costing(m.fs)
# Solver.solve(m, tee=True)

design_variables = [
    m.fs.flue_gas_in.pressure,
    m.fs.steam_sweep_feed.pressure,
    m.fs.RPB.ads.Tx,
    m.fs.RPB.des.Tx,
    m.fs.RPB.w_rpm,
    m.fs.RPB.ads.theta,
    m.fs.RPB.L,
    m.fs.RPB.D,
]

m.fs.RPB.ads.inlet_properties[0.0].flow_mol.fix()
m.fs.RPB.ads.Tx[0].setlb(273+25)

@m.fs.Objective()
def min_LCOC(b):
    # return b.RPB.energy_requirement[0]
    return b.costing.LCOC

for v in design_variables:
    v.unfix()
m.fs.RPB.ads.CO2_capture.fix(0.9)

# Solver.solve(m, tee=True)
# import numpy as np
# for cap in np.linspace(0.9, 0.99)
solver_methods.NEOS_solver(m.fs)