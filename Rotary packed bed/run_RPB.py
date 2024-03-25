import RPB_model
import pandas as pd
from idaes.core.util import to_json, from_json
import numpy as np
import pyomo.environ as pyo

def Remove_Pressure_Drop(b):
    # b.ads.R_dP.fix(0.001)
    # b.des.R_dP.fix(0.001)
    
    b.ads.F_in.fix()
    b.ads.P_in.fix()
    b.ads.P_out.unfix()
    
    b.des.F_in.fix()
    b.des.P_in.fix()
    b.des.P_out.unfix()

if __name__ == '__main__':
    
    has_pressure_drop = True
    RPB = RPB_model.full_model_creation(lean_temp_connection=True,
                                        configuration='counter-current')
    RPB_model.add_ads_inlet_comp_constraint(RPB.ads)
    

    y_in = {'CO2': 0.002201,
            'N2': 0.940705,
            'H2O': 0.05}
    y_in = {'CO2': 0.002,
            'N2': 0.948,
            'H2O': 0.05}
    
    y_in = {'CO2': 0.0391,
            'N2': 0.8768,
            'H2O': 0.0841}
    
    
    # @RPB.Expression()
    # def obj(RPB):
    #     return RPB.energy_requirement
    # RPB.objective = pyo.Objective(expr=RPB.obj)
    
    # RPB.ads.P_in.bounds = (0.99, 1.5)
    RPB.ads.P_in.setub(1.5)
    RPB.ads.P_in.setlb(0.99)
    # RPB.ads.P.bounds = (0.99, 1.5)
    RPB.ads.P.setub(1.5)
    # RPB.ads.P.setlb(0.99)
    
    if not has_pressure_drop:
        Remove_Pressure_Drop(RPB)
        
    homotopy_points = np.linspace(0.2, 1, 5)
    RPB_model.init_routine_2(RPB)#, homotopy_points)

    # from_json(RPB, fname="low_co2_initialization_no_dP.json.gz", gz=True)
    

    
    # RPB_model.solve_model(RPB,)# optarg={'max_iter': 0})
    Results = RPB_model.report(RPB)
    # Results_Inlet_Loading = RPB_model.report_loading(RPB)
    
    # co2_comp_list = np.linspace(0.04, 0.004, 10)
    # n2_comp_list = np.linspace(0.87, 0.94, 10)
    # for i in range(len(co2_comp_list)):
    #     print(i)
    #     RPB.ads.y_in["CO2"].fix(co2_comp_list[i])
    #     RPB.ads.y_in["N2"].fix(n2_comp_list[i])
    #     RPB_model.solve_model(RPB)
        
    
    

        
    