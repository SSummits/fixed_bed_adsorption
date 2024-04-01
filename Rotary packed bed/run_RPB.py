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
    
def set_polishing_bounds(sections):
    for section in sections:
        section.P_in.setub(10)
        section.L.setub(40)
        if str(section) == 'ads':
            section.Tx.setlb(298)
            section.Tx.setub(368)
        elif str(section) == 'des':
            section.Tx.setlb(373)
            section.Tx.setub(433)
        for z in section.z:
            for o in section.o:
                section.dPdz[z,o].setub(5)
                section.dPdz[z,o].setlb(-5)
                section.C_tot.setub(250)
                section.P[z,o].setub(10)
                section.vel[z,o].setub(15)

if __name__ == '__main__':
    
    has_pressure_drop = False
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
    
    RPB_model.scale_model(RPB.ads, gas_flow_direction=1, mode="adsorption")
    RPB_model.scale_model(RPB.des, gas_flow_direction=-1, mode="desorption")
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
        
    homotopy_points = np.linspace(0.1, 1, 10)
    # RPB_model.init_routine_1(RPB, homotopy_points)
    

    from_json(RPB, fname="polishing step optimized solution 031824.json.gz", gz=True)
    
    set_polishing_bounds([RPB.ads, RPB.des])
    
    for z in RPB.ads.z:
        for o in RPB.ads.o:
            RPB.ads.dPdz[z,o].setlb(-20)

    
    # opt_res = RPB_model.solve_model(RPB)#, optarg={'max_iter': 0})
    Results = RPB_model.report(RPB)
    # Results_Inlet_Loading = RPB_model.report_loading(RPB)
    
    # co2_comp_list = np.linspace(0.04, 0.004, 10)
    # n2_comp_list = np.linspace(0.87, 0.94, 10)
    # for i in range(len(co2_comp_list)):
    #     print(i)
    #     RPB.ads.y_in["CO2"].fix(co2_comp_list[i])
    #     RPB.ads.y_in["N2"].fix(n2_comp_list[i])
    #     RPB_model.solve_model(RPB)
        
    
    

        
    