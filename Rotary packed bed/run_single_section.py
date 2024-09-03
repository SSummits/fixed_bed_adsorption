from RPB_model import *
import pandas as pd
from idaes.core.util import to_json, from_json
import numpy as np
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale

def split_report(blk):
    items = [
        blk.ads.L,
        blk.ads.D,
        blk.ads.w_rpm,
        blk.ads.theta,
        blk.des.theta,
        blk.ads.P_in,
        blk.ads.P_out,
        blk.ads.F_in,
        blk.ads.Tg_in,
        blk.ads.Tx,
        blk.des.P_in,
        blk.des.P_out,
        blk.des.F_in,
        blk.des.Tg_in,
        blk.des.Tx,
        blk.ads.CO2_capture,
    ]

    names = []
    values = []
    fixed = []
    lb = []
    ub = []
    docs = []
    for item in items:
        names.append(item.to_string())
        values.append(item())
        if item.ctype != pyo.Var:
            fixed.append("N/A")
            lb.append("N/A")
            ub.append("N/A")
        else:
            fixed.append(item.fixed)
            lb.append(item.lb)
            ub.append(item.ub)
        docs.append(item.doc)

    report_df = pd.DataFrame(
        data={
            "Value": values,
            "Doc": docs,
            "Fixed": fixed,
            "Lower Bound": lb,
            "Upper Bound": ub,
        },
        index=names,
    )

    indexed_items = [
        blk.ads.y_in,
        blk.ads.y_out,
    ]

    names = []
    values = []
    docs = []
    fixed = []
    lb = []
    ub = []
    for item in indexed_items:
        names += [item[k].to_string() for k in item.keys()]
        values += [item[k]() for k in item.keys()]
        docs += [item.doc for k in item.keys()]
        fixed += [item[k].fixed for k in item.keys()]
        lb += [item[k].lb for k in item.keys()]
        ub += [item[k].ub for k in item.keys()]

    report_indexed_df = pd.DataFrame(
        data={
            "Value": values,
            "Doc": docs,
            "Fixed": fixed,
            "Lower Bound": lb,
            "Upper Bound": ub,
        },
        index=names,
    )

    report_df = pd.concat([report_df, report_indexed_df])

    return report_df

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
                
def toggle_design_variables(blk, fix=False):
    variable_list = [
                     blk.ads.L,
                     blk.ads.w_rpm,
                     blk.ads.theta,
                     blk.ads.Tx,
                     blk.ads.P_in,
                     blk.des.P_in,
                     blk.des.P_out,
                     blk.des.Tx,
                     ]
    if fix:
        for v in variable_list:
            v.fix()
    else:
        for v in variable_list:
            v.unfix()

if __name__ == '__main__':
    
    has_pressure_drop = True
    RPB = pyo.ConcreteModel(())
    RPB.ads = RPB_model(mode="adsorption", gas_flow_direction=1, has_pressure_drop=has_pressure_drop)
    RPB.des = RPB_model(mode="desorption", gas_flow_direction=-1, has_pressure_drop=has_pressure_drop)
    set_polishing_bounds([RPB.ads, RPB.des])
    
    L = 8
    RPB.ads.L.fix(L)
    RPB.des.L.fix(L)
    
    # w_rpm = 1.5e-4
    # RPB.ads.w_rpm.fix(w_rpm)
    # RPB.des.w_rpm.fix(w_rpm)
    
    # theta_ads = 0.96
    # RPB.ads.theta.fix(theta_ads)
    # RPB.des.theta.fix(1 - theta_ads)
    
    # RPB.ads.P_in.fix(1.75)
    # RPB.ads.F_in.fix(2572.5)
    # RPB.ads.P_out.unfix()
    
    from_json(RPB.ads, fname='json_files/99_PCC_ads_init.json.gz')
    from_json(RPB.des, fname='json_files/temp_des.json.gz')
    # single_section_init(RPB.des)
    
    