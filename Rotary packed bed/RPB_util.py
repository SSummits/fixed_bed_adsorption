import re
import pandas as pd

def set_bounds(fs):
    fs.flue_gas_in.pressure.setlb(1e5)
    fs.flue_gas_in.pressure.setub(2.5e5)

    fs.steam_sweep_feed.pressure.setlb(1e5)
    fs.steam_sweep_feed.pressure.setub(2.5e5)

    fs.RPB.ads.Tx.setlb(25+273.15)
    fs.RPB.ads.Tx.setub(100+273.15)

    fs.RPB.des.Tx.setlb(100+273.15)
    fs.RPB.des.Tx.setub(170+273.15)

    fs.RPB.w_rpm.setlb(1e-5)

    fs.RPB.L.setub(40)

    # fs.RPB.D.setub(10)
    # fs.RPB.D.setlb(1)

    fs.RPB.ads.pressure.setub(2.5e5)
    fs.RPB.ads.pressure_inlet.setub(2.5e5)
    fs.RPB.des.pressure.setub(2.5e5)
    fs.RPB.des.pressure_inlet.setub(2.5e5)

    fs.RPB.ads.vel.setub(100)
    fs.RPB.des.vel.setub(100)

def make_results_table(fs):
    res = pd.DataFrame()

    RPB_res = fs.RPB.report_custom()
    for i in RPB_res.index:
        res.at[i, 'Value'] = RPB_res['Value'][i]

    res.at['LCOC', 'Value'] = fs.costing.LCOC()
    res.at['Total Plant Cost', 'Value'] = fs.costing.total_TPC()
    res.at['Electricity Cost', 'Value'] = fs.costing.variable_operating_costs[0,'electricity']()
    res.at['Diamine Cost', 'Value'] = fs.costing.variable_operating_costs[0,'diamine']()
    res.at['Tetraamine Cost', 'Value'] = fs.costing.variable_operating_costs[0,'tetraamine']()

    return res