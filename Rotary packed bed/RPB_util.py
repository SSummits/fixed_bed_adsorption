def set_bounds(fs):
    fs.flue_gas_in.pressure.setlb(1.01325e5)
    fs.flue_gas_in.pressure.setub(4e5)

    fs.steam_sweep_feed.pressure.setlb(1.01325e5)
    fs.steam_sweep_feed.pressure.setub(3e5)

    fs.RPB.ads.Tx.setlb(298)
    fs.RPB.ads.Tx.setub(368)

    fs.RPB.des.Tx.setlb(373)
    fs.RPB.des.Tx.setub(433)

    fs.RPB.w_rpm.setlb(1e-5)
    fs.RPB.w_rpm.setub(5)

    fs.RPB.L.setlb(0.01)
    fs.RPB.L.setub(40)

    # fs.RPB.D.setub(10)
    # fs.RPB.D.setlb(1)

    fs.RPB.ads.pressure.setub(4e5)
    fs.RPB.ads.pressure_inlet.setub(4e5)
    fs.RPB.des.pressure.setub(3e5)
    fs.RPB.des.pressure_inlet.setub(3e5)

    fs.RPB.ads.vel.setub(100)
    fs.RPB.des.vel.setub(100)

    fs.RPB.ads.vel.setlb(0.01)