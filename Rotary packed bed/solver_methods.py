from pyomo.environ import SolverManagerFactory

def NEOS_solver(blk):
    
    solver_manager = SolverManagerFactory('neos')
    results = solver_manager.solve(
        blk,
        tee=True,
        # keedfiles=True,
        solver="conopt",
        # tmpdir="temp",
        options={'outlev': 2,
                'workfactor': 3,},
        # add_options=["gams_model.optfile=1;"],
    )
    results.write()
    
    return results