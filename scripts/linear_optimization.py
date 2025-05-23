import pandas as pd
import geopandas as gpd
import numpy as np
import pulp
from snakemake.logging import logger
import sys
import time
import gurobipy

# Configure logging
logger.info("Starting linear optimization...")


def run_optimization(demand_points_path, potential_sites_path, distance_matrix_path, 
                     output_sites_path, num_facilities=snakemake.params.num_facilities):
    try:
        # Load data
        demand_points = gpd.read_file(demand_points_path)
        potential_sites = gpd.read_file(potential_sites_path)
        matrix = pd.read_csv(distance_matrix_path)

        # Set the ID column as the index if it exists in the dataframe
        if 'ID' in matrix.columns:
            matrix = matrix.set_index('ID')
            logger.info("Set ID column as index")
        
        # Get dimensions
        num_demand_points = matrix.shape[1]
        num_potential_sites = matrix.shape[0]
        
        I = range(num_demand_points)
        J = range(num_potential_sites)
        
        # Set population per demand point
        pop = demand_points['total_est_pop'].to_dict()
        
        # Identify existing and potential sites
        existing = [j for j in J if potential_sites.iloc[j]['status'] == "open"]
        new_sites = [j for j in J if potential_sites.iloc[j]['status'] == "potential"]
        
        # Create facility open/closed indicator
        y = {}
        for j in J:
            if j in existing:
                y[j] = 1  # existing facility is fixed as open
            else:
                y[j] = pulp.LpVariable(f"Facility_{j}", cat='Binary')
        
        # Create the model
        prob = pulp.LpProblem("P-Median_Problem", pulp.LpMinimize)

        # Set the solver
        solver = pulp.GUROBI(msg=True)

        # Assignment decision variables
        x = pulp.LpVariable.dicts("Assign", [(i, j) for i in I for j in J], cat='Binary')
        
        # Objective Function: Minimize total weighted distance
        prob += pulp.lpSum(pop[i]*(matrix.iloc[j, i]) * x[(i, j)] for i in I for j in J)
        
        # Constraint: Each demand point is assigned to exactly one facility
        for i in I:
            prob += pulp.lpSum(x[(i, j)] for j in J) == 1
        
        # Constraint: A demand point can only be assigned to an open facility
        for i in I:
            for j in J:
                if j in new_sites:
                    prob += x[(i, j)] <= y[j]
                else:
                    prob += x[(i, j)] <= 1
        
        # Constraint: Exactly p new facilities are opened among potential sites
        prob += pulp.lpSum(y[j] for j in new_sites) == num_facilities
        
        # Try to solve with Gurobi if available, otherwise use the default solver
        start_time = time.time()
        try:
            solver_status = prob.solve(solver)
        except:
            logger.warning("Gurobi solver not available. Using default solver.")
            solver_status = prob.solve()
        solve_time = time.time() - start_time
        
        # Check if the model was solved successfully
        if solver_status != 1:
            logger.error(f"Optimization failed with status: {pulp.LpStatus[prob.status]}")
            sys.exit(1)
            
        logger.info(f"Status: {pulp.LpStatus[prob.status]}")
        logger.info(f"Total Weighted Distance: {pulp.value(prob.objective)}")
        
        # Get new opened facilities
        opened_new = [j for j in new_sites if pulp.value(y[j]) == 1]
        opened_facilities = existing + opened_new
        
        # Create a dataframe of selected sites
        selected_sites = potential_sites.iloc[opened_facilities].copy()
        selected_sites.to_file(output_sites_path, driver="GPKG")

        optimality_gap = solver.model.MIPGap

        #save optimality gap to a file
        with open(snakemake.output.optimality_gap, "w") as f:
            f.write(f"{optimality_gap:.4f}")

        logger.info(f"Optimality gap: {optimality_gap}")
        logger.info(f"Solver time: {solve_time:.2f} seconds")
        return
        
    except Exception as e:
        logger.error(f"Optimization failed: {e}")
        raise

if __name__ == "__main__":
    # This section runs when executed directly, not needed when called from Snakemake
    run_optimization(
        demand_points_path=snakemake.input.demand_points,
        potential_sites_path=snakemake.input.potential_sites,
        distance_matrix_path=snakemake.input.distance_matrix,
        output_sites_path=snakemake.output.sites,
        num_facilities=snakemake.params.num_facilities,
    )
