import microbial_evolution
import time

def main():
    # Setup simulation parameters.
    params = {
        # seed.
        "seed" : int(time.time()), 
        
        # Each epoch is so many h.
        "epoch" : 1,
        
        # Number of epochs.
        "epochs" : 100,

        # Volume in L.
        "volume" : 1e-6,
        
        # Total Phosphorous concentration in nmol-P per L.
        "P_tot" : 5000,
        
        # Initial virus population size.
        "V_pop" : 1000, 
        
        # Initial host population size.
        "H_pop" : 1000, 

        # Maximum mass of a host cell as a multiple of its nmol-P value.
        "mass_max" : 1, 
        "nmol_P_max" : 1.6e-8, 

        # Nutrient affinity of species 0 and strain 0, in L per h per nmol-P.
        "alpha00" : 0.4,  

        # Maximum growth of species 0 and strain 0, per h.
        "mu_max00" : 0.16, 

        # Other losses per h.
        "other_loss_type" : 0,       # mortality (0), metabolic (1), or both (2).
        "other_loss_rate" : 1.4e-2, 
        
        # Virus mutation (species, strain) probability.
        "V_mutation_prob" : [0.001, 0.01], 

        # Virus mutation (species, strain) width. 10% of beta for species and strain.
        "V_mutation_std" : [1.5e-11, 1.5e-11], 
        
        # Host mutation (species, strain) probability.
        "H_mutation_prob" : [0.005, 0.05], 

        # Host mutation (species, strain) width. 10% of alpha for species and mu_max for strain.
        "H_mutation_std" : [0.04, 0.016], 
        
        # Number of viruses produced per infection.
        "burst_size" : 5, 

        # Effective adsorption coefficient of species 0 and strain 0, in L per h per individual.
        "beta00" : 1.5e-10, 

        # Infection efficiency.
        "memory" : 0.7, 

        # Virus decay rate per h.
        "decay_rate" : 1.4e-2, 
    }

    # Convert units to correspond to simulated volume and time, using normalized cell mass.
    params["P_tot"] =  params["P_tot"] * params["volume"] / params["nmol_P_max"]
    params["alpha00"] =  params["alpha00"] / params["volume"] * params["epoch"] * params["nmol_P_max"]
    params["mu_max00"] = params["mu_max00"] * params["epoch"]
    params["other_loss_rate"] = params["other_loss_rate"] * params["epoch"]
    params["beta00"] =  params["beta00"] / params["volume"] * params["epoch"]
    params["dacay_rate"] = params["decay_rate"] * params["epoch"]
    params["V_mutation_std"][0] = params["V_mutation_std"][0] / params["volume"] * params["epoch"]
    params["V_mutation_std"][1] = params["V_mutation_std"][1] / params["volume"] * params["epoch"]
    params["H_mutation_std"][0] = params["H_mutation_std"][0] * params["epoch"]
    params["H_mutation_std"][1] = params["H_mutation_std"][1] * params["epoch"]

    # Run simulation.
    microbial_evolution.run(params)

if __name__ == "__main__":
    main()




