import microbial_evolution

def main():
    params = {
        # Time step in h.
        "time_step" : 1,
        
        # Total number of time steps.
        "iterations" : 1000,

        # Volume in L.
        "volume" : 1e-6,
        
        # Total Phosphorous concentration in nmol-P per L (*)
        "P_tot" : 5000,
        
        # Initial virus population size.
        "V_init" : 2000, 
        
        # Initial host population size.
        "H_init" : 100, 

        # Maximum mass of a host cell as a multiple of its nmol-P value.
        "mass_max" : 1, 
        "nmol_p_max" : 1.6e-8, 
    
        # Nutrient affinity in L per h per nmol-P (*)
        "alpha00" : 0.4,  

        # Maximum growth per h (*)
        "mu_max00" : 0.16, 

        # Other losses per h.
        "other_loss_type" : 0, # mortality (0), metabolic (1), or both (2).
        "other_loss_rate" : 1.4e-2,  # (*)
        
        # Virus mutation (species, strain) probability.
        "V_mutation_prob" : [0.001, 0.01], 

        # Virus mutation (species, strain) width.
        "V_mutation_std" : [1.5e-11, 1.5e-11], 
        
        # Host mutation (species, strain) probability.
        "H_mutation_prob" : [0.001, 0.01], 

        # Host mutation (species, strain) width.
        "H_mutation_std" : [0.01, 0.01], 
        
        # Number of viruses produced per infection.
        "burst_size" : 5, 

        # Effective adsorption coefficient in L per h per individual.
        "beta00" : 1.5e-10, 

        # Infection efficiency.
        "memory" : 0.9, 

        # Virus decay rate per h.
        "decay_rate" : 1.4e-2, 
    }

    # Unit conversions.
    params["P_tot"] =  params["P_tot"] / params["volume"] / params["nmol_p_max"]
    params["alpha00"] =  params["alpha00"] / params["volume"] * params["time_step"] * params["nmol_p_max"]
    params["mu_max00"] = params["mu_max00"] * params["time_step"]
    params["other_loss_rate"] = params["other_loss_rate"] * params["time_step"]
    params["beta00"] =  params["beta00"] / params["volume"] * params["time_step"]
    params["dacay_rate"] = params["decay_rate"] * params["time_step"]

    microbial_evolution.run(params)

if __name__ == "__main__":
    main()




