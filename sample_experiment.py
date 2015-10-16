import microbial_evolution

def main():
    params = {
        # Time step in h.
        "time_step" : 1,
        
        # Total number of time steps.
        "iterations" : 1000,

        # Volume in L.
        "volume" : 1e-6,
        
        # Total Phosphorous concentration in nmol-P per L
        "P_tot" : 500,
        
        # Initial virus population size.
        "V_init" : 2000, 
        
        # Initial host population size.
        "H_init" : 1000, 

        # Maximum mass of a host cell as a multiple of its nmol-P value.
        "mass_max" : 1, 
        "nmol_p_max" : 1.6e-8, 
    
        # Nutrient affinity in L per h per individual.
        "alpha" : 6.4e-9,  

        # Maximum growth per h.
        "mu_max0" : 0.16, 

        # Other losses per h.
        "other_loss_rate" : 1.4e-2, 
        
        # Virus mutation probability.
        "V_mutation_prob" : 0.01, 

        # Virus mutation width.
        "V_mutation_std" : 1.5e-11, 
        
        # Host mutation probability.
        "H_mutation_prob" : 0.01, 

        # Host mutation width.
        "H_mutation_std" : 0.01, 
        
        # Number of viruses produced per infection.
        "burst_size" : 5, 

        # Effective adsorption coefficient in L per h per individual.
        "beta0" : 1.5e-10, 

        # Infection efficiency.
        "memory" : 0.9, 

        # Virus decay rate per h.
        "decay_rate" : 1.4e-2, 
    }

    microbial_evolution.run(params)

if __name__ == "__main__":
    main()




