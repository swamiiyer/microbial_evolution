import microbial_evolution
import time

def main():
    # Setup simulation parameters.
    params = {
        # Number of epochs in h.
        "epochs" : 1000,

        # Simulation volume in L.
        "volume" : 1e-6,
        
        # Total Phosphorous concentration in nmol-P per L.
        "P_tot" : 500,
        
        # Initial virus population size.
        "V_pop" : 0, 
        
        # Initial host population size.
        "H_pop" : 1000, 

        # nmol-P per host cell.
        "nmol_P_max" : 1.6e-8, 

        # Nutrient affinity of host strain 0 in L per h per nmol-P.
        "alpha0" : 0.4,  

        # Maximum growth rate of host strain 0 per h.
        "mu_max0" : 0.16, 

        # Mortality probability of host per h. 
        "mortality_rate" : 1.4e-2, 

        # Metabolic loss rate (in units of individuals) of host per h. 
        "metabolic_loss_rate" : 1.4e-2, 
        
        # Virus mutation probability.
        "V_mutation_prob" : 0.01, 

        # Virus mutation width (standard deviation); 10% of beta.
        "V_mutation_std" : 1.5e-11, 
        
        # Host mutation probability.
        "H_mutation_prob" : 0.05, 

        # Host mutation width (standard deviation); 10% of alpha.
        "H_mutation_std" : 0.016, 
        
        # Number of viruses produced per infection.
        "burst_size" : 5, 

        # Effective adsorption coefficient of virus strain 0 in L per h per individual.
        "beta0" : 1.5e-10, 

        # Infection efficiency.
        "memory" : 0.7, 

        # Virus decay rate per h.
        "decay_rate" : 1.4e-2, 
    }

    # Convert units to correspond to simulated volume and time, using normalized cell mass.
    params["P_tot"] =  params["P_tot"] * params["volume"] / params["nmol_P_max"]
    params["alpha0"] =  params["alpha0"] / params["volume"] * params["nmol_P_max"]
    params["beta0"] =  params["beta0"] / params["volume"]
    params["H_mutation_std"] = params["H_mutation_std"] / params["volume"] * params["nmol_P_max"]
    params["V_mutation_std"] = params["V_mutation_std"] / params["volume"]

    # Run simulation.
    microbial_evolution.run(params)

if __name__ == "__main__":
    main()




