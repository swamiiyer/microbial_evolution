import microbial_evolution, sys, time

def main():
    # Setup simulation parameters.
    params = {
        # Number of epochs in h.
        "epochs" : 300,

        # Simulation volume in L.
        "volume" : 1e-6,
        
        # Total Phosphorous concentration in nmol-P per L.
        "P_tot" : 500,
        
        # Initial virus population size.
        "V_pop" : 500, 
        
        # Initial host population size.
        "H_pop" : 50, 

        # nmol-P per host cell.
        "nmol_P_max" : 1.6e-8, 

        # Nutrient affinity of initial host strain in L per h per nmol-P.
        "alpha" : 0.4,  

        # Maximum growth rate of initial host strain per h.
        "mu_max" : 0.16, 

        # Mortality probability of host per h. 
        "mortality_rate" : 1.4e-2, 

        # Metabolic loss rate (in units of individuals) of host per h. 
        "metabolic_loss_rate" : 1.4e-2, 
        
        # Virus mutation probability.
        "V_mutation_prob" : 0.01, 

        # Host mutation probability.
        "H_mutation_prob" : 0.05, 
        
        # Number of viruses produced per infection.
        "burst_size" : 5, 

        # Effective adsorption coefficient of initial virus strain in L
        # per h per individual.
        "beta" : 1.5e-10, 

        # Infection efficiency.
        "memory" : 0.7, 

        # Virus decay rate per h.
        "decay_rate" : 1.4e-2, 
    }

    params["V_mutation_std"] = 0.1 * params["beta"]
    params["H_mutation_std"] = 0.1 * params["mu_max"]
    
    # Convert units to correspond to simulated volume and time, using
    # normalized cell mass.
    params["P_tot"] =  params["P_tot"] * params["volume"] / params["nmol_P_max"]
    params["alpha"] =  params["alpha"] / params["volume"] * params["nmol_P_max"]
    params["beta"] =  params["beta"] / params["volume"]
    params["H_mutation_std"] = params["H_mutation_std"] / params["volume"] * params["nmol_P_max"]
    params["V_mutation_std"] = params["V_mutation_std"] / params["volume"]

    params["V_bin_width"] = 0.2 * params["V_mutation_std"]
    params["H_bin_width"] = 0.2 * params["H_mutation_std"]
    
    # Run simulation.
    fname = sys.argv[0].split('/')[-1].replace('.py', '.pkl')
    microbial_evolution.run(params, fname)

if __name__ == "__main__":
    main()
