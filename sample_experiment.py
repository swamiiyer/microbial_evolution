import microbial_evolution, sys, time

def main():

    #
    # Define parameters.
    #
    
    params = {}

    # Seed for the random number generator.
    params["seed"] = int(time.time())
        
    # Number of epochs (h).
    params["epochs"] = 300

    # Simulation volume (L).
    params["volume"] = 1e-6
        
    # Total Phosphorous concentration (nmol-P per L).
    params["P_tot"] = 500
        
    # Initial virus population size.
    params["V_pop"] = 500
        
    # Adsorption coefficient of initial virus (L per h per individual).
    params["beta"] = 1.5e-10

    # Standard deviation of "beta".
    params["beta_std"] = 0.1 * params["beta"]
    
    # Number of viruses produced per infection.
    params["burst_size"] = 5

    # Infection efficiency.
    params["memory"] = 0.7

    # Virus mutation probability.
    params["V_mutation_prob"] = 0.01

    # Virus decay rate (per h).
    params["decay_rate"] = 1.4e-2

    # Initial host population size.
    params["H_pop"] = 50

    # Phosphorous concentration in a single host (nmol-P).
    params["nmol_P_max"] = 1.6e-8
    
    # Nutrient affinity of initial host (L per h per individual).
    params["alpha"] = 1.5e-9
    
    # Maximum growth rate of initial host (per h).
    params["mu_max"] = 0.16

    # Standard deviation of "mu_max".
    params["mu_max_std"] = 0.1 * params["mu_max"]

    # Host mutation probability.
    params["H_mutation_prob"] = 0.05
    
    # Host mortality rate (per h). 
    params["mortality_rate"] = 1.4e-2
    
    # Host metabolic loss rate (per h). 
    params["metabolic_loss_rate"] = 1.4e-2

    #
    # Express parameters in terms of simulated volume and normalized cell mass.
    #
    
    params["P_tot"]      *= params["volume"] / params["nmol_P_max"]
    params["beta"]       /= params["volume"]
    params["beta_std"]   /= params["volume"]
    params["alpha"]      /= params["volume"]
    params["mu_max_std"] /= params["volume"]

    #
    # Run simulation.
    #
    
    fname = sys.argv[0].split('/')[-1].replace('.py', '.pkl')
    microbial_evolution.run(params, fname)

if __name__ == "__main__":
    main()
