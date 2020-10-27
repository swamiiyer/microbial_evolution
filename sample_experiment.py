import math, microbial_evolution, sys, time

def main(args):
    #
    # Define parameters.
    #

    params = {}

    # Seed for the random number generator.
    params["seed"] = int(time.time())

    # Number of epochs (h).
    params["epochs"] = 720

    # Simulation volume (L).
    params["volume"] = 1e-6

    # Total Phosphorous concentration (nmol-P per L).
    params["P_tot"] = 71036

    # Phosphorous concentration in a single host (nmol-P).
    params["nmol_P_max"] = 8.3e-5

    # Washout rate (per h).
    params["washout_rate"] = 0.2

    # Initial host genotype.
    params["H_genotype"] = 0.1

    # Standard deviation for the mutations in host genotype.
    params["H_genotype_std"] = 0.005

    # Initial virus genotype.
    params["V_genotype"] = 0.1

    # Standard deviation for the mutations in virus genotype.
    params["V_genotype_std"] = 0.005

    # Initial host population size.
    params["H_pop"] = 46

    # Initial virus population size.
    params["V_pop"] = 810

    # Nutrient affinity of host (L per h per individual).
    params["alpha"] = 1.6e-7

    # Maximum growth rate of host (per h).
    params["mu_max"] = 0.738

    # Growth factor of host (a function of its genotype).
    params["growth_factor"] = "lambda h: 0.8 + 0.4 * math.exp(-20 * (h - 0.5) ** 2)"

    # Maximum adsorption rate of virus (L per h per individual).
    params["beta"] = 6.2e-11

    # Number of viruses produced per infection.
    params["burst_size"] = 10

    # Specificity of virus.
    params["specificity"] = 200

    # Host mutation probability.
    params["H_mutation_prob"] = 0.006

    # Virus mutation probability.
    params["V_mutation_prob"] = 0.006

    # Host mortality rate (per h).
    params["mortality_rate"] = 1.4e-2

    # Host metabolic loss rate (per h). 
    params["metabolic_loss_rate"] = 1.4e-2

    # Virus decay rate (per h).
    params["decay_rate"] = 1.4e-2

    #
    # Express parameters in terms of simulated volume and normalized cell mass.
    #

    params["P_tot"] *= params["volume"] / params["nmol_P_max"]
    params["alpha"] /= params["volume"]
    params["beta"] /= params["volume"]

    #
    # Run simulation.
    #

    microbial_evolution.run(params, args[1])


if __name__ == "__main__":
    main(sys.argv)
