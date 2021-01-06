import math, microbial_evolution, random, sys, time


def main(args):
    params = {}

    # Seed for the random number generator.
    params["seed"] = int(time.time())

    # Number of epochs (h).
    params["epochs"] = 720

    # Simulation volume (L).
    params["volume"] = 1e-6

    # Total Phosphorous concentration (nmol-P per L).
    params["pTotal"] = 71036

    # Phosphorous concentration in a single host (nmol-P).
    params["pMax"] = 8.3e-5

    # Washout rate (per h).
    params["washout"] = 0.2

    # Initial host genotype.
    params["hG0"] = 0.1

    # Standard deviation for the mutations in host genotype.
    params["hGSigma"] = 0.005

    # Initial virus genotype.
    params["vG0"] = 0.1

    # Standard deviation for the mutations in virus genotype.
    params["vGSigma"] = 0.005

    # Initial host population size.
    params["h0"] = 46

    # Initial virus population size.
    params["v0"] = 810

    # Maximum nutrient affinity of host (L per h per individual).
    params["alphaMax"] = 1.6e-7

    # Maximum growth rate of host (per h).
    params["muMax"] = 0.738

    # Maximum adsorption rate of virus (L per h per individual).
    params["betaMax"] = 6.2e-11

    # Number of viruses produced per infection.
    params["burst"] = 10

    # Host mutation probability.
    params["hMu"] = 0.006

    # Virus mutation probability.
    params["vMu"] = 0.006

    # Host mortality rate (per h).
    params["mortality"] = 1.4e-2

    # Host metabolic loss rate (per h). 
    params["metabolicLoss"] = 1.4e-2

    # Virus decay rate (per h).
    params["decay"] = 1.4e-2

    # Express parameters in terms of simulated volume and normalized cell mass.
    params["pTotal"] *= params["volume"] / params["pMax"]
    params["alphaMax"] /= params["volume"]
    params["betaMax"] /= params["volume"]

    # Nutrient affinity of host as a linear function of its genotype.
    def alphaLinear(hG):
        return params["alphaMax"] / params["hG0"] * hG

    # Nutrient affinity of host as a Gaussian function of its genotype.
    def alphaGaussian(hG):
        return params["alphaMax"] * math.exp(-(hG - params["hG0"]) ** 2)

    # Adsorption coefficient as a function of nutrient affinity, and host and virus genotypes.
    def betaTradeoff(alpha, hG, vG):
        return params["betaMax"] / params["alphaMax"] * alpha(hG) * math.exp(-(vG - hG) ** 2)

    # Adsorption coefficient as a random value.
    def betaRandom(alpha, hG, vG):
        return random.uniform(0, params["betaMax"] / params["alphaMax"] * alpha(hG))

    # Run simulation.
    microbial_evolution.run(params, alphaGaussian, betaTradeoff, args[1])


if __name__ == "__main__":
    main(sys.argv)
