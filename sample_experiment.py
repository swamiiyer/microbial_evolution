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

    # Initial host genotype (nutrient affinity).
    params["gAlpha0"] = 0.1

    # Standard deviation for the mutations in host genotype.
    params["hGSigma"] = 0.005

    # Initial virus genotype (memory).
    params["gMemory0"] = 0.05

    # Initial virus genotype (adsorption coefficient).
    params["gBeta0"] = 0.1

    # Standard deviation for the mutations in virus genotype.
    params["vGSigma"] = 0.005

    # Initial host population size.
    params["h0"] = 46

    # Initial virus population size.
    params["v0"] = 810

    # Proportionality constant for the nutrient affinity of host (L per h per individual).
    params["alpha"] = 1.6e-7

    # Maximum growth rate of host (per h).
    params["muMax"] = 0.738

    # Proportionality constant for the adsorption rate of virus (L per h per individual).
    params["beta"] = 6.2e-11

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

    # Parameters in terms of simulated volume and normalized cell mass.
    params["pTotal"] *= params["volume"] / params["pMax"]
    params["alpha"] /= params["volume"]
    params["beta"] /= params["volume"]

    # Returns nutrient affinity of the given host as a linear function of its genotype.
    def alpha(host):
        return params["alpha"] * host.gAlpha

    # Returns the compatibility score of the given host and virus as an exponential function of
    # their genotypes.
    def genotypicMatch(host, virus):
        if virus.gMemory == 0.0:
            return 1.0 if host.gAlpha == virus.gBeta else 0.0
        return math.exp(-(abs(host.gAlpha - virus.gBeta) / virus.gMemory))

    # Returns the compatibility score of the given host and virus as a random value.
    def randomMatch(host, virus):
        return random.random()

    # Returns the adsorption coefficient of the given virus as a function its genotype.
    def genotypicBeta(virus):
        if virus.gMemory == 0.0:
            return params["beta"] * virus.gBeta
        return params["beta"] * virus.gBeta / virus.gMemory

    # Returns the adsorption coefficient of the given virus as a random value.
    def randomBeta(virus):
        return random.uniform(0, params["beta"])

    # Run simulation.
    microbial_evolution.run(params, alpha, genotypicMatch, genotypicBeta, args[1])


if __name__ == "__main__":
    main(sys.argv)
