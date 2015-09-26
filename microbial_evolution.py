import json, random, sys

class Virus(object):
    """
    Represents a virus cell.
    """

    burst_size = None
    mutation_rate = None
    mutation_std = None
    decay_rate = None
    memory = None

    def __init__(self, strain, beta):
        self.strain = strain
        self.beta = beta

    def __str__(self):
        return "V(%d): beta = %e" %(self.strain, self.beta)

class Host(object):
    """
    Represents a host cell.
    """

    mutation_rate = None
    mutation_std = None
    lysis_loss_rate = None
    other_loss_rate = None
    alpha = None
    mass_max = None

    def __init__(self, strain, mass, mu_max):
        self.strain = strain
        self.mass = mass
        self.mu_max = mu_max

    def __str__(self):
        return "H(%d): mass = %e, mu_max = %e" %(self.strain, self.mass, \
                                                 self.mu_max)

def main(args):
    """
    Entry point.
    """
    
    if len(args) != 1:
        sys.exit("Usage: python microbial_evolution.py <params file>")

    # Parse parameters
    params = json.load((open(args[0], "r")))
    
    V_tot, H_tot = 0, 0 # total number of viruses and hosts
    nV, nH = 0, 0       # total number of virus and host strains

    P_tot = params["P_tot"]
    V_init = [Virus(0, params["beta0"]) for i in range(params["V_init"])]
    H_init = [Host(0, random.uniform(0.5 * params["mass_max"], 
                                     params["mass_max"]), params["mu_max0"]) \
              for i in range(params["H_init"])]
    population = V_init
    population.extend(H_init)
    random.shuffle(population)
    biomass = sum([host.mass for host in H_init])
    DIP = P_tot - biomass
    print P_tot, biomass, DIP + biomass

if __name__ == "__main__":
    main(sys.argv[1:])
