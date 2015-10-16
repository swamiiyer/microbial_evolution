import numpy, random, sys

class Virus(object):
    """
    Represents a virus cell.
    """

    def __init__(self, strain, beta):
        self.strain = strain
        self.beta = beta

    def __str__(self):
        return "V(%d): beta = %e" %(self.strain, self.beta)

class Host(object):
    """
    Represents a host cell.
    """

    def __init__(self, strain, mass, mu_max):
        self.strain = strain
        self.mass = mass
        self.mu_max = mu_max

    def __str__(self):
        return "H(%d): mass = %e, mu_max = %e" %(self.strain, self.mass, \
                                                 self.mu_max)

def run(params):
    """
    Entry point.
    """
    
    V_tot, H_tot = 0, 0 # total number of viruses and hosts
    nV, nH = 0, 0       # total number of virus and host strains

    P_tot = params["P_tot"]
    V_init = [Virus(0, params["beta0"]) for i in range(params["V_init"])]
    H_init = [Host(0, random.uniform(0.5 * params["mass_max"], 
                                     params["mass_max"]), params["mu_max0"]) \
              for i in range(params["H_init"])]

    # Biomass of simulated host population in nmol-P.
    biomass_sim = sum([host.mass for host in H_init]) * params["nmol_p_max"]
    biomass_conc = biomass_sim / params["volume"]
    DIP = P_tot - biomass_conc

    for t in range(params["iterations"]):
        # Growth
        H_new = []
        for j in H_init:
            mu_j = (params["alpha"] * DIP) / (1 + params["alpha"] * DIP / j.mu_max)            
            j.mass += j.mass * mu_j * params["time_step"]
 
            
            if j.mass >= params["mass_max"]:
                # Create daughter cells with possible mutations
                r = random.random()
                strain = j.strain
                mu_max = j.mu_max
                if r < params["H_mutation_prob"]:
                    mu_max = max(0, random.gauss(j.mu_max, params["H_mutation_std"]))
                    strain += 1 if mu_max - j.mu_max > params["H_mutation_std"] else -1
                  
                d1 = Host(strain, j.mass / 2, mu_max)
                d2 = Host(strain, j.mass / 2, mu_max)
                H_new.append(d1)
                H_new.append(d2)
            else:
                H_new.append(j)
        
            if random.random() > params["other_loss_rate"] * params["time_step"]:
                H_new.append(j)
            else:
                DIP += j.mass * params["nmol_p_max"] / params["volume"]

        H_init = H_new
        
        V_strains = {V.strain for V in V_init}
        for i in V_strains:
            H_strains = {H.strain for H in H_init if H.strain <= i}
            V_i = [V for V in V_init if V.strain == i]
            beta = numpy.average([V.beta for V in V_i])
            for j in H_strains:
                H_j = [H for H in H_init if H.strain == j]
                infected = min(int(params["memory"] ** abs(i - j) * beta * len(V_i) * len(H_j) * params["time_step"] / params["volume"]), len(H_j))

                print i, len(V_i), j, len(H_j)
                for x in random.sample(H_j, infected):
                    DIP += x.mass * params["nmol_p_max"] / params["volume"]
                    H_init.remove(x)

                for x in random.sample(V_i, infected):
                    for y in range(params["burst_size"]):
                        strain = i
                        beta = x.beta
                        r = random.random()
                        if r < params["V_mutation_prob"]:
                            beta = max(0, random.gauss(x.beta, params["V_mutation_std"]))
                            strain += 1 if beta - x.beta > params["V_mutation_std"] else -1
                        V_init.append(Virus(strain, beta))
                    V_init.remove(x)
        
        V_new = []
        for V in V_init:
            if random.random() > params["decay_rate"] * params["time_step"]:
                V_new.append(V)
        V_init = V_new
        print len(V_init), len(H_init)


    # # Biomass concentration of simulated host population in nmol-P per L.

    # print biomass_conc, DIP, DIP + biomass_conc
