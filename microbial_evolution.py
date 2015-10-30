import random, sys

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
        H_new = []
        for j in H_init:
            # Loss due to metabolism
            if params["other_loss_type"] == 1 or params["other_loss_type"] == 2:
                mass_loss = j.mass * params["other_loss_rate"] * params["time_step"]
                j.mass -= mass_loss
                DIP += mass_loss
            
            # Loss due to mortality
            if params["other_loss_type"] == 0 or params["other_loss_type"] == 2:
                if random.random() < params["other_loss_rate"] * params["time_step"]:
                    DIP += j.mass * params["nmol_p_max"] / params["volume"]
                    continue

            # Growth.
            mu_j = (params["alpha"] * DIP) / (1 + params["alpha"] * DIP / j.mu_max)            
            j.mass = j.mass + j.mass * mu_j * params["time_step"]
            print mu_j
            if j.mass >= params["mass_max"]:
                # Create daughter cells with possible mutations
                r = random.random()
                strain = j.strain
                mu_max = j.mu_max
                if r < params["H_mutation_prob"]:
                    mu_max = max(0, random.gauss(j.mu_max, params["H_mutation_std"]))
                    if abs(mu_max - j.mu_max) > params["H_mutation_std"]:
                        strain += 1 if mu_max - j.mu_max > 0 else -1
                d1 = Host(strain, j.mass / 2, mu_max)
                d2 = Host(strain, j.mass / 2, mu_max)
                H_new.append(d1)
                H_new.append(d2)
            else:
                H_new.append(j)

        V_new = []
        for i in V_init:
            # Loss.
            if random.random() < params["decay_rate"] * params["time_step"]:
                continue

            # Infection.
            infected = False
            for j in H_init:
                p = params["memory"] ** abs(i.strain - j.strain) * i.beta * params["time_step"] / params["volume"]
                if random.random() < p: # Virus infects host and multiplies
                    DIP += j.mass * params["nmol_p_max"] / params["volume"]
                    for k in range(params["burst_size"]):
                        strain = i.strain
                        beta = i.beta
                        if random.random() < params["V_mutation_prob"]:
                            beta = max(0, random.gauss(i.beta, params["V_mutation_std"]))
                            if abs(beta - i.beta) > params["V_mutation_std"]:
                                strain += 1 if beta - i.beta > 0 else -1
                        V_new.append(Virus(strain, beta))
                    infected = True
                    if j in H_new:
                        H_new.remove(j)
                    break
            
            # No infection.
            if not infected:
                V_new.append(i)

        H_init = H_new
        V_init = V_new

        V_strains = {v.strain for v in V_init}
        H_strains = {h.strain for h in H_init}
        
        print t, len(V_init), len(H_init), len(V_strains), len(H_strains)
