import random, sys

class Virus(object):
    """
    Represents a virus cell.
    """

    def __init__(self, species, strain, beta):
        self.species = species
        self.strain = strain
        self.beta = beta

    def __str__(self):
        return "V(%d): beta = %e" %(self.strain, self.beta)

class Host(object):
    """
    Represents a host cell.
    """

    def __init__(self, species, strain, mass, mu_max, alpha):
        self.species = species
        self.strain = strain
        self.mass = mass
        self.mu_max = mu_max
        self.alpha = alpha

    def __str__(self):
        return "H(%d): mass = %e, mu_max = %e" %(self.strain, self.mass, \
                                                 self.mu_max)

def run(params):
    """
    Entry point.
    """

    random.seed(params["seed"])
    
    P_tot = params["P_tot"]
    V_init = [Virus(0, 0, params["beta00"]) for i in range(params["V_init"])]
    H_init = [Host(0, 0, random.uniform(0.5 * params["mass_max"], 
                                        params["mass_max"]), params["mu_max00"], params["alpha00"]) \
              for i in range(params["H_init"])]

    # Biomass of simulated host population.
    biomass_sim = sum([host.mass for host in H_init])
    #biomass_conc = biomass_sim / params["volume"] (TBD)
    DIP = P_tot - biomass_sim

    print "seed = %d" %(params["seed"])

    for t in range(params["iterations"]):
        V_strains = len({v.strain for v in V_init})
        V_species = len({v.species for v in V_init})
        H_strains = len({h.strain for h in H_init})
        H_species = len({h.species for h in H_init})
                
        biomass_sim = sum([j.mass for j in H_init])
        print "%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d" %(DIP, biomass_sim, DIP + biomass_sim, len(V_init), len(H_init), V_strains, V_species, H_strains, H_species)

        H_new = []
        for j in H_init:
            # Loss due to metabolism
            if params["other_loss_type"] == 1 or params["other_loss_type"] == 2:
                mass_loss = j.mass * params["other_loss_rate"]
                j.mass -= mass_loss
                DIP += mass_loss
            
            # Loss due to mortality
            if params["other_loss_type"] == 0 or params["other_loss_type"] == 2:
                if random.random() < params["other_loss_rate"]:
                    DIP += j.mass
                    continue

            # Growth.
            mu_j = (j.alpha * DIP) / (1 + j.alpha * DIP / j.mu_max)            
            delta = j.mass * mu_j
            j.mass = j.mass + delta
            DIP -= delta
            if j.mass >= params["mass_max"]:
                # Create daughter cells with possible mutations
                r = random.random()
                species = j.species
                strain = j.strain
                mu_max = j.mu_max
                alpha = j.alpha

                # species level
                if r < params["H_mutation_prob"][0]:
                    alpha = max(0, random.gauss(j.alpha, params["H_mutation_std"][0]))

                    # mu_max has an inverse relationship with alpha
                    mu_max = mu_max - mu_max * (alpha - j.alpha) / j.alpha

                    if alpha - j.alpha > params["H_mutation_std"][0]:
                        species -= 1
                    else:
                        species += 1

                # strain level
                if r < params["H_mutation_prob"][1]:
                    mu_max = max(0, random.gauss(j.mu_max, params["H_mutation_std"][1]))
                    if abs(mu_max - j.mu_max) > params["H_mutation_std"][1]:
                        strain += 1 if mu_max - j.mu_max > 0 else -1

                d1 = Host(species, strain, j.mass / 2, mu_max, alpha)
                d2 = Host(species, strain, j.mass / 2, mu_max, alpha)
                H_new.append(d1)
                H_new.append(d2)
            else:
                H_new.append(j)

        V_new = []
        for i in V_init:
            # Loss.
            if random.random() < params["decay_rate"]:
                continue

            # Infection.
            infected = False
            for j in H_init:
                if i.species != j.species:
                    continue

                p =  i.beta * params["memory"] ** abs(i.strain - j.strain)
                if random.random() < p: # Virus infects host and multiplies
                    DIP += j.mass
                    for k in range(params["burst_size"]):
                        species = i.species
                        strain = i.strain
                        beta = i.beta

                        # Species level
                        if random.random() < params["V_mutation_prob"][0]:
                            beta = max(0, random.gauss(i.beta, params["V_mutation_std"][0]))
                            if abs(beta - i.beta) > params["V_mutation_std"][0]:
                                species += 1 if beta - i.beta > 0 else -1

                        # Strain level
                        if random.random() < params["V_mutation_prob"][1]:
                            beta = max(0, random.gauss(i.beta, params["V_mutation_std"][1]))
                            if abs(beta - i.beta) > params["V_mutation_std"][1]:
                                strain += 1 if beta - i.beta > 0 else -1

                        V_new.append(Virus(species, strain, beta))
                    infected = True
                    if j in H_new:
                        H_new.remove(j)
                    break
            
            # No infection.
            if not infected:
                V_new.append(i)

        H_init = H_new
        V_init = V_new

        
        
        
        
