import gzip, pickle, random, sys

def slot(a, b, x, n = 10):
    """
    Breaks the interval [a, b] into n slots and returns the slot number 
    [1, n] that x belongs to; returns -1 if x is not in the interval.
    """
    
    lo = 1
    hi = n
    dx = (b - a) / n
    while lo <= hi:
        mid = lo + (hi - lo) // 2
        mid_lo = mid * dx
        mid_hi = mid_lo + dx
        if x >= mid_lo and x < mid_hi:
            return mid
        elif x < mid_lo:
            hi = mid - 1
        else:
            lo = mid + 1
    return -1

class Virus(object):
    """
    Represents a virus cell.
    """

    def __init__(self, strain, beta):
        """ 
        Construct a virus given its strain id, and adsorption coefficient.
        """
        self.strain = strain
        self.beta = beta

    def __str__(self):
        """
        Return a string representation of this virus.
        """
        return "V[%d]: %e" %(self.strain, self.beta)

class Host(object):
    """
    Represents a host cell.
    """

    def __init__(self, strain, mass, mu_max):
        """
        Construct a host given its strain id, mass, and maximum growth rate.
        """
        self.strain = strain
        self.mass = mass
        self.mu_max = mu_max

    def __str__(self):
        """
        Return a string representation of this host.
        """
        return "H[%d]: %e, %e" %(self.strain, self.mass, self.mu_max)

def run(params):
    """
    Entry point.
    """

    P_tot = params["P_tot"]
    H_pop = [Host(slot(0.05, 0.5, params["mu_max0"], 100), \
                  random.uniform(0.5, 1.0), params["mu_max0"]) \
             for i in range(params["H_pop"])]
    V_pop = [Virus(slot(1e-5, 1e-3, params["beta0"], 100), params["beta0"]) \
             for i in range(params["V_pop"])]

    # Biomass of host population in units of individuals.
    biomass_sim = sum([host.mass for host in H_pop])

    # Dissolved Phosphorous in units of individuals.
    DIP = P_tot - biomass_sim

    output = gzip.open("results.pklz", "wb")
    pickle.dump(params, output)
    for t in range(params["epochs"]):
        print("Epoch %d..." %(t))
        pickle.dump((V_pop, H_pop), output)

        H_pop_new = []
        for j in H_pop:
            # Loss due to mortality
            if random.random() < params["mortality_rate"]:
                DIP += j.mass
                continue
            
            # Loss due to metabolism
            mass_loss = j.mass * params["metabolic_loss_rate"]
            j.mass -= mass_loss
            DIP += mass_loss

            # Growth.
            mu_j = (params["alpha0"] * DIP) / (1 + params["alpha0"] * DIP / j.mu_max)
            delta = j.mass * mu_j
            j.mass += delta
            DIP -= delta
            if j.mass >= 1.0:
                # Create daughter cells with possible mutations
                r = random.random()
                strain = j.strain
                mu_max = j.mu_max
                if r < params["H_mutation_prob"]:
                    mu_max = max(0, random.gauss(j.mu_max, params["H_mutation_std"]))
                    strain = slot(0.05, 0.5, mu_max, 100)
                    assert strain >= 1 and strain <= 100
                d1 = Host(strain, j.mass / 2, mu_max)
                d2 = Host(strain, j.mass / 2, mu_max)
                H_pop_new.append(d1)
                H_pop_new.append(d2)
            else:
                H_pop_new.append(j)

        V_pop_new = []
        for i in V_pop:
            # Loss.
            if random.random() < params["decay_rate"]:
                continue

            # Infection.
            infected = False
            for j in H_pop_new:
                p =  i.beta * len(H_pop_new) * params["memory"] ** abs(i.strain - j.strain)
                if random.random() < p: # Virus infects host and multiplies
                    infected = True
                    H_pop_new.remove(j)
                    DIP += j.mass
                    for k in range(params["burst_size"]):
                        strain = i.strain
                        beta = i.beta
                        if random.random() < params["V_mutation_prob"]:
                            beta = max(0, random.gauss(i.beta, params["V_mutation_std"]))
                            strain = slot(1e-5, 1e-3, beta, 100)
                            assert strain >= 1 and strain <= 100
                        V_pop_new.append(Virus(strain, beta))
                    break
            
            # No infection.
            if not infected:
                V_pop_new.append(i)

        H_pop = H_pop_new
        V_pop = V_pop_new

    output.close()
