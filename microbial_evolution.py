import gzip, pickle, random, sys

class Virus(object):
    """
    Represents a virus cell.
    """

    def __init__(self, ctime, beta):
        """ 
        Constructs a virus given its creation time and adsorption coefficient.
        """
        
        self.ctime = ctime
        self.beta = beta

    def __str__(self):
        """
        Returns a string representation of this virus.
        """
        
        return "V[%d]: %e" %(self.ctime, self.beta)

class Host(object):
    """
    Represents a host cell.
    """

    def __init__(self, ctime, mass, mu_max):
        """
        Constructs a host given its creation time, mass, and maximum 
        growth rate.
        """
        
        self.ctime = ctime
        self.mass = mass
        self.mu_max = mu_max

    def __str__(self):
        """
        Returns a string representation of this host.
        """
        
        return "H[%d]: %e, %e" %(self.ctime, self.mass, self.mu_max)

def run(params, fname):
    """
    Simulates microbial evolution using the specified parameters and saves
    the results in a file with the given name.
    """

    # File for saving simulation results.
    results = gzip.open(fname, "wb")

    # Save the parameter values.
    pickle.dump(params, results)

    # Initial virus and host populations.
    viruses = [Virus(0, params["beta"]) for i in range(params["V_pop"])]
    hosts = [Host(0, random.uniform(0.5, 1.0), params["mu_max"])
             for i in range(params["H_pop"])]

    # Biomass of host population in units of individuals.
    biomass_sim = sum([host.mass for host in hosts])

    # Dissolved Phosphorous in units of individuals.
    DIP = params["P_tot"] - biomass_sim

    for epoch in range(params["epochs"]):
        print("Epoch %d..." %(epoch))

        # Report and quit if virus/host population goes extinct.
        if len(viruses) == 0:
            print("Virus population extinct!")
            break
        if len(hosts) == 0:
            print("Host population extinct!")
            break

        # Save virus/host populations at the current epoch.
        pickle.dump(viruses, results)
        pickle.dump(hosts, results)

        # Host dynamics.
        temp_hosts = []
        for host in hosts:
            # Loss due to mortality.
            if random.random() < params["mortality_rate"]:
                DIP += host.mass
                continue

            # Loss due to metabolism.
            mass_loss = host.mass * params["metabolic_loss_rate"]
            host.mass -= mass_loss
            DIP += mass_loss

            # Growth.
            mu = params["alpha"] * DIP / \
                 (1 + params["alpha"] * DIP / host.mu_max)
            delta = host.mass * mu
            host.mass += delta
            DIP = max(0, DIP - delta)
            if host.mass > 1.0:
                # Host divides into two daughter cells, with possible mutation.
                mu_max = host.mu_max
                ctime = host.ctime
                if random.random() < params["H_mutation_prob"]:
                    mu_max = max(1e-5, random.gauss(host.mu_max,
                                                 params["mu_max_std"]))
                    ctime = epoch
                temp_hosts.append(Host(ctime, host.mass / 2, mu_max))
                temp_hosts.append(Host(ctime, host.mass / 2, mu_max))
            else:
                # Host does not divide.
                temp_hosts.append(host)

        hosts = temp_hosts
        
    # Write the simulation results to the file system.    
    results.close()
