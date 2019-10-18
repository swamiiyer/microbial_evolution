import gzip, pickle, random, sys

class Virus(object):
    """
    Represents a virus cell.
    """

    def __init__(self, stime, beta):
        """ 
        Constructs a virus given its time of speciation and adsorption 
        coefficient.
        """

        self.stime = stime
        self.beta = beta

    def __str__(self):
        """
        Returns a string representation of this virus.
        """
        
        return "V:%d,%e" %(self.stime, self.beta)

class Host(object):
    """
    Represents a host cell.
    """

    def __init__(self, btime, stime, mass, mu_max):
        """
        Constructs a host given its time of birth, time of speciation, 
        mass, and maximum growth rate.
        """
        
        self.btime = btime
        self.stime = stime
        self.mass = mass
        self.mu_max = mu_max

    def __str__(self):
        """
        Returns a string representation of this host.
        """
        
        return "H:%d,%d,%e,%e" %(self.btime, self.stime, self.mass, self.mu_max)

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
    hosts = [Host(0, 0, random.uniform(0.5, 1.0), params["mu_max"])
             for i in range(params["H_pop"])]

    # Biomass of host population in units of individuals.
    biomass_sim = sum([host.mass for host in hosts])

    # Dissolved Phosphorous in units of individuals.
    DIP = params["P_tot"] - biomass_sim

    # Evolutionary dynamics of hosts/viruses.
    for epoch in range(1, params["epochs"] + 1):
        print("Epoch %d..." %(epoch))

        # Report if virus/host population goes extinct.
        if len(viruses) == 0:
            print("Virus population extinct!")
        if len(hosts) == 0:
            print("Host population extinct!")

        # Save the virus and host populations.
        pickle.dump(viruses, results)
        pickle.dump(hosts, results)
            
        # Host dynamics.
        next_epoch_hosts = []
        cell_divisions = 0 # number of cell divisions in this epoch
        gtime = 0.0        # average host generation time in this epoch
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
            mu = params["alpha"] * DIP / (1 + params["alpha"] * DIP / host.mu_max)
            host.mass += mu
            DIP = max(0, DIP - mu)
            if host.mass > 1.0:
                # Host divides into two daughter cells, with possible mutation.
                cell_divisions += 1
                gtime += epoch - host.btime
                mu_max = host.mu_max
                stime = host.stime
                if random.random() < params["H_mutation_prob"]:
                    mu_max = max(params["epsilon"],
                                 random.gauss(host.mu_max, params["mu_max_std"]))
                    stime = epoch
                next_epoch_hosts.append(Host(epoch, stime, host.mass / 2, mu_max))
                next_epoch_hosts.append(Host(epoch, stime, host.mass / 2, mu_max))
            else:
                # Host does not divide.
                next_epoch_hosts.append(host)
        gtime = gtime / cell_divisions if cell_divisions > 0 else 0.0
        hosts = next_epoch_hosts

        # Virus dynamics.
        next_epoch_viruses = []
        interactions = []
        for virus in viruses:
            # Loss due to decay.
            if random.random() < params["decay_rate"]:
                continue
            
            # Infection.
            uninfected_hosts = []
            p = virus.beta * params["memory"] ** abs(virus.stime - host.stime)
            for host in hosts:
                if random.random() < p:
                    # virus infects host and multiplies.
                    interactions.append((virus, host))
                    DIP += host.mass
                    stime = virus.stime
                    beta = virus.beta
                    for i in range(params["burst_size"]):
                        if random.random() < params["V_mutation_prob"]:
                            beta = max(params["epsilon"],
                                       random.gauss(virus.beta, params["beta_std"]))
                            stime = epoch
                        next_epoch_viruses.append(Virus(stime, beta))
                else:
                    uninfected_hosts.append(host)
            if len(hosts) == len(uninfected_hosts):
                # virus failed to infect any host.
                next_epoch_viruses.append(virus)
            hosts = uninfected_hosts
        viruses = next_epoch_viruses
        
        # Save the virus-host interactions and the average host
        # generation time.
        pickle.dump(interactions, results)
        pickle.dump(gtime, results)

        # DEBUG
        print("    DIP = %.2f, hosts = %d, gtime = %.2f, viruses = %d, infections = %d"
              %(DIP, len(hosts), gtime, len(viruses), len(interactions)))
        
    # Write the simulation results to the file system.    
    results.close()
