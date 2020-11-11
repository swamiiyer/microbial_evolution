import lzma, math, pickle, random


# Represents a host cell.
class Host(object):
    # Constructs a host given its genotype and mass.
    def __init__(self, genotype, mass):
        self.genotype = genotype
        self.mass = mass

    # Returns a string representation of this host.
    def __str__(self):
        return "H:%e,%e" % (self.genotype, self.mass)


# Represents a virus cell.
class Virus(object):
    # Constructs a virus given its genotype.
    def __init__(self, genotype):
        self.genotype = genotype

    # Returns a string representation of this virus.
    def __str__(self):
        return "V:%e" % (self.genotype)


# Simulates microbial evolution using the specified parameters and saves the results in a file
# with the given name.
def run(params, fname):
    # Initial host and virus populations.
    hosts = [Host(params["H_genotype"], random.uniform(0.5, 1.0)) for i in range(params["H_pop"])]
    viruses = [Virus(params["V_genotype"]) for i in range(params["V_pop"])]

    # Dissolved Phosphorous in units of individuals.
    biomass = sum([host.mass for host in hosts])
    DIP0 = params["P_tot"] - biomass
    DIP = DIP0

    # Evolutionary dynamics of hosts/viruses.
    DIP_list, hosts_list, viruses_list, infections_list = [], [], [], []
    for epoch in range(params["epochs"]):
        print("Epoch %d..." % (epoch))

        # Report if host/virus population goes extinct.
        if len(hosts) == 0:
            print("Host population extinct!")
        if len(viruses) == 0:
            print("Virus population extinct!")

        # In and out flow of Phosphorous.
        DIP += params["washout_rate"] * (DIP0 - DIP)

        # Save the DIP and the host and virus populations.
        DIP_list.append(DIP)
        hosts_list.append(hosts)
        viruses_list.append(viruses)

        # Host dynamics.
        next_epoch_hosts = []
        for host in hosts:
            # Loss due to washout.
            if random.random() < params["washout_rate"]:
                continue

            # Loss due to mortality.
            if random.random() < params["mortality_rate"]:
                DIP += host.mass
                continue

            # Loss due to metabolism.
            mass_loss = host.mass * params["metabolic_loss_rate"]
            host.mass -= mass_loss
            DIP += mass_loss

            # Growth.
            alpha = params["alpha"] / params["H_genotype"] * host.genotype
            mu = alpha * DIP / (1 + alpha * DIP / params["mu_max"])
            mass_gain = mu * 1
            host.mass += mass_gain
            DIP = max(0, DIP - mass_gain)
            if host.mass > 1.0:
                # Host divides into two daughter cells, with possible mutation.
                genotype = host.genotype
                if random.random() < params["H_mutation_prob"]:
                    genotype = min(max(0.0, random.gauss(genotype, params["H_genotype_std"])), 1.0)
                next_epoch_hosts.append(Host(genotype, host.mass / 2))
                next_epoch_hosts.append(Host(genotype, host.mass / 2))
            else:
                # Host does not divide.
                next_epoch_hosts.append(host)
        hosts = next_epoch_hosts

        # Virus dynamics.
        next_epoch_viruses = []
        infections = []
        for virus in viruses:
            # Loss due to washout.
            if random.random() < params["washout_rate"]:
                continue

            # Loss due to decay.
            if random.random() < params["decay_rate"]:
                continue

            # Infection.
            uninfected_hosts = list(hosts)
            for host in hosts:
                # Probability that virus infects host.
                beta = params["beta"] / params["H_genotype"] * host.genotype * \
                       math.exp((virus.genotype - host.genotype) ** 2)
                if random.random() < beta:
                    # Virus infects host and multiplies.
                    uninfected_hosts.remove(host)
                    infections.append((virus, host))
                    DIP += host.mass
                    genotype = virus.genotype
                    for i in range(params["burst_size"]):
                        if random.random() < params["V_mutation_prob"]:
                            genotype = min(max(0.0, random.gauss(genotype,
                                                                 params["V_genotype_std"])), 1.0)
                        next_epoch_viruses.append(Virus(genotype))
                    break

            if len(hosts) == len(uninfected_hosts):
                # Virus failed to infect any host.
                next_epoch_viruses.append(virus)

            hosts = uninfected_hosts
        viruses = next_epoch_viruses

        # Save the virus-host infections.
        infections_list.append(infections)

        # DEBUG
        print("    DIP = %.2f, hosts = %d, viruses = %d, infections = %d"
              % (DIP, len(hosts), len(viruses), len(infections)))

    # Write the simulation results to the file system.
    fh = lzma.open(fname, "wb")
    pickle.dump(params, fh)
    pickle.dump(DIP_list, fh)
    pickle.dump(hosts_list, fh)
    pickle.dump(viruses_list, fh)
    pickle.dump(infections_list, fh)
    fh.close()
