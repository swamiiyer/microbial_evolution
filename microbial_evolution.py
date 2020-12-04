import lzma, math, pickle, random


# Represents a host cell.
class Host(object):
    # Constructs a host given its genotype and mass.
    def __init__(self, g, mass):
        self.g = g
        self.mass = mass

    # Returns a string representation of this host.
    def __str__(self):
        return "H:%e,%e" % (self.g, self.mass)


# Represents a virus cell.
class Virus(object):
    # Constructs a virus given its genotype.
    def __init__(self, g):
        self.g = g

    # Returns a string representation of this virus.
    def __str__(self):
        return "V:%e" % (self.g)


# Simulates microbial evolution using the specified parameters and saves the results in a file
# with the given name.
def run(params, alpha, fname):
    # Initial host and virus populations.
    hosts = [Host(params["hG0"], random.uniform(0.5, 1.0)) for i in range(params["h0"])]
    viruses = [Virus(params["vG0"]) for i in range(params["v0"])]

    # Dissolved Phosphorous in units of individuals.
    biomass = sum([host.mass for host in hosts])
    dip0 = params["pTotal"] - biomass
    dip = dip0

    # Evolutionary dynamics of hosts/viruses.
    dipList, hostList, virusList, infectionList = [], [], [], []
    for epoch in range(params["epochs"]):
        print("Epoch %d..." % (epoch))

        # Report if host/virus population goes extinct.
        if len(hosts) == 0:
            print("Host population extinct!")
        if len(viruses) == 0:
            print("Virus population extinct!")

        # In and out flow of Phosphorous.
        dip += params["washout"] * (dip0 - dip)

        # Save the DIP and the host and virus populations.
        dipList.append(dip)
        hostList.append(hosts)
        virusList.append(viruses)

        # Host dynamics.
        nextEpochHosts = []
        for host in hosts:
            # Loss due to washout.
            if random.random() < params["washout"]:
                continue

            # Loss due to mortality.
            if random.random() < params["mortality"]:
                dip += host.mass
                continue

            # Loss due to metabolism.
            massLoss = host.mass * params["metabolicLoss"]
            host.mass -= massLoss
            dip += massLoss

            # Growth.
            mu = alpha(host.g) * dip / (1 + alpha(host.g) * dip / params["muMax"])
            massGain = mu * 1
            host.mass += massGain
            dip = max(0, dip - massGain)
            if host.mass > 1.0:
                # Host divides into two daughter cells, with possible mutation.
                hG = host.g
                if random.random() < params["hMu"]:
                    hG = min(max(0.0, random.gauss(hG, params["hGSigma"])), 1.0)
                nextEpochHosts.append(Host(hG, host.mass / 2))
                nextEpochHosts.append(Host(hG, host.mass / 2))
            else:
                # Host does not divide.
                nextEpochHosts.append(host)
        hosts = nextEpochHosts

        # Virus dynamics.
        nextEpochViruses = []
        infections = []
        for virus in viruses:
            # Loss due to washout.
            if random.random() < params["washout"]:
                continue

            # Loss due to decay.
            if random.random() < params["decay"]:
                continue

            # Infection.
            uninfectedHosts = list(hosts)
            for host in hosts:
                # Probability that virus infects host.
                beta = params["betaMax"] / params["alphaMax"] * alpha(host.g) \
                       * math.exp(-(virus.g - host.g) ** 2)

                if random.random() < beta:
                    # Virus infects host and multiplies.
                    uninfectedHosts.remove(host)
                    infections.append((virus, host))
                    dip += host.mass
                    vG = virus.g
                    for i in range(params["burst"]):
                        if random.random() < params["vMu"]:
                            vG = min(max(0.0, random.gauss(vG, params["vGSigma"])), 1.0)
                        nextEpochViruses.append(Virus(vG))
                    break

            if len(hosts) == len(uninfectedHosts):
                # Virus failed to infect any host.
                nextEpochViruses.append(virus)

            hosts = uninfectedHosts
        viruses = nextEpochViruses

        # Save the virus-host infections.
        infectionList.append(infections)

        # DEBUG
        print("    DIP = %.2f, hosts = %d, viruses = %d, infections = %d"
              % (dip, len(hosts), len(viruses), len(infections)))

    # Write the simulation results to the file system.
    fh = lzma.open(fname, "wb")
    pickle.dump(params, fh)
    pickle.dump(dipList, fh)
    pickle.dump(hostList, fh)
    pickle.dump(virusList, fh)
    pickle.dump(infectionList, fh)
    fh.close()
