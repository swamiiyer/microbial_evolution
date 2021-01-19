import lzma, pickle, random


# Represents a host cell.
class Host(object):
    # Constructs a host given its genotype and mass.
    def __init__(self, gAlpha, mass):
        self.gAlpha = gAlpha
        self.mass = mass


# Represents a virus cell.
class Virus(object):
    # Constructs a virus given its genotype.
    def __init__(self, gMemory, gBeta):
        self.gMemory = gMemory
        self.gBeta = gBeta


# Simulates microbial evolution using the specified parameters and functions (alpha,
# compatibility, and beta), and saves the results in a file with the given name.
def run(params, alpha, compatibility, beta, fname):
    # Initial host and virus populations.
    hosts = [Host(params["gAlpha0"], random.uniform(0.5, 1.0)) for i in range(params["h0"])]
    viruses = [Virus(params["gMemory0"], params["gBeta0"]) for i in range(params["v0"])]

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
            mu = alpha(host) * dip / (1 + alpha(host) * dip / params["muMax"])
            massGain = mu * 1
            host.mass += massGain
            dip = max(0, dip - massGain)
            if host.mass > 1.0:
                # Host divides into two daughter cells, with possible mutation.
                gAlpha = host.gAlpha
                if random.random() < params["hMu"]:
                    gAlpha = min(max(0.0, random.gauss(gAlpha, params["hGSigma"])), 1.0)
                nextEpochHosts.append(Host(gAlpha, host.mass / 2))
                nextEpochHosts.append(Host(gAlpha, host.mass / 2))
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
                if random.random() < compatibility(host, virus):
                    # Virus compatible with host.
                    if random.random() < beta(virus):
                        # Virus infects host and multiplies.
                        uninfectedHosts.remove(host)
                        infections.append((virus, host))
                        dip += host.mass
                        gMemory = virus.gMemory
                        gBeta = virus.gBeta
                        for i in range(params["burst"]):
                            if random.random() < params["vMu"]:
                                gMemory = min(max(0.0, random.gauss(gMemory, params["vGSigma"])),
                                              1.0)
                                gBeta = min(max(0.0, random.gauss(gBeta, params["vGSigma"])), 1.0)
                            nextEpochViruses.append(Virus(gMemory, gBeta))
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
