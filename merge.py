import dill, glob, gzip, math, numpy, sys


# Entry point.
def main(args):
    # The directory containing the replicates (ie, .pkl files)
    dirname = args[1]

    # Number of bins.
    bins = int(args[2])

    # Get stats across all the replicates.
    min_host_genotype, max_host_genotype, min_virus_genotype, max_virus_genotype = stats(dirname)

    min_host_mass, max_host_mass = 0.5, 1.0

    # For binning.
    hbinwidth = (max_host_genotype - min_host_genotype) / bins
    hbinlist = numpy.arange(min_host_genotype, max_host_genotype, hbinwidth)
    mbinwidth = (max_host_mass - min_host_mass) / bins
    mbinlist = numpy.arange(min_host_mass, max_host_mass, mbinwidth)
    vbinwidth = (max_virus_genotype - min_virus_genotype) / bins
    vbinlist = numpy.arange(min_virus_genotype, max_virus_genotype, vbinwidth)
    assert len(hbinlist) == bins
    assert len(mbinlist) == bins
    assert len(vbinlist) == bins

    # List of .pkl files within dirname.
    pklfiles = glob.glob("%s/*.pkl" % (dirname))

    # Data structures.
    PARAMS = None
    DIP, HOST_COUNT, VIRUS_COUNT, INFECTION_COUNT = None, None, None, None
    INFECTION_MAP = None
    HOST_GENOTYPE_DIST, HOST_MASS_DIST, VIRUS_GENOTYPE_DIST = None, None, None

    # For each replicate...
    print("Merging all replicates...")
    for i, pklfile in enumerate(sorted(pklfiles)):
        # Skip over the summarized .pkl file.
        if pklfile.endswith("summary.pkl"):
            continue

        print("  Processing %s..." % (pklfile))

        # Open the .pkl file for the ith replicate.
        fh = gzip.open(pklfile, "rb")

        # Load the experimental parameters.
        PARAMS = dill.load(fh)

        # Initialize data structures.
        if i == 0:
            DIP = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            HOST_COUNT = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            VIRUS_COUNT = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            INFECTION_COUNT = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            INFECTION_MAP = [None] * PARAMS["epochs"]
            HOST_GENOTYPE_DIST = numpy.zeros((bins - 1, PARAMS["epochs"]))
            HOST_MASS_DIST = numpy.zeros((bins - 1, PARAMS["epochs"]))
            VIRUS_GENOTYPE_DIST = numpy.zeros((bins - 1, PARAMS["epochs"]))

        # For each epoch...
        for j in range(0, PARAMS["epochs"]):
            # Load data for the jth epoch.
            dip = dill.load(fh)
            hosts = dill.load(fh)
            viruses = dill.load(fh)
            infections = dill.load(fh)

            # The net DIP, number of hosts, number of viruses, and number of infections in the jth
            # epoch.
            DIP[i, j] = dip
            HOST_COUNT[i, j] = len(hosts)
            VIRUS_COUNT[i, j] = len(viruses)
            INFECTION_COUNT[i, j] = len(infections)

            # Accumulate host genotype distribution for the jth epoch.
            dist = [host.genotype for host in hosts]
            HOST_GENOTYPE_DIST[:, j] += numpy.histogram(dist, hbinlist, density=False)[0]

            # Accumulate host mass distribution for the jth epoch.
            dist = [host.mass for host in hosts]
            HOST_MASS_DIST[:, j] += numpy.histogram(dist, mbinlist, density=False)[0]

            # Accumulate virus genotype distribution for the jth epoch.
            dist = [virus.genotype for virus in viruses]
            VIRUS_GENOTYPE_DIST[:, j] += numpy.histogram(dist, vbinlist, density=False)[0]

            # Accumulate infection map for the jth epoch.
            INFECTING_VIRUS_GENOTYPE = []
            INFECTED_HOST_GENOTYPE = []
            for virus, host in infections:
                INFECTING_VIRUS_GENOTYPE.append(virus.genotype)
                INFECTED_HOST_GENOTYPE.append(host.genotype)
            hist = numpy.histogram2d(INFECTED_HOST_GENOTYPE, INFECTING_VIRUS_GENOTYPE,
                                     (hbinlist, vbinlist), normed=False)[0]
            if (i == 0):
                INFECTION_MAP[j] = hist
            else:
                INFECTION_MAP[j] += hist

    # Compute averages across all replicates.
    for i in range(len(INFECTION_MAP)):
        INFECTION_MAP[i] /= len(pklfiles)
    HOST_GENOTYPE_DIST /= len(pklfiles)
    HOST_MASS_DIST /= len(pklfiles)
    VIRUS_GENOTYPE_DIST /= len(pklfiles)

    # Write merged results to summary.pkl.
    summary = gzip.open("%s/summary.pkl" % (dirname), "wb")
    dill.dump(PARAMS, summary)
    dill.dump(bins, summary)
    dill.dump(min_host_genotype, summary)
    dill.dump(max_host_genotype, summary)
    dill.dump(min_host_mass, summary)
    dill.dump(max_host_mass, summary)
    dill.dump(min_virus_genotype, summary)
    dill.dump(max_virus_genotype, summary)
    dill.dump(numpy.average(DIP, 0), summary)
    dill.dump(numpy.average(HOST_COUNT, 0), summary)
    dill.dump(numpy.average(VIRUS_COUNT, 0), summary)
    dill.dump(numpy.average(INFECTION_COUNT, 0), summary)
    dill.dump(HOST_GENOTYPE_DIST, summary)
    dill.dump(HOST_MASS_DIST, summary)
    dill.dump(VIRUS_GENOTYPE_DIST, summary)
    dill.dump(INFECTION_MAP, summary)
    summary.close()


# Returns the minimum and maximum values of the host and virus genotype across all the .pkl files
# in the specified directory.
def stats(dirname):
    print("Computing stats across all replicates...")
    pklfiles = glob.glob("%s/*.pkl" % (dirname))
    params = None
    min_host_genotype, max_host_genotype = math.inf, -math.inf
    min_virus_genotype, max_virus_genotype = math.inf, -math.inf
    for i, pklfile in enumerate(sorted(pklfiles)):
        if pklfile.endswith("summary.pkl"):
            continue
        print("  Processing %s..." % (pklfile))
        fh = gzip.open(pklfile, "rb")
        params = dill.load(fh)
        for j in range(0, params["epochs"]):
            dip = dill.load(fh)
            hosts = dill.load(fh)
            viruses = dill.load(fh)
            infections = dill.load(fh)
            genotypes = [host.genotype for host in hosts]
            min_host_genotype = min(min_host_genotype, min(genotypes))
            max_host_genotype = max(max_host_genotype, max(genotypes))
            genotypes = [virus.genotype for virus in viruses]
            min_virus_genotype = min(min_virus_genotype, min(genotypes))
            max_virus_genotype = max(max_virus_genotype, max(genotypes))
    return min_host_genotype, max_host_genotype, min_virus_genotype, max_virus_genotype


if __name__ == "__main__":
    main(sys.argv)
