import glob, lzma, math, numpy, pickle, sys


# Entry point.
def main(args):
    # The directory containing the replicates (ie, .pkl files)
    dirname = args[1]

    # Number of bins.
    bins = int(args[2])

    # Get stats across all the replicates.
    minHG, maxHG, minVG, maxVG = stats(dirname)
    minHMass, maxHMass = 0.5, 1.0

    # For binning.
    hbinwidth = (maxHG - minHG) / bins
    hbinlist = numpy.arange(minHG, maxHG, hbinwidth)
    hbinlist = hbinlist[:bins] if len(hbinlist) > bins else hbinlist
    mbinwidth = (maxHMass - minHMass) / bins
    mbinlist = numpy.arange(minHMass, maxHMass, mbinwidth)
    mbinlist = mbinlist[:bins] if len(mbinlist) > bins else mbinlist
    vbinwidth = (maxVG - minVG) / bins
    vbinlist = numpy.arange(minVG, maxVG, vbinwidth)
    vbinlist = vbinlist[:bins] if len(vbinlist) > bins else vbinlist

    # List of .pkl files within dirname.
    pklfiles = glob.glob("%s/*.pkl" % (dirname))

    # Data structures.
    params = None
    dipVal, hostCount, virusCount, infectionCount = None, None, None, None
    infectionMap = None
    hostGDist, hostMassDist, virusGDist = None, None, None

    # For each replicate...
    print("Merging all replicates...")
    for i, pklfile in enumerate(sorted(pklfiles)):
        # Skip over the summarized .pkl file.
        if pklfile.endswith("summary.pkl"):
            continue

        print("  Processing %s..." % (pklfile))

        # Read the results from the .pkl file for the ith replicate.
        fh = lzma.open(pklfile, "rb")
        params = pickle.load(fh)
        dipList = pickle.load(fh)
        hostList = pickle.load(fh)
        virusList = pickle.load(fh)
        infectionList = pickle.load(fh)
        fh.close()

        # Initialize data structures.
        if i == 0:
            dipVal = numpy.zeros((len(pklfiles), params["epochs"]))
            hostCount = numpy.zeros((len(pklfiles), params["epochs"]))
            virusCount = numpy.zeros((len(pklfiles), params["epochs"]))
            infectionCount = numpy.zeros((len(pklfiles), params["epochs"]))
            infectionMap = [None] * params["epochs"]
            hostGDist = numpy.zeros((bins - 1, params["epochs"]))
            hostMassDist = numpy.zeros((bins - 1, params["epochs"]))
            virusGDist = numpy.zeros((bins - 1, params["epochs"]))

        # For each epoch...
        for j in range(0, params["epochs"]):
            # Load data for the jth epoch.
            dip = dipList[j]
            hosts = hostList[j]
            viruses = virusList[j]
            infections = infectionList[j]

            # The net DIP, number of hosts, number of viruses, and number of infections in the jth
            # epoch.
            dipVal[i, j] = dip
            hostCount[i, j] = len(hosts)
            virusCount[i, j] = len(viruses)
            infectionCount[i, j] = len(infections)

            # Accumulate host genotype distribution for the jth epoch.
            dist = [host.g for host in hosts]
            hostGDist[:, j] += numpy.histogram(dist, hbinlist, density=False)[0]

            # Accumulate host mass distribution for the jth epoch.
            dist = [host.mass for host in hosts]
            hostMassDist[:, j] += numpy.histogram(dist, mbinlist, density=False)[0]

            # Accumulate virus genotype distribution for the jth epoch.
            dist = [virus.g for virus in viruses]
            virusGDist[:, j] += numpy.histogram(dist, vbinlist, density=False)[0]

            # Accumulate infection map for the jth epoch.
            infectingVirusG = []
            infectedHostG = []
            for virus, host in infections:
                infectingVirusG.append(virus.g)
                infectedHostG.append(host.g)
            hist = numpy.histogram2d(infectedHostG, infectingVirusG, (hbinlist, vbinlist),
                                     normed=False)[0]
            if (i == 0):
                infectionMap[j] = hist
            else:
                infectionMap[j] += hist

    # Compute averages across all replicates.
    for i in range(len(infectionMap)):
        infectionMap[i] /= len(pklfiles)
    hostGDist /= len(pklfiles)
    hostMassDist /= len(pklfiles)
    virusGDist /= len(pklfiles)

    # Write merged results to summary.pkl.
    summary = lzma.open("%s/summary.pkl" % (dirname), "wb")
    pickle.dump(params, summary)
    pickle.dump(bins, summary)
    pickle.dump(minHG, summary)
    pickle.dump(maxHG, summary)
    pickle.dump(minHMass, summary)
    pickle.dump(maxHMass, summary)
    pickle.dump(minVG, summary)
    pickle.dump(maxVG, summary)
    pickle.dump(numpy.average(dipVal, 0), summary)
    pickle.dump(numpy.average(hostCount, 0), summary)
    pickle.dump(numpy.average(virusCount, 0), summary)
    pickle.dump(numpy.average(infectionCount, 0), summary)
    pickle.dump(hostGDist, summary)
    pickle.dump(hostMassDist, summary)
    pickle.dump(virusGDist, summary)
    pickle.dump(infectionMap, summary)
    summary.close()


# Returns the minimum and maximum values of the host and virus genotype across all the .pkl files
# in the specified directory.
def stats(dirname):
    print("Computing stats across all replicates...")
    pklfiles = glob.glob("%s/*.pkl" % (dirname))
    params = None
    minHG, maxHG = math.inf, -math.inf
    minVG, maxVG = math.inf, -math.inf
    for i, pklfile in enumerate(sorted(pklfiles)):
        if pklfile.endswith("summary.pkl"):
            continue
        print("  Processing %s..." % (pklfile))
        fh = lzma.open(pklfile, "rb")
        params = pickle.load(fh)
        dipList = pickle.load(fh)
        hostList = pickle.load(fh)
        virusList = pickle.load(fh)
        infectionList = pickle.load(fh)
        fh.close()
        for j in range(0, params["epochs"]):
            dip = dipList[j]
            hosts = hostList[j]
            viruses = virusList[j]
            infections = infectionList[j]
            genotypes = [host.g for host in hosts]
            minHG = safeMin([minHG, safeMin(genotypes)])
            maxHG = safeMax([maxHG, safeMax(genotypes)])
            genotypes = [virus.g for virus in viruses]
            minVG = safeMin([minVG, safeMin(genotypes)])
            maxVG = safeMax([maxVG, safeMax(genotypes)])
    return minHG, maxHG, minVG, maxVG


# Returns the smallest value in a or 0.
def safeMin(a):
    return min(a) if len(a) > 0 else 0.0


# Returns the largest value in a or 0.
def safeMax(a):
    return max(a) if len(a) > 0 else 0.0


if __name__ == "__main__":
    main(sys.argv)
