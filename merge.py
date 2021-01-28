import glob, lzma, math, numpy, pickle, sys


# Entry point.
def main(args):
    # The directory containing the replicates (ie, .pkl files)
    dirname = args[1]

    # Number of bins.
    bins = int(args[2])

    # Get stats across all the replicates.
    minAlpha, maxAlpha, minMemory, maxMemory, minBeta, maxBeta = stats(dirname)
    minMass, maxMass = 0.5, 1.0

    # For binning.
    alphaBinWidth = (maxAlpha - minAlpha) / bins
    alphaBinList = numpy.arange(minAlpha, maxAlpha, alphaBinWidth)
    alphaBinList = alphaBinList[:bins] if len(alphaBinList) > bins else alphaBinList
    massBinWidth = (maxMass - minMass) / bins
    massBinList = numpy.arange(minMass, maxMass, massBinWidth)
    massBinList = massBinList[:bins] if len(massBinList) > bins else massBinList
    memoryBinWidth = (maxMemory - minMemory) / bins
    memoryBinList = numpy.arange(minMemory, maxMemory, memoryBinWidth)
    memoryBinList = memoryBinList[:bins] if len(memoryBinList) > bins else memoryBinList
    betaBinWidth = (maxBeta - minBeta) / bins
    betaBinList = numpy.arange(minBeta, maxBeta, betaBinWidth)
    betaBinList = betaBinList[:bins] if len(betaBinList) > bins else betaBinList

    # List of .pkl files within dirname.
    pklfiles = glob.glob("%s/*.pkl" % (dirname))

    # Data structures.
    params = None
    dipVal, hostCount, virusCount, infectionCount = None, None, None, None
    infectionMap = None
    alpha, mass, memory, beta = None, None, None, None
    alphaDist, massDist, memoryDist, betaDist = None, None, None, None

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
            alpha = numpy.zeros((len(pklfiles), params["epochs"]))
            mass = numpy.zeros((len(pklfiles), params["epochs"]))
            memory = numpy.zeros((len(pklfiles), params["epochs"]))
            beta = numpy.zeros((len(pklfiles), params["epochs"]))
            alphaDist = numpy.zeros((bins - 1, params["epochs"]))
            massDist = numpy.zeros((bins - 1, params["epochs"]))
            memoryDist = numpy.zeros((bins - 1, params["epochs"]))
            betaDist = numpy.zeros((bins - 1, params["epochs"]))
            infectionMap = [None] * params["epochs"]

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

            # Host genotype for the jth epoch.
            dist = [host.gAlpha for host in hosts]
            alpha[i, j] = numpy.average(dist)
            alphaDist[:, j] += numpy.histogram(dist, alphaBinList, density=False)[0]

            # Host mass for the jth epoch.
            dist = [host.mass for host in hosts]
            mass[i, j] = numpy.average(dist)
            massDist[:, j] += numpy.histogram(dist, massBinList, density=False)[0]

            # Virus genotype for the jth epoch.
            dist = [virus.gMemory for virus in viruses]
            memory[i, j] = numpy.average(dist)
            memoryDist[:, j] += numpy.histogram(dist, memoryBinList, density=False)[0]
            dist = [virus.gBeta for virus in viruses]
            beta[i, j] = numpy.average(dist)
            betaDist[:, j] += numpy.histogram(dist, betaBinList, density=False)[0]

            # Accumulate infection map for the jth epoch.
            infectingVirusG = []
            infectedHostG = []
            for virus, host in infections:
                infectingVirusG.append(virus.gBeta)
                infectedHostG.append(host.gAlpha)
            hist = numpy.histogram2d(infectedHostG, infectingVirusG, (alphaBinList, betaBinList),
                                     normed=False)[0]
            if (i == 0):
                infectionMap[j] = hist
            else:
                infectionMap[j] += hist

    # Compute averages across all replicates.
    for i in range(len(infectionMap)):
        infectionMap[i] /= len(pklfiles)
    alphaDist /= len(pklfiles)
    massDist /= len(pklfiles)
    betaDist /= len(pklfiles)

    # Write merged results to summary.pkl.
    summary = lzma.open("%s/summary.pkl" % (dirname), "wb")
    pickle.dump(params, summary)
    pickle.dump(bins, summary)
    pickle.dump(minAlpha, summary)
    pickle.dump(maxAlpha, summary)
    pickle.dump(minMass, summary)
    pickle.dump(maxMass, summary)
    pickle.dump(minMemory, summary)
    pickle.dump(maxMemory, summary)
    pickle.dump(minBeta, summary)
    pickle.dump(maxBeta, summary)
    pickle.dump(numpy.average(dipVal, 0), summary)
    pickle.dump(numpy.std(dipVal, 0), summary)
    pickle.dump(numpy.average(hostCount, 0), summary)
    pickle.dump(numpy.std(hostCount, 0), summary)
    pickle.dump(numpy.average(virusCount, 0), summary)
    pickle.dump(numpy.std(virusCount, 0), summary)
    pickle.dump(numpy.average(infectionCount, 0), summary)
    pickle.dump(numpy.std(infectionCount, 0), summary)
    pickle.dump(numpy.average(alpha, 0), summary)
    pickle.dump(numpy.std(alpha, 0), summary)
    pickle.dump(numpy.average(memory, 0), summary)
    pickle.dump(numpy.std(memory, 0), summary)
    pickle.dump(numpy.average(beta, 0), summary)
    pickle.dump(numpy.std(beta, 0), summary)
    pickle.dump(alphaDist, summary)
    pickle.dump(massDist, summary)
    pickle.dump(memoryDist, summary)
    pickle.dump(betaDist, summary)
    pickle.dump(infectionMap, summary)
    summary.close()


# Returns the minimum and maximum values of the host and virus genes across all the .pkl files
# in the specified directory.
def stats(dirname):
    print("Computing stats across all replicates...")
    pklfiles = glob.glob("%s/*.pkl" % (dirname))
    params = None
    minAlpha, maxAlpha = math.inf, -math.inf
    minMemory, maxMemory = math.inf, -math.inf
    minBeta, maxBeta = math.inf, -math.inf
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
            alpha = [host.gAlpha for host in hosts]
            minAlpha = safeMin([minAlpha, safeMin(alpha)])
            maxAlpha = safeMax([maxAlpha, safeMax(alpha)])
            memory = [virus.gBeta for virus in viruses]
            minMemory = safeMin([minBeta, safeMin(memory)])
            maxMemory = safeMax([maxBeta, safeMax(memory)])
            beta = [virus.gBeta for virus in viruses]
            minBeta = safeMin([minBeta, safeMin(beta)])
            maxBeta = safeMax([maxBeta, safeMax(beta)])
    return minAlpha, maxAlpha, minMemory, maxMemory, minBeta, maxBeta


# Returns the smallest value in a or 0.
def safeMin(a):
    return min(a) if len(a) > 0 else 0.0


# Returns the largest value in a or 0.
def safeMax(a):
    return max(a) if len(a) > 0 else 0.0


if __name__ == "__main__":
    main(sys.argv)
