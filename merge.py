import dill, glob, gzip, numpy, sys


def main(args):
    dirname = args[1]
    bins = int(args[2])

    pklfiles = glob.glob("%s/*.pkl" % (dirname))

    PARAMS = None
    DIP, HOST_COUNT, VIRUS_COUNT, INFECTION_COUNT = None, None, None, None
    INFECTION_MAP = None
    HOST_GTYPE_DIST, HOST_MASS_DIST, VIRUS_GTYPE_DIST = None, None, None

    # For each replicate...
    for i, pklfile in enumerate(sorted(pklfiles)):
        # Skip over the summarized .pkl file.
        if pklfile.endswith("summary.pkl"):
            continue

        print("Processing %s..." % (pklfile))

        # Open the .pkl file for the replicate.
        fh = gzip.open(pklfile, "rb")

        # Load the experimental parameters.
        PARAMS = dill.load(fh)

        # Initialize the data structures.
        if i == 0:
            DIP = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            HOST_COUNT = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            VIRUS_COUNT = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            INFECTION_COUNT = numpy.zeros((len(pklfiles), PARAMS["epochs"]))
            INFECTION_MAP = [None] * PARAMS["epochs"]
            HOST_GTYPE_DIST = numpy.zeros((bins, PARAMS["epochs"]))
            HOST_MASS_DIST = numpy.zeros((bins, PARAMS["epochs"]))
            VIRUS_GTYPE_DIST = numpy.zeros((bins, PARAMS["epochs"]))

        # For each epoch...
        for j in range(0, PARAMS["epochs"]):

            dip = dill.load(fh)
            hosts = dill.load(fh)
            viruses = dill.load(fh)
            infections = dill.load(fh)

            DIP[i, j] = dip
            HOST_COUNT[i, j] = len(hosts)
            VIRUS_COUNT[i, j] = len(viruses)
            INFECTION_COUNT[i, j] = len(infections)
            dist = [host.genotype for host in hosts]
            HOST_GTYPE_DIST[:, j] += numpy.histogram(dist, bins, density=False)[0]
            dist = [host.mass for host in hosts]
            HOST_MASS_DIST[:, j] += numpy.histogram(dist, bins, density=False)[0]
            dist = [virus.genotype for virus in viruses]
            VIRUS_GTYPE_DIST[:, j] += numpy.histogram(dist, bins, density=False)[0]

            INFECTING_VIRUS_GTYPE = []
            INFECTED_HOST_GTYPE = []
            for virus, host in infections:
                INFECTING_VIRUS_GTYPE.append(virus.genotype)
                INFECTED_HOST_GTYPE.append(host.genotype)
            hist = numpy.histogram2d(INFECTED_HOST_GTYPE, INFECTING_VIRUS_GTYPE, (bins, bins),
                                     normed=False)[0]
            if (i == 0):
                INFECTION_MAP[j] = hist
            else:
                INFECTION_MAP[j] += hist

    for i in range(len(INFECTION_MAP)):
        INFECTION_MAP[i] /= len(pklfiles)
    HOST_GTYPE_DIST /= len(pklfiles)
    HOST_MASS_DIST /= len(pklfiles)
    VIRUS_GTYPE_DIST /= len(pklfiles)

    summary = gzip.open("%s/summary.pkl" % (dirname), "wb")
    dill.dump(PARAMS, summary)
    dill.dump(numpy.average(DIP, 0), summary)
    dill.dump(numpy.average(HOST_COUNT, 0), summary)
    dill.dump(numpy.average(VIRUS_COUNT, 0), summary)
    dill.dump(numpy.average(INFECTION_COUNT, 0), summary)
    dill.dump(HOST_GTYPE_DIST, summary)
    dill.dump(HOST_MASS_DIST, summary)
    dill.dump(VIRUS_GTYPE_DIST, summary)
    dill.dump(INFECTION_MAP, summary)
    summary.close()


if __name__ == "__main__":
    main(sys.argv)
