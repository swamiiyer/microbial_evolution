import dill, glob, gzip, numpy, sys


def main(args):
    dirname = args[1]
    bins = int(args[2])
    pklfiles = glob.glob("%s/*.pkl" % (dirname))
    params, DIP, HOSTS, VIRUSES, INTERACTIONS = None, None, None, None, None
    HOST_GTYPE, HOST_MASS, VIRUS_GTYPE = None, None, None
    for i, pklfile in enumerate(sorted(pklfiles)):
        if pklfile.endswith("summary.pkl"):
            continue
        print("Processing %s..." % (pklfile))
        fh = gzip.open(pklfile, "rb")
        params = dill.load(fh)
        if i == 0:
            DIP = numpy.zeros((len(pklfiles), params["epochs"]))
            HOSTS = numpy.zeros((len(pklfiles), params["epochs"]))
            VIRUSES = numpy.zeros((len(pklfiles), params["epochs"]))
            INTERACTIONS = numpy.zeros((len(pklfiles), params["epochs"]))
            HOST_GTYPE = numpy.zeros((bins, params["epochs"]))
            HOST_MASS = numpy.zeros((bins, params["epochs"]))
            VIRUS_GTYPE = numpy.zeros((bins, params["epochs"]))
        for j in range(0, params["epochs"]):
            dip = dill.load(fh)
            hosts = dill.load(fh)
            viruses = dill.load(fh)
            interactions = dill.load(fh)
            DIP[i, j] = dip
            HOSTS[i, j] = len(hosts)
            VIRUSES[i, j] = len(viruses)
            INTERACTIONS[i, j] = len(interactions)
            dist = [host.genotype for host in hosts]
            HOST_GTYPE[:, j] += numpy.histogram(dist, bins, density=False)[0]
            dist = [host.mass for host in hosts]
            HOST_MASS[:, j] += numpy.histogram(dist, bins, density=False)[0]
            dist = [virus.genotype for virus in viruses]
            VIRUS_GTYPE[:, j] += numpy.histogram(dist, bins, density=False)[0]
    HOST_GTYPE /= len(pklfiles)
    HOST_MASS /= len(pklfiles)
    VIRUS_GTYPE /= len(pklfiles)

    summary = gzip.open("%s/summary.pkl" % (dirname), "wb")
    dill.dump(params, summary)
    dill.dump(numpy.average(DIP, 0), summary)
    dill.dump(numpy.average(HOSTS, 0), summary)
    dill.dump(numpy.average(VIRUSES, 0), summary)
    dill.dump(numpy.average(INTERACTIONS, 0), summary)
    dill.dump(HOST_GTYPE, summary)
    dill.dump(HOST_MASS, summary)
    dill.dump(VIRUS_GTYPE, summary)
    summary.close()


if __name__ == "__main__":
    main(sys.argv)
