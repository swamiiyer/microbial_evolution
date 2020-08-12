import dill, glob, gzip, numpy, sys


def main(args):
    dirname = args[1]
    pklfiles = glob.glob("%s/*.pkl" % (dirname))
    params, DIP_MAT, HOSTS_MAT, VIRUSES_MAT, INTERACTIONS_MAT = None, None, None, None, None
    for i, pklfile in enumerate(pklfiles):
        fh = gzip.open(pklfile, "rb")
        params = dill.load(fh)
        if i == 0:
            DIP_MAT = numpy.zeros((len(pklfiles), params["epochs"]))
            HOSTS_MAT = numpy.zeros((len(pklfiles), params["epochs"]))
            VIRUSES_MAT = numpy.zeros((len(pklfiles), params["epochs"]))
            INTERACTIONS_MAT = numpy.zeros((len(pklfiles), params["epochs"]))
        for j in range(0, params["epochs"]):
            DIP = dill.load(fh)
            hosts = dill.load(fh)
            viruses = dill.load(fh)
            interactions = dill.load(fh)
            DIP_MAT[i, j] = DIP
            HOSTS_MAT[i, j] = len(hosts)
            VIRUSES_MAT[i, j] = len(viruses)
            INTERACTIONS_MAT[i, j] = len(interactions)

    summary = gzip.open("%s/summary.pkl" % (dirname), "wb")
    dill.dump(params, summary)
    dill.dump(numpy.average(DIP_MAT, 0), summary)
    dill.dump(numpy.average(HOSTS_MAT, 0), summary)
    dill.dump(numpy.average(VIRUSES_MAT, 0), summary)
    dill.dump(numpy.average(INTERACTIONS_MAT, 0), summary)
    summary.close()


if __name__ == "__main__":
    main(sys.argv)
