import dill, gzip, numpy, pylab, sys


def main(args):
    fname = args[1]

    fh = gzip.open(fname, "rb")
    params = dill.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T = range(params["epochs"])
    DIP = dill.load(fh)
    H = dill.load(fh)
    V = dill.load(fh)
    I = dill.load(fh)
    HOST_GTYPE = dill.load(fh)
    HOST_MASS = dill.load(fh)
    VIRUS_GTYPE = dill.load(fh)
    INFECTION_MAP = dill.load(fh)
    fh.close()

    figure = pylab.figure(figsize=(8, 6), dpi=500)
    size = 6
    pylab.rcParams["axes.titlesize"] = size
    pylab.rcParams["axes.labelsize"] = size
    pylab.rcParams["xtick.labelsize"] = size
    pylab.rcParams["ytick.labelsize"] = size
    pylab.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(4, 1, 1)
    axes.set_ylabel("DIP")
    axes.plot(T, DIP, "k-", alpha=0.6)

    axes = figure.add_subplot(4, 1, 2)
    axes.set_ylabel("# of hosts")
    axes.plot(T, H, "b-", alpha=0.6)

    axes = figure.add_subplot(4, 1, 3)
    axes.set_ylabel("# of viruses")
    axes.plot(T, V, "r-", alpha=0.6)

    axes = figure.add_subplot(4, 1, 4)
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("# of infections")
    axes.plot(T, I, "g-", alpha=0.6)

    pylab.savefig("figure1.pdf", format="pdf", bbox_inches="tight")

    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_GTYPE, cmap=pylab.cm.Reds, origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host genotype distribution")
    bins = HOST_GTYPE.shape[0]
    ticks = numpy.arange(0, bins, 10)
    pylab.yticks(ticks, ticks / bins)
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("figure2.pdf", format="pdf", bbox_inches="tight")

    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_MASS, cmap=pylab.cm.Reds, origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host mass distribution")
    bins = HOST_MASS.shape[0]
    ticks = numpy.arange(0, bins, 10)
    pylab.yticks(ticks, numpy.arange(0.5, 1.0, 0.05))
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("figure3.pdf", format="pdf", bbox_inches="tight")

    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(VIRUS_GTYPE, cmap=pylab.cm.Reds, origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Virus genotype distribution")
    bins = VIRUS_GTYPE.shape[0]
    ticks = numpy.arange(0, bins, 10)
    pylab.yticks(ticks, ticks / bins)
    cb = pylab.colorbar(img)
    cb.set_label("# of viruses")
    pylab.savefig("figure4.pdf", format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv)
