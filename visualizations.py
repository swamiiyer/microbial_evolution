import dill, gzip, sys
import matplotlib.pyplot as plt


def main():
    fname = sys.argv[1]

    fh = gzip.open(fname, "rb")
    params = dill.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T, DIP, H, V, I = [], [], [], [], []
    for epoch in range(params["epochs"]):
        dip = dill.load(fh)
        hosts = dill.load(fh)
        viruses = dill.load(fh)
        interactions = dill.load(fh)
        T.append(epoch)
        DIP.append(dip)
        H.append(len(hosts))
        V.append(len(viruses))
        I.append(len(interactions))
    fh.close()

    figure = plt.figure(figsize=(8, 6), dpi=500)
    size = 6
    plt.rcParams["axes.titlesize"] = size
    plt.rcParams["axes.labelsize"] = size
    plt.rcParams["xtick.labelsize"] = size
    plt.rcParams["ytick.labelsize"] = size
    plt.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(4, 1, 1)
    axes.set_ylabel("DIP")
    axes.plot(T, DIP, "k-", alpha = 0.6)

    axes = figure.add_subplot(4, 1, 2)
    axes.set_ylabel("# of hosts")
    axes.plot(T, H, "b-", alpha = 0.6)

    axes = figure.add_subplot(4, 1, 3)
    axes.set_ylabel("# of viruses")
    axes.plot(T, V, "r-", alpha = 0.6)

    axes = figure.add_subplot(4, 1, 4)
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("# of infections")
    axes.plot(T, I, "g-", alpha = 0.6)

    plt.savefig("%s" % (fname.replace("pkl", "pdf")), format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
