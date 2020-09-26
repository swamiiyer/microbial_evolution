# visualization code for summary.pkl
import dill, gzip, numpy, pylab, sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pylab import *

def main(args):
    fname = args[1]
    # first read the summary file and extract the details
    fh = gzip.open(fname, "rb")
    params = dill.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T = range(params["epochs"])
    bins = dill.load(fh)
    min_Host_gen = dill.load(fh)
    max_Host_gen = dill.load(fh)
    min_Virus_gen = dill.load(fh)
    max_Virus_gen = dill.load(fh)
    DIP = dill.load(fh)
    H = dill.load(fh)
    V = dill.load(fh)
    I = dill.load(fh)
    HOST_GTYPE = dill.load(fh)
    HOST_MASS = dill.load(fh)
    VIRUS_GTYPE = dill.load(fh)
    INFECTION_MAP = dill.load(fh)
    fh.close()

    # find the max number of interaction
    max_Interaction=0
    for j in range(len(INFECTION_MAP[1])):    
        if max(INFECTION_MAP[1][j]) > max_Interaction:
            max_Interaction = max(INFECTION_MAP[1][j])

    # Figure 1: 1.1. DIP vs time, 1.2. host abundance vs time, 1.3. virus abundance vs time, 1.4. infection count vs time
    figure = pylab.figure(figsize=(8, 6), dpi=500)
    size = 6
    pylab.rcParams["axes.titlesize"] = size
    pylab.rcParams["axes.labelsize"] = size
    pylab.rcParams["xtick.labelsize"] = size
    pylab.rcParams["ytick.labelsize"] = size
    pylab.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(4, 1, 1)  # Figure 1.1. DIP vs time
    axes.set_ylabel("DIP")
    axes.plot(T, DIP, "k-", alpha=0.6)

    axes = figure.add_subplot(4, 1, 2)  # Figure 1.2. host abundance vs time
    axes.set_ylabel("# of hosts")
    axes.plot(T, H, "b-", alpha=0.6)

    axes = figure.add_subplot(4, 1, 3)  # Figure 1.3. virus abundance vs time
    axes.set_ylabel("# of viruses")
    axes.plot(T, V, "r-", alpha=0.6)

    axes = figure.add_subplot(4, 1, 4)  # Figure 1.4. infection count vs time
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("# of infections")
    axes.plot(T, I, "g-", alpha=0.6)

    pylab.savefig("figure1.pdf", format="pdf", bbox_inches="tight") 

    # Figure 2. host genotype distribution vs time
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_GTYPE, cmap=pylab.cm.jet, origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host genotype distribution")
    ticks = numpy.arange(0, bins, 10)
    pylab.yticks(ticks, ticks / bins)
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("figure2.pdf", format="pdf", bbox_inches="tight")

    # Figure 3. host mass distribution vs time
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_MASS, cmap=pylab.cm.jet, origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host mass distribution")
    ticks = numpy.arange(0, bins, 10)
    pylab.yticks(ticks, numpy.arange(0.5, 1.0, 0.05))
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("figure3.pdf", format="pdf", bbox_inches="tight")

    # Figure 4. virus genotype distribution vs time
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(VIRUS_GTYPE, cmap=pylab.cm.jet, origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Virus genotype distribution")
    ticks = numpy.arange(0, bins, 10)
    pylab.yticks(ticks, ticks / bins)
    cb = pylab.colorbar(img)
    cb.set_label("# of viruses")
    pylab.savefig("figure4.pdf", format="pdf", bbox_inches="tight")

    # Figure 5. time evolution of infection map as a movie
    fig, (ax0) = plt.subplots(1,1, figsize=(8, 6))
    ax0 = plt.gca()
    # updating the colorbar in each iteration
    ax0_divider = make_axes_locatable(ax0)
    cax0 = ax0_divider.append_axes("right", size="2%", pad="2%")

    def updateHist(frame):
        t=frame        
        print(" Epoch : %d" %t)
        
        # clear axis  
        ax0.clear()    
        cax0.cla()
                 
        norm = matplotlib.colors.Normalize(vmin=0, vmax=max_Interaction)   #fixed color bar for all iterations
        im = ax0.imshow(INFECTION_MAP[t], extent=(0, 1, 0, 1), cmap=pylab.cm.jet, origin='lower', norm=norm, animated=True)                 
        fig.colorbar(im,cax=cax0)

        ax0.set_ylabel("Virus, genotype (beta)")    
        ax0.set_xlabel("Host, genotype (mu_max)")    
        ax0.set_title('Epoch= {}'.format(t) )
        
        return fig, 

    simulation = animation.FuncAnimation(fig, updateHist, frames=len(T), blit=True)
    #simulation.save('Interaction.gif', writer='PillowWriter', fps=1)
    simulation.save('Interaction.mp4', fps=1, dpi=200)    #mp4 format

if __name__ == "__main__":
    main(sys.argv)
