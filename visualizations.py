# visualization code for summary.pkl
import pickle, lzma, numpy, pylab, sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pylab import *
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm
from itertools import chain 

def main(args):
    fname = args[1]
    # first read the summary file and extract the details
    fh = lzma.open(fname, "rb")
    params = pickle.load(fh)
    print("Parameters:")
    for param in params.keys():
        print(f"  {param} : {params[param]}")
    T = range(params["epochs"])
    bins = pickle.load(fh)
    min_H_gen = pickle.load(fh)
    max_H_gen = pickle.load(fh)
    min_H_mass = pickle.load(fh)  #new merge
    max_H_mass = pickle.load(fh)  #new merge
    min_V_gen = pickle.load(fh)
    max_V_gen = pickle.load(fh)
    DIP = pickle.load(fh)
    H = pickle.load(fh)
    V = pickle.load(fh)
    I = pickle.load(fh)
    HOST_GTYPE = pickle.load(fh)
    HOST_MASS = pickle.load(fh)
    VIRUS_GTYPE = pickle.load(fh)
    INFECTION_MAP = pickle.load(fh)
    fh.close()

    # 2 float numbers of min and max values
    min_Host_gen = float('%.2f' %(min_H_gen))
    max_Host_gen = float('%.2f' %(max_H_gen))
    min_Host_mass = float('%.2f' %(min_H_mass))
    max_Host_mass = float('%.2f' %(max_H_mass))
    min_Virus_gen = float('%.2f' %(min_V_gen))
    max_Virus_gen = float('%.2f' %(max_V_gen))
    
    # For binning, from merge.py
    hbinwidth = (max_Host_gen - min_Host_gen) / bins
    hbinlist = numpy.arange(min_Host_gen, max_Host_gen, hbinwidth)
    hbinlist = hbinlist[:bins] if len(hbinlist) > bins else hbinlist
    mbinwidth = (max_Host_mass - min_Host_mass) / bins
    mbinlist = numpy.arange(min_Host_mass, max_Host_mass, mbinwidth)
    mbinlist = mbinlist[:bins] if len(mbinlist) > bins else mbinlist
    vbinwidth = (max_Virus_gen - min_Virus_gen) / bins
    vbinlist = numpy.arange(min_Virus_gen, max_Virus_gen, vbinwidth)
    vbinlist = vbinlist[:bins] if len(vbinlist) > bins else vbinlist
    
    # find the max number of infections
    max_infections=0
    for v in INFECTION_MAP:
        max_infections = max(max_infections, v.max())        
        
    # frame speed for movies (figures 5 and 7): only show every 10th epoch
    fr_spd = 10
    fr = int(len(T) / fr_spd)    

    # Figure 1: 1.1. DIP vs time, 1.2. host abundance vs time, 1.3. virus abundance vs time, 1.4. infection count vs time
    print("Generating figure1.pdf...")
    figure = pylab.figure(figsize=(8, 6), dpi=500)
    size = 10
    pylab.rcParams["axes.titlesize"] = size
    pylab.rcParams["axes.labelsize"] = size
    pylab.rcParams["xtick.labelsize"] = size
    pylab.rcParams["ytick.labelsize"] = size
    pylab.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(4, 1, 1)  # Figure 1.1. DIP vs time
    axes.set_ylabel("DIP")
    axes.plot(T, DIP, "k-", alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(4, 1, 2)  # Figure 1.2. host abundance vs time
    axes.set_ylabel("# of hosts")
    axes.plot(T, H, "b-", alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(4, 1, 3)  # Figure 1.3. virus abundance vs time
    axes.set_ylabel("# of viruses")
    axes.plot(T, V, "r-", alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(4, 1, 4)  # Figure 1.4. infection count vs time
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("# of infections")
    axes.plot(T, I, "g-", alpha=0.6)

    pylab.savefig("figure1.pdf", format="pdf", bbox_inches="tight") 

    # Figure 2. host genotype distribution vs time
    print("Generating Host_Genotype.pdf...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_GTYPE, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host genotype")
    ticks = numpy.arange(0, bins, 10)
    show_tricks2 = hbinlist[0:-1:10]
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks2])
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("Host_Genotype.pdf", format="pdf", bbox_inches="tight")

    # Figure 3. host mass distribution vs time
    print("Generating Host_Mass.pdf...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_MASS, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host mass")
    ticks = numpy.arange(0, bins, 10)
    show_tricks3 = mbinlist[0:-1:10]
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks3])            # different method for figure 3 - ticks 
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("Host_Mass.pdf", format="pdf", bbox_inches="tight")

    # Figure 4. virus genotype distribution vs time
    print("Generating Virus_Genotype.pdf...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(VIRUS_GTYPE, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Virus genotype")
    ticks = numpy.arange(0, bins, 10)
    show_tricks4 = vbinlist[0:-1:10]
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks4])
    cb = pylab.colorbar(img)
    cb.set_label("# of viruses")
    pylab.savefig("Virus_Genotype.pdf", format="pdf", bbox_inches="tight")

    # Figure 5. time evolution of infection map as a movie
    print("Generating infections.mp4...")
    fig, (ax0) = plt.subplots(1,1)
    ax0 = plt.gca()
    # updating the colorbar in each iteration
    ax0_divider = make_axes_locatable(ax0)
    cax0 = ax0_divider.append_axes("right", size="2%", pad="2%")

    def updateHist(frame):
        t=frame * fr_spd
        print("  Epoch: %d" %t)
        
        # clear axis  
        ax0.clear()    
        cax0.cla()
                 
        norm_Log=colors.SymLogNorm(linthresh = 0.01, linscale=1.0, vmin=0, vmax=max_infections, clip=False)     #fixed colorbar
        im = ax0.imshow(INFECTION_MAP[t], extent=(min_Virus_gen, max_Virus_gen, min_Host_gen, max_Host_gen), norm=norm_Log , cmap='jet', origin='lower', aspect='auto')
        fig.colorbar(im,cax=cax0)
        cax0.set_ylabel('# of infections', rotation=270, fontsize= 8)
        #cax0.set_yticklabels(fontsize=10)

        ax0.set_ylabel("Host genotype")
        ax0.set_xlabel("Virus genotype")
        ax0.set_title(f"Epoch {t}")
        
        return fig, 

    simulation = animation.FuncAnimation(fig, updateHist, frames=fr, blit=True)
    simulation.save('infections.mp4', fps=5, dpi=200)

    # calculating Beta from summary.pkl
    Beta =[]
    for t in T:
        beta1 = []
        #print(t)
        for row in range(len(INFECTION_MAP[t][:])):
            beta2 = []
            for column in range(len(INFECTION_MAP[t][row][:])):
                host_GType_inf = INFECTION_MAP[t][row][column]
                virus_GType_inf = INFECTION_MAP[t][row][column]
                beta = (params["betaMax"] / params["hG0"] * host_GType_inf * math.exp(- (virus_GType_inf - host_GType_inf)** 2))
                beta2.append(beta)                
            beta1.append(beta2)            
        Beta.append(beta1)
    # """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # find the max value of Beta for binnings
    max_Beta = 0
    min_Beta = 0
    for t in T:
        for b in Beta[t]:
            max_Beta = max(max_Beta, max(b))
    print(max_Beta)

    bins_beta = 5 * bins                    # different bins for Beta
    Bbinwidth = (max_Beta - min_Beta) / (bins_beta)
    Bbinlist = numpy.arange(min_Beta, max_Beta, Bbinwidth)
    Bbinlist = Bbinlist[:bins] if len(Bbinlist) > bins_beta else Bbinlist

    # Figure 6. histogram of beta in each iteration
    print("Generating hist of beta.mp4...")
    fig1, (ax0) = plt.subplots(1, 1)
    def updateHist_B(frame):
        t=frame
        print("  Epoch: %d" %t)

        # clear axis  
        ax0.clear()    
        if len(list(chain.from_iterable(Beta[t]))) == 0:      # in case of extinction of virus 
            ax0.plot([])
        else:
            ax0.hist(list(chain.from_iterable(Beta[t])), bins=Bbinlist, density=False)
            ax0.set_yscale('log')           # using log for y-axis for better visualization

        ax0.set_ylabel("Abundency- Beta (log)")
        ax0.set_xlabel("bins- Beta")
        ax0.set_title(f"Epoch {t}")
        #ax0.set_ylim(0.0, 10)

        return fig1, 

    simulation1 = animation.FuncAnimation(fig1, updateHist_B, frames=len(T), blit=True)
    simulation1.save('hist_Beta.mp4', fps=1, dpi=200)
    
    # Figure 7. time evolution of Beta as a movie
    print("Generating Beta.mp4...")
    fig, (ax0) = plt.subplots(1,1)
    ax0 = plt.gca()
    # updating the colorbar in each iteration
    ax0_divider = make_axes_locatable(ax0)
    cax0 = ax0_divider.append_axes("right", size="2%", pad="2%")

    def updateHist(frame):
        
        t=frame * fr_spd       
        print("  Epoch: %d" %t)

        # clear axis  
        ax0.clear()    
        cax0.cla()
        
        norm_Log = colors.LogNorm(vmin=0.00001, vmax=max_Beta)  # vmin should be > 0.
        im = ax0.imshow(Beta[t], extent=(min_Virus_gen, max_Virus_gen, min_Host_gen, max_Host_gen), norm=norm_Log , cmap='jet', origin='lower', aspect='auto')  
        fig.colorbar(im,cax=cax0)
        cax0.set_ylabel('Beta', rotation=270, fontsize= 8)

        ax0.set_ylabel("Host genotype")
        ax0.set_xlabel("Virus genotype")
        ax0.set_title(f"Epoch {t}")

        return fig, 

    simulation = animation.FuncAnimation(fig, updateHist, frames=fr, blit=True)
    simulation.save('Beta.mp4', fps=5, dpi=200)
    
    
if __name__ == "__main__":
    main(sys.argv)
