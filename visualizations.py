# visualization code for summary.pkl
import pickle, lzma, numpy, pylab, sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pylab import *
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm
from itertools import chain 
from scipy.stats import entropy

def main(args):
    fname = args[1]
    # read the summary file and extract the details
    fh = lzma.open(fname, "rb")
    params = pickle.load(fh)
    print("Parameters:")
    for param in params.keys():
        print(f"  {param} : {params[param]}")
    T = range(params["epochs"])
    bins = pickle.load(fh)
    # min & max 
    minAlpha = pickle.load(fh)      
    maxAlpha = pickle.load(fh)     
    minMass = pickle.load(fh)      
    maxMass = pickle.load(fh)       
    minMemory = pickle.load(fh)     
    maxMemory = pickle.load(fh)     
    minBeta = pickle.load(fh)      
    maxBeta = pickle.load(fh)      
    # DIP
    DIP = pickle.load(fh)           # average of dipVal
    DIP_std = pickle.load(fh)       # std of dipVal
    # hostCount
    H = pickle.load(fh)             # average of hostCount
    H_std = pickle.load(fh)         # std of hostCount
    # virusCount
    V = pickle.load(fh)             # average of virusCount
    V_std = pickle.load(fh)         # std of virusCount
    # infectionCount
    I = pickle.load(fh)             # average of infectionCount
    I_std = pickle.load(fh)         # std of infectionCount

    alpha = pickle.load(fh)         # new merge
    alpha_std = pickle.load(fh)     # new merge
    memory = pickle.load(fh)        # new merge
    memory_std = pickle.load(fh)    # new merge
    beta = pickle.load(fh)          # new merge
    beta_std = pickle.load(fh)      # new merge

    # Distributions
    alphaDist = pickle.load(fh)     # host genotype
    HOST_MASS = pickle.load(fh)     # massDist
    memoryDist = pickle.load(fh)    # new merge
    betaDist = pickle.load(fh)      # virus genotype (beta)
    INFECTION_MAP = pickle.load(fh) 
    fh.close()


    # float numbers of min and max values
    min_Alpha = float('%.2f' %(minAlpha))       
    max_Alpha = float('%.2f' %(maxAlpha))       
    min_Mass = float('%.2f' %(minMass))    
    max_Mass = float('%.2f' %(maxMass))    
    min_Memory = float('%.2f' %(minMemory))  # new merge
    max_Memory = float('%.2f' %(maxMemory))  # new merge
    min_Beta = float('%.2f' %(minBeta))  
    max_Beta = float('%.2f' %(maxBeta))  
    
    # binning
    alphaBinWidth = (max_Alpha - min_Alpha) / bins
    alphaBinList = numpy.arange(min_Alpha, max_Alpha, alphaBinWidth)
    alphaBinList = alphaBinList[:bins] if len(alphaBinList) > bins else alphaBinList
    massBinWidth = (max_Mass - min_Mass) / bins
    massBinList = numpy.arange(min_Mass, max_Mass, massBinWidth)
    massBinList = massBinList[:bins] if len(massBinList) > bins else massBinList
    memoryBinWidth = (max_Memory - min_Memory) / bins
    memoryBinList = numpy.arange(min_Memory, max_Memory, memoryBinWidth)
    memoryBinList = memoryBinList[:bins] if len(memoryBinList) > bins else memoryBinList
    betaBinWidth = (max_Beta - min_Beta) / bins
    betaBinList = numpy.arange(min_Beta, max_Beta, betaBinWidth)
    betaBinList = betaBinList[:bins] if len(betaBinList) > bins else betaBinList



    # find the max number of infections
    max_infections=0
    for v in INFECTION_MAP:
        max_infections = max(max_infections, v.max())        
        
    # frame speed for movies (figures 7): only show every 10th epoch
    fr_spd = 10
    fr = int(len(T) / fr_spd)    

    # Figure 1: 1.1. DIP vs time, 1.2. host abundance vs time, 1.3. virus abundance vs time, 1.4. infection count vs time
    print("Generating figure1.pdf/jpeg...")
    figure = pylab.figure(figsize=(8, 6), dpi=500)
    size = 10
    pylab.rcParams["axes.titlesize"] = size
    pylab.rcParams["axes.labelsize"] = size
    pylab.rcParams["xtick.labelsize"] = size
    pylab.rcParams["ytick.labelsize"] = size
    pylab.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(4, 1, 1)  # Figure 1.1. DIP vs time
    axes.set_ylabel("DIP")
    axes.errorbar(T, DIP, yerr = DIP_std,  color='black', ecolor='lightgray', alpha=0.6)   # ls = '-' is a default
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(4, 1, 2)  # Figure 1.2. host abundance vs time
    axes.set_ylabel("# of hosts")
    axes.errorbar(T, H, yerr = H_std, color='blue', ecolor ='lightblue', alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(4, 1, 3)  # Figure 1.3. virus abundance vs time
    axes.set_ylabel("# of viruses")
    axes.errorbar(T, V, yerr = V_std, color='red', ecolor='mistyrose', alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(4, 1, 4)  # Figure 1.4. infection count vs time
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("# of infections")
    axes.errorbar(T, I, yerr = I_std, color='green', ecolor='lightgreen',  alpha=0.6)

    pylab.savefig("figure1.pdf", format="pdf", bbox_inches="tight") 
    pylab.savefig("figure1.jpeg", format="jpeg", bbox_inches="tight")


    # Figure 2: mean alpha, mean memory, mean beta w/ w.o std
    print("Generating genes.pdf/jpeg...")
    figure = pylab.figure(figsize=(10, 6), dpi=500)
    size = 8
    pylab.rcParams["axes.titlesize"] = size
    pylab.rcParams["axes.labelsize"] = size
    pylab.rcParams["xtick.labelsize"] = size
    pylab.rcParams["ytick.labelsize"] = size
    pylab.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(3, 2, 1)  
    axes.set_ylabel("Alpha")
    axes.plot(T, alpha, "b-", alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(3, 2, 3) 
    axes.set_ylabel("Memory")
    axes.plot(T, memory, "r-", alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(3, 2, 5)  
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("Beta")
    axes.plot(T, beta, "g-", alpha=0.6)

    # with error bars
    axes = figure.add_subplot(3, 2, 2)  
    axes.set_ylabel("Alpha")
    axes.errorbar(T, alpha, yerr = alpha_std,  color='blue', ecolor='lightblue', alpha=0.6)   
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(3, 2, 4)  
    axes.set_ylabel("Memory")
    axes.errorbar(T, memory, yerr = memory_std, color='red', ecolor='mistyrose', alpha=0.6)
    axes.set_xticks([])      # avoiding the mixing of xticks with another subplot

    axes = figure.add_subplot(3, 2, 6)  
    axes.set_xlabel("Time (h)")
    axes.set_ylabel("Beta")
    axes.errorbar(T, beta, yerr = beta_std, color='green', ecolor='lightgreen', alpha=0.6)

    pylab.savefig("genes.pdf", format="pdf", bbox_inches="tight") 
    pylab.savefig("genes.jpeg", format="jpeg", bbox_inches="tight")

    # Figure 3. host genotype (alpha) distribution vs time
    print("Generating Host_Genotype Alpha.pdf/jpeg...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(alphaDist, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("alpha")
    ticks = numpy.arange(0, bins, 10)
    show_tricks2 = alphaBinList[0:-1:10]      
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks2])
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("alphaDist.pdf", format="pdf", bbox_inches="tight")
    pylab.savefig("alphaDist.jpeg", format="jpeg", bbox_inches="tight")

    # Figure 4. host mass distribution vs time
    print("Generating Host_Mass.pdf/jpeg...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(HOST_MASS, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Host mass")
    ticks = numpy.arange(0, bins, 10)
    show_tricks3 = massBinList[0:-1:10]       
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks3])           
    cb = pylab.colorbar(img)
    cb.set_label("# of hosts")
    pylab.savefig("Host_Mass.pdf", format="pdf", bbox_inches="tight")
    pylab.savefig("Host_Mass.jpeg", format="jpeg", bbox_inches="tight")

    # Figure 5. virus genotype (memory) distribution vs time
    print("Generating Virus_Genotype memory.pdf/jpeg...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(memoryDist, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("Virus memory")
    ticks = numpy.arange(0, bins, 10)
    show_tricks4 = memoryBinList[0:-1:10]        # new merge
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks4])
    cb = pylab.colorbar(img)
    cb.set_label("# of viruses")
    pylab.savefig("memoryDist.pdf", format="pdf", bbox_inches="tight")
    pylab.savefig("memoryDist.jpegf", format="jpeg", bbox_inches="tight")

    # Figure 6. virus genotype (beta) distribution vs time
    print("Generating Virus_Genotype beta.pdf/jpeg...")
    pylab.figure(figsize=(8, 6), dpi=500)
    img = pylab.imshow(betaDist, norm=colors.SymLogNorm(linthresh = 1) , cmap='jet', origin='lower')
    pylab.axis("tight")
    img.axes.set_xlabel("Time (h)")
    img.axes.set_ylabel("beta")
    ticks = numpy.arange(0, bins, 10)
    show_tricks4 = betaBinList[0:-1:10]       
    pylab.yticks(ticks, [float('%.2f' %(show_item)) for show_item in show_tricks4])
    cb = pylab.colorbar(img)
    cb.set_label("# of viruses")
    pylab.savefig("betaDist.pdf", format="pdf", bbox_inches="tight")
    pylab.savefig("betaDist.jpeg", format="jpeg", bbox_inches="tight")

    # Figure 7. time evolution of infection map as a movie
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
        im = ax0.imshow(INFECTION_MAP[t], extent=(min_Beta, max_Beta, min_Alpha, max_Alpha), norm=norm_Log , cmap='jet', origin='lower', aspect='auto')
        fig.colorbar(im,cax=cax0)
        cax0.set_ylabel('# of infections', rotation=270, fontsize= 8)

        ax0.set_ylabel("Host genotype (Alpha)")
        ax0.set_xlabel("Virus genotype (Beta)")
        ax0.set_title(f"Epoch {t}")
        
        return fig, 
    
    simulation = animation.FuncAnimation(fig, updateHist, frames=fr, blit=True)
    simulation.save('infections.mp4', fps=5, dpi=200)  
    
    # Figure 7. Entropy of Host (alpha) & Virus (beta, Memory) Genotypes
    # calculate Entropy using python package
    EntrP_alpha = entropy(alphaDist, base=2)
    EntrP_beta = entropy(betaDist, base=2)
    EntrPy_memory = entropy(memoryDist, base=2)
    
    # plotting parts
    print("Generating Entropy.pdf...")
    figure = pylab.figure(figsize=(8, 6), dpi=500)
    size = 10
    pylab.rcParams["axes.titlesize"] = size
    pylab.rcParams["axes.labelsize"] = size
    pylab.rcParams["xtick.labelsize"] = size
    pylab.rcParams["ytick.labelsize"] = size
    pylab.rcParams["legend.fontsize"] = size

    axes = figure.add_subplot(3, 1, 1)  # Figure 1.1. DIP vs time
    axes.set_ylabel("Entropy (alpha)")
    axes.plot(T, EntrP_alpha, "b-", alpha=0.6)
    axes.set_xticks([])

    axes = figure.add_subplot(3, 1, 2)  # Figure 1.3. virus abundance vs time
    axes.set_ylabel("Entropy (beta)")
    axes.plot(T, EntrP_beta, "r-", alpha=0.6)
    axes.set_xticks([])

    axes = figure.add_subplot(3, 1, 3)  # Figure 1.3. virus abundance vs time
    axes.set_ylabel("Entropy (memory)")
    axes.plot(T, EntrPy_memory, "r-", alpha=0.6)
    axes.set_xlabel("Time (h)")

    pylab.savefig("Entropy.pdf", format="pdf", bbox_inches="tight")
    pylab.savefig("Entropy.jpeg", format="jpeg", bbox_inches="tight")
    
if __name__ == "__main__":
    main(sys.argv)
