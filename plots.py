import numpy, pandas, sys
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

def main(args):
    data = pandas.read_table(args[0], header = None, skiprows = 1)

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(1, figsize = (7, 4.5), dpi = 500)
    plt.xlabel(r"time $t$", fontproperties = font_prop)
    plt.ylabel(r"mass units", 
                 fontproperties = font_prop)
    T = range(len(data))
    DIP, biomass_sim, ptot_sim = list(data[0]), list(data[1]), list(data[2])

    plt.semilogy(T, DIP, "b-", linewidth = 2, alpha = 0.6, label = r"DIP")
    plt.semilogy(T, biomass_sim, "g-", linewidth = 2, alpha = 0.6, label = r"biomass (simulated)")
    plt.semilogy(T, ptot_sim, "r-", linewidth = 2, alpha = 0.6, label = r"ptot (simulated)")
    plt.xlim(0, max(T))
    plt.legend(loc = 'upper right', prop = font_prop)
    plt.savefig("figure1.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(1)


    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(2, figsize = (7, 4.5), dpi = 500)
    plt.xlabel(r"time $t$", fontproperties = font_prop)
    plt.ylabel(r"population", 
                 fontproperties = font_prop)
    V_pop, H_pop = list(data[3]), list(data[4])
    plt.semilogy(T, V_pop, "g-", linewidth = 2, alpha = 0.6, label = r"Virus")
    plt.semilogy(T, H_pop, "b-", linewidth = 2, alpha = 0.6, label = r"Host")
    plt.xlim(0, max(T))
    plt.legend(loc = 'upper right', prop = font_prop)
    plt.savefig("figure2.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(2)

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(3, figsize = (7, 4.5), dpi = 500)
    plt.xlabel(r"time $t$", fontproperties = font_prop)
    plt.ylabel(r"richness", 
                 fontproperties = font_prop)
    V_strains, V_species, H_strains, H_species = list(data[5]), list(data[6]), list(data[7]), list(data[8]), 
    plt.plot(T, V_strains, "g-", linewidth = 2, alpha = 0.6, label = r"Virus strains")
    plt.plot(T, V_species, "b-", linewidth = 2, alpha = 0.6, label = r"Virus species")
    plt.plot(T, H_strains, "r-", linewidth = 2, alpha = 0.6, label = r"Host strains")
    plt.plot(T, H_species, "m-", linewidth = 2, alpha = 0.6, label = r"Host species")
    plt.xlim(0, max(T))
    plt.legend(loc = 'upper right', prop = font_prop)
    plt.savefig("figure3.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(2)

if __name__ == "__main__":
    main(sys.argv[1:])
