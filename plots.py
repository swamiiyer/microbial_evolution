import numpy, gzip, pickle, sys
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.cm as cm

def symmetric(sorted_streams, stream_bounds):
    """Symmetric baseline"""
    lb, ub = numpy.min(stream_bounds[:,0,:],axis=0), numpy.max(stream_bounds[:,1,:],axis=0)
    return .5*(lb+ub)

def pos_only(sorted_streams, stream_bounds):
    """Lumps will only be positive"""
    lb, ub = numpy.min(stream_bounds[:,0,:],axis=0), numpy.max(stream_bounds[:,1,:],axis=0)
    return lb

def zero(sorted_streams, stream_bounds):
    """Zero baseline"""
    return numpy.zeros(stream_bounds.shape[2])

def min_weighted_wiggles(sorted_streams, stream_bounds):
    """Baseline recommended by Byron and Wattenberg"""
    lb, ub = numpy.min(stream_bounds[:,0,:],axis=0), numpy.max(stream_bounds[:,1,:],axis=0)
    weight = ub-lb
    
    sorted_streams = numpy.abs(sorted_streams)
    for i in xrange(len(sorted_streams)):
        sorted_streams[i,:] *= (-1)**i
    cusum_f = numpy.vstack((numpy.zeros(sorted_streams.shape[1]),
                        numpy.cumsum(sorted_streams[:-1,:], axis=0)))
    f_prime = numpy.diff(sorted_streams, axis=1)
    cusum_f_prime = numpy.diff(cusum_f, axis=1)
    g_prime = numpy.hstack(([0],-numpy.sum((f_prime*.5  + cusum_f_prime)*sorted_streams[:,1:],axis=0) / weight[1:]))
    g_prime[numpy.where(weight==0)] = 0
    g = numpy.cumsum(g_prime)
    return g

def stacked_graph(streams, cmap=plt.cm.bone, color_seq='linear', baseline_fn=min_weighted_wiggles):
    """
    Produces stacked graphs using matplotlib.
    
    Reference: 'Stacked graphs- geometry & aesthetics' by Byron and Wattenberg
    http://www.leebyron.com/else/streamgraph/download.php?file=stackedgraphs_byron_wattenberg.pdf
    
    Parameters:
      - streams: A list of time-series of positive values. Each element must be of the same length.
      - cmap: A matplotlib color map. Defaults to 'bone'.
      - colo_seq: 'linear' or 'random'.
      - baseline_fn: Current options are symmetric, pos_only, zero and min_weighted_wiggles.
    """
    # Sort by onset times
    onset_times = [numpy.where(numpy.abs(stream)>0)[0][0] for stream in streams]
    order = numpy.argsort(onset_times)
    streams = numpy.asarray(streams)
    sorted_streams = streams[order]
    
    t = numpy.arange(streams.shape[1])
    
    # Establish bounds
    stream_bounds = [numpy.vstack((numpy.zeros(streams.shape[1]), sorted_streams[0])),
                    numpy.vstack((-sorted_streams[1], (numpy.zeros(streams.shape[1]))))]
    side = -1
    for stream in sorted_streams[2:]:
        side *= -1
        if side==1:
            stream_bounds.append(numpy.vstack((stream_bounds[-2][1], stream_bounds[-2][1]+stream)))
        else:
            stream_bounds.append(numpy.vstack((stream_bounds[-2][0]-stream, stream_bounds[-2][0])))
            
    stream_bounds = numpy.array(stream_bounds)
    
    # Compute baseline
    baseline = baseline_fn(sorted_streams, stream_bounds)
    
    # Choose colors
    t_poly = numpy.hstack((t,t[::-1]))
    if color_seq=='linear':
        colors = numpy.linspace(0,1,streams.shape[1])
    elif color_seq=='random':
        colors = numpy.random.random(size=streams.shape[1])
    else:
        raise ValueError, 'Color sequence %s unrecognized'%color_seq
    
    # Plot    
    plt.axis('off')        
    for i in xrange(len(stream_bounds)):
        bound = stream_bounds[i]
        color = cmap(colors[i])
        plt.fill(t_poly, numpy.hstack((bound[0]-baseline,(bound[1]-baseline)[::-1])), facecolor=color, linewidth=0.,edgecolor='none')
        
def main(args):
    fh = gzip.open("results.pklz", "rb")
    params = pickle.load(fh)
    T = range(params["epochs"])
    biomass_sim, DIP = [], []
    V_size, H_size = [], []
    V_sp, V_st, H_sp, H_st = [], [], [], []
    V_sp_min, V_sp_max, V_st_min, V_st_max, H_sp_min, H_sp_max, H_st_min, H_st_max = 0, 0, 0, 0, 0, 0, 0, 0
    for t in T:
        print "pass 1, epoch %d..." %(t)
        V_pop, H_pop = pickle.load(fh)
        V_size.append(len(V_pop))
        H_size.append(len(H_pop))
        mass, sp, st = 0.0, set(), set()
        for host in H_pop:
            mass += host.mass
            sp.add(host.species)
            st.add(host.strain)
            H_sp_min = host.species if host.species < H_sp_min else H_sp_min
            H_sp_max = host.species if host.species > H_sp_max else H_sp_max
            H_st_min = host.strain if host.strain < H_st_min else H_st_min
            H_st_max = host.strain if host.strain > H_st_max else H_st_max                
        H_sp.append(len(sp))
        H_st.append(len(st))
        biomass_sim.append(mass)
        DIP.append(params["P_tot"] - mass)
        sp, st = set(), set()
        for virus in V_pop:
            sp.add(virus.species)
            st.add(virus.strain)
            V_sp_min = virus.species if virus.species < V_sp_min else V_sp_min
            V_sp_max = virus.species if virus.species > V_sp_max else V_sp_max
            V_st_min = virus.strain if virus.strain < V_st_min else V_st_min
            V_st_max = virus.strain if virus.strain > V_st_max else V_st_max
        V_sp.append(len(sp))
        V_st.append(len(st))
    fh.close()

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(1, figsize = (7, 4.5), dpi = 500)
    plt.xlabel(r"time $t$", fontproperties = font_prop)
    plt.ylabel(r"normalized mass $m_N$", 
                 fontproperties = font_prop)
    plt.semilogy(T, DIP, "b-", linewidth = 2, alpha = 0.6, label = r"DIP")
    plt.semilogy(T, biomass_sim, "g-", linewidth = 2, alpha = 0.6, label = r"Biomass")
    plt.xlim(0, max(T))
    plt.legend(loc = 'upper right', prop = font_prop)
    plt.savefig("figure1.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(1)

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(2, figsize = (7, 4.5), dpi = 500)
    plt.xlabel(r"time $t$", fontproperties = font_prop)
    plt.ylabel(r"total population $N$", 
                 fontproperties = font_prop)
    plt.semilogy(T, V_size, "g-", linewidth = 2, alpha = 0.6, label = r"Virus ($N_V$)")
    plt.semilogy(T, H_size, "b-", linewidth = 2, alpha = 0.6, label = r"Host ($N_H$)")
    plt.xlim(0, max(T))
    plt.legend(loc = 'upper right', prop = font_prop)
    plt.savefig("figure2.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(2)

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(3, figsize = (7, 4.5), dpi = 500)
    plt.xlabel(r"time $t$", fontproperties = font_prop)
    plt.ylabel(r"richness $n$", 
                 fontproperties = font_prop)
    plt.plot(T, V_sp, "b-", linewidth = 2, alpha = 0.6, label = r"Virus ($n^V_{sp}$)")
    plt.plot(T, V_st, "g-", linewidth = 2, alpha = 0.6, label = r"Virus ($n^V_{st}$)")
    plt.plot(T, H_sp, "m-", linewidth = 2, alpha = 0.6, label = r"Host ($n^H_{sp}$)")
    plt.plot(T, H_st, "r-", linewidth = 2, alpha = 0.6, label = r"Host ($n^H_{st}$)")
    plt.xlim(0, max(T))
    plt.ylim(0, max(max(V_st), max(V_sp), max(H_st), max(H_sp)))
    plt.legend(loc = 'upper right', prop = font_prop)
    plt.savefig("figure3.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(3)


    # font_prop = font_manager.FontProperties(size = 12)
    # plt.figure(4, figsize = (7, 4.5), dpi = 500)
    # fh = open("results.pkl", "r")
    # pickle.load(fh) # params
    # dsets = numpy.zeros((H_sp_max - H_sp_min + 1, params["epochs"]))
    # for t in T:
    #     print "pass 2, epoch %d..." %(t)
    #     V_pop, H_pop = pickle.load(fh)
    #     for i, sp in enumerate(range(H_sp_min, H_sp_max + 1)):
    #         count = sum([1 for host in H_pop if host.species == sp])
    #         dsets[i, t] = count
    # stacked_graph(dsets.tolist(), baseline_fn = min_weighted_wiggles, color_seq='random')
    # plt.savefig("figure4.pdf", format = "pdf", bbox_inches = "tight")
    # plt.close(4)

    fh = gzip.open("results.pklz", "rb")
    pickle.load(fh) # params
    H_dsets = numpy.zeros((H_sp_max - H_sp_min + 1, params["epochs"]))
    V_dsets = numpy.zeros((V_sp_max - V_sp_min + 1, params["epochs"]))
    H_survivors = []
    for t in T:
        print "pass 2, epoch %d..." %(t)
        V_pop, H_pop = pickle.load(fh)
        for i, sp in enumerate(range(H_sp_min, H_sp_max + 1)):
            count = sum([1 for host in H_pop if host.species == sp])
            H_dsets[i, t] = count
            if t == params["epochs"] - 1:
                freq = numpy.zeros(H_st_max - H_st_min + 1)
                for j, st in enumerate(range(H_st_min, H_st_max + 1)):
                    count = sum([1 for host in H_pop if host.species == sp and host.strain == st]) 
                    freq[j] = count
                H_survivors.append(freq)

        for i, sp in enumerate(range(V_sp_min, V_sp_max + 1)):
            count = sum([1 for virus in V_pop if virus.species == sp])
            V_dsets[i, t] = count

    print len(H_survivors), H_survivors[0].shape

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(4, figsize = (7, 4.5), dpi = 500)
    im = plt.imshow(H_dsets, interpolation = "nearest", origin = "l", 
                    cmap = cm.RdYlGn, extent = [0, params["epochs"], H_sp_min, H_sp_max])
    plt.colorbar(im)
    plt.savefig("figure4.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(4)

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(5, figsize = (7, 4.5), dpi = 500)
    im = plt.imshow(V_dsets, interpolation = "nearest", origin = "l", 
                    cmap = cm.RdYlGn, extent = [0, params["epochs"], V_sp_min, V_sp_max])
    plt.colorbar(im)
    plt.savefig("figure5.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(5)

    font_prop = font_manager.FontProperties(size = 12)
    plt.figure(6, figsize = (7, 4.5), dpi = 500)
    im = plt.imshow(H_survivors, interpolation = "nearest", origin = "l", 
                    cmap = cm.RdYlGn, extent = [H_st_min, H_st_max, H_sp_min, H_sp_max])
    plt.colorbar(im)
    plt.savefig("figure6.pdf", format = "pdf", bbox_inches = "tight")
    plt.close(5)

if __name__ == "__main__":
    main(sys.argv[1:])
