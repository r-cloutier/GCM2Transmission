from imports import *

global path2exotransmit
path2exotransmit = '/Users/ryancloutier/Research/Exo_Transmit'

def get_species_opacity(chem_str):
    # get opacity data
    try:
        f = open('%s/Opac/opac%s.dat'%(path2exotransmit, chem_str), 'r')
        g = f.readlines()
        f.close()
    except IOError:
        raise IOError('%s is not a valid ExoTransmit chemical.'%chem_str)

    # get T and P grids
    Ts = np.array([float(i) for i in g[0].split(' ')[:-1]])
    Ps = np.array([float(i) for i in g[1].split(' ')[:-1]])  # pascals
    
    # read in opacity at any T and P
    Nwl, wli = 4616, -1
    wls, opac = np.zeros(Nwl), np.zeros((Nwl, Ts.size, Ps.size))
    for i in range(2,len(g)):

        if len(g[i]) == 15:
            wli += 1
            wls[wli] = float(g[i])*1e6
            
        elif len(g[i]) in [404,405]:
            pi = np.where(Ps == 10**float(g[i].split(' ')[0]))[0][0]
            opac[wli,:,pi] = np.array([float(j) for j in g[i].split(' ')[1:-1]])

        else:
            print len(g[i])
            raise ValueError('Should not be stuck here.')
        
    return wls, opac, Ts, Ps
