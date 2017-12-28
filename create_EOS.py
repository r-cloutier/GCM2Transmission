from imports import *


def setup_Earthlike_EOS(P_cell, T_cell, H2O_cell):
    '''
    Create a custom EOS file based on the approximate terrestrial atmospheric 
    abundances and variable water content for each cell that a ray is 
    transmitted through with a fixed TP-profile.
    '''
    assert P_cell.size == T_cell.size
    assert T_cell.size == H2O_cell.size

    # get template EOS file
    fname = '/Users/ryancloutier/Research/Exo_Transmit/EOS/eos_1Xsolar_gas.dat'
    f = open(fname, 'r')
    g = f.readlines()
    f.close()
    col_names = np.ascontiguousarray(g[0].split('\t\t')[:-1])

    # create unique EOS file for each cell
    for i in range(P_cell.size):
        
        for j in range(2,len(g)):

            # add abundances
            if len(g[j]) == 495:
                
                # add fixed chemical abundances
                Xs = np.ascontiguousarray(g[j].split('\t')[:-1])
                Xs = _set_abundances(Xs, col_names)
                
                # add water
                T = float(g[j].split('\t')[0])  # K
                P = float(g[j].split('\t')[1])  # Pa
                Xs[16] = _get_water_abundance(T, P, T_ray[i], P_ray[i],
                                              H2O_ray[i])

                # check abundances
                if Xs.sum() > 1:
                    raise ValueError("Cumulative atmospheric abundances "+ \
                                     "exceed 1.")

            # Convert to exotransmit format
            
            


def _set_abundances(Xs, col_names):
    '''
    Set abundances to Earth-like with the exception of water which is done 
    separately.
    '''
    assert col_names.size == 38
    Xout = np.zeros(38)
    for i in range(38):

        if col_names[i] in ['T','P'] :
            Xout[i] = float(Xs[i])
        elif col_names[i] == 'CH4':
            Xout[i] = 1.7e-6        
        elif col_names[i] == 'CO2':
            Xout[i] = 400e-6
        elif col_names[i] == 'H2':
            Xout[i] = .55e-6
        elif col_names[i] == 'He':
            Xout[i] = 5.24e-6
        elif col_names[i] == 'N2':
            Xout[i] = .7808
        elif col_names[i] == 'O2':
            Xout[i] = .2095
        else:
            Xout[i] = 0.

    return Xout


def _get_water_abundance(T, P, T_ray, P_ray, H2O_ray):
    '''
    
    '''
