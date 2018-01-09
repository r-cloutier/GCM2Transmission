from imports import *


def setup_Earthlike_EOS(Parr, Tarr, X_H2Oarr):
    '''
    Create a custom EOS file based on the approximate terrestrial atmospheric 
    abundances and variable water content for a vertical column.
    '''
    # get template EOS file
    fname = '/Users/ryancloutier/Research/Exo_Transmit/EOS/eos_1Xsolar_gas.dat'
    f = open(fname, 'r')
    g = f.readlines()
    f.close()
    col_names = np.ascontiguousarray(g[0].split('\t\t')[:-1])

    # create a unique EOS
    for j in range(2,len(g)):

        # add abundances
        if len(g[j]) == 495:

            # add fixed chemical abundances
            Xs = np.ascontiguousarray(g[j].split('\t')[:-1])
            Xs = _set_abundances(Xs, col_names)
                
            # add water
            T = float(g[j].split('\t')[0])  # K
            P = float(g[j].split('\t')[1])  # Pa
            Xs[16] = _set_water_abundance(X_H2O_cell, Xs[2:])

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


def _set_water_abundance(Xs, X_H2O_cell, P_cell, T_cell):
    '''
    Set water abundance in EOS file at the expense of the other chemical 
    constituents.
    '''
    assert Xs.size == 36
    
    
