from imports import *
import rvs

global path2exotransmit, EOSfile, TPfile, userfile, chemfile
EOSfile, TPfile, userfile, chemfile = 'eos_GCM.dat', 't_p_GCM.dat', \
                                      'userInput.in', 'selectChem.in'
path2exotransmit = '/Users/ryancloutier/Research/Exo_Transmit'


def setup_Earthlike_EOS(Parr, X_H2Oarr):
    '''
    Create a custom EOS file based on the approximate terrestrial atmospheric 
    abundances and variable water content for a vertical column.
    '''
    Parr, X_H2Oarr = np.ascontiguousarray(Parr), np.ascontiguousarray(X_H2Oarr)
    assert Parr.size == X_H2Oarr.size

    # get template EOS file
    fname = '%s/EOS/eos_1Xsolar_gas.dat'%path2exotransmit
    f = open(fname, 'r')
    g = f.readlines()
    f.close()
    col_names = np.ascontiguousarray(g[0].split('\t\t')[:-1])
    
    # create a unique EOS
    h = g[:4]
    for j in range(4,len(g)):

        # add abundances
        if len(g[j]) == 495:

            # add fixed chemical abundances
            Xs = np.ascontiguousarray(g[j].split('\t')[:-1])
            Xs = _set_abundances(Xs, col_names, Parr, X_H2Oarr)
            print Xs[2:].sum()
            if not np.isclose(Xs[2:].sum(), 1, 1e-5):
                raise ValueError("Cumulative atmospheric abundances "+ \
                                 "exceed 1.")

            # convert to exotransmit format
            htmp = ''
            for i in range(Xs.size):
                htmp += '%.6e\t'%Xs[i]
            h.append('%s\n'%htmp)
                
        # add 
        else:
            h.append(g[j])

    # write file
    f = open('%s/EOS/%s'%(path2exotransmit, EOSfile), 'w')
    f.write(''.join(h))
    f.close()


def _set_abundances(Xs, col_names, Parr, X_H2Oarr):
    '''
    Set abundances to Earth-like with variable water content.
    '''
    assert Xs.size == 38
    assert Xs.size == col_names.size

    # add water at this pressure level
    fill_values = float(X_H2Oarr[Parr == Parr.min()]), \
                  float(X_H2Oarr[Parr == Parr.max()])
    fint = interp1d(Parr, X_H2Oarr, bounds_error=False, fill_value=fill_values)
    P = float(Xs[1])
    X_H2O = float(fint(P))

    # get scaling factor for non-water species
    Earth_Xs = np.array([1.7e-6, 400e-6, .55e-6, 5.24e-6, .7808, .2095])
    scale = (1. - X_H2O) / Earth_Xs.sum()
    Earth_Xs *= scale
    
    Xout = np.zeros(38)
    for i in range(38):
        
        if col_names[i] in ['T','P'] :
            Xout[i] = float(Xs[i])
        elif col_names[i] == 'CH4':
            Xout[i] = Earth_Xs[0]
        elif col_names[i] == 'CO2':
            Xout[i] = Earth_Xs[1]
        elif col_names[i] == 'H2':
            Xout[i] = Earth_Xs[2]
        elif col_names[i] == 'He':
            Xout[i] = Earth_Xs[3]
        elif col_names[i] == 'N2':
            Xout[i] = Earth_Xs[4]
        elif col_names[i] == 'O2':
            Xout[i] = Earth_Xs[5]
        elif col_names[i] == 'H2O':
            Xout[i] = X_H2O
        else:
            Xout[i] = 0.

    return Xout    

    
def setup_TP_file(Parr, Tarr):
    '''
    Create a TP-profile for ExoTransmit.

    Parameters
    ----------
    `Parr': array-like

    '''
    # get template EOS file
    fname = '%s/T_P/t_p_300K.dat'%path2exotransmit
    f = open(fname, 'r')
    g = f.readlines()
    f.close()

    # create interpolation function
    Parr, Tarr = np.ascontiguousarray(Parr), np.ascontiguousarray(Tarr) 
    assert Parr.size == Tarr.size
    fill_values = float(Tarr[Parr == Parr.min()]), \
                  float(Tarr[Parr == Parr.max()])
    fint = interp1d(Parr, Tarr, bounds_error=False, fill_value=fill_values)

    # create new file with temperatures interpolated along the pressure
    h = []
    h.append(g[0])
    for i in range(1, len(g)):
        ls = np.ascontiguousarray(g[i].split(' '))
        inds = np.where(ls != '')[0]
        P = float(ls[inds[1]])  # pascals
        T = fint(P)
        ls[inds[2]] = '%.7e'%T
        h.append(' '.join(ls))

    # write file
    f = open('%s/T_P/%s'%(path2exotransmit, TPfile), 'w')
    f.write(''.join(h))
    f.close()


def create_input_file(g, rp, Rs, outfile='GCM.dat'):
    '''
    Create ExoTransmit input file.

    Parameters
    ----------
    `outfile': str
        Name of the output file containing the transmission spectrum
    `g': scalar
        Planetary surface gravity in m/s^2
    `rp': scalar
        Planetary radius in Earth radii
    `Rs': scalar
        Stellar radius in Solar radii

    '''
    # get template
    fname = '%s/userInput_template.in'%path2exotransmit
    f = open(fname, 'r')
    h = f.readlines()
    f.close()

    h[3] = '%s\n'%path2exotransmit
    h[5] = '/T_P/%s\n'%TPfile
    h[7] = '/EOS/%s\n'%EOSfile
    h[9] = '/Spectra/%s\n'%outfile
    h[11] = '%.2f\n'%g
    h[13] = '%.2e\n'%rvs.Rearth2m(rp)
    h[15] = '%.2e\n'%rvs.Rsun2m(Rs)

    # write file
    f = open('%s/%s'%(path2exotransmit, userfile), 'w')
    f.write(''.join(h))
    f.close()


def create_chem_file():
    '''
    Create ExoTransmit chemical file.
    '''
    # get template
    fname = '%s/selectChem_template.in'%path2exotransmit
    f = open(fname, 'r')
    g = f.readlines()
    f.close()

    h, chems2keep = [], ['CH4','CO2','H2','H2O','He','N2','O2']
    h.append(g[0])
    for i in range(1,len(g)-3):
        ls = np.ascontiguousarray(g[i].split(' '))
        inds = np.where(ls != '')[0]
        ls[-1] = '1\n' if ls[0] in chems2keep else '0\n'
        h.append(' '.join(ls))

    for i in range(len(g)-3,len(g)):
        h.append(g[i])
        
    # write file
    f = open('%s/%s'%(path2exotransmit, chemfile), 'w')
    f.write(''.join(h))
    f.close()
