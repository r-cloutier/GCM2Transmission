from imports import *
import rvs
import create_exotransmit as exo

global mp, kb, rp, g, Rs, path2exotransmit, path2plasim
mp, kb = 1.6726219e-27, 1.38064852e-23
rp, g, Rs = 1., 9.8, .2
path2exotransmit = '/Users/ryancloutier/Research/Exo_Transmit'
path2plasim = '/Users/ryancloutier/Research/PlaSim'


#main('TIDAL1.0.001', outpathprefix='GCM/terminator', outfile='Tidallylocked.dat')
def main(simname, t=39, outpathprefix='GCM/terminator', outfile='GCMtidallylocked.dat'):
    '''
    Compute the transmission spectrum from a GCM by setting up and running 
    ExoTransmit for each vertical atmospheric column in the GCM and combining 
    them using a weighted-mean.

    Parameters
    ----------
    `t': scalar
        Time in GCM units at which to compute the transmission spectrum

    '''
    # get GCM data
    time, lon, lat, _, Ps, P, T,_,_ = \
                            _get_GCM_data('plasim_samples/%s.nc'%simname)
    if t not in time:
        raise ValueError('t not in GCM time array.')

    # compute atmospheric depth
    depth = _P2h(P, Ps, T)

    # setup and run exotransit from each vertical column
    _, Nh, Nlat, Nlon = T.shape
    clean, t = True, int(t)
    mass = np.zeros((Nh-1, Nlat-1, Nlon-1))
    for i in range(Nlat-1):
        for j in range(Nlon-1):
            
            # check that transmission occurs through this column
            #if _is_transmission(lat[i], lon[j], depth[time==t,:,i,j].max()):
            if _at_terminator(lat[i], lon[j]):
            
                # get column mass for weighting coefficients
                #mass[:,i,j] = _get_mass(P, T, depth, lat, lon, t, i, j)
                mass[:,i,j] = 1.
                
                # compute transmission spectrum
		outpath = '%s/%s'%(outpathprefix, outfile.replace('.dat',''))
                _setup_exotransmit(simname, t, i, j,
                                   outfile='%s/%s_%i_%i.dat'%(outpath,outfile.replace('.dat',''),i,j))
		sys.exit('stop')
                _run_exotransmit(clean)
                clean = False
                
    # compute the mass-coefficient for each column
    mass = np.sum(mass, 0)
    coeffs = mass / mass.sum()
    assert coeffs.sum() == 1
    hdu = fits.PrimaryHDU(coeffs)
    hdu.writeto('%s/Spectra/%s/coefficients_%s.fits'%(path2exotransmit,outpath,simname[:-4]), overwrite=True)
    
    # compute the master transmission spectrum
    # ie: send rays at fixed (y,z) and add up the mass weighted transmission
    # spectra
    wl, spectrum, hdr = _coadd_spectra(coeffs, '%s/%s'%(outpath, outfile[:-4]))
    np.savetxt('%s/Spectra/%s/%s'%(path2exotransmit, outpath, outfile),
               np.array([wl, spectrum]).T, fmt='%.6e', delimiter='\t',
               header=hdr)

    

def _setup_exotransmit(simname, tindex, latindex, lonindex,
                       outfile='default.dat'):
    '''
    Setup ExoTransmit files to compute the transmission spectrum for a single 
    vertical column from the GCM.

    Parameters
    ----------
    `t': scalar
        Time in GCM units at which to compute the transmission spectrum
    `latitude': scalar
        Latitude of the vertical column in degrees from the equator
    `longitude': scalar
        Longitude of the vertical column in degrees from the substellar point

    '''
    # get GCM data
    _,_,_,_,_,P,T,X_H2O,clouds = _get_GCM_data('plasim_samples/%s.nc'%simname)
    Ntime, Nh, Nlat, Nlon = T.shape

    # create exotransmit EOS file for this column
    exo.setup_Earthlike_EOS(P, X_H2O[tindex,:,latindex,lonindex])
    #P_layer = P[8]
    #exo.setup_singlelayer_Earthlike_EOS(P, X_H2O[tindex,:,latindex,lonindex],
    #                                    P_layer)

    # create exotransmit TP profile for this column
    exo.setup_TP_file(P, T[tindex,:,latindex,lonindex])

    # get cloudtop pressure from the cloud-weighted-mean pressure 
    #cloud_col = clouds[tindex,:,latindex,lonindex] + 0
    #cloud_col = cloud_col/cloud_col.sum() if cloud_col.sum() > 0 else 0.
    #cloudP = np.sum(P*cloud_col)
    cloudP = 5e3   # TEMP    

    # create exotransmit input files
    exo.create_input_file(g, rp, Rs, cloudP, outfile)
    exo.create_chem_file()

    
def _run_exotransmit(clean=False):
    '''
    Compute the transmission spectrum using ExoTransmit.
    '''
    os.chdir(path2exotransmit)
    if clean:
        os.system('make clean')
        os.system('make')
    os.system('./Exo_Transmit')
    os.chdir(path2plasim)
    

def _get_GCM_data(fname):
    '''
    Read-in temporal, spatial, and thermodynamic variables from the netcdf 
    output of the GCM.
    '''
    data = nc.Dataset(fname, 'r').variables
    time, lon, lat, lev = data['time'][:], data['lon'][:], data['lat'][:], \
                          data['lev'][:]
    T, Ps, X_H2O = data['ta'][:], data['ps'][:], data['hus'][:]
    clouds = data['cl'][:]
    P, Ps = lev*10, Ps*10  # hPa to Pa
    h = _P2h(P, Ps, T)
    return time, lon, lat, h, Ps, P, T, X_H2O, clouds


def _is_transmission(lat, lon, H):
    '''
    Check that transmission occurs at least partially through a vertical column 
    with a specified latitude and longitude.

    Parameters
    ----------
    `lat': scalar
        Latitude of the vertical column in degrees from the equator
    `lon': scalar
         Longitude of the vertical column in degrees from the substellar point
    `H': scalar
        Height of the GCM atmosphere in Earth radii

    '''
    phi, theta = _geo2sphere(lat, lon)
    y, z = _sphere2yz(rp+H, phi, theta)
    return np.sqrt(y*y + z*z) > rp


def _at_terminator(lat, lon):
    '''
    Check if this column is at the terminator (i.e. x=0)
    '''
    phi, theta = _geo2sphere(lat, lon)
    x = np.sin(theta)*np.cos(phi)
    return np.isclose(x, 0, atol=1e-4)
    
    
def _P2h(P, Ps, T, mu=28.97):
    '''
    Convert the atmospheric pressure to a depth in Earth radii.

    Parameters
    ----------
    `P': array-like
        Atmospheric pressure as a function of depth in Pascals
    `T': array-like
        4D atmospheric temperature in kelvin

    '''
    # reshape arrays
    Ntime, NP, Nlat, Nlon = T.shape
    P4 = np.repeat(np.repeat(np.repeat(P[np.newaxis,:], Ntime,
                                       axis=0)[:,:,np.newaxis], Nlat,
                             axis=2)[:,:,:,np.newaxis], Nlon, axis=3)
    assert P4.shape == T.shape
    Ps4 = np.repeat(Ps[:,np.newaxis,:,:], NP, axis=1)    
    assert Ps4.shape == T.shape

    # compute depth vs pressure
    H = kb*T / (g*mu*mp)
    return rvs.m2Rearth(H * np.log(Ps4 / P4))


def _geo2sphere(lat, lon):
    '''
    Convert latitude and lonitude into spherical angles.
    '''
    phi = np.deg2rad(lon)
    theta = np.deg2rad(90. - lat)
    return phi, theta


def _geo2cart(r, lat, lon):
    phi, theta = _geo2sphere(lat, lon)
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)
    return x, y, z


def _sphere2cart2D(r, theta):
    '''
    Convert radial/azimuthal coordinates to cartesian x,y.
    '''
    x, y = r*np.cos(theta), r*np.sin(theta)
    return x, y 


def _sphere2yz(r, phi, theta):
    '''
    Convert radial/azimuthal coordinates to cartesian x,y.
    '''
    y, z = r*np.sin(theta)*np.sin(phi), r*np.cos(theta)
    return y, z


def _cart2geo3D(x, y, z):
    '''
    Convert 3D cartesian coordinates to longitudes, latitudes, and heights.
    '''
    r, phi, theta = _cart2sphere3D(x, y, z)
    h   = r
    lon = np.rad2deg(phi)
    lat = np.rad2deg(np.pi/2 - theta)
    return h, lon, lat


def _cart2sphere3D(x, y, z):
    '''
    Convert 3D cartesian coordinates to spherical coordinates.
    '''
    r     = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi   = np.arctan(y / x)
    return r, phi, theta


def _get_3d_cartesian_map(depth, lon, lat, tindex):
    depth3d = np.mean(np.array([depth[tindex,:-1,:-1,:-1],
                                depth[tindex,1:,1:,1:]]), 0)
    Nh, Nlat, Nlon = depth3d.shape
    latcell = np.mean(np.array([lat[:-1],  lat[1:]]), 0)
    loncell = np.mean(np.array([lon[:-1],  lon[1:]]), 0)
    lon3d = np.stack([np.stack([loncell for i in range(Nlat)])
                      for i in range(Nh)])
    lat3d = np.stack([np.stack([latcell for i in range(Nlon)]).T
                      for i in range(Nh)])

    # convert to spherical coords
    r = rp + depth3d
    phi = np.deg2rad(lon3d)
    theta = np.deg2rad(90. - lat3d)

    # get cartesian map
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)
    return x, y, z
    


def _get_masscoeff_grid_FAIL(tindex, Ps, P):
    '''
    Compute the grid of dP/g (proportional to mass) and normalize to get the 
    mass-weighted coefficients as a function of lat, lon, and depth.
    '''
    Ntime, Nlat, Nlon = Ps.shape
    NP = P.size
    Pgrid3d = np.stack([np.stack([P for i in range(Nlat)], 0)
                        for i in range(Nlon)], 0).T
    Psgrid3d = np.stack([Ps[tindex] for i in range(NP)])
    P = np.insert(Pgrid3d, 0, Psgrid3d, axis=1)
    dP = np.diff(Pgrid, axis=0)
    return dP/g


def _get_mass(P, T, depth, lat, lon, t, i, j, folder='VolumeCalculations'):
    Pcell = P[:-1] - P[1:]
    Tcell = T[t,:-1,i,j] - T[t,1:,i,j]
    rhocell = Pcell / Tcell
    x1, y1, z1  = _geo2cart(depth[t,:-1,i,j]-depth[t,1:,i,j],
                            lat[i+1], lon[j+1])
    x2, y2, z2  = _geo2cart(depth[t,:-1,i+1,j+1]-depth[t,1:,i+1,j+1],
                            lat[i+1], lon[j+1])
    Vcell = abs(x1-x2) * abs(y1-y2) * abs(z1-z2)
    return rhocell * Vcell

    
def _get_masscoeff_grid(P, Ps, T, depth, lat, lon, tindex, N=1e7):
    '''
    Compute the grid of unnormalized cell masses assuming an ideal gas and 
    renormalizing to get the mass-weighted coefficients as a function of lat, 
    lon, and depth.
    '''
    Ntime, Nh, Nlat, Nlon = T.shape
    assert P.size == Nh
    Pcell = np.array([np.mean([P[i], P[i+1]]) for i in range(P.size-1)])
    P3d = np.stack([np.stack([Pcell for i in range(Nlat-1)], 0)
                    for i in range(Nlon-1)], 0).T
    T3d = np.mean(np.array([T[tindex,:-1,:-1,:-1], T[tindex,1:,1:,1:]]), 0)
    assert T3d.shape == P3d.shape
    
    # compute cell masses modolo a constant
    rho3d = P3d / T3d
    #V3d = _compute_cell_volumes(depth, lat, lon, tindex, N)
    V3d = _read_cell_volumes(depth)
    mass3d = rho3d * V3d

    # compute mass-weighted coefficients
    coeffs = mass3d / mass3d.sum()
    return coeffs


def _read_cell_volumes(depth, folder='VolumeCalculations'):
    Ntime, Nh, Nlat, Nlon = depth.shape
    V3d = np.zeros((Nh-1, Nlat-1, Nlon-1))
    for i in range(Nh-1):
        d = fits.open('%s/V3d_TIDAL1.0.001.nc_depth%.2d'%(folder, i))[0].data
        V3d[i] = d.reshape(Nlat-1, Nlon-1)
    return V3d


def _compute_cell_volumes(depth, lat, lon, tindex, N=1e7):
    '''
    Compute the 3D spatial map of cell volumes.
    '''
    Ntime, Nh, Nlat, Nlon = depth.shape
    assert lat.size == Nlat
    assert lon.size == Nlon

    V3d = np.zeros((Nh-1, Nlat-1, Nlon-1))
    for i in range(Nh-1):
        for j in range(Nlat-1):
            for k in range(Nlon-1):
                cell_bnds1 = depth[tindex,i,j,k], lat[j], lon[k]
                cell_bnds2 = depth[tindex,i+1,j+1,k+1], lat[j+1], lon[k+1]
                V3d[i,j,k] = _compute_cell_V(cell_bnds1, cell_bnds2,
                                             depth.max(), N)

    return V3d


def _compute_cell_V(cell_bnds1, cell_bnds2, H, N=1e7):
    '''
    Compute the approximate cell volume by drawing spherical coordinates 
    within the planet's atmosphere (of known volume) and computing the fraction
    of points that lie within the cell bounded by the specified depths, 
    latitudes, and longitudes. 
    '''    
    # get total atmosphere volume
    assert rp > H
    V = 4*np.pi/3 * ((rp+H)**3 - rp**3)

    # draw random locations in the atmosphere
    N = int(N)
    hs   = np.random.uniform(0, H, N)
    lats = np.random.uniform(-90, 90, N)
    lons = np.random.uniform(0, 360, N)

    # read boundaries of the cell whose volume is being measured
    h1, lat1, lon1 = cell_bnds1
    h2, lat2, lon2 = cell_bnds2
    hcell, latcell, loncell = np.append(h1,h2), np.append(lat1,lat2), \
                              np.append(lon1,lon2)

    # compute cell volume from fraction of draws within the cell
    in_cell = (hs >= hcell.min()) & (hs <= hcell.max()) & \
              (lats >= latcell.min()) & (lats <= latcell.max()) & \
              (lons >= loncell.min()) & (lons <= loncell.max())
    Vcell = V * in_cell.sum() / float(N)
    return Vcell


def _coadd_spectra(coeffs, prefix):
    '''Get the GCM spectra and compute the weighted mean.'''
    # initialize spectrum array
    fs = glob.glob('%s/Spectra/%s_*.dat'%(path2exotransmit, prefix))
    wl, spectrum = np.loadtxt(fs[0], skiprows=2).T
    spectrum = np.zeros(spectrum.size)

    # get coefficients only where a transmission spectru was calculated
    Nlat, Nlon = coeffs.shape
    coeffsv2 = np.zeros_like(coeffs)
    for i in range(Nlat):
        for j in range(Nlon):
	    fname = '%s/Spectra/%s_%i_%i.dat'%(path2exotransmit, prefix, i, j)
	    try:
		_,spec = np.loadtxt(fname, skiprows=2).T
		coeffsv2[i,j] = coeffs[i,j]
	    except IOError:
		pass
    coeffsv2 = coeffsv2 / coeffsv2.sum()
    assert coeffsv2.sum() == 1    

    # coadd mass-weighted spectra
    nspectra = 0.
    for i in range(Nlat):
        for j in range(Nlon):
	    try:		
		fname = '%s/Spectra/%s_%i_%i.dat'%(path2exotransmit, prefix, i, j)
                _,spec = np.loadtxt(fname, skiprows=2).T
		assert np.all(np.isfinite(spec))
                spectrum += coeffsv2[i,j] * spec
            except IOError:
                pass
    
    # get header
    fs = glob.glob('%s/Spectra/%s_*.dat'%(path2exotransmit, prefix))
    f = open(fs[0], 'r')
    g = f.readlines()
    f.close()
    hdr = ''.join(g[:2])[:-1]
    return wl, spectrum, hdr
