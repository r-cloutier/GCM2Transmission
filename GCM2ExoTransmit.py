from imports import *
import rvs
import create_exotransmit as exo

global mp, kb, rp, g, Rs
mp, kb = 1.6726219e-27, 1.38064852e-23
rp, g = 1., 9.8
Rs = .2


def main(t=39):
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
    time, lon, lat, _, Ps, P, T,_ = _get_GCM_data(
        'plasim_samples/TIDAL1.0.001.nc')
    if t not in time:
        raise ValueError('t not in GCM time array.')

    # compute atmospheric depth
    depth = _P2h(P, Ps, T)

    # setup and run exotransit from each vertical column
    _,_, Nlat, Nlon = T.shape
    for i in range(Nlat):
        for j in range(Nlon):
            
            # check that transmission occurs through this column
            if _is_transmission(lat[i], lon[j], depth[time==t,:,i,j].max()):

                _setup_exotransmit(t, lat[i], lon[j],
                                   outfile='output_%i_%i.dat'%(i,j))
                clean = True if (i==0) & (j==0) else False
                _run_exotransmit(clean)

    # compute the mass-coefficient map
    coeffs = _get_masscoeff_grid(P, T[time == t])
    
    # compute the master transmission spectrum



def _setup_exotransmit(t, latitude, longitude, outfile='default.dat'):
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
    time, lon, lat, h, Ps, P, T, X_H2O = _get_GCM_data(
        'plasim_samples/TIDAL1.0.001.nc')
    Ntime, Nh, Nlat, Nlon = T.shape

    # create exotransmit EOS file for this column
    tindex    = time == t
    lattindex = lat == latitude
    lonindex  = lon == longitude
    exo.setup_Earthlike_EOS(P, X_H2O[tindex,:,latindex,lonindex])

    # create exotransmit TP profile for this column
    exo.setup_TP_file(P, T[tindex,:,latindex,longindex])

    # create exotransmit input files
    exo.create_input_file(g, rp, Rs, outfile)
    exo.create_chem_file()

    
def _run_exotransmit(clean=False):
    '''
    Compute the transmission spectrum using ExoTransmit.
    '''
    if clean:
        os.system('make clean')
        os.system('make')
    os.system('./Exo_Transmit')            


def _get_GCM_data(fname):
    '''
    Read-in temporal, spatial, and thermodynamic variables from the netcdf 
    output of the GCM.
    '''
    data = nc.Dataset(fname, 'r').variables
    time, lon, lat, lev = data['time'][:], data['lon'][:], data['lat'][:], \
                          data['lev'][:]
    T, Ps, X_H2O = data['ta'][:], data['ps'][:], data['hus'][:]
    P, Ps = lev*10, Ps*10  # hPa to Pa 
    h = _P2h(P, Ps, T)
    return time, lon, lat, h, Ps, P, T, X_H2O


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


def _get_masscoeff_grid(P, Ps, T, tindex):
    '''
    Compute the grid of unnormalized cell masses assuming an ideal gas and 
    renormalizing to get the mass-weighted coefficients as a function of lat, 
    lon, and depth.
    '''
    Ntime, Nh, Nlat, Nlon = T.shape
    assert P.size == Nh
    P3d = np.stack([np.stack([P for i in range(Nlat)], 0)
                    for i in range(Nlon)], 0).T

    # compute cell masses modolo a constant
    rho3d = P3d / T[tindex]
    depth = _P2h(P, Ps, T)[tindex]
    x = 2*depth*np.tan(dlon/2.)
    y = 2*depth*np.tan(dlat/2.)
    z = np.diff(depth) #??
    V3d = x*y*z
    mass3d = rho3d * V3d


def _compute_cell_volume(cell_bnds1, cell_bnds2, rp, H, N=1e5):
    '''
    Compute the approximate cell volume by drawing spherical coordinates 
    within the planet's atmosphere (of known volume) and computing the fraction
    of points that lie within the cell bounded by the specified depths, 
    latitudes, and longitudes. 
    '''    
    # get atmosphere volume
    V = 4*np.pi/3 * ((rp+H)**3 - rp**3)

    # draw random locations in the atmosphere
    N = int(N)
    hs   = np.random.uniform(rp, rp+H, N)
    lats = np.random.uniform(-90, 90, N)
    lons = np.random.uniform(0, 360, N)

    # read cell boundaries
    h1, lat1, lon1 = cell_bnds1
    h2, lat2, lon2 = cell_bnds2
    hcell, latcell, loncell = np.append(h1,h2), np.append(lat1,lat2), \
                              np.append(lon1,lon2)

    # compute cell volume
    incell = (hs >= hcell.min()) & (hs <= hcell.max()) & \
             (lats >= latcell.min()) & (lats <= latcell.max()) & \
             (lons >= loncell.min()) & (lons <= loncell.max())
    Vcell = V * incell.sum() / incell.size
    print 'fractional volume = %.6e'%(Vcell / V)
    return Vcell
