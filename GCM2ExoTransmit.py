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
    

def _create_interpolators(time, h, lon, lat, P, T, X_H2O):
    '''
    Initialize interpolation functions to get P, T, and X_H2O.
    ##Interpolate along the ray path to get GCM varaiables of interest. 
    '''
    # form grid of spatial and temporal points
    Ntime, NP, Nlat, Nlon = T.shape
    Npnts = np.product(T.shape)
    points = np.zeros((Npnts, 4))
    points[:,0] = np.ascontiguousarray(list(time) * NP * Nlat * Nlon)    
    points[:,1] = h.reshape(Npnts)
    points[:,2] = np.ascontiguousarray(list(lat) * Ntime * NP * Nlon)
    points[:,3] = np.ascontiguousarray(list(lon) * Ntime * NP * Nlat)
    
    # interpolate to get pressures
    P4 = np.repeat(np.repeat(np.repeat(P[np.newaxis,:], Ntime,
                                       axis=0)[:,:,np.newaxis], Nlat,
                             axis=2)[:,:,:,np.newaxis], Nlon, axis=3)
    assert P4.shape == T.shape
    Pint = lint(points, P4.reshape(Npnts))
    
    # interpolate to get temperatures
    Tint = lint(points, T.reshape(Npnts))

    # interpolate to get specific humidity
    Xint = lint(points, X_H2O.reshape(Npnts))

    return Pint, Tint, Xint
