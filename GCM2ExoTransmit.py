from imports import *
import rvs

global mp, kb, rp, g
mp, kb = 1.6726219e-27, 1.38064852e-23
rp, g = rvs.Rearth2m(1), 9.8


def main():
    '''
    Compute the transmission spectrum from a GCM output file.
    '''
    # get data
    time, lon, lat, h, Ps, P, T, X_H2O = _get_GCM_data(
        'plasim_samples/TIDAL1.0.001.nc')
    Ntime, NP, Nlat, Nlon = T.shape

    # compute mass weighted spectrum for each cell along each light ray
    ##nrays = Nlon * NP
    x_ray = np.linspace(+rp+h.max(), -rp-h.max(), Nlat)
    R_rays, Theta_rays = np.linspace(rp, rp+h.max(), NP), \
                         np.arange(0, 2*np.pi, 2*np.pi/Nlon)
    for i in range(NP):
        for j in range(Nlon):

            y_ray, z_ray = _sphere2cart2D(R_rays[i], Theta_rays[j])

            h_ray, lon_ray, lat_ray = _cart2geo3D(x_ray, y_ray, z_ray)

            _write_EOS_file(h, lon, lat, Ps, P, T, X_H2O,
                            h_ray, lon_ray, lat_ray)
            
            
            
            

def _get_GCM_data(fname):
    '''
    Read-in temporal, spatial, and thermodynamic variables from the netcdf 
    output of the GCM.
    '''
    data = nc.Dataset(fname, 'r').variables
    time, lon, lat, lev = data['time'][:], data['lon'][:], data['lat'][:], \
                          data['lev'][:]
    T, Ps, X_H2O = data['ta'][:], data['ps'][:], data['hus'][:]
    P, Ps = lev*10, Ps*10
    h = _P2h(P, Ps, T)
    return time, lon, lat, h, Ps, P, T, X_H2O


def _P2h(P, Ps, T, mu=28.97):
    '''
    Convert the atmospheric pressure to a depth.

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
    return H * np.log(Ps4 / P4)


def _sphere2cart2D(r, theta):
    '''
    Convert radial/azimuthal coordinates to cartesian x,y.
    '''
    x, y = r*np.cos(theta), r*np.sin(theta)
    return x, y 


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


def _write_EOS_file(h, lon, lat, Ps, P, T, X_H2O,
                    h_ray, lon_ray, lat_ray):
    '''
    Write the ExoTransmit EOS file based on the GCM output and the ray path.
    '''
    # get GCM output along the ray
    _get_GCM_variables_along_ray(h, lon, lat, Ps, P, T, X_H2O,
                                 h_ray, lon_ray, lat_ray)

    # write EOS file
    

def _get_GCM_variables_along_ray(h, lon, lat, Ps, P, T, X_H2O,
                                 h_ray, lon_ray, lat_ray):
    '''
    Interpolate along the ray path to get GCM varaiables of interest. 
    '''
    
    
