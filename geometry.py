from imports import *
import netCDF4 as nc


global r0, heff, G, Mp, Rp 
r0, heff, G, Mp, Rp = 1., 1e-2, 6.67408e-11, 5.972e24, 6.371e6 

def _lng2phi(lng):
    return np.deg2rad(lng)

def _lat2theta(lat):
    return np.deg2rad(90.-lat)

def _sphere2cart(r, phi, theta):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

'''
def example_plasim_output_TMP():
    Nlng, Nlat, Nr, Ntime = 65, 33, 11, 1
    lng_grid = np.arange(0, 360, 360./Nlng)
    lat_grid = np.linspace(90, -90, Nlat)
    r_grid   = np.linspace(r0, r0*(1+heff), Nr)
    t_grid   = np.arange(Ntime)
    P_grid4d   = abs(np.random.normal(0, 1, size=(Ntime, Nr, Nlat, Nlng)))
    T_grid4d   = abs(np.random.normal(0, 1, size=(Ntime, Nr, Nlat, Nlng)))
    H2O_grid4d = abs(np.random.normal(0, 1, size=(Ntime, Nr, Nlat, Nlng)))
    return P_grid4d, T_grid4d, H2O_grid4d, t_grid, r_grid, lat_grid, lng_grid
'''

def get_plasim_data():
    data = nc.Dataset("plasim_samples/TIDAL1.0.001.nc","r")
    t_grid = data.variables['time'][:]
    T_grid4d = data.variables['ta'][:]
    Ntime, Nz, Nlat, Nlng = T.shape
    lng_grid = np.arange(0, 360, 360/(Nlng+1.))
    lat_grid = np.linspace(90, -90, Nlat+1)
    P_grid = data.variables['lev'][:]
    
    Pbnds_grid4d = _get_boundary_pressures(P_grid, data.variables['ps'][:])

    T_grid4d = data.variables['ta'][:]
    H2O_grid4d = abs(np.random.normal(0, 1, size=(Ntime, Nr, Nlat, Nlng)))#TEMP
    return P, Psurf, T, H2O_grid4d, t_grid, lng_grid, lat_grid, r_grid


def _get_boundary_pressures(P_grid, Psurf):
    Nz = P_grid.size
    Ntime, Nlat, Nlng = Psurf.shape
    P_gridtmp = np.append(0, P_grid)
    P_grid4d = np.stack([np.stack([np.stack([P_gridtmp for i in range(Nlng+1)],
                                            1) for i in range(Nlat+1)], 1)
                         for i in range(Ntime)], 0)
    Pbnds_grid4d = np.zeros((Ntime, Nz+1, Nlat+1, Nlng+1))
    for i in range(Nz+1):
        if i < Nz:
            Pbnds_grid4d[:,i] = np.mean([P_grid4d[:,i], P_grid4d[:,i+1]],
                                          axis=0)
        else:
            Pbnds_grid4d[:,i] = Psurf
            
    return Pbnds_grid4d


def get_ray_H20():
    T_grid4d, H2O_grid4d, Psurf_grid3d, t_grid, P_grid, lat_grid, lng_grid = \
                                                        get_plasim_output()
    
    Ntime, Nz, Nlat, Nlng = T_grid4d.shape
    assert lng_grid.size == Nlng
    assert lat_grid.size == Nlat
    assert P_grid.size == Nz
    assert t_grid.size == Ntime

    # get spherical coordinates for all grid cells (3 spatial dimensions)
    phi_grid4d = np.stack([np.stack([np.stack([_lng2phi(lng_grid)
                                               for i in range(Nlat)], 0)
                                     for i in range(Nz)], 0)
                           for i in range(Ntime)], 0)
    theta_grid4d = np.stack([np.stack([np.stack([_lat2theta(lat_grid)
                                                 for i in range(Nlng)], 0).T
                                       for i in range(Nz)], 0)
                             for i in range(Ntime)], 0)
    P_grid4d = np.stack([np.stack([np.stack([P_grid
                                             for i in range(Nlat)], 0)
                                   for i in range(Nlng)], 0).T
                         for i in range(Ntime)], 0)
    t_grid4d = np.stack([np.stack([np.stack([t_grid
                                             for i in range(Nlng)], 1)
                                   for i in range(Nlat)], 1)
                         for i in range(Nz)], 1)
    assert P_grid4d.shape == phi_grid4d.shape
    assert phi_grid4d.shape == theta_grid4d.shape
    assert theta_grid4d.shape == P_grid4d.shape
    assert P_grid4d.shape == t_grid4d.shape

    # get cartesian coordinates for all grid cells
    x_grid4d, y_grid4d, z_grid4d = _sphere2cart(P_grid4d, phi_grid4d,
                                                theta_grid4d)

    # get mass-weights for each cell
    mass_grid4d = _get_cell_mass(P_grid4d, Psurf_grid3d)

    # send rays through the atmosphere along x at fixed y and z
    N_azi = 2*(Nlat-1)
    Nrays = N_azi*(Nz-1)
    projected_planet_angles = np.arange(0, 2*np.pi, 2*np.pi/(N_azi+1))
    P_ray, T_ray, H2O_ray, rayind = np.zeros(Nrays), np.zeros(Nrays), \
                                    np.zeros(Nrays), 0
    
    for i in range(Nz-1):
        for j in range(N_azi):

            # Get cell y and z position
            y1, z1 = P_grid[i]*np.cos(projected_planet_angles[j]), \
                     P_grid[i]*np.sin(projected_planet_angles[j])
            y2, z2 = P_grid[i+1]*np.cos(projected_planet_angles[j+1]), \
                     P_grid[i+1]*np.sin(projected_planet_angles[j+1])
            ys, zs = [y1,y2], [z1,z2]
                        
            # get cells along the ray path
            g = (y_grid4d >= np.min(ys)) & (y_grid4d <= np.max(ys)) & \
                (z_grid4d >= np.min(zs)) & (z_grid4d <= np.max(zs))

            if g.sum() > 0:
                # get mass-weighted P,T,H2O along the ray
                mass_coeffs = mass_grid4d / mass_grid4d[g].sum()
                P_ray[rayind] = np.sum(mass_coeffs[g] * P_grid4d[g])
                T_ray[rayind] = np.sum(mass_coeffs[g] * T_grid4d[g])
                H2O_ray[rayind] = np.sum(mass_coeffs[g] * H2O_grid4d[g])
            else:
                T_ray[rayind] = -99

            rayind += 1

    g = (T_ray != -99)
    return P_ray[g], T_ray[g], H2O_ray[g]


def get_plasim_output():
    data = nc.Dataset("plasim_samples/TIDAL1.0.001.nc","r")
    t_grid = data.variables['time'][:]
    P_grid = data.variables['lev'][:]
    lat_grid = data.variables['lat'][:]
    lng_grid = data.variables['lon'][:]
    T_grid4d = data.variables['ta'][:]
    H2O_grid4d = np.random.normal(0,1,T_grid4d.shape)
    Psurf_grid3d = data.variables['ps'][:]
    return T_grid4d, H2O_grid4d, Psurf_grid3d, \
        t_grid, P_grid, lat_grid, lng_grid
    

def _get_cell_mass(P_grid4d, Psurf_grid3d):
    P_grid4d_wPsurf = np.insert(P_grid4d[:,::-1],0,Psurf_grid3d,axis=1)[:,::-1]
    dP = np.diff(P_grid4d_wPsurf, axis=1)
    g = G*Mp/Rp**2
    return dP/g
    

#dP/g to get vertically in mass in the cell
#P, T is given at the middle of the cell
