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


def get_ray_H20(P_grid4d, T_grid4d, H2O_grid4d, t_grid, P_grid,
                lat_grid, lng_grid):

    Ntime, Nr, Nlat, Nlng = P_grid4d.shape
    assert lng_grid.size == Nlng
    assert lat_grid.size == Nlat
    assert r_grid.size == Nr
    assert t_grid.size == Ntime

    # get spherical coordinates for all grid cells (3 spatial dimensions)
    phi_grid4d = np.stack([np.stack([np.stack([_lng2phi(lng_grid)
                                               for i in range(Nlat)], 0)
                                     for i in range(Nr)], 0)
                           for i in range(Ntime)], 0)
    theta_grid4d = np.stack([np.stack([np.stack([_lat2theta(lat_grid)
                                                 for i in range(Nlng)], 0).T
                                       for i in range(Nr)], 0)
                             for i in range(Ntime)], 0)
    P_grid4d = np.stack([np.stack([np.stack([P_grid
                                             for i in range(Nlat)], 0)
                                   for i in range(Nlng)], 0).T
                         for i in range(Ntime)], 0)
    t_grid4d = np.stack([np.stack([np.stack([t_grid
                                             for i in range(Nlng)], 1)
                                   for i in range(Nlat)], 1)
                         for i in range(Nr)], 1)
    assert P_grid4d.shape == phi_grid4d.shape
    assert phi_grid4d.shape == theta_grid4d.shape
    assert theta_grid4d.shape == r_grid4d.shape
    assert r_grid4d.shape == t_grid4d.shape

    # get cartesian coordinates for all grid cells
    x_grid4d, y_grid4d, z_grid4d = _sphere2cart(r_grid4d, phi_grid4d,
                                                theta_grid4d)

    # send rays through the atmosphere along x at fixed y and z
    N_azi = 2*(Nlat-1)
    Nrays = N_azi*(Nr-1)
    projected_planet_angles = np.arange(0, 2*np.pi, 2*np.pi/(N_azi+1))
    P_ray, T_ray, H2O_ray, ind = np.zeros(Nrays), np.zeros(Nrays), \
                                 np.zeros(Nrays), 0
    mass_grid4d = np.ones(P_grid4d.shape) # TEMP
    
    for i in range(Nr-1):
        for j in range(N_azi):

            # Get cell y and z position
            y1, z1 = r_grid[i]*np.cos(projected_planet_angles[j]), \
                     r_grid[i]*np.sin(projected_planet_angles[j])
            y2, z2 = r_grid[i+1]*np.cos(projected_planet_angles[j+1]), \
                     r_grid[i+1]*np.sin(projected_planet_angles[j+1])
            ys, zs = [y1,y2], [z1,z2]
            
            # get cells along the ray path
            g = (y_grid4d >= np.min(ys)) & (y_grid4d <= np.max(ys)) & \
                (z_grid4d >= np.min(zs)) & (z_grid4d <= np.max(zs))

            # get normalized mass-weights for each cell along the ray path
            # is mass or density information even given??
            #mass_grid4d
            #mass_grid4d /= mass_grid[g].sum()

            # get mass-weighted P,T,H2O along the ray
            mass_grid4d /= g[g].size  # TEMP
            P_ray[ind] = np.sum(mass_grid4d[g] * P_grid4d[g])
            T_ray[ind] = np.sum(mass_grid4d[g] * T_grid4d[g])
            #H2O_ray[ind] = np.sum(mass_grid4d[g] * H2O_grid4d[g])
            ind += 1
            
    return P_ray, T_ray, H2O_ray


def _get_cell_mass(P_grid4d, Psurf_grid4d):
    P_grid = P_grid4d[0,:,0,0]
    dP = np.diff(P_grid)
    g = G*mp/rp**2
    return dP/g
    

#dP/g to get vertically in mass in the cell
#P, T is given at the middle of the cell
