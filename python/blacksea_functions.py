def calc_uv_track(lon,lat,sla):

    import numpy as np
    import numexpr as ne
    '''
        Returns the geostrophic velocity components
        Gradient approx. by finite differences
    '''
    g = 9.81                # gravity
    omega = 2*np.pi/(86400) # earth rotation

    # Define great-circle distance for vectors of
    # longitude and latitude
    
    def haversine_dist(lon1, lat1, lon2, lat2):
        '''
        Haversine formula to calculate distance between two lon/lat points
        Uses mean earth radius in metres (from ROMS scalars.h) = 6371315.0
        Input:  lon1, lat1, lon2, lat2
        Return: distance (m)
        '''
        lon1, lat1, lon2, lat2 = lon1.copy(), lat1.copy(), lon2.copy(), lat2.copy()
        dlat = np.deg2rad(lat2 - lat1)
        dlon = np.deg2rad(lon2 - lon1)
        np.deg2rad(lat1, out=lat1)
        np.deg2rad(lat2, out=lat2)
        a = ne.evaluate('sin(0.5 * dlon) * sin(0.5 * dlon)')
        a = ne.evaluate('a * cos(lat1) * cos(lat2)')
        a = ne.evaluate('a + (sin(0.5 * dlat) * sin(0.5 * dlat))')
        c = ne.evaluate('2 * arctan2(sqrt(a), sqrt(1 - a))')
        return ne.evaluate('6371315.0 * c') # Return the distance in meters

    # Compute distance
    distSLA=haversine_dist(lon[0:-1],lat[0:-1],lon[1:],lat[1:])
    # Convert degrees to radian
    lon, lat = lon.copy(), lat.copy()
    np.deg2rad(lat, out=lat)
    np.deg2rad(lon, out=lon)
    # Compute f and geostrophic velocity
    fcoriolis = 2*omega*np.sin(0.5*(lat[0:-1]+lat[1:]))
    geostr_vel = g*np.diff(sla)/(fcoriolis*distSLA)
    # Compute angle (used for plotting)
    alpha = np.arctan2(lat[-1]-lat[0],lon[-1]-lon[0])
    # Compute velocity components
    ugeostr = -geostr_vel*np.sin(alpha)
    vgeostr = geostr_vel*np.cos(alpha)

    #print 'alpha = '+str(np.rad2deg(abs(alpha)))+' degrees'
    return geostr_vel,ugeostr,vgeostr,distSLA,alpha


def uvpmask(mask):
        '''
        Get mask at u, v, psi points
        '''
        Mp, Lp = mask.shape
        M = Mp - 1
        L = Lp - 1
        u_mask = mask[:,:L] * mask[:,1:Lp]
        v_mask = mask[:M] * mask[1:Mp]
        return u_mask, v_mask

def half_interp(hone, htwo):
    import numpy as np
    import numexpr as ne
    
    '''
    Speed up frequent operations
    '''
    return ne.evaluate('0.5 * (hone + htwo)')

def distLonLat(lon1, lat1, lon2, lat2):
    import numpy as np
    import numexpr as ne
    '''
    Haversine formula to calculate distance between one point and another
    Uses mean earth radius in metres (from scalars.h) = 6371315.0
    '''
    lon1, lat1, lon2, lat2 = lon1.copy(), lat1.copy(), lon2.copy(), lat2.copy()
    dlat = np.deg2rad(lat2 - lat1)
    dlon = np.deg2rad(lon2 - lon1)
    np.deg2rad(lat1, out=lat1)
    np.deg2rad(lat2, out=lat2)
    a = ne.evaluate('sin(0.5 * dlon) * sin(0.5 * dlon)')
    a = ne.evaluate('a * cos(lat1) * cos(lat2)')
    a = ne.evaluate('a + (sin(0.5 * dlat) * sin(0.5 * dlat))')
    c = ne.evaluate('2 * arctan2(sqrt(a), sqrt(1 - a))')
    # Return the distance
    return ne.evaluate('6371315.0 * c')

def u2rho_2d(uu_in):
    import numpy as np
    import numexpr as ne
    '''
    Convert a 2D field at u points to a field at rho points
    Checked against Jeroen's u2rho.m
    '''
    def uu2ur(uu_in, Mp, Lp):
        L = Lp - 1
        Lm = L  - 1
        u_out = np.zeros((Mp, Lp))
        u_out[:, 1:L] = half_interp(uu_in[:, 0:Lm], uu_in[:, 1:L])
        u_out[:, 0] = u_out[:, 1]
        u_out[:, L] = u_out[:, Lm]
        return (np.squeeze(u_out))
    Mshp, Lshp = uu_in.shape
    return uu2ur(uu_in, Mshp, Lshp + 1)


def v2rho_2d(vv_in):
    import numpy as np
    import numexpr as ne
    # Convert a 2D field at v points to a field at rho points
    def vv2vr(vv_in, Mp, Lp):
        M = Mp - 1
        Mm = M  - 1
        v_out = np.zeros((Mp, Lp))
        v_out[1:M] = half_interp(vv_in[:Mm], vv_in[1:M])
        v_out[0] = v_out[1]
        v_out[M] = v_out[Mm]
        return (np.squeeze(v_out))
    Mshp, Lshp = vv_in.shape
    return vv2vr(vv_in, Mshp + 1, Lshp)

def getSurfGeostrVel(f, zeta, pm, pn, umask, vmask):
    import numpy as np
    import numexpr as ne
    '''
    Returns u and v geostrophic velocity at
    surface from ROMS variables f, zeta, pm, pn...
    Note: output at rho points
    Adapted from IRD surf_geostr_vel.m function
    by Evan Mason
    Changed to gv2, 14 May 07...
    '''
    def gv2(f, zeta, pm, pn, umask, vmask): # Pierrick's version
        gof = ne.evaluate('9.81 / f')
        ugv = -gof * v2rho_2d(vmask * (zeta[1:] - zeta[:-1]) \
                   * half_interp(pn[1:], pn[:-1]))
        vgv =  gof * u2rho_2d(umask * (zeta[:, 1:] - zeta[:, :-1]) \
                   * half_interp(pm[:, 1:], pm[:, :-1]))
        return ugv, vgv
    ugv, vgv = gv2(f, zeta, pm, pn, umask, vmask)
    return ugv, vgv

def create_netcdf_eke(netcdffile,nlat,nlon,valex,time2write,\
                      lon2write,lat2write,u2write,\
                      v2write,EKE2write):

    import netCDF4 as netcdf
    from netCDF4 import Dataset
    import time as ttime

    ''' Create a NetCDF file that will contain
    the geostrophic velocity and the EKE
    '''
    rootgrp = Dataset(netcdffile,'w', format='NETCDF4')
    #print rootgrp.file_format

    rootgrp.history = 'Created ' + ttime.ctime(ttime.time())
    rootgrp.author = 'ctroupin'
    rootgrp.institute = 'SOCIB'
     
    # Create dimensions
    time = rootgrp.createDimension('time', None)     # to have record dimension
    lat = rootgrp.createDimension('lat', nlat)
    lon = rootgrp.createDimension('lon', nlon)

    # Create variables
    time = rootgrp.createVariable('time','f4',('time',))
    latitudes = rootgrp.createVariable('latitude','f4',('lat',))
    longitudes = rootgrp.createVariable('longitude','f4',('lon',))
    u = rootgrp.createVariable('u','f4',('time','lat','lon',))
    v = rootgrp.createVariable('v','f4',('time','lat','lon',))
    EKE = rootgrp.createVariable('EKE','f4',('time','lat','lon',))

    # Put attributes
    time.long_name = 'Time' ;
    time.standard_name = 'Time' ;
    time.units = 'days since 1950-01-01 00:00:00' ;

    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    u.units = 'meters per second'
    u.long_name = 'Zonal component of geostrophical velocity'
    u.missing_value = valex
    v.units = 'meters per second'
    v.long_name = 'Meridional component of geostrophical velocity'
    v.missing_value = valex
    EKE.units = 'cm^2/s^2'
    EKE.long_name = 'Eddy Kinetic Energy'
    EKE.missing_value = valex

    # write variables
    time[:] = time2write
    latitudes[:]=lat2write
    longitudes[:] = lon2write
    u[0,:,:] = u2write
    v[0,:,:] = v2write
    EKE[0,:,:] = EKE2write
    
    # close NetCDF
    rootgrp.close()
