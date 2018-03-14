#!/usr/bin/env python
# python

import sys
import numpy as np
from netCDF4 import Dataset
import geo2rot
import geom

def point(ncd_name,lon,lat):
    ##########################
    # Inputs
    # ncd = path to netcdf file
    # lon = geographical longitude
    # lat = geographical latitude
    # Outputs
    # dq = convective heat flux [W m-2]
    #############################
    # constants
    cps = 1000 # J kg-1 K-1 heat capacity of the air
    # read grid variables
    ncd = Dataset(ncd_name,'r')
    rlon = ncd.variables['rlon'][:]
    rlat = ncd.variables['rlat'][:]
    vcoord = ncd.variables['vcoord']
    # indexes from lon and lat
    [lon_r,lat_r] = geo2rot.g2r(lon,lat)
    idx = np.abs(rlon-(lon_r)).argmin()
    idy = np.abs(rlat-(lat_r)).argmin()
    # read variables            
    u = ncd.variables['U'][0,-1,:,:]
    v = ncd.variables['V'][0,-1,:,:]
    t = ncd.variables['T'][0,-1,:,:]
    w = ncd.variables['W'][0,-1,:,:]
    tup = ncd.variables['T'][0,-2,:,:]
    rho = ncd.variables['RHO'][0,-1,:,:]
    rhoup = ncd.variables['RHO'][0,-2,:,:]
    # additional variables
    lat = ncd.variables['lat'][:,:]
    lon = ncd.variables['lon'][:,:]
    # cell size
    dx1 = abs(geom.dist(lon[idy,idx],lat[idy,idx],lon[idy,idx-1],lat[idy,idx-1]))
    dx2 = abs(geom.dist(lon[idy,idx],lat[idy,idx],lon[idy,idx+1],lat[idy,idx+1]))
    dx = (dx1 + dx2)/2
    dy1 = abs(geom.dist(lon[idy,idx],lat[idy,idx],lon[idy-1,idx],lat[idy-1,idx]))
    dy2 = abs(geom.dist(lon[idy,idx],lat[idy,idx],lon[idy+1,idx],lat[idy+1,idx]))
    dy = (dy1 + dy2)/2
    # convert to horizontal
    h = vcoord[-2] - vcoord[-1]
    Ax1 = dx1 * h
    Ax2 = dx2 * h
    Ay1 = dy1 * h
    Ay2 = dy2 * h
    # advective flux
    q_in_x = Ax1 * cps * ((rho[idy,idx]+rho[idy,idx-1])/2) *   \
             u[idy,idx] * ((t[idy,idx]+t[idy,idx-1])/2)
    q_out_x = Ax2 * cps * ((rho[idy,idx]+rho[idy,idx+1])/2) *  \
             u[idy,idx+1] * ((t[idy,idx]+t[idy,idx+1])/2)
    q_in_y = Ay1 * cps * ((rho[idy,idx]+rho[idy-1,idx])/2) *   \
             v[idy,idx] * ((t[idy,idx]+t[idy-1,idx])/2)
    q_out_y = Ay2 * cps * ((rho[idy,idx]+rho[idy+1,idx])/2) *  \
             v[idy+1,idx] * ((t[idy,idx]+t[idy+1,idx])/2)
    Dq_tot = q_in_x - q_out_x + q_in_y - q_out_y # in W
    Dq_vol = Dq_tot / (dx * dy * h) # volumetric (W/m3)
    Dq_flux = Dq_vol * h  # W / m2 as other turbulent (unresolved) fluxes
    # vertical component not considered
    #q_out_z = cps * ((rhoup[idy,idx]+rho[idy,idx])/2) *         \
    #         w[idy+1,idx] * ((tup[idy,idx]+t[idy,idx])/2)
    return Dq_flux

def spatial(nc_name):
    ##########################
    # Inputs
    # ncd = path to netcdf file
    # lon = geographical longitude
    # lat = geographical latitude
    # Outputs
    # dq = convective heat flux [W m-2]
    # Hypotesis
    # 1. t can be interpolated linearly
    # 2. vertical advective fluxes neglected (need also w)
    # constant (to be edited event.)
    #############################
    # constants
    cps = 718 # J kg-1 K-1 at constant volume
    ds = 280 # m
    # read grid variables
    nc = Dataset(nc_name,'r')
    vcoord = nc.variables['vcoord']
    # read variables        
    u = nc.variables['U'][0,-1,:,:]
    v = nc.variables['V'][0,-1,:,:]
    #w = nc.variables['W'][0,-1,:,:]
    t = nc.variables['T'][0,-1,:,:]
    tup = nc.variables['T'][0,-1,:,:]
    rho = nc.variables['RHO'][0,-1,:,:]
    rhoup = nc.variables['RHO'][0,-2,:,:]
    # define output variable
    Dq = np.zeros((np.shape(u)[0],np.shape(u)[1]))
    # convert to horizontal
    h = vcoord[-2] - vcoord[-1]
    # advective flux
    for j in range(0, np.shape(u)[0]-1):
        for k in range(0, np.shape(u)[1]-1):
            # advection
            q_in_x = h/ds * cps * (rho[j,k]+rho[j,k-1])**0.5 *   \
                     u[j,k] * (t[j,k]+t[j,k-1])**0.5
            q_out_x = h/ds * cps * (rho[j,k]+rho[j,k+1])**0.5 *  \
                     u[j,k+1] * (t[j,k]+t[j,k+1])**0.5
            q_in_y = h/ds * cps * (rho[j,k]+rho[j-1,k])**0.5 *   \
                     v[j,k] * (t[j,k]+t[j-1,k])**0.5
            q_out_y = h/ds * cps * (rho[j,k]+rho[j+1,k])**0.5 *  \
                     v[j+1,k] * (t[j,k]+t[j+1,k])**0.5
            #q_out_z = cps * (rhoup[j,k]+rho[j,k])**0.5 *          \
            #         w[j,k] * t[j,k]
            Dq[j,k] = q_in_x - q_out_x + q_in_y - q_out_y # - q_out_z


    return Dq
=======
    # advective flux
    q_in_x = rho*cps*u[idy,idx]*t[idy,idx-1]
    q_out_x = rho*cps*u[idy,idx+1]*t[idy,idx+1]
    q_in_y = rho*cps*v[idy,idx]*t[idy-1,idx]
    q_out_y = rho*cps*v[idy+1,idx]*t[idy+1,idx]
    Dq = q_in_x - q_out_x + q_in_y - q_out_y
    # convert to horizontal
    h = vcoord[-2] - vcoord[-1]
    Av = h * ds
    Ah = ds**2
    dh = Av/Ah
    Dq_hor = Dq * dh
    return Dq_hor

def spatial(ncd_name):
    ##########################
    # Inputs
    # ncd = path to netcdf file
    # Outputs
    # dq = convective heat flux [W m-2]
    # Hypotesis
    # - vertical advective fluxes neglected (need also w)
    # constant (to be edited event.)
    #############################
    # constants
    rho = 1 # kg m-3
    cps = 1000 # J kg-1 K-1
    ds = 280 # m
    # read grid variables
    ncd = Dataset(ncd_name,'r')
    vcoord = ncd.variables['vcoord']
    # read variables            
    u = ncd.variables['U'][0,-1,:,:]
    v = ncd.variables['V'][0,-1,:,:]
    t = ncd.variables['T'][0,-1,:,:]
    # define output variable
    Dq = np.zeros((np.shape(u)[0]-1,np.shape(u)[1]-1))
    # advective flux
    for j in range(0, np.shape(u)[0]-1):
        for k in range(0, np.shape(u)[1]-1):
            q_in_x = rho*cps*u[j,k]*t[j,k-1]
            q_out_x = rho*cps*u[j,k+1]*t[j,k+1]
            q_in_y = rho*cps*v[j,k]*t[j-1,k]
            q_out_y = rho*cps*v[j+1,k]*t[j+1,k]
            Dq[j,k] = q_in_x - q_out_x + q_in_y - q_out_y


    # convert to horizontal
    h = vcoord[-2] - vcoord[-1]
    Av = h * ds
    Ah = ds**2
    dh = Av/Ah
    Dq_hor = Dq * dh
    return Dq_hor

