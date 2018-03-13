#!/usr/bin/env python
# python

import numpy as np
from netCDF4 import Dataset
import geo2rot

def point(ncd_name,lon,lat):
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
    rho = 1 # kg m-3
    cps = 1000 # J kg-1 K-1
    ds = 280 # m
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
    # advective flux
    q_in_x = rho*cps*u[idy,idx]*(t[idy,idx]+t[idy,idx-1])**0.5
    q_out_x = rho*cps*u[idy,idx+1]*(t[idy,idx]+t[idy,idx+1])**0.5
    q_in_y = rho*cps*v[idy,idx]*(t[idy,idx]+t[idy-1,idx])**0.5
    q_out_y = rho*cps*v[idy+1,idx]*(t[idy,idx]+t[idy+1,idx])**0.5
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
            q_in_x = rho*cps*u[j,k]*(t[j,k]+t[j,k-1])**0.5
            q_out_x = rho*cps*u[j,k+1]*(t[j,k]+t[j,k+1])**0.5
            q_in_y = rho*cps*v[j,k]*(t[j,k]+t[j-1,k])**0.5
            q_out_y = rho*cps*v[j+1,k]*(t[j,k]+t[j+1,k])**0.5
            Dq[j,k] = q_in_x - q_out_x + q_in_y - q_out_y


    # convert to horizontal
    h = vcoord[-2] - vcoord[-1]
    Av = h * ds
    Ah = ds**2
    dh = Av/Ah
    Dq_hor = Dq * dh
    return Dq_hor

