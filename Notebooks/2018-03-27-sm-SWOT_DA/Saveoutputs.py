################################
# SSH save function            #
################################ 
import numpy as np 
import netCDF4 as nc
from datetime import datetime
from datetime import timedelta
import os

from params import *

def Save_present_time_outputs(state_vectors_names,n_ens,date0):
    """
    NAME 
        Saveoutputs 

    DESCRIPTION
        Save a netcdf file of an ensemble of SSH fields

        Args:  
            state_vectors_names(string): path and name of the propagated state_vectors file (one-ensemble member) 
            n_ens (integer): number of ensemble member 
            date0 (datetime object): current date
            
        Internal args: 
            name_exp (string): name of the experiment

        Returns:  
            None 

    """
    #date0=datetime(2012,6,14)
    name_exp=name_exp_save
    path=path_save 
    if os.path.isdir(path)==False: 
        cmd1='mkdir '+path
        os.system(cmd1)
    name=path+name_exp+'_y'+str(date0.year)+'m'+str(date0.month).zfill(2)+'d'+str(date0.day).zfill(2)+'h'+str(date0.hour).zfill(2)+'_SSHdegrad.nc'
    cmd2='cp '+state_vectors_names+' '+name  
    os.system(cmd2) 

    
    return 
        

    
"""
    ncout = nc.Dataset(path+name_exp+'_y'+str(date0.year)+'m'+str(date0.month)+'d'+str(date0.day)+'h'+str(date0.hour)+'_SSHdegrad.nc', 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('x', lon2d.shape[0])
    ncout.createDimension('y', lat2d.shape[1])
    ncout.createDimension('member', n_ens)
    nclon = ncout.createVariable('lon', 'f', ('x','y',))
    nclat = ncout.createVariable('lat', 'f', ('x','y',)) 
    ncens = ncout.createVariable('ens', 'f', ('member',)) 
    nchei = ncout.createVariable('SSH', 'f', ('member','x','y',)) 
    nclat[:,:] = lat2d 
    nclon[:,:] = lon2d  
    ncens[:] = range(n_ens) 
    nchei[:,:,:] = field0 
    ncout.close()
"""