##########################################
# State vectors initialization functions #
##########################################
import netCDF4 as nc
import numpy as np
import os


def Initialization(function,n_ens,*args):
    """
    NAME 
        Initialization 

    DESCRIPTION
        Global initialization function

        Args:
            function (function): initialization function 
            n_ensemble (integer): number of ensemble members

        Returns:
            function(n_ens,*args): output of *function*   
    """            
            
    return function(n_ens,*args)


def NATL60state(n_ens=1):
    """
    NAME 
        NATL60state

    DESCRIPTION 
        Copy the SSH initial field in *file_name_init_SSH_field* at *path_init_SSH_field* as the current state_vectors0 files

        Args: 
            n_ens (integer): number of ensemble members

        Internal args:
            file_name_init_SSH_field (string): name of the initialization file needed to create the inital SSH field(s) 
            path_init_SSH_field (string): path of the initialization file needed to create the inital SSH field(s)

        Returns: 
            state_vectors0_names(array of strings): arrays of path and names of the state_vectors0 files 

    """

    # Initial SSH field file name
    file_name_init_SSH_field='NATL60OSMO-CJM165_y2012m10d01h00.1h_SSHdegrad.nc'
    # Initial SSH field path
    path_init_SSH_field='/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/OSMOSIS/QGSWOSMO-IniNATL60_1h/'+file_name_init_SSH_field


    if n_ens>1:
        print('Warning: NATL60state only works for one-member-ensemble') 
    fid = nc.Dataset(path_init_SSH_field)
    lon=np.array(fid.variables["nav_lon"][:])
    lat=np.array(fid.variables["nav_lat"][:]) 
    multiplefields=np.array(fid.variables["degraded_sossheig"][:,:]) 
      
    state_vectors0_names='TMP_DA/state_vectors0.nc'
    ncout = nc.Dataset(state_vectors0_names, 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('x', lon.shape[0])
    ncout.createDimension('y', lat.shape[1])
    ncout.createDimension('member', n_ens)
    ncens = ncout.createVariable('ens', 'd', ('member',))  
    nclon = ncout.createVariable('nav_lon', 'f', ('x','y',))
    nclat = ncout.createVariable('nav_lat', 'f', ('x','y',)) 
    nclat[:,:] = lat 
    nclon[:,:] = lon 
    nchei = ncout.createVariable('degraded_sossheig', 'f', ('member','x','y',))
    ncens[:] = range(n_ens) 
    for i_ens in range(n_ens):  
        nchei[i_ens,:,:] = multiplefields[0,:,:] 
    ncout.close()
         
    
    return state_vectors0_names


def Reinitialization(statevectors0_names,statevectors_names):
    """
    NAME 
        Reinitialization 

    DESCRIPTION
        Re-initialize the state vectors to loop 

        Args:
            statevectors0_names
            statevectors_names
    
    """            
    cmd1='cp '+statevectors_names+' '+statevectors0_names
    os.system(cmd1)
            
    return 