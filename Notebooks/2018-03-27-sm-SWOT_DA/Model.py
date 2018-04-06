################################
#       Model functions        #
################################ 
import sys,os,shutil
sys.path.insert(0,'/home/jovyan/work/Boost-SWOT/2018/Notebooks/2018-03-27-sm-SWOT_DA/') 
import numpy as np 
import netCDF4 as nc

# Model specific libraries 
from modelQG_libraries import * 
 


def EnsembleModel(function,state_vectors0_names,n_ens,time0,state_vectors_names,*args):
    """
    NAME 
        EnsembleModel 

    DESCRIPTION
        Global model function

        Args: 
            function (function): model function
            state_vectors0_names (list of strings): ensemble of path and names of the current state_vectors0 files (all members) 
            n_ens (integer): number of ensemble member 
            time0 (integer [in seconds]): time of propagation

        Returns: 
            state_vectors_names (string): ensemble of path and names of the propagated state_vectors files (all members)

    """    
    ncout = nc.Dataset(state_vectors_names, 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('member', n_ens)
    ncens = ncout.createVariable('ens', 'd', ('member',))
    ncens[:] = range(n_ens) 
        
    for i_ens in range(n_ens):  
        function(state_vectors0_names,time0,i_ens,ncout,*args) 
         
    ncout.close()  
        
    return state_vectors_names
        


def QG_SW(state_vectors0_names,time0,i_ens,ncout): 
    """
    NAME 
        QG_SW

    DESCRIPTION 
        Function calling the quasi-geostrophic shallow water model (C. Ubelmann)

        Args: 
            state_vectors0_name(string): path and name of the current state_vectors0 file (one-ensemble member)   
            time0 (float [in seconds]): time of propagation
            i_ens (integer): ensemble member index

        Internal args:
            c (2D array of dimensions [n_lon,n_lat]): Rossby first baroclinic phase speed  
            dtmodel (integer [in seconds]): model time step 
            deltat (integer [in seconds]): output time step 

        Returns: 
            None

    """      
    
    fid_deg = nc.Dataset(state_vectors0_names)
    lon2d=np.array(fid_deg.variables["nav_lon"][:,:]) 
    lat2d=np.array(fid_deg.variables["nav_lat"][:,:])   
    current_SSH=np.array(fid_deg.variables["degraded_sossheig"][i_ens,:,:]) 
    
    # Model parameters: Internal args
    c=current_SSH*0. + 2.5 # in m/s 
    dtmodel=600 # propagator time step
    deltat=time0
     
    
    # Model propagation
    propagated_SSH, trash = qgsw.qgsw(Hi=current_SSH, c=c, lon=lon2d, lat=lat2d, tint=time0, dtout=deltat, dt=dtmodel,rappel=None,snu=0.) 
     
     
    
    if i_ens==0:
        ncout.createDimension('x', lon2d.shape[0])
        ncout.createDimension('y', lat2d.shape[1])   
        nclon = ncout.createVariable('nav_lon', 'f', ('x','y',))
        nclat = ncout.createVariable('nav_lat', 'f', ('x','y',))  
        nchei = ncout.createVariable('degraded_sossheig', 'f', ('member','x','y',)) 
        nclat[:,:] = lat2d 
        nclon[:,:] = lon2d   
        nchei[0,:,:] = propagated_SSH[-1,:,:] 
    else:    
        ncout.variables["degraded_sossheig"][i_ens,:,:]= propagated_SSH[-1,:,:] 
        
        
    
    return 