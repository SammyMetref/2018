################################
#       Analysis functions        #
################################ 
import sys,os,shutil 
import numpy as np 
import netCDF4 as nc 
import matplotlib.pylab as plt

from Observations import *
 
from params import * 
 

############################
# Global analysis function #
############################

def Analysis(function,state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,obsop,*args):
    """
    NAME 
        Analysis 

    DESCRIPTION
        Global analysis function

        Args: 
            function (function): model function
            state_vectors_names (list of strings): ensemble of path and names of the current state_vectors files (all members) 
            n_ens (integer): number of ensemble member  

        Returns: 
            analyzed_vectors_names (string): ensemble of path and names of the analyzed_vectors files (all members)

    """   
 
    analyzed_vectors_names=function(state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,obsop,*args) 
          
        
    return analyzed_vectors_names
        

######################
# Analysis functions #
######################

def NoAnalysis(state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,obsop): 
    """
    NAME 
        Replacement

    DESCRIPTION 
        Function replacing state vectors by observations

        Args: 
            state_vectors_names (string): path and name of the current state_vectors file    
            obs_file (string): path and name of the current observation file 
            n_ens (integer): number of ensemble members
 

        Returns: 
            analyzed_vectors_names (string): ensemble of path and names of the analyzed_vectors files (all members)

    """      
    
    
    
    analyzed_vectors_names=tmp_DA_path+'state_analyzed.nc'
    
    cmd1="cp "+state_vectors_names+" "+analyzed_vectors_names 
        
    
    return analyzed_vectors_names

def DirectInsertion(state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,obsop): 
    """
    NAME 
        DirectInsertion

    DESCRIPTION 
        Function replacing state vectors by observations

        Args: 
            state_vectors_names(string): path and name of the current state_vectors file    
            obs_file (string): path and name of the current observation file   
            tmp_DA_path (string): temporary directory path
            sosie_path (string): sosie directory path
            name_sosie_output (string): name of sosie output 
            name_sosie_map (string): name of sosie map (2nd sosie output)
            n_ens (integer): number of ensemble members
 

        Returns: 
            analyzed_vectors_names (string): ensemble of path and names of the analyzed_vectors files (all members)

    """      
    
    analyzed_vectors_names=tmp_DA_path+'state_analyzed.nc'
     
    
    n_var=1                               # To be moved to Exp1_params.py
    name_var=["degraded_sossheig"]        # To be moved to Exp1_params.py
    
    # Retrieve forecast as ensemble of state vectors
    [forecast,lon2d,lat2d] = EnsState_vectorize(state_vectors_names,n_ens,n_var,name_var) 
    
    # Retrieve observation as ensemble of state vectors (H^{-1})
    obs_inv_file=obs_file[:-3]+"_inv.nc"
    [obs_inv,lon2d,lat2d] = EnsState_vectorize(obs_inv_file,1,n_var,["ssh"],name_lon="nav_lon",name_lat="nav_lat")#,name_var) 
    # Note: The inverse observations must be created offline (2018-03-27-sm-pretreatment-SWOTinputs)
      
    # Analysis default definition 
    analysis = np.zeros_like(forecast)
    analysis[:,:]=forecast[:,:]
     
          
    for i_ens in range(n_ens):
        analysis[i_ens,obs_inv[0,:]>-50] = obs_inv[0,obs_inv[0,:]>-50]  
         
    
    EnsStateVector_save(analysis,lon2d,lat2d,n_ens,n_var,name_var,analyzed_vectors_names) 
    
      
    
    return analyzed_vectors_names


def Nudging(state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,obsop): 
    """
    NAME 
        Replacement 

    DESCRIPTION 
        Function replacing state vectors by observations

        Args: 
            state_vectors_names(string): path and name of the current state_vectors file    
            obs_file (string): path and name of the current observation file   
            tmp_DA_path (string): temporary directory path
            sosie_path (string): sosie directory path
            name_sosie_output (string): name of sosie output 
            name_sosie_map (string): name of sosie map (2nd sosie output)
            n_ens (integer): number of ensemble members
 

        Returns: 
            analyzed_vectors_names (string): ensemble of path and names of the analyzed_vectors files (all members)

    """      
     
    
    analyzed_vectors_names=tmp_DA_path+'state_analyzed.nc'
     
    
    n_var=1                               # To be moved to Exp1_params.py
    name_var=["degraded_sossheig"]        # To be moved to Exp1_params.py
    
    # Retrieve forecast as ensemble of state vectors
    [forecast,lon2d,lat2d] = EnsState_vectorize(state_vectors_names,n_ens,n_var,name_var) 
    
    # Retrieve observation as ensemble of state vectors (H^{-1})
    obs_inv_file=obs_file[:-3]+"_inv.nc"
    [obs_inv,lon2d,lat2d] = EnsState_vectorize(obs_inv_file,1,n_var,["ssh"],name_lon="nav_lon",name_lat="nav_lat")#,name_var) 
    # Note: The inverse observations must be created offline (2018-03-27-sm-pretreatment-SWOTinputs)
     
    # Analysis default definition 
    analysis = np.zeros_like(forecast)
    analysis[:,:]=forecast[:,:]
      
    for i_ens in range(n_ens):
        analysis[i_ens,obs_inv[0,:]>-50] = obs_inv[0,obs_inv[0,:]>-50]  
         
    
    EnsStateVector_save(analysis,lon2d,lat2d,n_ens,n_var,name_var,analyzed_vectors_names) 
    
      
    
    return analyzed_vectors_names


####################################
# Other analysis related-functions #
####################################

def Innovation(state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,obsop): 
    """
    NAME 
        Innovation

    DESCRIPTION 
        Function compute innovation using ObsOperator

        Args: 
            state_vectors_names(string): path and name of the current state_vectors file    
            obs_file (string): path and name of the current obs_file file    
            n_ens (integer): number of ensemble members
 

        Returns: 
            analyzed_vectors_names (string): ensemble of path and names of the analyzed_vectors files (all members)

    """     
                    
    state_proj_name=ObsOperator(obsop,state_vectors_names,obs_file,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens)
    
    
    fid_deg1 = nc.Dataset(obs_file)
    lon2d=np.array(fid_deg1.variables["lon"][:,:]) 
    lat2d=np.array(fid_deg1.variables["lat"][:,:])   
    obs=np.array(fid_deg1.variables["ssh_model"][:,:])  
    
    n_obs_var=1                 # To be moved to Exp1_params.py
    name_obs_var=["ssh"]        # To be moved to Exp1_params.py
    
    fid_deg = nc.Dataset(state_proj_name[0])
    n_tot_obs=0
    for i_var in range(n_obs_var):
        n_tot_obs=n_tot_obs+np.shape(fid_deg.variables[name_obs_var[i_var]])[1]*np.shape(fid_deg.variables[name_obs_var[i_var]])[2]
    innov=np.zeros([n_ens,n_tot_obs],) 
     
    for i_ens in range(n_ens):
        fid_deg = nc.Dataset(state_proj_name[i_ens])
        lon2d=np.array(fid_deg.variables["lon"][:,:]) 
        lat2d=np.array(fid_deg.variables["lat"][:,:])   
        i_innov=0
        for i_var in range(n_obs_var):                                               
            for i_lon in range(np.shape(fid_deg.variables["lon"][:,:])[0]):
                for j_lat in range(np.shape(fid_deg.variables["lat"][:,:])[1]):
                    innov[i_ens,i_innov]=np.array(fid_deg.variables[name_obs_var[i_var]][0,i_lon,j_lat]) - obs[i_lon,j_lat]
                    i_innov=i_innov+1
    
    innov[(innov<=-9900)]=float('Inf') 
       
    
    return innov


def EnsState_vectorize(state_vectors_names,n_ens,n_var,name_var,name_lon="nav_lon",name_lat="nav_lat"):
    """
    NAME 
        EnsState_vectorize

    DESCRIPTION 
        Create a 2D-array of [ensemble members, total variables]

        Args: 
            state_vectors_names (string): path and name of the current state_vectors file        
            n_ens (integer): number of ensemble members
            n_var (integer): number of state variables
            name_var (list of strings): names of state variables
 

        Returns: 
            state_vectors (2D-array): 'n_ens' vectors of the vectorize ensemble member states

    """       
    
    fid_deg = nc.Dataset(state_vectors_names)
    lon=np.array(fid_deg.variables[name_lon][:,:])
    lat=np.array(fid_deg.variables[name_lat][:,:]) 
    n_tot=n_var*np.shape(lon)[0]*np.shape(lat)[1]
    state_vectors=np.zeros([n_ens,n_tot],)
     
    
    i_tot=0
    for i_var in range(n_var):  
        for i_lon in range(np.shape(lon)[0]):
            for j_lat in range(np.shape(lat)[1]):
                state_vectors[:,i_tot]=np.array(fid_deg.variables[name_var[i_var]][:,i_lon,j_lat]) 
                i_tot=i_tot+1
                
    return [state_vectors, lon, lat]
            

def EnsStateVector_save(state_vectors,lon,lat,n_ens,n_var,name_var,state_vectors_names):
    """
    NAME 
        EnsStateVector_save

    DESCRIPTION 
        Save a netcdf file from a 2D-array of [ensemble members, total variables]

        Args: 
            state_vectors_names (string): path and name of the current state_vectors file   
            lon (array): longitud of state_vectors 
            lat (array): latitud of state_vectors
            n_ens (integer): number of ensemble members
            n_var (integer): number of state variables
            name_var (list of strings): names of state variables
            state_vectors_names (string): name of 'EnsStateVector_save' saved file
 

        Returns: 
            state_vectors (2D-array): 'n_ens' vectors of the vectorize ensemble member states

    """      


    ncout = nc.Dataset(state_vectors_names, 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('member', n_ens)
    ncens = ncout.createVariable('ens', 'd', ('member',))
    ncens[:] = range(n_ens) 
    ncout.createDimension('x', lon.shape[0])
    ncout.createDimension('y', lat.shape[1])   
    nclon = ncout.createVariable('nav_lon', 'f', ('x','y',))
    nclat = ncout.createVariable('nav_lat', 'f', ('x','y',))  
    nclat[:,:] = lat
    nclon[:,:] = lon   
    i_tot=0
    for i_var in range(n_var):
        nchei = ncout.createVariable(name_var[i_var], 'f', ('member','x','y',)) 
        for i_lon in range(np.shape(lon)[0]):
            for j_lat in range(np.shape(lat)[1]):
                nchei[:,i_lon,j_lat] = state_vectors[:,i_tot]
                i_tot=i_tot+1
                
    if False: 
        plt.figure()
        plt.pcolormesh(lon,lat,nchei[0,:,:])
        plt.colorbar()
        plt.show()
    
    ncout.close()
                
    return

