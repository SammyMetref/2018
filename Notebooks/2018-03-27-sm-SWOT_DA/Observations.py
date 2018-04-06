################################
# SSH initialization functions #
################################
import netCDF4 as nc
import numpy as np
import os.path 


def ObservationCheck(time1,deltat,observationpath,observationprefixe,ncycle,pass_names,initobsdate):
    """
    NAME 
        ObservationCheck 

    DESCRIPTION
        ObservationCheck verifies if observations are available in [time1-deltat,time1+deltat]

        Args:
            time1 (datetime object): current time 
            deltat (datetime.timedelta object): data assimilation cycle time  
            observationpath (string): path to observations
            observationprefixe (string): prefixe of observations
            ncycle (int): cycle number
            initobsdate (datetime object): initial observation time (yyyy,mm,dd,01,01,01 of the NATL60 input of SWOT simulator)

        Returns:
            observations_available (boolean): Informs if observations are available in [time1-deltat,time1+deltat]
            file_obs (1D string array): the list of file names where observations where encountered
             
    """          
    obspathpref=observationpath+observationprefixe+'_c'+str(ncycle).zfill(2) 
    
    time1mdeltat_seconds_fromfirstobs=(time1-deltat/2-initobsdate).total_seconds()
    time1pdeltat_seconds_fromfirstobs=(time1+deltat/2-initobsdate).total_seconds()
     
    
    obs_encountered=0 
    file_obs=np.array([])
    for pass_name in pass_names:
        file=obspathpref+'_'+pass_name+'.nc'  
        if os.path.isfile(file): 
            fid = nc.Dataset(file)
            obstimes=np.array(fid.variables["time_sec"][:]) 
            check_obs_in_pass=(time1mdeltat_seconds_fromfirstobs < obstimes[0]) and (obstimes[0] < time1pdeltat_seconds_fromfirstobs) or (time1mdeltat_seconds_fromfirstobs < obstimes[-1]) and (obstimes[-1] < time1pdeltat_seconds_fromfirstobs)  
            
            if check_obs_in_pass: 
                file_obs=np.append(file_obs,file)
                obs_encountered=obs_encountered + 1    


    if obs_encountered!=0: 
        return [True,file_obs]   
         
    return [False,None]

def GetObsTime(observation_path,observation_name):
    """
    NAME 
        GetObsTime 

    DESCRIPTION
        GetObsTime retrive the time of *observation_name* at *observation_path*

        Args: 
            observation_path (string): path to observation
            observation_name (string): name of observation

        Returns:
            obs_time (float, in seconds): observation time 
    """    
    
    file=observation_path+observation_name 
    fid = nc.Dataset(file)
    obstimes=np.array(fid.variables["time_sec"][:])
    obstime=obstimes[0]
    
    return obstime


def ObsOperator(function,state_vectors_names,observation_name,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,*args):
    """
    NAME 
        ObsOperator 

    DESCRIPTION
        Global observation operator function

        Args:
            function (function): observation operator function 
            state_vectors_names (array of strings): path and names of state_vectors
            observation_name (array of strings): path and name of observation
            n_ensemble (integer): number of ensemble members

        Returns:
            function(state_vectors_names,n_ens,*args): output of *function*  
    """            
            
    return function(state_vectors_names,observation_name,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens,*args)


def ObsOperator_SOSIE_SWOT(state_vectors_names,observation_name,tmp_DA_path,sosie_path,name_sosie_output,name_sosie_map,n_ens=1):
    """
    NAME 
        ObsOperator_SOSIE_SWOT

    DESCRIPTION 
        Call SOSIE to create the SWOT observation projection of the model states (from *state_vectors_names* files)
        
        Args: 
            state_vectors_names (array of strings): path and names of state_vectors
            observation_name (array of strings): path and name of observation
            n_ens (integer): number of ensemble members 

        Returns: 
            state_projections_names (array of strings): path and names of model state projection onto the SWOT observation plane 

    """ 
     
    
    state_projection_names=[]
    for i_ens in range(n_ens):
        
        # Erasing previous state_projection
        name_output="state_projections_"+str(i_ens).zfill(2)+".nc"
        cmd1="rm "+tmp_DA_path+name_output    
        os.system(cmd1)
        
        # Seperating state_vectors by ensemble member
        stringiens=str(i_ens).zfill(2)
        state_vector_iens_name=state_vectors_names[:-4]+"_"+stringiens+".nc" 
        
        cmd2="ncks -d member,"+stringiens+","+stringiens+" "+state_vectors_names+" "+state_vector_iens_name 
        os.system(cmd2) 
        
        file=state_vector_iens_name 
        if os.path.isfile(file)==False:
            print('Error: No state_vector'+str(i_ens).zfill(2)+' file in '+tmp_DA_path)  
        
        #
        # SOSIE reads source.nc and target.nc and writes output in TMP_DA/ directory (details in namelist)
        #
        
        # Copying state_vector_iens and observation as source and target for SOSIE
        cmd3="cp "+state_vector_iens_name+" "+tmp_DA_path+"source.nc"  
        os.system(cmd3)  
        cmd4="cp "+observation_name+" "+tmp_DA_path+"target.nc" 
        os.system(cmd4)  
               
        # Running SOSIE 
        cmd5=sosie_path+"sosie.x -f "+sosie_path+"namelist1"
        os.system(cmd5)  
            
        # Renaming SOSIE output as the i_ens state_projection 
        if os.path.isfile(tmp_DA_path+name_sosie_output)==False:
            print('Error: SOSIE failed to work')
        cmd6="mv "+tmp_DA_path+name_sosie_output+" "+tmp_DA_path+name_output  
        os.system(cmd6)  
        
        state_projection_names=np.append(state_projection_names,tmp_DA_path+name_output) 
            
    # Erasing SOSIE map 
    cmd7="rm "+tmp_DA_path+name_sosie_map    
    os.system(cmd7)
    
    return state_projection_names


             