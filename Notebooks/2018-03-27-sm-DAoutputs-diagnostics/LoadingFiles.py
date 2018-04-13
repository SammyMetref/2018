import sys,os,shutil
import numpy as np
import matplotlib.pylab as plt
import time
import netCDF4 as nc 
from datetime import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta 
import scipy
from scipy import signal 
import pdb
from math import pi
from scipy.fftpack import fft


####################
# Loading Free run #
####################

def LoadEnsforecast(init_date0,final_date0,deltat0,ntime0,nens0) :
    present_date0=init_date0
    # Ensforecast SSH
    name_assim_time_step='hourly'
    # First time step
    itime=0
    year1=str(present_date0.year)
    month1=str(present_date0.month).zfill(2)
    day1=str(present_date0.day).zfill(2)
    hour1=str(present_date0.hour).zfill(2)
    name_exp='QGSWOSMO_IniNATL60_'+name_assim_time_step+'EnsForecast_'+str(nens0).zfill(2)+'ens'
    file_deg='/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/DA_OSMOSIS/'+name_exp+'/'+name_exp+'_y'+year1+'m'+month1+'d'+day1+'h'+hour1+'_SSHdegrad.nc'
    fid_deg = nc.Dataset(file_deg)
    lon2d0=np.array(fid_deg.variables["nav_lon"][:,:]) 
    lat2d0=np.array(fid_deg.variables["nav_lat"][:,:])  
    SSH_Ensforecast_init=np.array(fid_deg.variables["degraded_sossheig"][:,:,:]) 
    SSH_Ensforecast0=np.zeros([ntime0,np.shape(SSH_Ensforecast_init)[0],np.shape(SSH_Ensforecast_init)[1],np.shape(SSH_Ensforecast_init)[2]],)
    SSH_Ensforecast0[itime,:,:,:]=SSH_Ensforecast_init
    present_date0=present_date0+deltat0
    itime=itime+1   

    # Time loop
    while present_date0<final_date0: 
        # OSMOSIS: SSH from NATL60, daily outputs, degradated grid (3x3)
        year1=str(present_date0.year)
        month1=str(present_date0.month).zfill(2)
        day1=str(present_date0.day).zfill(2)
        hour1=str(present_date0.hour).zfill(2)
        file_deg='/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/DA_OSMOSIS/'+name_exp+'/'+name_exp+'_y'+year1+'m'+month1+'d'+day1+'h'+hour1+'_SSHdegrad.nc'
        fid_deg = nc.Dataset(file_deg)
        lon2d0=np.array(fid_deg.variables["nav_lon"][:,:]) 
        lat2d0=np.array(fid_deg.variables["nav_lat"][:,:])    
        SSH_Ensforecast0[itime,:,:,:]=np.array(fid_deg.variables["degraded_sossheig"][:,:,:]) 
        present_date0=present_date0+deltat0
        itime=itime+1   
    
    return SSH_Ensforecast0,lon2d0,lat2d0
    



####################
# Loading DA runs  #
####################

def LoadEnsSSHDAprod(init_date0,final_date0,deltat0,ntime0,nameDAprods0) :
    
    SSH_DAprods0=[]

    for i_DAprod in range(np.shape(nameDAprods0)[0]):  
        present_date0=init_date0  

        # Da product SSH
        nameDAprod=nameDAprods0[i_DAprod]  
        # First time step
        itime=0
        year1=str(present_date0.year)
        month1=str(present_date0.month).zfill(2)
        day1=str(present_date0.day).zfill(2)
        hour1=str(present_date0.hour).zfill(2)
        name_exp='QGSWOSMO_Ini'+nameDAprod
        file_deg='/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/DA_OSMOSIS/'+name_exp+'/'+name_exp+'_y'+year1+'m'+month1+'d'+day1+'h'+hour1+'_SSHdegrad.nc'
        fid_deg = nc.Dataset(file_deg)
        lon2d0=np.array(fid_deg.variables["nav_lon"][:,:]) 
        lat2d0=np.array(fid_deg.variables["nav_lat"][:,:])  
        SSH_DAprod_init=np.mean(np.array(fid_deg.variables["degraded_sossheig"][:,:,:]) ,0)
        SSH_DAprod=np.zeros([ntime0,np.shape(SSH_DAprod_init)[0],np.shape(SSH_DAprod_init)[1]],)
        SSH_DAprod[itime,:,:]=SSH_DAprod_init
        present_date0=present_date0+deltat0
        itime=itime+1  

        # Time loop
        while present_date0<final_date0:   
            year1=str(present_date0.year)
            month1=str(present_date0.month).zfill(2)
            day1=str(present_date0.day).zfill(2)
            hour1=str(present_date0.hour).zfill(2)
            file_deg='/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/DA_OSMOSIS/'+name_exp+'/'+name_exp+'_y'+year1+'m'+month1+'d'+day1+'h'+hour1+'_SSHdegrad.nc'
            fid_deg = nc.Dataset(file_deg) 
            lon2d0=np.array(fid_deg.variables["nav_lon"][:,:]) 
            lat2d0=np.array(fid_deg.variables["nav_lat"][:,:])  
            SSH_DAprod[itime,:,:]=np.mean(np.array(fid_deg.variables["degraded_sossheig"][:,:,:]) ,0) 
            present_date0=present_date0+deltat0
            itime=itime+1       

        SSH_DAprods0.append(SSH_DAprod) 

    return SSH_DAprods0,lon2d0,lat2d0



##################################
# Loading NATL60 run (reference) #
##################################

def LoadingNATL60SSH(final_date0,deltat0,ntime0):

    init_date1=datetime(2012,10,1,0)  
    deltatmonth=relativedelta(months=1)
    present_date1=init_date1 
    final_date1=final_date0

    ntime1=(final_date1-present_date1).days*24+1 

    # Warning: Script made for 1h outputs !!!
    # OSMOSIS: SSH from NATL60, daily outputs, degradated grid (3x3)

    file_deg='/mnt/meom/MODEL_SET/NATL60/NATL60-CJM165-S/1h/OSMOSIS/NATL60OSMO-CJM165_y2012m10.1h_SSHdegrad.nc'
    fid_deg = nc.Dataset(file_deg)
    lon2d=np.array(fid_deg.variables["nav_lon"][:,:]) 
    lat2d=np.array(fid_deg.variables["nav_lat"][:,:])  
    SSH_NATL60_degrad1=np.zeros([ntime1,np.shape(lon2d)[0],np.shape(lat2d)[1]],)
    itime=0 
    while present_date1 < final_date1:
        SSHindex=np.arange(int( (min(final_date1,present_date1+deltatmonth)-present_date1).days*24+ (min(final_date1,present_date1+deltatmonth)-present_date1).seconds/3600)) 
        year1=str(present_date1.year)
        month1=str(present_date1.month).zfill(2)
        day1=str(present_date1.day).zfill(2) 
        file_deg='/mnt/meom/MODEL_SET/NATL60/NATL60-CJM165-S/1h/OSMOSIS/NATL60OSMO-CJM165_y'+year1+'m'+month1+'.1h_SSHdegrad.nc'
        fid_deg = nc.Dataset(file_deg)
        lon2d=np.array(fid_deg.variables["nav_lon"][:,:]) 
        lat2d=np.array(fid_deg.variables["nav_lat"][:,:])     
        SSH_NATL60_degrad1[itime:itime+SSHindex[-1]+1,:,:]=np.array(fid_deg.variables["degraded_sossheig"][SSHindex,:,:]) 
        present_date1=present_date1+deltatmonth 
        itime=itime+SSHindex[-1]+1     
 
        
        
    SSH_NATL60_degrad0=SSH_NATL60_degrad1[1:,:,:]    
 
    return SSH_NATL60_degrad0

    