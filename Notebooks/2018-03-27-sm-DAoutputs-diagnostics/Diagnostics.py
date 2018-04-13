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


######################
#   SSH field plots  #
######################

def PlotSSH(SSH_DAprods0,nameDAprods0,SSH_Ensforecast0,lon2d0,lat2d0,ntime0):
    timet0=0
    timet1=int(ntime0/2) 
    timet2=ntime0-2 
    
    for i_DAprod in range(np.shape(nameDAprods0)[0]): 

        plt.figure(figsize=(20, 10))

        SSH_DAprod=SSH_DAprods0[i_DAprod]

        meanSSH_DAprod= SSH_DAprod 

        plt.subplot(131)
        plt.pcolormesh(lon2d0,lat2d0,meanSSH_DAprod[timet0,:,:])
        plt.colorbar(extend='both', fraction=0.042, pad=0.04)
        plt.title('ensemble mean SSH(%sd)'%(timet0+1));

        plt.subplot(132)
        plt.pcolormesh(lon2d0,lat2d0,meanSSH_DAprod[timet1,:,:])
        plt.colorbar(extend='both', fraction=0.042, pad=0.04)
        plt.title('ensemble mean SSH(%sd)'%(timet1+1));

        plt.subplot(133)
        plt.pcolormesh(lon2d0,lat2d0,meanSSH_DAprod[timet2,:,:])
        plt.colorbar(extend='both', fraction=0.042, pad=0.04)
        plt.title('ensemble mean SSH(%sd)'%(timet2+1));
 
    plt.figure(figsize=(20, 10))

    meanSSH_Ensforecast=np.mean(SSH_Ensforecast0[:,:,:,:],1) 

    plt.subplot(131)
    plt.pcolormesh(lon2d0,lat2d0,meanSSH_Ensforecast[timet0,:,:])
    plt.colorbar(extend='both', fraction=0.042, pad=0.04)
    plt.title('ensemble mean SSH(%sd)'%(timet0+1));

    plt.subplot(132)
    plt.pcolormesh(lon2d0,lat2d0,meanSSH_Ensforecast[timet1,:,:])
    plt.colorbar(extend='both', fraction=0.042, pad=0.04)
    plt.title('ensemble mean SSH(%sd)'%(timet1+1));

    plt.subplot(133)
    plt.pcolormesh(lon2d0,lat2d0,meanSSH_Ensforecast[timet2,:,:])
    plt.colorbar(extend='both', fraction=0.042, pad=0.04)
    plt.title('ensemble mean SSH(%sd)'%(timet2+1));
    
    return


#################################
#   SSH difference field plots  #
#################################

def PlotPt2PtDiff(SSH_DAprods0,nameDAprods0,SSH_Ensforecast0,SSH_NATL60_degrad0,lon2d0,lat2d0,ntime0):
    timet0=0
    timet1=int(ntime0/2) 
    timet2=ntime0-2 

    max_range=0.02 
    plt.figure(figsize=(20, 10))

    plt.subplot(131)
    plt.pcolormesh(lon2d0,lat2d0,np.mean(SSH_Ensforecast0[timet0,:,:,:],0)-SSH_NATL60_degrad0[timet0,:,:],cmap=plt.cm.get_cmap('bwr'))
    plt.colorbar(extend='both', fraction=0.042, pad=0.04)
    plt.title('DA_SSH-NATL60_SSH(%sd)'%(timet0+1));
    plt.clim(-max_range,max_range)

    plt.subplot(132)
    plt.pcolormesh(lon2d0,lat2d0,np.mean(SSH_Ensforecast0[timet1,:,:,:],0)-SSH_NATL60_degrad0[timet1,:,:],cmap=plt.cm.get_cmap('bwr'))
    plt.colorbar(extend='both', fraction=0.042, pad=0.04)
    plt.title('DA_SSH-NATL60_SSH(%sd)'%(timet1+1));
    plt.clim(-max_range,max_range)

    plt.subplot(133)
    plt.pcolormesh(lon2d0,lat2d0,np.mean(SSH_Ensforecast0[timet2,:,:,:],0)-SSH_NATL60_degrad0[timet2,:,:],cmap=plt.cm.get_cmap('bwr'))
    plt.colorbar(extend='both', fraction=0.042, pad=0.04)
    plt.title('DA_SSH-NATL60_SSH(%sd)'%(timet2+1));
    plt.clim(-max_range,max_range)



    for i_DAprod in range(np.shape(nameDAprods0)[0]): 

        SSH_DAprod=SSH_DAprods0[i_DAprod]

        plt.figure(figsize=(20, 10))

        plt.subplot(131)
        plt.pcolormesh(lon2d0,lat2d0,SSH_DAprod[timet0,:,:]-SSH_NATL60_degrad0[timet0,:,:],cmap=plt.cm.get_cmap('bwr'))
        plt.colorbar(extend='both', fraction=0.042, pad=0.04)
        plt.title('DA_SSH-NATL60_SSH(%sd)'%(timet0+1));
        plt.clim(-max_range,max_range)

        plt.subplot(132)
        plt.pcolormesh(lon2d0,lat2d0,SSH_DAprod[timet1,:,:]-SSH_NATL60_degrad0[timet1,:,:],cmap=plt.cm.get_cmap('bwr'))
        plt.colorbar(extend='both', fraction=0.042, pad=0.04)
        plt.title('DA_SSH-NATL60_SSH(%sd)'%(timet1+1));
        plt.clim(-max_range,max_range)

        plt.subplot(133)
        plt.pcolormesh(lon2d0,lat2d0,SSH_DAprod[timet2,:,:]-SSH_NATL60_degrad0[timet2,:,:],cmap=plt.cm.get_cmap('bwr'))
        plt.colorbar(extend='both', fraction=0.042, pad=0.04)
        plt.title('DA_SSH-NATL60_SSH(%sd)'%(timet2+1));
        plt.clim(-max_range,max_range)
        
    return

########################
# SSH RMSE computation #
########################

def ComputingRMSE(SSH_NATL60_degrad0,SSH_DAprods0,nameDAprods0,SSH_Ensforecast0,ntime1):
    ntime0=ntime1-1
    RMSE_DA=np.zeros(ntime0,)  
    RMSE_Free0=np.zeros(ntime0,)
    nRMSE_DA=np.zeros(ntime0,)  
    nRMSE_Free0=np.zeros(ntime0,)


    print(ntime0)

    nlon=np.shape(SSH_NATL60_degrad0)[1]
    nlat=np.shape(SSH_NATL60_degrad0)[2]
    RMSE_DAs0=[]
    nRMSE_DAs0=[]
    for i_DAprod in range(np.shape(nameDAprods0)[0]):  
        SSH_DAprod=SSH_DAprods0[i_DAprod]
        RMSE_DA=np.zeros(ntime0,)  
        for itime in range(ntime0): 
            RMSE_DA[itime]=np.sqrt(np.sum(np.sum(np.square(SSH_DAprod[itime,:,:]-SSH_NATL60_degrad0[itime,:,:])))/nlon/nlat)
            RMSE_Free0[itime]=np.sqrt(np.sum(np.sum(np.square(np.mean(SSH_Ensforecast0[itime,:,:,:],0)-SSH_NATL60_degrad0[itime,:,:])))/nlon/nlat)

            nRMSE_DA[itime]=np.sqrt(np.sum(np.sum(np.square(SSH_DAprod[itime,:,:]-SSH_NATL60_degrad0[itime,:,:])))/nlon/nlat) / ( np.max(np.max(SSH_NATL60_degrad0[itime,:,:]))-np.min(np.min(SSH_NATL60_degrad0[itime,:,:])) )
            nRMSE_Free0[itime]=np.sqrt(np.sum(np.sum(np.square(np.mean(SSH_Ensforecast0[itime,:,:,:],0)-SSH_NATL60_degrad0[itime,:,:])))/nlon/nlat) / ( np.max(np.max(SSH_NATL60_degrad0[itime,:,:]))-np.min(np.min(SSH_NATL60_degrad0[itime,:,:])) )

        RMSE_DAs0.append(RMSE_DA)
        nRMSE_DAs0.append(nRMSE_DA)
        
    return RMSE_DAs0,nRMSE_DAs0,RMSE_Free0,nRMSE_Free0
    

#####################
#   SSH RMSE plots  #
#####################
def PlotRMSE(RMSE_DAs0,nRMSE_DAs0,RMSE_Free0,nRMSE_Free0,nameDAprods0):
    spinup=2*24
    plt.figure(figsize=(15, 5)) 
    for i_DAprod in range(np.shape(nameDAprods0)[0]):
        RMSE_DA=RMSE_DAs0[i_DAprod]
        plt.plot(RMSE_DA[spinup:], label=nameDAprods0[i_DAprod])    
    plt.plot(RMSE_Free0[spinup:],'k--', label='Free')
    plt.legend()
    plt.xlabel('Time (hours)')
    plt.ylabel('RMSE')
    plt.title('SSH')

    plt.figure(figsize=(15, 5)) 
    for i_DAprod in range(np.shape(nameDAprods0)[0]):
        nRMSE_DA=nRMSE_DAs0[i_DAprod]
        plt.plot(nRMSE_DA[spinup:]*100., label=nameDAprods0[i_DAprod])  
    plt.plot(nRMSE_Free0[spinup:]*100.,'k--', label='Free') 
    plt.legend()
    plt.xlabel('Time (hours)')
    plt.ylabel('Normalized RMSE (%)') 
    plt.title('SSH')

    return

######################
# Spectral coherence #
######################
def psd1d(hh=None,dx=1.,tap=0.05, detrend=True):


  hh=hh-np.mean(hh)
  nx=np.shape(hh)[0]


  if detrend:
    hh=scipy.signal.detrend(hh)

  if tap>0:  
    ntaper = np.int(tap * nx + 0.5)
    taper = np.zeros(nx)+1.
    taper[:ntaper]=np.cos(np.arange(ntaper)/(ntaper-1.)*pi/2+3*pi/2)
    taper[-ntaper:] = np.cos(-np.arange(-ntaper+1,1)/(ntaper-1.)*pi/2+3*pi/2)
    hh=hh*taper

  ss=fft(hh)

  ff=np.arange(1,nx/2-1)/(nx*dx)

  PSD=2*dx/(nx)*np.abs(ss[1:nx/2-1])**2


  return ff, PSD


def cpsd1d(hh1=None,hh2=None,dx=1.,tap=0.05, detrend=True):


  hh1=hh1-np.mean(hh1)
  hh2=hh2-np.mean(hh2)
  nx=np.shape(hh1)[0]


  if detrend:
    hh1=scipy.signal.detrend(hh1)
    hh2=scipy.signal.detrend(hh2)

  if tap>0:
    ntaper = np.int(tap * nx + 0.5)
    taper = np.zeros(nx)+1.
    taper[:ntaper]=np.cos(np.arange(ntaper)/(ntaper-1.)*pi/2+3*pi/2)
    taper[-ntaper:] = np.cos(-np.arange(-ntaper+1,1)/(ntaper-1.)*pi/2+3*pi/2)
    hh1=hh1*taper
    hh2=hh2*taper

  ss1=fft(hh1)
  ss2=fft(hh2) 
      
  ff=np.arange(1,nx/2-1)/(nx*dx) 

  C=np.cos(np.angle(ss1[1:int(nx/2)-1])-np.angle(ss2[1:int(nx/2)-1]))

  return ff, C