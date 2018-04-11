##########################################
##########################################
##     Experiment 1 - Parameters        ##
########################################## 
##########################################
from datetime import datetime
from datetime import timedelta 

 
################################
#   Initialization parameters  #
################################
n_ensemble=1                                                            # ensemble size 
from initialization_functions import NATL60state as initialization      # initialization function
  
################################
#       Time parameters        #
################################
init_date=datetime(2012,10,1,0)                                         # initial date(yyyy,mm,dd,hh)
present_date=init_date                                                  # present date(yyyy,mm,dd,hh)
final_date=datetime(2012,11,1,0)                                        # final date (yyyy,mm,dd,hh)
assimilation_time_step=timedelta(hours=6)                                # assimilation time step   
name_assim_time_step='hourly-6h'

################################
#       Model parameters       #
################################ 
from Model import QG_SW as model                                        # Model function 

################################
#    Observation parameters    #
################################
obs_path="/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/OBS/swot292-OSMOSIS-1h/"
                                                                        # observation path
obs_prefixe="swot292-OSMOSIS"                                           # observation prefixe 
from Observations import ObsOperator_SOSIE_SWOT as obsoperator          # observation operator


################################
#      Analysis parameters     #
################################
 
name_analysis='DirectInsertion'
from Analysis import DirectInsertion as analysisop                
    
################################
#      Outputs parameters      #
################################
saveoutputs=True                                                        # save all outputs
name_exp_save='QGSWOSMO_IniNATL60_'+name_assim_time_step+name_analysis+'_'+str(n_ensemble).zfill(2)+'ens'
path_save='/home/metrefs/Boost-SWOT/2018/Notebooks/2018-03-27-sm-SWOT_DA/OUTPUTS/'+name_exp_save+'/' 
 