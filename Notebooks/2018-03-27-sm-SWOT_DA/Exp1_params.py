##########################################
##########################################
##     Experiment 1 - Parameters        ##
########################################## 
name_experiment = 'Exp1'
##########################################
from datetime import datetime
from datetime import timedelta 

 
################################
#   Initialization parameters  #
################################
n_ensemble = 1                                                             # ensemble size 
# NATL60state, NoisyNATL60state
from initialization_functions import NoisyNATL60state as initialization  # initialization function
sigma_initnoise = 0.01                                                   # standard-deviation for initial noise (in NoisyNATL60)
name_assim_init = 'NoisyNATL60state'
  
################################
#       Time parameters        #
################################
init_date = datetime(2012,10,1,0)                                         # initial date(yyyy,mm,dd,hh)
present_date = init_date                                                  # present date(yyyy,mm,dd,hh)
final_date = datetime(2012,11,1,0)                                        # final date (yyyy,mm,dd,hh)
propagation_time_step = timedelta(hours=1)                                # assimilation time step 
assimilation_time_step = timedelta(hours=3)                               # assimilation time step   
name_assim_time_step = 'hourly-3h'

################################
#       Model parameters       #
################################ 
from Model import QG_SW as model                                        # Model function 

################################
#    Observation parameters    #
################################
obs_path = "/mnt/meom/workdir/metrefs/Boost-SWOT/2018/DATA/OBS/swot292-OSMOSIS-1h/"
                                                                        # observation path
obs_prefixe = "swot292-OSMOSIS"                                           # observation prefixe 
from Observations import ObsOperator_SOSIE_SWOT as obsoperator          # observation operator


################################
#      Analysis parameters     #
################################
 
name_analysis = 'Nudging'
from Analysis import Nudging as analysisop    

# Nudging
nudging_coef = 0.6                                                      # nudging towards obs. coeff: no (=0) or full (=1) nudging
    
################################
#      Outputs parameters      #
################################
saveoutputs = True                                                        # save all outputs
name_exp_save = name_experiment+'_QGSWOSMO_Ini'+name_assim_init+'_'+name_assim_time_step+'_'+name_analysis+'_'+str(n_ensemble).zfill(2)+'ens'
path_save = '/home/metrefs/Boost-SWOT/2018/Notebooks/2018-03-27-sm-SWOT_DA/OUTPUTS/'+name_exp_save+'/' 
 