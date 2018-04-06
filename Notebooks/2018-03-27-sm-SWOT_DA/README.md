# Assimilation of SWOT data 

## Description

### 1. Experiment setting
#### 1.1. Initialization settings 

###### Definitions

- <span style="color:blue">n_ensemble</span>: number of ensemble members, *<span style="color:gray">integer</span>*
  
###### Imports

- Initialization from initialization_functions: global initialization function calling *initialization*
- initialization from initialization_functions:
    * <span style="color:blue">NATL60state</span>: initialization of a single SSH field (n_ensemble==1) using an internally prescribed SSH field of dimensions

#### 1.2.  Time settings

###### Definitions

- <span style="color:blue">init_date</span>: experimental starting date, *<span style="color:gray">datetime object (yyyy,mm,dd,hh)</span>*
- present_date: current date starting at *init_date* and ending at *final_date* during time loop, *<span style="color:gray">datetime object (yyyy,mm,dd,hh)</span>*
- <span style="color:blue">final_date</span>: experimental ending date, *<span style="color:gray">datetime object (yyyy,mm,dd,hh)</span>*
- <span style="color:blue">assimilation_time_step</span>: assimilation cycle time step, *<span style="color:gray">timedelta object (in hours)</span>* 
- cycle_time: SWOT cycle time step, *<span style="color:gray">timedelta object (in hours)</span>* 

#### 1.3.  Model settings 

###### Definitions
 
###### Imports
- EnsembleModel from Model: global model function propagating an ensemble of state vectors over a time period by calling *model* *n_ensemble*-times 

- model from Model:   
    * <span style="color:blue">QG_SW</span>: quasi-geostrophic shallow-water model propagating an SSH field over a time period 
        
 

#### 1.4.  Observation settings

###### Definitions

- <span style="color:blue">obs_path</span>: observation directory path, *<span style="color:gray">string</span>*  

###### Imports

- ObservationCheck from Observation: function checking if observations are available at assimilation time 


#### 1.5.  Analysis settings

###### Parameters

- a 

###### Functions

- a

#### 1.6.  Outputs settings

###### Definitions

- <span style="color:blue"> outputs</span>: command the "plot" or "save" of state vectors, *<span style="color:gray"> string </span>*

###### Functions

- Save_present_time_outputs from Saveoutputs: save an ensemble state vectors at a definite time in a file

### 2. Initialization

###### Definitions
  
- state_vectors0: ensemble of initial state vectors produced by *Initialization* function, *<span style="color:gray">tuple (3D array of dimensions [n_ensemble,n_lon, n_lat])</span>*
- lon: longitud field produced by *Initialization* function, *<span style="color:gray">2D array of dimensions [n_lon, n_lat]</span>*
- lat: latitud field produced by *Initialization* function, *<span style="color:gray">2D array of dimensions [n_lon, n_lat]</span>*

### 3. Time loop

- i_cycle: temporary cycle counter, *<span style="color:gray">integer</span>* 
- present_date: date of the current cycle, *<span style="color:gray">datetime object (yyyy,mm,dd,hh)</span>*


#### 3.1. Propagation 

###### Definitions
- state_vectors: ensemble of *present_date* state vectors produced by *model* function, *<span style="color:gray">tuple (3D array of dimensions [n_ensemble,n_lon,n_lat])</span>*

#### 3.2. Observation check
#### 3.3. Observation retrieval
#### 3.4. Analysis 
#### 3.5. Outputs 
### 4. Final outputs
 
