##########################################################



# Loading libraries
source("include_library.R")

# Loading initial settings
source("initial_setting.R")
source("read_rawdata.R")
source("cut_region.R")

# Loading functions
source("R/CommonGrid.R")
source("R/ComputePOD.R")
source("R/FitKriging.R")
source("R/GPdev.ind.R")
source("R/PredNewCase.R")
sourceCpp("R/Cpp/bigInnerProd.cpp")
sourceCpp("R/Cpp/bigInnerProd2.cpp")
sourceCpp("R/Cpp/computeDev.cpp")
sourceCpp("R/Cpp/computeS.cpp")
sourceCpp("R/Cpp/computeR.cpp") # add on 07/27/2016
sourceCpp("R/Cpp/GPpred.cpp")

### read design data ###
design.df <- read.csv(paste0(working_dir, "/", design_file_name), header = FALSE)
design.df <- design.df[Case.Index,]
num.case <- nrow(design.df) # set number of sample cases

design.df <- design.df[,-1]
design.df <- matrix(design.df, ncol = 1)
  
### Step1. Common Grid ###
Golden_Case <- CommonGrid(num.case = num.case, 
                          Number.NN = 10, output_path = working_dir)

### Step2. Compute POD ###
ComputePOD(num.case = num.case, 
           energy.cutoff = POD_Mode_Energy.cutoff, 
           golden_case = Golden_Case, 
           CommonGrid_filepath = working_dir, 
           POD_maxnum.cutoff = POD_Mode_number.cutoff,
           output_path = paste0(working_dir,"/POD_results"))

### Step3. Fit Kriging Model ###
FitKriging(num.case = num.case, des = design.df, 
           ind_assumption.fg = TRUE, 
           POD_filepath = paste0(working_dir,"/POD_results"), 
           output_path = paste0(working_dir,"/Kriging_results"))

### Step4. Prediction ###
newCase <- matrix(c(10.5), ncol = 1)
PredNewCase(num.case = num.case, newCase = newCase, 
            golden_case = Golden_Case, 
            des = design.df, 
            CommonGrid_filepath = working_dir, 
            POD_filepath = paste0(working_dir,"/POD_results"), 
            Kriging_filepath = paste0(working_dir,"/Kriging_results"), 
            output_path = paste0(working_dir,"/Prediction_results"))

