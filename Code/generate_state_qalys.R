# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters state utilities

generate_state_qalys <- function(input_parameters, 
                                 treatment_names = treatment_names, 
                                 state_names = state_names,
                                 initial_age = initial_age,
                                 final_age = final_age,
                                 starting_age = starting_age,
                                 gender = gender,
                                 sensitivity = NULL) {
  
  n_treatments <- length(treatment_names)
  n_samples <- dim(input_parameters)[1]
  n_states <- length(state_names)
  
  utilities_primary <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_primary")
  utilities_revision <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_revision")
  
  #un_utilities <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_unadjusted")
  
 
  # Construct array of state utilities with default value 0
  state_utilities <- array(dim = c(n_samples, n_treatments, n_states), dimnames = list(NULL, treatment_names,state_names))
  event_disutilities <- array(dim = c(n_samples, n_treatments, n_states), dimnames = list(NULL, treatment_names,state_names))
  # State utilities are a mixture of state event  utilities, and disutilities of transient events
  multi_disutilities = rlnorm(n_samples, log(4.5), log(1.5))
  utilities_revision_row_index <- which(utilities_revision[, "age"] == ini_age &
                                          utilities_revision[, "gender"] == gender)
  
  utilities_primary_row_index <- which(utilities_primary[, "treatment"] == "Cemented" &
                                         utilities_primary[, "age"] == ini_age &
                                         utilities_primary[, "gender"] == gender)
  
  disutilities_mean_Cem = as.numeric(utilities_revision[utilities_revision_row_index,3])- as.numeric(utilities_primary[utilities_primary_row_index,7])
  disutilities_UL_Cem = as.numeric(utilities_revision[utilities_revision_row_index,5])- as.numeric(utilities_primary[utilities_primary_row_index,9])
  disutilities_LL_Cem = as.numeric(utilities_revision[utilities_revision_row_index,4])- as.numeric(utilities_primary[utilities_primary_row_index,8])
  
  norm_disutilities_revision_Cem <- rnorm(n_samples, mean = disutilities_mean_Cem,  
                             sd = (disutilities_UL_Cem-disutilities_LL_Cem)/2*1.96)
  
  utilities_primary_row_index <- which(utilities_primary[, "treatment"] == "Uncemented" &
                                         utilities_primary[, "age"] == ini_age &
                                         utilities_primary[, "gender"] == gender)
  
  disutilities_mean_Unc = as.numeric(utilities_revision[utilities_revision_row_index,3])- as.numeric(utilities_primary[utilities_primary_row_index,7])
  disutilities_UL_Unc = as.numeric(utilities_revision[utilities_revision_row_index,5])- as.numeric(utilities_primary[utilities_primary_row_index,9])
  disutilities_LL_Unc = as.numeric(utilities_revision[utilities_revision_row_index,4])- as.numeric(utilities_primary[utilities_primary_row_index,8])
  
  norm_disutilities_revision_Unc <- rnorm(n_samples, mean = disutilities_mean_Unc,  
                                          sd = (disutilities_UL_Unc-disutilities_LL_Unc)/2*1.96)
  
  utilities_primary_row_index <- which(utilities_primary[, "treatment"] == "Hybrid" &
                                         utilities_primary[, "age"] == ini_age &
                                         utilities_primary[, "gender"] == gender)
  
  disutilities_mean_Hyb = as.numeric(utilities_revision[utilities_revision_row_index,3])- as.numeric(utilities_primary[utilities_primary_row_index,7])
  disutilities_UL_Hyb = as.numeric(utilities_revision[utilities_revision_row_index,5])- as.numeric(utilities_primary[utilities_primary_row_index,9])
  disutilities_LL_Hyb = as.numeric(utilities_revision[utilities_revision_row_index,4])- as.numeric(utilities_primary[utilities_primary_row_index,8])
  
  norm_disutilities_revision_Hyb <- rnorm(n_samples, mean = disutilities_mean_Hyb,  
                                          sd = (disutilities_UL_Hyb-disutilities_LL_Hyb)/2*1.96)
  norm_disutilities_revision <- cbind(
    norm_disutilities_revision_Cem,
    norm_disutilities_revision_Unc,
    norm_disutilities_revision_Hyb
  )
  colnames(norm_disutilities_revision) <- treatment_names
  
  norm_disutilities_rerevision <- rnorm(n_samples, mean = as.numeric(utilities_revision[utilities_revision_row_index,15]),
                                        sd = (as.numeric(utilities_revision[utilities_revision_row_index,17])-as.numeric(utilities_revision[utilities_revision_row_index,16]))/2*1.96)
  
  
  #if(!is.null(sensitivity)) {
  #if(sensitivity == "un_utilities") {
  #norm_disutilities <- rnorm(n_samples, mean = as.numeric(un_utilities[utility_row_index,15]),  sqrt((as.numeric(un_utilities[utility_row_index,18])^2+(as.numeric(un_utilities[utility_row_index,19])^2))))
  #}
  #}
  for(treatment_name in treatment_names) {
    event_disutilities[ ,treatment_name , "State Post THR <2 years"] <- norm_disutilities_revision[,treatment_name] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<2", treatment_name)])))* multi_disutilities/2
    event_disutilities[ ,treatment_name , "State Post THR >=2 years < 10 years"] <- norm_disutilities_revision[,treatment_name] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_2-10", treatment_name)])))* multi_disutilities/2
    event_disutilities[ ,treatment_name , "State Post THR >=10 years"] <- norm_disutilities_revision[,treatment_name] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))* multi_disutilities/2
    
  }
  
  
  event_disutilities[ , , "State Early revision"] <- norm_disutilities_rerevision * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))* multi_disutilities/2
  event_disutilities[ , , "State middle revision"] <- norm_disutilities_rerevision * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))* multi_disutilities/2
  event_disutilities[ , , "State late revision"] <- norm_disutilities_rerevision * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))* multi_disutilities/2
  event_disutilities[ , , "State second revision"] <- norm_disutilities_rerevision * (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))* multi_disutilities/2
  
  
  state_utilities[,"Cemented" , "State Post THR <2 years"] <- input_parameters[ ,"qalys_primary_Cemented"] + event_disutilities[,"Cemented" , "State Post THR <2 years"]
  state_utilities[,"Uncemented" , "State Post THR <2 years"] <- input_parameters[ ,"qalys_primary_Uncemented"] + event_disutilities[,"Uncemented" , "State Post THR <2 years"]
  state_utilities[,"Hybrid" , "State Post THR <2 years"] <- input_parameters[ ,"qalys_primary_Hybrid"] + event_disutilities[,"Hybrid" , "State Post THR <2 years"]
  
  state_utilities[,"Cemented" , "State Post THR >=2 years < 10 years"] = input_parameters[ ,"qalys_primary_Cemented"] + event_disutilities[, "Cemented", "State Post THR >=2 years < 10 years"]
  state_utilities[, "Uncemented", "State Post THR >=2 years < 10 years"] = input_parameters[ ,"qalys_primary_Uncemented"] + event_disutilities[, "Uncemented", "State Post THR >=2 years < 10 years"]
  state_utilities[,"Hybrid" , "State Post THR >=2 years < 10 years"] = input_parameters[ ,"qalys_primary_Hybrid"] + event_disutilities[, "Hybrid", "State Post THR >=2 years < 10 years"]
  
  state_utilities[, "Cemented", "State Post THR >=10 years"] = input_parameters[ ,"qalys_primary_Cemented"] + event_disutilities[, "Cemented", "State Post THR >=10 years"]
  state_utilities[, "Uncemented", "State Post THR >=10 years"] = input_parameters[ ,"qalys_primary_Uncemented"] + event_disutilities[, "Uncemented", "State Post THR >=10 years"]
  state_utilities[, "Hybrid", "State Post THR >=10 years"] = input_parameters[ ,"qalys_primary_Hybrid"] + event_disutilities[, "Hybrid", "State Post THR >=10 years"]
  
  state_utilities[, , "State Early revision"] = rep(input_parameters[ ,"qalys_State Early revision"], n= n_treatments) + event_disutilities[, , "State Early revision"]
  state_utilities[, , "State middle revision"] = rep(input_parameters[ ,"qalys_State middle revision"], n= n_treatments) + event_disutilities[, , "State middle revision"]
  state_utilities[, , "State late revision"] = rep(input_parameters[ ,"qalys_State late revision"], n= n_treatments) + event_disutilities[, , "State late revision"]
  state_utilities[, , "State second revision"] = rep(input_parameters[ ,"qalys_State second revision"], n= n_treatments) + event_disutilities[, , "State second revision"]
  state_utilities[, , "State Death"] = rep(input_parameters[ ,"qalys_State Death"], n= n_treatments)
  
  return(state_utilities)
} # End function
