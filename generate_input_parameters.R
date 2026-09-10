generate_input_parameters <- function(n_samples, treatment_names = treatment_names, 
                                      state_names = state_names,
                                      initial_age = initial_age,
                                      final_age = final_age,
                                      starting_age = starting_age,
                                      gender = gender,
                                      sensitivity = NULL) {
  
  n_treatments <- length(treatment_names)
  
  
  parameter_names <- c(paste0("log_rate_1st_revision_<2", treatment_names),
                       paste0("log_rate_1st_revision_2-10", treatment_names),
                       paste0("log_rate_1st_revision_>10", treatment_names),
                       "log_rate_2nd_revision_early", "log_rate_2nd_revision_middle", 
                       "log_rate_2nd_revision_late", "log_rate_higher_revision",
                       paste0("cost_primary_", treatment_names),
                       "cost_revision",
                       "cost_rerevision",
                       "qalys_primary_Cemented","qalys_primary_Uncemented","qalys_primary_Hybrid",
                       "qalys_State Early revision","qalys_State middle revision", "qalys_State late revision", 
                       "qalys_State second revision", "qalys_State Death")
  n_parameters <- length(parameter_names)
  input_parameters <- array(dim = c(n_samples, n_parameters), dimnames = list(NULL, parameter_names))
  
  
  # data: female 55-64 years old group
  lifetables <- read_excel(paste0(data_directory, "/KNIPS Main input data.xlsx"), sheet = "uk_lifetables")
  lograte_revision <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "revision_log_rate")
  primary_costs <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "primary_unsuccess_costs")
  revision_costs <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "revision_unsuccess_costs")
  rerevision_costs <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "rerevision_unsuccess_costs")
  
  utilities_primary <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_primary")
  utilities_revision <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_revision")
  un_utilities <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_unadjusted")
  log_rate_1st_revision <- read_excel(paste0(data_directory,"/", paste0(gender, "-", initial_age, "-", "rate.xlsx")))
  
  
  # Impute implants to be the average over all implants
  
  input_parameters[ , "log_rate_1st_revision_<2Cemented"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cemented"],
                                                                                               sd = ((UL1[treatment == "Cemented"]-LL1[treatment == "Cemented"])/(2*1.96))))    
  input_parameters[ , "log_rate_1st_revision_<2Uncemented"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Uncemented"],
                                                                                                 sd = ((UL1[treatment == "Uncemented"]-LL1[treatment == "Uncemented"])/(2*1.96))))                                                                                      
  input_parameters[ , "log_rate_1st_revision_<2Hybrid"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Hybrid"],
                                                                                             sd = ((UL1[treatment == "Hybrid"]-LL1[treatment == "Hybrid"])/(2*1.96))))                                                                                       
  
  
  
  input_parameters[ , "log_rate_1st_revision_2-10Cemented"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cemented"],
                                                                                                 sd = ((UL2[treatment == "Cemented"]-LL2[treatment == "Cemented"])/(2*1.96))))
  input_parameters[ , "log_rate_1st_revision_2-10Uncemented"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Uncemented"],
                                                                                                   sd = ((UL2[treatment == "Uncemented"]-LL2[treatment == "Uncemented"])/(2*1.96))))                                                                                       
  input_parameters[ , "log_rate_1st_revision_2-10Hybrid"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Hybrid"],
                                                                                               sd = ((UL2[treatment == "Hybrid"]-LL2[treatment == "Hybrid"])/(2*1.96))))                                                                                     
  
  
  input_parameters[ , "log_rate_1st_revision_>10Cemented"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cemented"],
                                                                                                sd = ((UL3[treatment == "Cemented"]-LL3[treatment == "Cemented"])/(2*1.96))))
  input_parameters[ , "log_rate_1st_revision_>10Uncemented"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Uncemented"],
                                                                                                  sd = ((UL3[treatment == "Uncemented"]-LL3[treatment == "Uncemented"])/(2*1.96))))                                                                                      
  input_parameters[ , "log_rate_1st_revision_>10Hybrid"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Hybrid"],
                                                                                              sd = ((UL3[treatment == "Hybrid"]-LL3[treatment == "Hybrid"])/(2*1.96))))                                                                                      
  
  
  
  input_parameters[ , "log_rate_2nd_revision_early"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "early_revision"],
                                                                                     sd = ((UL[parameter == "early_revision"]-LL[parameter == "early_revision"])/(2*1.96))))
  input_parameters[ , "log_rate_2nd_revision_middle"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "middle_revision"],
                                                                                      sd = ((UL[parameter == "middle_revision"]-LL[parameter == "middle_revision"])/(2*1.96))))
  input_parameters[ , "log_rate_2nd_revision_late"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "late_revision"],
                                                                                    sd = ((UL[parameter == "late_revision"]-LL[parameter == "late_revision"])/(2*1.96))))
  
  input_parameters[ , "log_rate_higher_revision"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "third_revision"],
                                                                                  sd = ((UL[parameter == "third_revision"]-LL[parameter == "third_revision"])/(2*1.96))))
  ##################################
  ####### costs############
  generate_lognormal <- function(mean_val, sd_val, n_samples) {
    mean_val <- max(mean_val, 0.01)
    sd_val <- max(sd_val, 0.01)
    
    log_mean <- log(mean_val^2 / sqrt(sd_val^2 + mean_val^2))
    log_sd <- sqrt(log(1 + (sd_val^2 / mean_val^2)))
    
    return(rlnorm(n_samples, meanlog = log_mean, sdlog = log_sd))
  }
  
  for(treatment_name in treatment_names) {
    row_index <- which(primary_costs[, "treatment"] == treatment_name & 
                         primary_costs[, "age"] == ini_age & 
                         primary_costs[, "gender"] == gender)
    
    mean_val <- as.numeric(primary_costs[row_index, "mean"])
    sd_val <- as.numeric(primary_costs[row_index, "SE"])
    
    input_parameters[, paste0("cost_primary_", treatment_name)] <- 
      generate_lognormal(mean_val, sd_val, n_samples)
  }
  

  
  
  ##########################################
  #############revision costs###############
  row_index <- which(revision_costs[, "age"] == ini_age & 
                       revision_costs[, "gender"] == gender)
  
  mean_val <- as.numeric(revision_costs[row_index, "mean"])
  sd_val <- as.numeric(revision_costs[row_index, "SE"])
  
  input_parameters[, "cost_revision"] <- 
    generate_lognormal(mean_val, sd_val, n_samples)
  
  
  
  row_index <- which(rerevision_costs[, "age"] == ini_age & 
                       rerevision_costs[, "gender"] == gender)
  
  mean_val <- as.numeric(rerevision_costs[row_index, "mean"])
  sd_val <- as.numeric(rerevision_costs[row_index, "SE"])
  
  input_parameters[, "cost_rerevision"] <- 
    generate_lognormal(mean_val, sd_val, n_samples)
  
  
  
  #################################
  ######### states qalys###########
  #########primary impants#########
  utilities_primary_row_index <- which(utilities_primary[, "treatment"] == "Cemented" &
                                         utilities_primary[, "age"] == ini_age &
                                         utilities_primary[, "gender"] == gender)
  
  input_parameters[ ,"qalys_primary_Cemented"] <- rnorm(n_samples, mean = as.numeric(utilities_primary[utilities_primary_row_index,7]), 
                                                        sd = (as.numeric(utilities_primary[utilities_primary_row_index,9])-as.numeric(utilities_primary[utilities_primary_row_index,8]))/(2*1.96))
  
  utilities_primary_row_index <- which(utilities_primary[, "treatment"] == "Uncemented" &
                                         utilities_primary[, "age"] == ini_age &
                                         utilities_primary[, "gender"] == gender)
  input_parameters[ ,"qalys_primary_Uncemented"] <- rnorm(n_samples, mean = as.numeric(utilities_primary[utilities_primary_row_index,7]), 
                                                          sd = (as.numeric(utilities_primary[utilities_primary_row_index,9])-as.numeric(utilities_primary[utilities_primary_row_index,8]))/(2*1.96))
  
  utilities_primary_row_index <- which(utilities_primary[, "treatment"] == "Hybrid" &
                                         utilities_primary[, "age"] == ini_age &
                                         utilities_primary[, "gender"] == gender)
  input_parameters[ ,"qalys_primary_Hybrid"] <- rnorm(n_samples, mean = as.numeric(utilities_primary[utilities_primary_row_index,7]), 
                                                      sd = (as.numeric(utilities_primary[utilities_primary_row_index,9])-as.numeric(utilities_primary[utilities_primary_row_index,8]))/(2*1.96))
  
  #########revision states#########
  utilities_revision_row_index <- which(utilities_revision[, "age"] == ini_age &
                                          utilities_revision[, "gender"] == gender)
  
  input_parameters[ ,"qalys_State Early revision"] <- input_parameters[ ,"qalys_State middle revision"] <-  input_parameters[ ,"qalys_State late revision"] <- 
    rep(rnorm(n_samples, mean = as.numeric(utilities_revision[utilities_revision_row_index,6]), sd = (as.numeric(utilities_revision[utilities_revision_row_index,8])-as.numeric(utilities_revision[utilities_revision_row_index,7]))/(2*1.96)), n= n_treatments)
  
  input_parameters[ ,"qalys_State second revision"] <- rep(rnorm(n_samples, mean = as.numeric(utilities_revision[utilities_revision_row_index,12]), 
                                                                 sd = (as.numeric(utilities_revision[utilities_revision_row_index,14])-as.numeric(utilities_revision[utilities_revision_row_index,13]))/(2*1.96)), n= n_treatments)
  
  input_parameters[ ,"qalys_State Death"]<- rep(rep(0, each = n_samples), n= n_treatments)
  
  
  return(input_parameters)
} # End function

