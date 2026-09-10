# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Script conduct cost-effectiveness analysis of atrial fibrillation using example Markov model

# Load BCEA library to help analyse and visualise results

library(BCEA)
library(readxl)
library(ggplot2)
library(reshape)
library(dplyr)
library(RColorBrewer)

set.seed(14142234)


source('code/generate_input_parameters.R')
source('code/generate_transition_matrices.R')
source('code/convert_transition_matrices_to_df.R')
source('code/generate_state_qalys.R')
source('code/generate_state_costs.R')
source('code/generate_net_benefit_lapply.R')


data_directory = "C:/Users/yx18392/Desktop/HIPPY/semi Markov data"


# Define global simulation parameters
n_samples <- 1000

# Define global model structure parameters
n_states <- 8
state_names <- paste("State", c("Post THR <2 years", "Post THR >=2 years < 10 years", "Post THR >=10 years", 
                                "Early revision",  "middle revision", "late revision", "second revision", "Death"))

treatment_names <- c("Cemented", "Uncemented", "Hybrid")
n_treatments <- length(treatment_names)

event_names <- c("early primary", "middle primary", "late primary", "early revision", "middle revision", "late revision", "reresion")


#age_range <- "0-69"
age_range <- "70-100"

initial_age <- substring(age_range, 1, 2)
final_age <- substring(age_range, 4, 6)
if (final_age == "9") {final_age <- 69}

if (initial_age == "0-") {ini_age <- 60}else{
  ini_age <- as.numeric(initial_age)}


if(initial_age == "0-") {starting_age <- 60 }else{
  tarting_age <- 70}
# Specify the gender
gender <- "female"

# Define global scenario parameters
  n_cycles <- 100- as.numeric(ini_age)

# Generate the input parameters
# This will be converted into transition matrix, state costs, and state utilities
input_parameters <- generate_input_parameters(n_samples,
                                              treatment_names = treatment_names, 
                                              state_names = state_names,
                                              initial_age = initial_age,
                                              final_age = final_age,
                                              starting_age = starting_age,
                                              gender = gender)

# Run the Markov model to get the model outputs
model_outputs <- generate_net_benefit(input_parameters,
                                      treatment_names = treatment_names, 
                                      state_names = state_names,
                                      initial_age = initial_age,
                                      final_age = final_age,
                                      starting_age = starting_age,
                                      gender = gender,
                                      n_cycles = n_cycles,
                                      lambda = 20000)

####################################################################################################
## Quick manual check of resutls ###################################################################
####################################################################################################


# Expected net benefit at ?20,000
#with(model_outputs, colMeans(20000 * total_qalys - total_costs))

#rowMeans(model_outputs$incremental_net_benefit)
##################################################################################################################
## Now use BCEA to analyse the outputs   #########################################################################
##################################################################################################################
# Create a bcea object for the knips model
# Note costs and QALYs need to be transposed in this example for BCEA to run
knips_bcea <- bcea(e = t(total_qalys), c = t(total_costs), ref = 1, interventions = treatment_names) 
# Get summary statistics
#summary(knips_bcea, wtp = 20000)

# Plot a CEAC
setwd('C:/Users/yx18392/Desktop')

my_colors <- c("red", "blue", "green")  # 3个治疗方案用3种颜色

png(file=paste0(gender,"-", age_range,"_multice.png"))
knips_multi_ce <- multi.ce(knips_bcea)
ceac.plot(knips_multi_ce, 
          graph = "ggplot2",
          line = list(color = my_colors),  # 直接指定颜色向量
          pos = c(1,1))
dev.off()


dev.off()

format_results <- function(x, n_digits = 2) {
  paste0(format(mean(x), nsmall = n_digits, digits = n_digits), " (",
         format(quantile(x, prob = 0.025), nsmall = n_digits, digits = n_digits), " ,",
         format(quantile(x, prob = 0.975), nsmall = n_digits, digits = n_digits), ")")
}

results_matrix <- matrix(nrow = n_treatments, ncol = 4)
rownames(results_matrix) <- treatment_names
for(treatment in treatment_names) {
  results_matrix[treatment, ] <- c(format_results(total_qalys[treatment,]),
                                   format_results(total_costs[treatment,]),
                                   format_results(net_benefit[treatment,]),
                                   format_results(incremental_net_benefit[treatment,]))
}

format_results <- function(x, n_digits = 0) {  # 改为0位小数
  paste0(format(round(mean(x), n_digits), nsmall = n_digits, scientific = FALSE), " (",
         format(round(quantile(x, prob = 0.025), n_digits), nsmall = n_digits, scientific = FALSE), " ,",
         format(round(quantile(x, prob = 0.975), n_digits), nsmall = n_digits, scientific = FALSE), ")")
}

results_matrix <- matrix(nrow = n_treatments, ncol = 4)
rownames(results_matrix) <- treatment_names
for(treatment in treatment_names) {
  results_matrix[treatment, ] <- c(format_results(total_qalys[treatment,], n_digits = 2),
                                   format_results(total_costs[treatment,], n_digits = 0),
                                   format_results(net_benefit[treatment,], n_digits = 0),
                                   format_results(incremental_net_benefit[treatment,], n_digits = 0))
}
results_matrix = cbind(treatment_names,results_matrix, round(ICER,2))
colnames(results_matrix)  = c("Treatment_names", "QALYs", "Costs", "Net_benefit", "Inc.Net_benefit", "ICER")

for (treatment in treatment_names){
  if(results_matrix[treatment,6]<0) {
    if (mean(total_qalys[treatment,])-mean(total_qalys[1,])>0){results_matrix[treatment,6] = "Dominant"}
    if (mean(total_qalys[treatment,])-mean(total_qalys[1,])<0){results_matrix[treatment,6] = "Dominated"} 
  }
  
}

icer_table = data.frame(results_matrix)
View(icer_table)
write.csv(icer_table, file = paste0(gender,"-", age_range,"icer_results.csv"))


##################################################################################
## VoI analysis - takes a few minutes ############################################
##################################################################################

evpi_table <- matrix(nrow = 11, ncol = 2)
rownames(evpi_table) <- c("Total", "1st revision probabilities_<2", "1st revision probabilities_2-10",
                          "1st revision probabilities_>10","2nd and higher revision probabilities",
                          "primary utilities", "revision utilities","successful revision costs","unsuccessful revision costs",
                          "successful primary costs","unsuccessful primary costs")
colnames(evpi_table) <- c("Per person", "Population")

# Assume implant decision remains relevant for at least 10 years
# And that there are 160000 THR per year
# https://www.njrcentre.org.uk/njrcentre/Patients/Joint-replacement-statistics#:~:text=In%20England%20and%20Wales%20there,the%20practice%20is%20growing%20rapidly.
# This is regardless of age and gender but assume EVPI same for each category
technology_horizon <- 10
discounted_population_size <- sum((1/1.035)^(0:(technology_horizon - 1))) * 0.249142143355744 * 116089 

# Total EVPI
evpi_table["Total", c("Per person", "Population")] <-  knips_bcea$evi[201] * c(1, discounted_population_size)

colnames(input_parameters) <- gsub("<", "lt", colnames(input_parameters))
colnames(input_parameters) <- gsub(">", "gt", colnames(input_parameters))
colnames(input_parameters) <- gsub("-", "_", colnames(input_parameters))
# Log rate of first revision (same EVPPI as if calculating probability)
# Use GP instead of GAM as 4 parameters
evppi_gp_1st_revision_1 <- evppi(param_idx = c("log_rate_1st_revision_lt2Cemented",
                                   "log_rate_1st_revision_lt2Uncemented",
                                   "log_rate_1st_revision_lt2Hybrid"),
                                 input = input_parameters, he = knips_bcea, method = 'gp')
evpi_table["1st revision probabilities_<2", c("Per person", "Population")] <- evppi_gp_1st_revision_1$evppi[201] * c(1, discounted_population_size)

evppi_gp_1st_revision_2 <- evppi(param_idx = c("log_rate_1st_revision_2_10Cemented",
                                                "log_rate_1st_revision_2_10Uncemented",
                                                "log_rate_1st_revision_2_10Hybrid"),
                                  input = input_parameters, he = knips_bcea, method = 'gp')
evpi_table["1st revision probabilities_2-10", c("Per person", "Population")] <- evppi_gp_1st_revision_2$evppi[201] * c(1, discounted_population_size)

evppi_gp_1st_revision_3 <- evppi(param_idx = c("log_rate_1st_revision_gt10Cemented",
                                               "log_rate_1st_revision_gt10Uncemented",
                                               "log_rate_1st_revision_gt10Hybrid"),
                                 input = input_parameters, he = knips_bcea, method = 'gp')
evpi_table["1st revision probabilities_>10", c("Per person", "Population")] <- evppi_gp_1st_revision_3$evppi[201] * c(1, discounted_population_size)


# Log rate 2nd and higher revision probabilities
evppi_gam_2nd_and_higher_revision <- evppi(param_idx =  
                                             c("log_rate_2nd_revision_early", "log_rate_2nd_revision_middle", 
                                               "log_rate_2nd_revision_late", "log_rate_higher_revision"), 
                                           input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["2nd and higher revision probabilities", c("Per person", "Population")] <- evppi_gam_2nd_and_higher_revision$evppi[201] * c(1, discounted_population_size)

# Utilities
evppi_gam_primary_utilities <- evppi(param_idx =  
                               c( "qalys_primary_Cemented", 
                                  "qalys_primary_Uncemented", 
                                  "qalys_primary_Hybrid"), 
                             input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["primary utilities", c("Per person", "Population")] <- evppi_gam_primary_utilities$evppi[201] * c(1, discounted_population_size)

colnames(input_parameters) <- make.names(colnames(input_parameters))
evppi_gam_revision_utilities <- evppi(param_idx =  
                               c( "qalys_State.Early.revision", 
                                  "qalys_State.middle.revision", 
                                  "qalys_State.late.revision", 
                                  "qalys_State.second.revision"), 
                             input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["revision utilities", c("Per person", "Population")] <- evppi_gam_revision_utilities$evppi[201] * c(1, discounted_population_size)

evppi_gam_suc_revision_costs <- evppi(param_idx =  
                                   c( "cost_revision_success",
                                      "cost_rerevision_success"), 
                                 input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["successful revision costs", c("Per person", "Population")] <- evppi_gam_suc_revision_costs$evppi[201] * c(1, discounted_population_size)

evppi_gam_unsuc_revision_costs <- evppi(param_idx =  
                                        c( "cost_revision_unsuccess",
                                           "cost_rerevision_unsuccess"), 
                                      input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["unsuccessful revision costs", c("Per person", "Population")] <- evppi_gam_unsuc_revision_costs$evppi[201] * c(1, discounted_population_size)

evppi_gam_primary_costs <- evppi(param_idx =  
                                        c( "cost_primary_success_Cemented", 
                                           "cost_primary_success_Uncemented",
                                           "cost_primary_success_Hybrid"), 
                                      input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["successful primary costs", c("Per person", "Population")] <- evppi_gam_primary_costs$evppi[201] * c(1, discounted_population_size)

evppi_gam_revision_costs <- evppi(param_idx =  
                                   c( "cost_primary_unsuccess_Cemented", 
                                      "cost_primary_unsuccess_Uncemented",
                                      "cost_primary_unsuccess_Hybrid"), 
                                 input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["unsuccessful primary costs", c("Per person", "Population")] <- evppi_gam_revision_costs$evppi[201] * c(1, discounted_population_size)


write.csv(evpi_table, file = paste0(gender,"-", age_range,"evppi.csv"))

View(evpi_table)

