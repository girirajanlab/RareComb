#' Compare the observed frequencies of combinations with their expected frequencies under the assumption of independence within a single group
#'
#' @description
#' This function takes a Boolean dataframe as input and compares the observed frequency of combinations
#' that meet the criteria specified by the users with their corresponding expectation derived
#' under the assumption of independence between the constituent elements of each combination
#'
#' @usage
#' compare_expected_vs_observed(boolean_input_df, combo_length, min_indv_threshold,
#'                              max_freq_threshold, input_format,
#'                              pval_filter_threshold, adj_pval_type)
#'
#' @param boolean_input_df An input Boolean dataframe with multiple input variables
#' @param combo_length The length of the combinations specified by the user
#' @param min_indv_threshold Minimum number of instances that support the combination
#' @param max_freq_threshold Maximum fraction of the cohort size that could support a combination (i.e., filter out highly frequent events)
#' @param pval_filter_threshold Optional | p-value cut-off to use for multiple testing adjustment (Default = 0.05)
#' @param input_format Optional | Naming convention used for input variables (Default = 'Input_')
#' @param adj_pval_type Optional | Type of multiple testing corrections to use (Default = 'BH'; Alternative option = 'bonferroni')
#'
#' @return A dataframe with the list of multiple-testing adjusted statistically significant combinations
#'         along with quantitative measures (frequencies, p-values etc) that support the findings.
#'
#' @examples
#'
#'    compare_expected_vs_observed(boolean_input_df, 2, 10, 0.25, 0.05,
#'                                 input_format = 'Input_', adj_pval_type = 'BH')
#'
#'
#' @author
#' Vijay Kumar Pounraja
#'
#' @importFrom methods as
#' @importFrom stats binom.test
#' @importFrom stats p.adjust
#' @importFrom arules apriori
#' @importFrom arules DATAFRAME
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_count
#' @importFrom stringr str_sort
#' @importFrom tidyr separate
#' @import pwr
#' @export

compare_expected_vs_observed <- function(boolean_input_df, combo_length, min_indv_threshold, max_freq_threshold, input_format = 'Input_',
                                         pval_filter_threshold = 0.05, adj_pval_type = 'BH') {

  ###################################################
  # Section 1: Format the input data and parameters #
  ###################################################
  # Identify all the input variables
  input_colname_list <- colnames(boolean_input_df)[grepl(paste0("^" ,input_format), colnames(boolean_input_df))]

  number_of_indv <- dim(boolean_input_df)[1]
  max_instances <- round(number_of_indv * max_freq_threshold)

  apriori_input_df <- boolean_input_df[,input_colname_list]
  apriori_input_df <- as.data.frame(sapply(apriori_input_df, as.numeric))
  apriori_input_df <- apriori_input_df[,((colSums(apriori_input_df) >= min_indv_threshold & colSums(apriori_input_df) < max_instances))]
  apriori_input_df <- apriori_input_df[,stringr::str_sort(colnames(apriori_input_df), numeric = TRUE)]
  apriori_input_df <- as.data.frame(sapply(apriori_input_df, as.factor))

  sel_input_colname_list <- colnames(apriori_input_df)[grepl(paste0("^" ,input_format), colnames(apriori_input_df))]

  #####################################################################################
  # Section 2: Calculate the frequencies for all the combinations based on user input #
  ###########################################################################################################################
  # APRIORI (Combo): Generate frequent itemset of a given size in which mutations/events and phenotypes co-occur among them #
  ###########################################################################################################################
  support_threshold = min_indv_threshold/dim(apriori_input_df)[1]
  include_output_ind <- "N"

  freqitems_df <- run_apriori_freqitems(apriori_input_df, combo_length, support_threshold, sel_input_colname_list)

  ##################################################################################
  # APRIORI (Individual): Generate frequencies of event for each individual entity #
  ##################################################################################
  support_threshold = min_indv_threshold/dim(apriori_input_df)[1]
  include_output_ind <- "N"

  freqitems_size1_df <- run_apriori_freqitems(apriori_input_df, 1 , support_threshold, sel_input_colname_list)

  diff_colnames <- "Y"
  all_freqitems_df <- custom_left_join(freqitems_df, freqitems_size1_df, combo_length = combo_length, diff_colnames = diff_colnames)
  colnames(all_freqitems_df)[(combo_length + 1):dim(all_freqitems_df)[2]] <- c("Obs_Count_Combo", paste0("Obs_Count_I", 1:combo_length))

  exp_prob_calc_string <- paste0("(all_freqitems_df$Obs_Count_I", 1:combo_length, "/number_of_indv)", collapse = " * ")
  all_freqitems_df$Exp_Prob_Combo <- eval(parse(text=exp_prob_calc_string))
  all_freqitems_df$Obs_Prob_Combo <- all_freqitems_df$Obs_Count_Combo/number_of_indv

  all_freqitems_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                           x = all_freqitems_df$Obs_Count_Combo, n = number_of_indv,
                                                           p = all_freqitems_df$Exp_Prob_Combo,
                                                           alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

  all_freqitems_df$Adj_Pval_bonf <- round(p.adjust(all_freqitems_df$pvalue_more , "bonferroni"), 3)
  all_freqitems_df$Adj_Pval_BH <- round(p.adjust(all_freqitems_df$pvalue_more , "BH"), 3)

  if (adj_pval_type == 'BH') {
    all_sig_freqitems_df <- subset(all_freqitems_df, all_freqitems_df[["Adj_Pval_BH"]] < pval_filter_threshold)
  } else if (adj_pval_type == 'bonferroni') {
    all_sig_freqitems_df <- subset(all_freqitems_df, all_freqitems_df[["Adj_Pval_bonf"]] < pval_filter_threshold)
  }

  return(all_sig_freqitems_df)

}

