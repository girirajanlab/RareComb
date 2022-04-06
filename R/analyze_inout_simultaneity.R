#' Analyze relationships between rare events among multiple input and output variables
#'
#' @description
#' This function takes a Boolean dataframe as input and analyzes the
#' relationship between input and output variables for the combinations
#' that that include at least a single output variable andmeet all the
#' input criteria specified by the user.
#'
#' @usage
#' analyze_in_out_simultaneity(boolean_input_mult_df, combo_length, min_output_count,
#'                             max_output_count, min_indv_threshold, max_freq_threshold,
#'                             input_format, output_format, pval_filter_threshold,
#'                             adj_pval_type)
#'
#' @param boolean_input_mult_df An input Boolean dataframe with multiple input and outcome variables
#' @param combo_length The length of the combinations specified by the user
#' @param min_output_count Minimum number of output variables present in the combination
#' @param max_output_count Maximum number of output variables present in the combination
#' @param min_indv_threshold Minimum number of instances that support the combination
#' @param max_freq_threshold Maximum fraction of the cohort size that could support a combination (i.e., filter out highly frequent events)
#' @param input_format Optional | Naming convention used for input variables (Default = 'Input_')
#' @param output_format Optional | Naming convention used for output variables (Default = 'Output_')
#' @param pval_filter_threshold Optional | p-value cut-off to use to identify significant combinations (Default = 0.05)
#' @param adj_pval_type Optional | Type of multiple testing corrections to use (Default = 'BH'; Alternative option = 'bonferroni')
#'
#' @return A dataframe with the list of multiple-testing adjusted statistically significant combinations
#'         along with quantitative measures (frequencies, p-values etc) that support the findings.
#'
#' @examples
#'     analyze_in_out_simultaneity(boolean_input_mult_df, 3, 1, 2, 5, 0.25,
#'                                 input_format = 'Input_', output_format = 'Output_',
#'                                 pval_filter_threshold = 0.05, adj_pval_type = 'BH')
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

analyze_in_out_simultaneity <- function(boolean_input_mult_df, combo_length, min_output_count, max_output_count, min_indv_threshold,
                                        max_freq_threshold, input_format = 'Input_', output_format = 'Output_',
                                        pval_filter_threshold = 0.05, adj_pval_type = 'BH') {

  ###################################################
  # Section 1: Format the input data and parameters #
  ###################################################
  # Identify all the input and output variables
  input_colname_list <- colnames(boolean_input_mult_df)[grepl(paste0("^" ,input_format), colnames(boolean_input_mult_df))]
  output_colname_list <- colnames(boolean_input_mult_df)[grepl(paste0("^" ,output_format), colnames(boolean_input_mult_df))]

  number_of_indv <- dim(boolean_input_mult_df)[1]
  max_instances <- round(number_of_indv * max_freq_threshold)

  apriori_input_df <- boolean_input_mult_df[,input_colname_list]
  apriori_input_df <- as.data.frame(sapply(apriori_input_df, as.numeric))
  apriori_input_df <- apriori_input_df[,((colSums(apriori_input_df) >= min_indv_threshold & colSums(apriori_input_df) < max_instances))]
  apriori_input_df <- apriori_input_df[,str_sort(colnames(apriori_input_df), numeric = TRUE)]

  apriori_output_df <- boolean_input_mult_df[,output_colname_list]
  apriori_output_df <- as.data.frame(sapply(apriori_output_df, as.numeric))
  apriori_output_df <- apriori_output_df[,((colSums(apriori_output_df) >= min_indv_threshold & colSums(apriori_output_df) < max_instances))]
  apriori_output_df <- apriori_output_df[,str_sort(colnames(apriori_output_df), numeric = TRUE)]

  apriori_input_only_df <- apriori_input_df
  apriori_input_only_df <- as.data.frame(sapply(apriori_input_only_df, as.factor))

  apriori_input_df <- cbind(apriori_input_df, apriori_output_df)
  apriori_input_df <- as.data.frame(sapply(apriori_input_df, as.factor))

  sel_input_colname_list <- colnames(apriori_input_df)[grepl(paste0("^" ,input_format), colnames(apriori_input_df))]
  sel_output_colname_list <- colnames(apriori_output_df)[grepl(paste0("^" ,output_format), colnames(apriori_output_df))]

  #####################################################################################
  # Section 2: Calculate the frequencies for all the combinations based on user input #
  ###########################################################################################################################
  # APRIORI (Combo): Generate frequent itemset of a given size in which mutations/events and phenotypes co-occur among them #
  ###########################################################################################################################
  support_threshold = min_indv_threshold/dim(apriori_input_df)[1]

  in_out_freqitems_df <- run_apriori_rules_inout_simult(apriori_input_df, combo_length, support_threshold, sel_input_colname_list, output_colname_list = sel_output_colname_list)

  ##################################################################################
  # APRIORI (Individual): Generate frequencies of event for each individual entity #
  ##################################################################################
  support_threshold = min_indv_threshold/dim(apriori_input_df)[1]
  include_output_ind <- "N"

  in_out_freqitems_size1_df <- run_apriori_freqitems(apriori_input_df, 1 , support_threshold, c(sel_input_colname_list, sel_output_colname_list),
                                                      include_output_ind = include_output_ind)

  diff_colnames <- "Y"
  all_in_out_freqitems_df <- custom_left_join(in_out_freqitems_df, in_out_freqitems_size1_df, combo_length = combo_length, diff_colnames = diff_colnames)
  colnames(all_in_out_freqitems_df)[(combo_length + 1):dim(all_in_out_freqitems_df)[2]] <- c("Obs_Count_Combo", paste0("Case_Obs_Count_I", 1:combo_length))

  all_in_out_freqitems_df$Output_Count <- apply(apply(all_in_out_freqitems_df[,1:combo_length],1,stringr::str_count, output_format),2,sum)

  all_in_out_freqitems_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] <= max_output_count &
                                                             all_in_out_freqitems_df[["Output_Count"]] >= min_output_count)

  exp_prob_calc_string <- paste0("(all_in_out_freqitems_df$Case_Obs_Count_I", 1:combo_length, "/number_of_indv)", collapse = " * ")
  all_in_out_freqitems_df$Exp_Prob_Combo <- eval(parse(text=exp_prob_calc_string))
  all_in_out_freqitems_df$Obs_Prob_Combo <- all_in_out_freqitems_df$Obs_Count_Combo/number_of_indv

  ####################################################################################################
  # Section 3: Calculate the simultaneity rates exclusively for the input variables where applicable #
  ####################################################################################################
  if (combo_length == 2) {
    all_in_out_freqitems_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] == 1)
    all_in_out_freqitems_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                             x = all_in_out_freqitems_df$Obs_Count_Combo, n = number_of_indv,
                                                             p = all_in_out_freqitems_df$Exp_Prob_Combo,
                                                             alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))
  }


  if (combo_length == 3) {
    sel_in_out_freqitems_singleton_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] %in% c(0,2,3))
    sel_in_out_freqitems_pairs_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] == 1)

    singleton_count <- dim(sel_in_out_freqitems_singleton_df)[1]
    pairs_count <- dim(sel_in_out_freqitems_pairs_df)[1]

    freqitems_singleton_df <- data.frame()
    freqitems_pairs_df <- data.frame()

    # Process combinations for which input-only p-values does NOT have to be calculated
    if (singleton_count > 0) {
      freqitems_singleton_df <- sel_in_out_freqitems_singleton_df
      freqitems_singleton_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                         x = freqitems_singleton_df$Obs_Count_Combo, n = number_of_indv,
                                                                         p = freqitems_singleton_df$Exp_Prob_Combo,
                                                                         alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))
      freqitems_singleton_df$input_only_pvalue_more <- NA
    }

    # Process combinations for which input-only p-values have to be calculated
    if (pairs_count > 0) {
      unique_input_list <- unique(c(unique(sel_in_out_freqitems_pairs_df$Item_1), unique(sel_in_out_freqitems_pairs_df$Item_2)))
      apriori_input_sel_input_df <- apriori_input_only_df[,unique_input_list]
      apriori_input_sel_input_df <- apriori_input_sel_input_df[,str_sort(names(apriori_input_sel_input_df), numeric = TRUE)]
      cont_input_colname_list <- colnames(apriori_input_sel_input_df)

      support_threshold = min_indv_threshold/dim(apriori_input_only_df)[1]
      include_output_ind = "N"
      cont_freqitems_df <- run_apriori_freqitems(apriori_input_sel_input_df, 2 , support_threshold, cont_input_colname_list,
                                                         include_output_ind = include_output_ind)

      diff_colnames <- "Y"
      cont_combo_length <- 2
      in_out_cont_freqitems_df <- custom_left_join(cont_freqitems_df, in_out_freqitems_size1_df, combo_length = cont_combo_length, diff_colnames = diff_colnames)
      colnames(in_out_cont_freqitems_df)[(cont_combo_length + 1):dim(in_out_cont_freqitems_df)[2]] <- c("Cont_Obs_Count_Combo", paste0("Cont_Obs_Count_I", 1:cont_combo_length))

      in_out_w_cont_freqitems_df <- dplyr::left_join(sel_in_out_freqitems_pairs_df, in_out_cont_freqitems_df, by = c("Item_1" = "Item_1", "Item_2" = "Item_2"))

      exp_prob_calc_string <- paste0("(in_out_w_cont_freqitems_df$Cont_Obs_Count_I", 1:(combo_length - 1), "/number_of_indv)", collapse = " * ")
      in_out_w_cont_freqitems_df$Cont_Exp_Prob_Combo <- eval(parse(text=exp_prob_calc_string))
      in_out_w_cont_freqitems_df$Cont_Obs_Prob_Combo <- in_out_w_cont_freqitems_df$Cont_Obs_Count_Combo/number_of_indv

      in_out_w_cont_freqitems_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                     x = in_out_w_cont_freqitems_df$Obs_Count_Combo, n = number_of_indv,
                                                                     p = in_out_w_cont_freqitems_df$Exp_Prob_Combo,
                                                                     alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

      in_out_w_cont_freqitems_df$input_only_pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                               x = in_out_w_cont_freqitems_df$Cont_Obs_Count_Combo, n = number_of_indv,
                                                                               p = in_out_w_cont_freqitems_df$Cont_Exp_Prob_Combo,
                                                                               alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))
      freqitems_pairs_df <- in_out_w_cont_freqitems_df[,c(1:dim(sel_in_out_freqitems_pairs_df)[2], (dim(in_out_w_cont_freqitems_df)[2] - 1),
                                                                            dim(in_out_w_cont_freqitems_df)[2])]
    }

    # Merge the results from singletons and pairs
    all_in_out_freqitems_df <- dplyr::bind_rows(freqitems_singleton_df, freqitems_pairs_df)

  }


  if (combo_length == 4) {
    sel_in_out_freqitems_singleton_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] %in% c(0,3,4))
    sel_in_out_freqitems_pairs_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] == 2)
    sel_in_out_freqitems_triplets_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Output_Count"]] == 1)

    singleton_count <- dim(sel_in_out_freqitems_singleton_df)[1]
    pairs_count <- dim(sel_in_out_freqitems_pairs_df)[1]
    triplets_count <- dim(sel_in_out_freqitems_triplets_df)[1]

    freqitems_singleton_df <- data.frame()
    freqitems_pairs_df <- data.frame()
    freqitems_triplets_df <- data.frame()

    if (singleton_count > 0) {
      freqitems_singleton_df <- sel_in_out_freqitems_singleton_df
      freqitems_singleton_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                         x = freqitems_singleton_df$Obs_Count_Combo, n = number_of_indv,
                                                                         p = freqitems_singleton_df$Exp_Prob_Combo,
                                                                         alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))
      freqitems_singleton_df$input_only_pvalue_more <- NA
    }

    if (pairs_count > 0) {
      unique_input_list <- unique(c(unique(sel_in_out_freqitems_pairs_df$Item_1), unique(sel_in_out_freqitems_pairs_df$Item_2)))
      apriori_input_sel_in_df <- apriori_input_only_df[,unique_input_list]
      apriori_input_sel_in_df <- apriori_input_sel_in_df[,str_sort(names(apriori_input_sel_in_df), numeric = TRUE)]
      cont_input_colname_list <- colnames(apriori_input_sel_in_df)

      support_threshold = min_indv_threshold/dim(apriori_input_only_df)[1]
      include_output_ind = "N"
      cont_freqitems_df <- run_apriori_freqitems(apriori_input_sel_input_df, 2 , support_threshold, cont_input_colname_list,
                                                 include_output_ind = include_output_ind)

      diff_colnames <- "Y"
      cont_combo_length <- 2
      in_out_cont_freqitems_df <- custom_left_join(cont_freqitems_df, in_out_freqitems_size1_df, combo_length = cont_combo_length, diff_colnames = diff_colnames)
      colnames(in_out_cont_freqitems_df)[(cont_combo_length + 1):dim(in_out_cont_freqitems_df)[2]] <- c("Cont_Obs_Count_Combo", paste0("Cont_Obs_Count_I", 1:cont_combo_length))

      diff_colnames <- "N"
      in_out_w_cont_freqitems_df <- custom_left_join(sel_in_out_freqitems_pairs_df, in_out_cont_freqitems_df, combo_length = cont_combo_length, diff_colnames = diff_colnames)

      in_out_w_cont_freqitems_df$Cont_Exp_Prob_Combo <- (in_out_w_cont_freqitems_df$Cont_Obs_Count_I1/number_of_indv *
                                                           in_out_w_cont_freqitems_df$Cont_Obs_Count_I2/number_of_indv)
      in_out_w_cont_freqitems_df$Cont_Obs_Prob_Combo <- in_out_w_cont_freqitems_df$Cont_Obs_Count_Combo/number_of_indv

      in_out_w_cont_freqitems_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                     x = in_out_w_cont_freqitems_df$Obs_Count_Combo, n = number_of_indv,
                                                                     p = in_out_w_cont_freqitems_df$Exp_Prob_Combo,
                                                                     alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

      in_out_w_cont_freqitems_df$input_only_pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                               x = in_out_w_cont_freqitems_df$Cont_Obs_Count_Combo, n = number_of_indv,
                                                                               p = in_out_w_cont_freqitems_df$Cont_Exp_Prob_Combo,
                                                                               alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

      freqitems_pairs_df <- in_out_w_cont_freqitems_df[,c(1:dim(sel_in_out_freqitems_pairs_df)[2],(dim(in_out_w_cont_freqitems_df)[2] - 1),
                                                                            dim(in_out_w_cont_freqitems_df)[2])]
    }

    if (triplets_count > 0) {
      unique_input_list <- unique(c(unique(sel_in_out_freqitems_triplets_df$Item_1), unique(sel_in_out_freqitems_triplets_df$Item_2),
                                   unique(sel_in_out_freqitems_triplets_df$Item_3)))
      apriori_input_sel_in_df <- apriori_input_only_df[,unique_input_list]
      apriori_input_sel_in_df <- apriori_input_sel_in_df[,str_sort(names(apriori_input_sel_in_df), numeric = TRUE)]
      cont_input_colname_list <- colnames(apriori_input_sel_in_df)

      support_threshold = min_indv_threshold/dim(apriori_input_only_df)[1]
      include_output_ind = "N"
      cont_freqitems_df <- run_apriori_freqitems(apriori_input_sel_input_df, 3 , support_threshold, cont_input_colname_list,
                                                 include_output_ind = include_output_ind)

      diff_colnames <- "Y"
      cont_combo_length <- 3
      in_out_cont_freqitems_df <- custom_left_join(cont_freqitems_df, in_out_freqitems_size1_df, combo_length = cont_combo_length, diff_colnames = diff_colnames)
      colnames(in_out_cont_freqitems_df)[(cont_combo_length + 1):dim(in_out_cont_freqitems_df)[2]] <- c("Cont_Obs_Count_Combo", paste0("Cont_Obs_Count_I", 1:cont_combo_length))

      diff_colnames <- "N"
      in_out_w_cont_freqitems_df <- custom_left_join(sel_in_out_freqitems_triplets_df, in_out_cont_freqitems_df, combo_length = cont_combo_length, diff_colnames = diff_colnames)

      in_out_w_cont_freqitems_df$Cont_Exp_Prob_Combo <- (in_out_w_cont_freqitems_df$Cont_Obs_Count_I1/number_of_indv *
                                                           in_out_w_cont_freqitems_df$Cont_Obs_Count_I2/number_of_indv *
                                                           in_out_w_cont_freqitems_df$Cont_Obs_Count_I3/number_of_indv)
      in_out_w_cont_freqitems_df$Cont_Obs_Prob_Combo <- in_out_w_cont_freqitems_df$Cont_Obs_Count_Combo/number_of_indv

      in_out_w_cont_freqitems_df$pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                     x = in_out_w_cont_freqitems_df$Obs_Count_Combo, n = number_of_indv,
                                                                     p = in_out_w_cont_freqitems_df$Exp_Prob_Combo,
                                                                     alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

      in_out_w_cont_freqitems_df$input_only_pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                               x = in_out_w_cont_freqitems_df$Cont_Obs_Count_Combo, n = number_of_indv,
                                                                               p = in_out_w_cont_freqitems_df$Cont_Exp_Prob_Combo,
                                                                               alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

      freqitems_triplets_df <- in_out_w_cont_freqitems_df[,c(1:dim(sel_in_out_freqitems_triplets_df)[2],(dim(in_out_w_cont_freqitems_df)[2] - 1),
                                                                               dim(in_out_w_cont_freqitems_df)[2])]
    }

    all_in_out_freqitems_df <- dplyr::bind_rows(freqitems_singleton_df, freqitems_pairs_df, freqitems_triplets_df)

  }

  #all_in_out_freqitems_df <- subset(all_in_out_freqitems_df, is.na(input_only_pvalue_more) == TRUE |
  #                                    input_only_pvalue_more > 1/10^30)
  all_in_out_freqitems_df$Adj_Pval_bonf <- round(p.adjust(all_in_out_freqitems_df$pvalue_more , "bonferroni"), 3)
  all_in_out_freqitems_df$Adj_Pval_BH <- round(p.adjust(all_in_out_freqitems_df$pvalue_more , "BH"), 3)
  all_in_out_freqitems_df <- all_in_out_freqitems_df[order(-all_in_out_freqitems_df$Output_Count,
                                                           all_in_out_freqitems_df$pvalue_more),]

  all_sig_in_out_freqitems_df <- subset(all_in_out_freqitems_df, all_in_out_freqitems_df[["Adj_Pval_BH"]] < 0.01)

  return(all_sig_in_out_freqitems_df)

}
