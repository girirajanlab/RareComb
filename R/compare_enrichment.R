#' Compare the enrichment in combinations of input variables between the binary outcomes (case/control)
#'
#' @description
#' This function takes a Boolean dataframe as input and quantifies the enrichment in the observed frequency of combinations
#' that meet the criteria specified by the users compared to their corresponding expectation derived under the assumption
#' of independence between the constituent elements of each combination. The function then reports the multiple-testing
#' adjusted significant combinations in which enrichment is observed in cases but not in controls.
#'
#' @usage
#' compare_enrichment(boolean_input_df, combo_length, min_indv_threshold,
#'                    max_freq_threshold, input_format, output_format,
#'                    pval_filter_threshold, adj_pval_type, min_power_threshold,
#'                    sample_names_ind)
#'
#' @param boolean_input_df An input Boolean dataframe with multiple input and a single binary outcome variable
#' @param combo_length The length of the combinations specified by the user
#' @param min_indv_threshold Minimum number of instances that support the combination
#' @param max_freq_threshold Maximum fraction of the cohort size that could support a combination (i.e., filter out highly frequent events)
#' @param input_format Optional | Naming convention used for input variables (Default = 'Input_')
#' @param output_format Optional | Naming convention used for output variables (Default = 'Output_')
#' @param pval_filter_threshold Optional | p-value cut-off to use to identify significant combinations in cases (Default = 0.05)
#' @param adj_pval_type Optional | Type of multiple testing corrections to use (Default = 'BH'; Alternative option = 'bonferroni')
#' @param min_power_threshold Optional | Minimum statistical power (at 5\% sig.threshold) required for significant combinations to be returned in the results (Default = 0.7)
#' @param sample_names_ind Optional | Indicator to specify if the output should includes row names that support each significant combination (Default = 'N'; Alternative option = 'Y')
#'
#' @return A dataframe with the list of multiple-testing adjusted statistically significant combinations
#'         along with quantitative measures (frequencies, p-values etc) that support the findings.
#'
#' @examples
#'
#'
#'    compare_enrichment(boolean_input_df, 3, 5, 0.25, input_format = 'Input_',
#'                       output_format = 'Output_', adj_pval_type = 'bonferroni',
#'                       sample_names_ind = 'N')
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
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_count
#' @importFrom stringr str_sort
#' @importFrom tidyr separate
#' @import pwr
#' @import sqldf
#' @export

compare_enrichment <- function(boolean_input_df, combo_length, min_indv_threshold, max_freq_threshold,
                               input_format = 'Input_', output_format = 'Output_', pval_filter_threshold = 0.05,
                               adj_pval_type = 'BH', min_power_threshold = 0.7, sample_names_ind = 'N') {

  # Identify all the input and output variables
  input_colname_list <- colnames(boolean_input_df)[grepl(paste0("^" ,input_format), colnames(boolean_input_df))]
  output_column <- colnames(boolean_input_df)[grepl(paste0("^" ,output_format), colnames(boolean_input_df))]

  cases_sample_list <- subset(boolean_input_df, boolean_input_df[[output_column]] == 1, select = c("Sample_Name"))
  controls_sample_list <- subset(boolean_input_df, boolean_input_df[[output_column]] == 0, select = c("Sample_Name"))

  tmp_apriori_input_cases_df <- subset(boolean_input_df, boolean_input_df[[output_column]] == 1, select = input_colname_list)
  tmp_apriori_input_controls_df <- subset(boolean_input_df, boolean_input_df[[output_column]] == 0, select = input_colname_list)
  tmp_apriori_input_cases_df <- as.data.frame(sapply(tmp_apriori_input_cases_df, as.numeric))
  tmp_apriori_input_controls_df <- as.data.frame(sapply(tmp_apriori_input_controls_df, as.numeric))

  number_of_cases <- dim(tmp_apriori_input_cases_df)[1]
  max_instances <- round(number_of_cases * max_freq_threshold)

  apriori_input_cases_df <- tmp_apriori_input_cases_df[,((colSums(tmp_apriori_input_cases_df) >= min_indv_threshold
                                                          & colSums(tmp_apriori_input_controls_df) >= 1 &
                                                            colSums(tmp_apriori_input_cases_df) < max_instances))]
  apriori_input_cases_df <- apriori_input_cases_df[, stringr::str_sort(colnames(apriori_input_cases_df), numeric = TRUE)]
  apriori_input_cases_df <- as.data.frame(sapply(apriori_input_cases_df, factor))

  apriori_input_controls_df <- tmp_apriori_input_controls_df[,((colSums(tmp_apriori_input_cases_df) >= min_indv_threshold
                                                                & colSums(tmp_apriori_input_controls_df) >= 1 &
                                                                  colSums(tmp_apriori_input_cases_df) < max_instances))]
  apriori_input_controls_df <- apriori_input_controls_df[, stringr::str_sort(colnames(apriori_input_controls_df), numeric = TRUE)]
  apriori_input_controls_df <- as.data.frame(sapply(apriori_input_controls_df, factor))

  sel_input_colname_list <- colnames(apriori_input_cases_df)[grepl(paste0("^" ,input_format), colnames(apriori_input_cases_df))]

  ############################
  # CASES / SEVERE Phenotype #
  ############################################################################################################
  # APRIORI (Combo): Generate frequent itemset of a given size in which mutations/events co-occur among them #
  ############################################################################################################
  support_threshold = min_indv_threshold/dim(apriori_input_cases_df)[1]
  include_output_ind <- "N"

  case_freqitems_df <- run_apriori_freqitems(apriori_input_cases_df, combo_length , support_threshold, sel_input_colname_list,
                                                     include_output_ind = include_output_ind)
  colnames(case_freqitems_df)[dim(case_freqitems_df)[2]] <- 'Case_Obs_Count_Combo'

  unique_items_string <- paste0("unique(c(", paste0("as.character(case_freqitems_df$Item_", 1:combo_length, ")", collapse = ", "), "))")
  uniq_combo_items <- eval(parse(text=unique_items_string))

  print(paste0('Number of initial combinations identified for cases: ', dim(case_freqitems_df)[1]))
  print(paste0('Number of unique items in cases: ', length(uniq_combo_items)))

  ##################################################################################
  # APRIORI (Individual): Generate frequencies of event for each individual entity #
  ##################################################################################
  support_threshold = min_indv_threshold/dim(apriori_input_cases_df)[1]
  include_output_ind <- "N"

  case_freqitems_size1_df <- run_apriori_freqitems(apriori_input_cases_df, 1 , support_threshold, uniq_combo_items,
                                             include_output_ind = include_output_ind)
  colnames(case_freqitems_size1_df)[dim(case_freqitems_size1_df)[2]] <- 'Obs_Count_A'

  diff_colnames <- "Y"
  all_case_freqitems_df <- custom_left_join(case_freqitems_df, case_freqitems_size1_df, combo_length = combo_length, diff_colnames = diff_colnames)
  colnames(all_case_freqitems_df)[(dim(all_case_freqitems_df)[2] - (combo_length - 1)):dim(all_case_freqitems_df)[2]] <- paste0("Case_Obs_Count_I", 1:combo_length)

  exp_prob_calc_string <- paste0("(all_case_freqitems_df$Case_Obs_Count_I", 1:combo_length, "/number_of_cases)", collapse = " * ")
  all_case_freqitems_df$Case_Exp_Prob_Combo <- eval(parse(text=exp_prob_calc_string))

  all_case_freqitems_df$Case_Obs_Prob_Combo <- all_case_freqitems_df$Case_Obs_Count_Combo/number_of_cases
  all_case_freqitems_df$Case_pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                              x = all_case_freqitems_df$Case_Obs_Count_Combo, n = number_of_cases,
                                                              p = all_case_freqitems_df$Case_Exp_Prob_Combo,
                                                              alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

  #############################
  # CONTROLS / MILD Phenotype #
  ############################################################################################################
  # APRIORI (Combo): Generate frequent itemset of a given size in which mutations/events co-occur among them #
  ############################################################################################################
  apriori_input_controls_df <- apriori_input_controls_df[,uniq_combo_items]
  apriori_input_controls_df <- apriori_input_controls_df[,stringr::str_sort(colnames(apriori_input_controls_df), numeric = TRUE)]
  apriori_input_controls_df <- as.data.frame(sapply(apriori_input_controls_df, factor))
  number_of_controls <- dim(apriori_input_controls_df)[1]
  print(paste0('Number of controls: ', number_of_controls))

  support_threshold = 2/number_of_controls
  include_output_ind <- "N"

  cont_freqitems_df <- run_apriori_freqitems(apriori_input_controls_df, combo_length , support_threshold, uniq_combo_items,
                                             include_output_ind = include_output_ind)
  colnames(cont_freqitems_df)[dim(cont_freqitems_df)[2]] <- 'Temp_Obs_Count_Combo'

  print(paste0('Number of combinations with support of at least 2 in controls: ', dim(cont_freqitems_df)[1]))

  ##################################################################################
  # APRIORI (Individual): Generate frequencies of event for each individual entity #
  ##################################################################################
  support_threshold = 1/dim(apriori_input_controls_df)[1]
  include_output_ind <- "N"

  cont_freqitems_size1_df <- run_apriori_freqitems(apriori_input_controls_df, 1 , support_threshold, uniq_combo_items,
                                                   include_output_ind = include_output_ind)
  colnames(cont_freqitems_size1_df)[dim(cont_freqitems_size1_df)[2]] <- 'Obs_Count_A'

  diff_colnames <- "N"
  case_cont_freqitems_df <- custom_left_join(all_case_freqitems_df, cont_freqitems_df, combo_length = combo_length, diff_colnames = diff_colnames)
  case_cont_freqitems_df$Temp_Obs_Count_Combo[is.na(case_cont_freqitems_df$Temp_Obs_Count_Combo)] <- 0

  diff_colnames <- "Y"
  all_case_cont_freqitems_df <- custom_left_join(case_cont_freqitems_df, cont_freqitems_size1_df, combo_length = combo_length, diff_colnames = diff_colnames)
  colnames(all_case_cont_freqitems_df)[(dim(all_case_cont_freqitems_df)[2] - (combo_length - 1)):dim(all_case_cont_freqitems_df)[2]] <- paste0("Cont_Obs_Count_I", 1:combo_length)

  exp_prob_calc_string <- paste0("(all_case_cont_freqitems_df$Cont_Obs_Count_I", 1:combo_length, "/number_of_controls)", collapse = " * ")
  all_case_cont_freqitems_df$Cont_Exp_Prob_Combo <- eval(parse(text=exp_prob_calc_string))

  all_case_cont_freqitems_df$Temp_Control_pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                           x = all_case_cont_freqitems_df$Temp_Obs_Count_Combo,
                                                                           n = number_of_controls,
                                                                           p = all_case_cont_freqitems_df$Cont_Exp_Prob_Combo,
                                                                           alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

  sel_case_cont_freqitems_df <- all_case_cont_freqitems_df

  print(paste0('Number of combinations considered for multiple testing correction: ', dim(sel_case_cont_freqitems_df)[1]))

  sel_case_cont_freqitems_df$Case_Adj_Pval_bonf <- round(p.adjust(sel_case_cont_freqitems_df$Case_pvalue_more , "bonferroni"), 3)
  sel_case_cont_freqitems_df$Case_Adj_Pval_BH <- round(p.adjust(sel_case_cont_freqitems_df$Case_pvalue_more , "BH"), 3)


  if (adj_pval_type == 'BH') {
    all_sig_case_cont_freqitems_df <- subset(sel_case_cont_freqitems_df, sel_case_cont_freqitems_df[["Case_Adj_Pval_BH"]] < pval_filter_threshold &
                                                                         sel_case_cont_freqitems_df[["Temp_Control_pvalue_more"]] > pval_filter_threshold)
  } else if (adj_pval_type == 'bonferroni') {
    all_sig_case_cont_freqitems_df <- subset(sel_case_cont_freqitems_df, sel_case_cont_freqitems_df[["Case_Adj_Pval_bonf"]] < pval_filter_threshold &
                                                                         sel_case_cont_freqitems_df[["Temp_Control_pvalue_more"]] > pval_filter_threshold)
  }

  multtest_sig_comb_count <- dim(all_sig_case_cont_freqitems_df)[1]
  print(paste0('Number of combinations that are significant after multiple testing correction: ', multtest_sig_comb_count))


  #################################################################################################
  # Check if there is at least a single significant combination after multiple testing correction #
  #################################################################################################
  if (multtest_sig_comb_count > 0) {

    ###################################################
    # Check if zero frequency cases exist in controls #
    ###################################################
    cont_combos_w_zero_freq_df <- subset(all_sig_case_cont_freqitems_df, all_sig_case_cont_freqitems_df[["Temp_Obs_Count_Combo"]] == 0)
    zero_freq_combo_count <- dim(cont_combos_w_zero_freq_df)[1]

    print(paste0('Number of combinations with support less than 2 in controls: ', zero_freq_combo_count))

    if (zero_freq_combo_count > 0) {

      ####################################################################################
      # Select the items for which frequency of combinations in controls must be refined #
      ####################################################################################
      unique_items_string <- paste0("unique(c(", paste0("as.character(cont_combos_w_zero_freq_df$Item_", 1:combo_length, ")", collapse = ", "), "))")
      case_uniq_sig_items <- eval(parse(text=unique_items_string))

      ##############################
      # REFINE CONTROL FREQUENCIES #
      #######################################################################################################
      # APRIORI (Combo): Refine the frequency of events that cooccur either once or never in control groups #
      #######################################################################################################
      support_threshold = 0
      include_output_ind <- "N"
      refine_freqitems_df <- run_apriori_freqitems(apriori_input_controls_df, combo_length, support_threshold, case_uniq_sig_items,
                                                       include_output_ind = include_output_ind)
      colnames(refine_freqitems_df)[dim(refine_freqitems_df)[2]] <- 'Cont_Ref_Count_Combo'
      table(refine_freqitems_df$Cont_Ref_Count_Combo)
      sel_refine_freqitems_df <- subset(refine_freqitems_df, refine_freqitems_df[["Cont_Ref_Count_Combo"]] <= 1)

      diff_colnames <- "N"
      ref_sig_case_cont_freqitems_df <- custom_left_join(all_sig_case_cont_freqitems_df, sel_refine_freqitems_df, combo_length = combo_length, diff_colnames = diff_colnames)
      ref_sig_case_cont_freqitems_df$Cont_Ref_Count_Combo[is.na(ref_sig_case_cont_freqitems_df$Cont_Ref_Count_Combo)] <- 0
      ref_sig_case_cont_freqitems_df$Cont_Obs_Count_Combo <- ref_sig_case_cont_freqitems_df$Temp_Obs_Count_Combo + ref_sig_case_cont_freqitems_df$Cont_Ref_Count_Combo

      ref_sig_case_cont_freqitems_df$Control_pvalue_more <- as.numeric(mapply(function(x,p,n, ...){binom.test(x, n, p, ...)$p.value},
                                                                              x = ref_sig_case_cont_freqitems_df$Cont_Obs_Count_Combo,
                                                                              n = number_of_controls,
                                                                              p = ref_sig_case_cont_freqitems_df$Cont_Exp_Prob_Combo,
                                                                              alternative = "greater", conf.level = 0.95, SIMPLIFY = FALSE))

      cols_to_drop <- c("Temp_Obs_Count_Combo", "Cont_Ref_Count_Combo", "Temp_Control_pvalue_more")
      significant_case_cont_freqitems_df <- ref_sig_case_cont_freqitems_df[ , !(names(ref_sig_case_cont_freqitems_df) %in% cols_to_drop)]

    } else {
      colnames(all_sig_case_cont_freqitems_df)[colnames(all_sig_case_cont_freqitems_df) == 'Temp_Obs_Count_Combo'] <- 'Cont_Obs_Count_Combo'
      colnames(all_sig_case_cont_freqitems_df)[colnames(all_sig_case_cont_freqitems_df) == 'Temp_Control_pvalue_more'] <- 'Control_pvalue_more'
      significant_case_cont_freqitems_df <- all_sig_case_cont_freqitems_df
    }

    significant_case_cont_freqitems_df$Cont_Obs_Prob_Combo <- significant_case_cont_freqitems_df$Cont_Obs_Count_Combo/number_of_controls
    significant_case_cont_freqitems_df$Case_Exp_Count_Combo <- round((significant_case_cont_freqitems_df$Case_Exp_Prob_Combo * number_of_cases), 2)
    significant_case_cont_freqitems_df$Cont_Exp_Count_Combo <- round((significant_case_cont_freqitems_df$Cont_Exp_Prob_Combo * number_of_controls),2)
    significant_case_cont_freqitems_df$Effect_Size <- as.numeric(mapply(function(p1, p2, ...){ES.h(p1, p2, ...)},
                                                                        p1 = significant_case_cont_freqitems_df$Case_Obs_Count_Combo/number_of_cases,
                                                                        p2 = significant_case_cont_freqitems_df$Cont_Obs_Count_Combo/number_of_controls ))
    significant_case_cont_freqitems_df$Power_One_Pct <- round(as.numeric(mapply(function(e,n1,n2, ...){pwr.2p2n.test(e, n1, n2, ...)$power},
                                                                                e = significant_case_cont_freqitems_df$Effect_Size,
                                                                                n1 = number_of_cases, n2 = number_of_controls, sig.level = 0.01)), 3)
    significant_case_cont_freqitems_df$Power_Five_Pct <- round(as.numeric(mapply(function(e,n1,n2, ...){pwr.2p2n.test(e, n1, n2, ...)$power},
                                                                                 e = significant_case_cont_freqitems_df$Effect_Size,
                                                                                 n1 = number_of_cases, n2 = number_of_controls, sig.level = 0.05)), 3)

    select_cols_string <- paste0("significant_case_cont_freqitems_df[order(-significant_case_cont_freqitems_df$Power_Five_Pct, -significant_case_cont_freqitems_df$Effect_Size), c(",
                                 paste0("'Item_", 1:combo_length, "'", collapse = ","), ",", paste0("'Case_Obs_Count_I", 1:combo_length, "'", collapse = ","), ",",
                                 "'Case_Exp_Prob_Combo', 'Case_Obs_Prob_Combo','Case_Exp_Count_Combo', 'Case_Obs_Count_Combo', 'Case_pvalue_more'", ",",
                                 paste0("'Cont_Obs_Count_I", 1:combo_length, "'", collapse = ","), ",",
                                 "'Cont_Exp_Prob_Combo', 'Cont_Obs_Prob_Combo', 'Cont_Exp_Count_Combo','Cont_Obs_Count_Combo', 'Control_pvalue_more', 'Case_Adj_Pval_BH', 'Case_Adj_Pval_bonf', 'Effect_Size', 'Power_One_Pct', 'Power_Five_Pct'", ")]")
    output_sig_case_cont_freqitems_df <- eval(parse(text = select_cols_string))

    output_sig_case_cont_freqitems_df <- subset(output_sig_case_cont_freqitems_df, output_sig_case_cont_freqitems_df[["Power_Five_Pct"]] >= min_power_threshold)


    if (dim(output_sig_case_cont_freqitems_df)[1] > 0) {
      if (sample_names_ind == 'N') {
        # Return the output without adding sample names
        output_sig_case_cont_freqitems_df
      } else if (sample_names_ind == 'Y') {

        apriori_input_cases_w_samples_df <- cbind(apriori_input_cases_df, cases_sample_list)
        apriori_input_controls_w_samples_df <- cbind(apriori_input_controls_df, controls_sample_list)

        input_by_sample_cases_df <- reshape2::melt(apriori_input_cases_w_samples_df, id = "Sample_Name")
        input_by_sample_cases_df <- input_by_sample_cases_df[input_by_sample_cases_df$value > 0, !(names(input_by_sample_cases_df) %in% c("value"))]
        input_by_sample_cases_df$Sample_Type <- 'Case_Samples'

        input_by_sample_controls_df <- reshape2::melt(apriori_input_controls_w_samples_df, id = "Sample_Name")
        input_by_sample_controls_df <- input_by_sample_controls_df[input_by_sample_controls_df$value > 0, !(names(input_by_sample_controls_df) %in% c("value"))]
        input_by_sample_controls_df$Sample_Type <- 'Control_Samples'

        input_by_sample_df <- rbind(input_by_sample_cases_df, input_by_sample_controls_df)
        input_by_sample_df[] <- lapply(input_by_sample_df, as.character)

        if (combo_length == 2) {
           sig_combos_df <- subset(output_sig_case_cont_freqitems_df, select = c("Item_1", "Item_2"))
           find_samples_step1_df <- dplyr::left_join(sig_combos_df, input_by_sample_df, by = c("Item_1" = "variable"))
           find_samples_step2_df <- dplyr::inner_join(find_samples_step1_df, input_by_sample_df,
                                                      by = c("Item_2" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
           find_samples_step3_df <- sqldf::sqldf("select Item_1, Item_2, Sample_Type, GROUP_CONCAT(Sample_Name) as Sample_List
                                          from find_samples_step2_df group by Item_1, Item_2, Sample_Type")
           find_samples_step4_df <- reshape2::dcast(find_samples_step3_df, Item_1 + Item_2 ~ Sample_Type, value.var = "Sample_List")
           out_sig_case_cont_freqitems_w_samples_df <- dplyr::left_join(output_sig_case_cont_freqitems_df, find_samples_step4_df,
                                                                        by = c("Item_1", "Item_2"))
        }

        if (combo_length == 3) {
          sig_combos_df <- subset(output_sig_case_cont_freqitems_df, select = c("Item_1", "Item_2", "Item_3"))
          find_samples_step1_df <- dplyr::left_join(sig_combos_df, input_by_sample_df, by = c("Item_1" = "variable"))
          find_samples_step2_df <- dplyr::inner_join(find_samples_step1_df, input_by_sample_df,
                                                     by = c("Item_2" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step3_df <- dplyr::inner_join(find_samples_step2_df, input_by_sample_df,
                                                     by = c("Item_3" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step4_df <- sqldf::sqldf("select Item_1, Item_2, Item_3, Sample_Type, GROUP_CONCAT(Sample_Name) as Sample_List
                                          from find_samples_step3_df group by Item_1, Item_2, Item_3, Sample_Type")
          find_samples_step5_df <- reshape2::dcast(find_samples_step4_df, Item_1 + Item_2 + Item_3 ~ Sample_Type, value.var = "Sample_List")
          out_sig_case_cont_freqitems_w_samples_df <- dplyr::left_join(output_sig_case_cont_freqitems_df, find_samples_step5_df,
                                                                       by = c("Item_1", "Item_2", "Item_3"))
        }

        if (combo_length == 4) {
          sig_combos_df <- subset(output_sig_case_cont_freqitems_df, select = c("Item_1", "Item_2", "Item_3", "Item_4"))
          find_samples_step1_df <- dplyr::left_join(sig_combos_df, input_by_sample_df, by = c("Item_1" = "variable"))
          find_samples_step2_df <- dplyr::inner_join(find_samples_step1_df, input_by_sample_df,
                                                     by = c("Item_2" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step3_df <- dplyr::inner_join(find_samples_step2_df, input_by_sample_df,
                                                     by = c("Item_3" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step4_df <- dplyr::inner_join(find_samples_step3_df, input_by_sample_df,
                                                     by = c("Item_4" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step5_df <- sqldf::sqldf("select Item_1, Item_2, Item_3, Item_4, Sample_Type, GROUP_CONCAT(Sample_Name) as Sample_List
                                          from find_samples_step4_df group by Item_1, Item_2, Item_3, Item_4, Sample_Type")
          find_samples_step6_df <- reshape2::dcast(find_samples_step5_df, Item_1 + Item_2 + Item_3 + Item_4 ~ Sample_Type, value.var = "Sample_List")
          out_sig_case_cont_freqitems_w_samples_df <- dplyr::left_join(output_sig_case_cont_freqitems_df, find_samples_step6_df,
                                                                       by = c("Item_1", "Item_2", "Item_3", "Item_4"))
        }

        if (combo_length == 5) {
          sig_combos_df <- subset(output_sig_case_cont_freqitems_df, select = c("Item_1", "Item_2", "Item_3", "Item_4", "Item_5"))
          find_samples_step1_df <- dplyr::left_join(sig_combos_df, input_by_sample_df, by = c("Item_1" = "variable"))
          find_samples_step2_df <- dplyr::inner_join(find_samples_step1_df, input_by_sample_df,
                                                     by = c("Item_2" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step3_df <- dplyr::inner_join(find_samples_step2_df, input_by_sample_df,
                                                     by = c("Item_3" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step4_df <- dplyr::inner_join(find_samples_step3_df, input_by_sample_df,
                                                     by = c("Item_4" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step5_df <- dplyr::inner_join(find_samples_step4_df, input_by_sample_df,
                                                     by = c("Item_5" = "variable", "Sample_Name" = "Sample_Name", "Sample_Type" = "Sample_Type"))
          find_samples_step6_df <- sqldf::sqldf("select Item_1, Item_2, Item_3, Item_4, Item_5, Sample_Type, GROUP_CONCAT(Sample_Name) as Sample_List
                                          from find_samples_step5_df group by Item_1, Item_2, Item_3, Item_4, Item_5, Sample_Type")
          find_samples_step7_df <- reshape2::dcast(find_samples_step6_df, Item_1 + Item_2 + Item_3 + Item_4 + Item_5 ~ Sample_Type, value.var = "Sample_List")
          out_sig_case_cont_freqitems_w_samples_df <- dplyr::left_join(output_sig_case_cont_freqitems_df, find_samples_step7_df,
                                                                       by = c("Item_1", "Item_2", "Item_3", "Item_4", "Item_5"))
        }

        # Return all significant combinations with corresponding sample names as the output
        out_sig_case_cont_freqitems_w_samples_df
    }


   } else {
     warning("No significant combinations that meet the specified power threshold")
     warning("Returning ONLY the non-significant combinations")
     # Return all non-significant combinations as the output
     sel_case_cont_freqitems_df
   }


  } else {
    warning("No significant combinations were found after multiple testing correction")
    warning("Returning ONLY the non-significant combinations")
    # Return all non-significant combinations as the output
    sel_case_cont_freqitems_df
  } # End if (multtest_sig_comb_count > 0)


} # End function
