#' Generate rules using the apriori algorithm
#'
#' @description
#' This function takes in a factorized Boolean matrix and generate rules that meet
#' all the user provided criteria while allowing the outcome variables to be
#' part of either LHS or RHS of the rules but restricting the input variables
#' to the LHS of the rules.
#'
#' @param apriori_input_df An input factorized Boolean dataframe with multiple input and outcome variables
#' @param combo_length The length of the combinations specified by the user
#' @param support_threshold Minimum support value calculated based on the minimum absolute observed frequency threshold specified by the user
#' @param input_colname_list A list of column names that identify the input variables
#' @param output_colname_list Optional | A list of column names that identify the outcome variables
#'
#' @return A list of rules that meet all the constraints supplied to the apriori algorithm
#'
#' @details
#' This is a function leveraged by few of the four main methods available to the users.
#'
#' @author
#' Vijay Kumar Pounraja
#'
#' @importFrom arules apriori
#' @importFrom arules DATAFRAME
#' @importFrom stringr str_replace_all
#' @importFrom tidyr separate
#' @export

run_apriori_rules_inout_simult <- function(apriori_input_df, combo_length, support_threshold, input_colname_list, output_colname_list = output_colname_list) {

  lhs_length <- combo_length - 1
  trans_input <- as(apriori_input_df, "transactions")
  all_apriori_freqitems <- arules::apriori(trans_input, parameter = list(supp = support_threshold, target = "rules", confidence = 0,
                                                                         minlen = combo_length, maxlen = combo_length, maxtime = 0),
                                             appearance = list(lhs=paste0(input_colname_list, "=1"), both=paste0(output_colname_list, "=1")))

  all_apriori_freqitems_df <- arules::DATAFRAME(all_apriori_freqitems)
  all_apriori_freqitems_df$LHS <- stringr::str_replace_all(all_apriori_freqitems_df$LHS, "[{}]", "")
  all_apriori_freqitems_df$LHS <- stringr::str_replace_all(all_apriori_freqitems_df$LHS, "=1", "")
  all_apriori_freqitems_df$RHS <- stringr::str_replace_all(all_apriori_freqitems_df$RHS, "[{}]", "")
  all_apriori_freqitems_df$RHS <- stringr::str_replace_all(all_apriori_freqitems_df$RHS, "=1", "")
  all_apriori_freqitems_df <- tidyr::separate(data = all_apriori_freqitems_df, col = "LHS", into = c(paste0("Item_", 1:lhs_length)), sep = ",")
  colnames(all_apriori_freqitems_df)[combo_length] <- paste0('Item_', combo_length)
  all_apriori_freqitems_df <- subset(all_apriori_freqitems_df, select = c(paste0('Item_', seq(1:combo_length)), 'count'))

  return(all_apriori_freqitems_df)
}


