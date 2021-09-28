#' Generate frequent items using the apriori algorithm
#'
#' @description
#' This function takes in a factorized Boolean matrix and generate frequent itemsets
#' that meet all the user provided criteria provided by the calling function.
#'
#' @param apriori_input_df An input factorized Boolean dataframe with multiple input and outcome variables
#' @param combo_length The length of the combinations specified by the user
#' @param support_threshold Minimum support value calculated based on the minimum absolute observed frequency threshold specified by the user
#' @param input_colname_list A list of column names that identify the input variables
#' @param confidence_threshold Minimum confidence threshold specified by the user
#' @param include_output_ind Specifies if the outcome variables must also be made part of the analysis using the algorithm
#' @param output_colname_list A list of column names that identify the outcome variables
#'
#' @return A list of frequent item sets that meet all the constraints supplied to the apriori algorithm
#'
#' @details
#' This is a function leveraged by few of the four main methods available to the users.
#'
#' @author
#' Vijay Kumar Pounraja
#'
#' @importFrom arules apriori
#' @importFrom arules DATAFRAME
#' @import arules
#' @importFrom stringr str_replace_all
#' @importFrom tidyr separate
#' @export

run_apriori_freqitems <- function(apriori_input_df, combo_length, support_threshold, input_colname_list,
                                  confidence_threshold = confidence_threshold, include_output_ind = include_output_ind,
                                  output_colname_list = output_colname_list) {

    trans_input <- as(apriori_input_df, "transactions")
    if (missing(include_output_ind) || include_output_ind == 'N') {
      all_apriori_freqitems <- arules::apriori(trans_input, parameter = list(supp = support_threshold, target = "frequent itemsets",
                                                                             minlen = combo_length, maxlen = combo_length, maxtime = 0),
                                               appearance = list(items=c(paste0(input_colname_list, "=1"))))
    } else if (include_output_ind == 'Y') {
      all_apriori_freqitems <- arules::apriori(trans_input, parameter = list(supp = support_threshold, target = "frequent itemsets",
                                                                             minlen = combo_length, maxlen = combo_length, maxtime = 0),
                                               appearance = list(items=c(paste0(input_colname_list, "=1"), paste0(output_colname_list, "=1"))))
      all_apriori_freqitems <- subset(all_apriori_freqitems, subset = items %in% c(paste0(output_colname_list, "=1")))
    }

    all_apriori_freqitems_df <- arules::DATAFRAME(all_apriori_freqitems)
    all_apriori_freqitems_df$items <- stringr::str_replace_all(all_apriori_freqitems_df$items, "[{}]", "")
    all_apriori_freqitems_df$items <- stringr::str_replace_all(all_apriori_freqitems_df$items, "=1", "")
    all_apriori_freqitems_df <- tidyr::separate(data = all_apriori_freqitems_df, col = items, into = c(paste0("Item_", 1:combo_length)), sep = ",")
    all_apriori_freqitems_df <- subset(all_apriori_freqitems_df, select = c(paste0('Item_', seq(1:combo_length)), 'count'))
    colnames(all_apriori_freqitems_df)[dim(all_apriori_freqitems_df)[2]] <- 'Obs_Count_Combo'

  return(all_apriori_freqitems_df)
}
