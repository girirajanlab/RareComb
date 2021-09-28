#' Perform successive left joins to fetch information about the constituent elements of the combinations
#'
#' @description
#' Fetching the frequency of multiple individual elements that make up the combinations of varying length
#' and hence varying variable names or to join two similar data frames using identical variable names
#' necessitates this function that supplements and joins data based on the length of the combinations.
#'
#' @param left_df The data frame with information about the combinations
#' @param right_df The data frame with information either about the combinations or their constituent elements
#' @param combo_length The length of the combinations specified by the user used to determine the number of successive joins to attempt
#' @param diff_colnames Indicator that specifies if the joins are to be made based on same or different column names
#'
#' @return An output dataframe with the results of the join operation
#'
#' @author
#' Vijay Kumar Pounraja
#'
#' @importFrom dplyr left_join
#' @importFrom stringr str_replace
#' @import magrittr
#' @export

custom_left_join <- function(left_df, right_df, combo_length = combo_length, diff_colnames = diff_colnames){
  if (diff_colnames == "Y") {
    join_string <- paste0("left_join(., right_df, by = c('Item_", 1:combo_length, "' = 'Item_1'))", collapse = " %>% ")
    join_string <- stringr::str_replace(join_string, "\\.", "left_df")
    join_results_df <- eval(parse(text=join_string))
  } else if (diff_colnames == "N") {
    join_string <- paste0("left_join(left_df, right_df, by = c(", paste0("'", 'Item_', (1:combo_length), "'", collapse = ","), "))")
    join_results_df <- eval(parse(text=join_string))
  }

  return(join_results_df)

}
