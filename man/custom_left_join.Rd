% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_left_join.R
\name{custom_left_join}
\alias{custom_left_join}
\title{Perform successive left joins to fetch information about the constituent elements of the combinations}
\usage{
custom_left_join(
  left_df,
  right_df,
  combo_length = combo_length,
  diff_colnames = diff_colnames
)
}
\arguments{
\item{left_df}{The data frame with information about the combinations}

\item{right_df}{The data frame with information either about the combinations or their constituent elements}

\item{combo_length}{The length of the combinations specified by the user used to determine the number of successive joins to attempt}

\item{diff_colnames}{Indicator that specifies if the joins are to be made based on same or different column names}
}
\value{
An output dataframe with the results of the join operation
}
\description{
Fetching the frequency of multiple individual elements that make up the combinations of varying length
and hence varying variable names or to join two similar data frames using identical variable names
necessitates this function that supplements and joins data based on the length of the combinations.
}
\author{
Vijay Kumar Pounraja
}
