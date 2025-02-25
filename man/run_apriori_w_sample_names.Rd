% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_apriori_w_sample_names.R
\name{run_apriori_w_sample_names}
\alias{run_apriori_w_sample_names}
\title{Generate frequent items along with the names of supporting observations using the apriori algorithm}
\usage{
run_apriori_w_sample_names(
  apriori_input_df,
  combo_length,
  support_threshold,
  input_colname_list,
  input_sample_list,
  confidence_threshold = confidence_threshold,
  include_output_ind = include_output_ind,
  output_colname_list = output_colname_list
)
}
\arguments{
\item{apriori_input_df}{An input factorized Boolean dataframe with multiple input and outcome variables}

\item{combo_length}{The length of the combinations specified by the user}

\item{support_threshold}{Minimum support value calculated based on the minimum absolute observed frequency threshold specified by the user}

\item{input_colname_list}{A list of column names that identify the input variables}

\item{input_sample_list}{A list of row names that identify the samples/observations}

\item{confidence_threshold}{Minimum confidence threshold specified by the user}

\item{include_output_ind}{Specifies if the outcome variables must also be made part of the analysis using the algorithm}

\item{output_colname_list}{A list of column names that identify the outcome variables}
}
\value{
A list of frequent item sets that meet all the constraints supplied to the apriori algorithm
}
\description{
This function takes in a factorized Boolean matrix and generate frequent item sets
that meet all the user provided criteria provided by the calling function. This
function includes in it's output the identifiers of observations that support
each significant combination.
}
\details{
This is a function leveraged by few of the four main methods available to the users.
}
\author{
Vijay Kumar Pounraja
}
