% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_enrichment_modifiers.R
\name{compare_enrichment_modifiers}
\alias{compare_enrichment_modifiers}
\title{Compare the enrichment in combinations of input variables between the binary outcomes (case/control)}
\usage{
compare_enrichment_modifiers(boolean_input_df, combo_length, min_indv_threshold,
                             max_freq_threshold, primary_input_entities, input_format,
                             output_format, pval_filter_threshold, adj_pval_type,
                             min_power_threshold, sample_names_ind)
}
\arguments{
\item{boolean_input_df}{An input Boolean dataframe with multiple input and a single binary outcome variable}

\item{combo_length}{The length of the combinations specified by the user}

\item{min_indv_threshold}{Minimum number of instances that support the combination}

\item{max_freq_threshold}{Maximum fraction of the cohort size that could support a combination (i.e., filter out highly frequent events)}

\item{primary_input_entities}{List of variables that MUST be part of the combinations identified by the method}

\item{input_format}{Optional | Naming convention used for input variables (Default = 'Input_')}

\item{output_format}{Optional | Naming convention used for output variables (Default = 'Output_')}

\item{pval_filter_threshold}{Optional | p-value cut-off to use to identify significant combinations in cases (Default = 0.05)}

\item{adj_pval_type}{Optional | Type of multiple testing corrections to use (Default = 'BH'; Alternative option = 'bonferroni')}

\item{min_power_threshold}{Optional | Minimum statistical power (at 5\% sig.threshold) required for significant combinations to be returned in the results (Default = 0.7)}

\item{sample_names_ind}{Optional | Indicator to specify if the output should includes row names that support each significant combination (Default = 'N'; Alternative option = 'Y')}
}
\value{
A dataframe with the list of multiple-testing adjusted statistically significant combinations
        along with quantitative measures (frequencies, p-values etc) that support the findings.
}
\description{
This function takes a Boolean dataframe as input and quantifies the enrichment in the observed frequency of combinations
that include at least one of the input variables supplied by the user as well as meet other user-specified criteria
compared to their corresponding expectation derived under the assumption of independence between the constituent elements
of each combination. The function then reports the combinations in which enrichment is observed in cases but not in controls.
}
\examples{


   compare_enrichment_modifiers(boolean_input_df, 2, 4, 0.25, input_format = 'Input_',
                                output_format = 'Output_', primary_input_entities = input_list,
                                adj_pval_type = 'bonferroni', sample_names_ind = 'N')

}
\author{
Vijay Kumar Pounraja
}
