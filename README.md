# RareComb 
*RareComb* is a combinatorial framework that couples the apriori algorithm with binomial tests to systematically analyze patterns of rare event combinations between groups of interest to identify specific combinations that significantly associate with phenotypes. This generalizable, modular and extensible framework does not depend on apriori knowledge and can detect rare patterns from high-dimensional genetic data and generate interpretable results, making it readily useful for analyzing cohorts of all size ranges and providing a structured approach to dissect the genetic basis of complex disorders.

# Citation

[A general framework for identifying rare variant combinations in complex disorders](https://www.biorxiv.org/content/10.1101/2021.10.01.462832v1)

Vijay Kumar Pounraja, Santhosh Girirajan



# Software Requirements
*RareComb* is an easy-to-use R package with five built-in user-facing functions that take a sparse Boolean dataframe as their input along with multiple input parameters to constrain the execution of these functions and the results they generate. The R package can be installed using the ***'install.packages('RareComb')'*** command and loaded into memory using the ***'library(RareComb)'*** command. *RareComb* depends on the following R packages,

> 1) arules
> 2) dplyr
> 3) pwr
> 4) reshape2
> 5) sqldf
> 6) stats
> 7) stringr
> 8) tidyr

**The five user-facing functions supported by the package along with their descriptions are as follows,**

> **1) analyze_in_out_simultaneity *(Comorbidity Analysis)*** : Analyze the relationship between input and output variables for all combinations that include at least a single output variable and meet all the input criteria specified by the user.

> **2) compare_enrichment *(Enrichment in cases + Non-enrichment in controls)*** : Quantify the enrichment in the observed frequency of cooccurring rare events in combinations that meet all the input criteria specified by the user compared to their corresponding expectation derived under the assumption of independence between the constituent elements within each combination. The function then reports the multiple-testing adjusted significant combinations in which enrichment is observed in cases but not in controls.

> **3) compare_enrichment_depletion *(Enrichment in cases + Depletion in controls)*** : Quantify the enrichment in the observed frequency of cooccurring rare events in combinations that meet all the input criteria specified by the user compared to their corresponding expectation derived under the assumption of independence between the constituent elements within each combination. The function then reports the multiple-testing adjusted significant combinations in which enrichment is observed in cases and depletion is observed in controls.

> **4) compare_enrichment_modifiers *(Must include one of the user-supplied input variables in the significant combinations)*** : Quantify the enrichment in the observed frequency of combinations that include at least one of the input variables supplied by the user as well as meet other user-specified criteria compared to their corresponding expectation derived under the assumption of independence between the constituent elements of each combination. The function then reports the combinations in which enrichment is observed in cases but not in controls.

> **5) compare_expected_vs_observed *(Compare observed frequency with the expected frequency within a single group)*** : Compare the observed frequency of combinations
that meet all the user-specified criteria with their corresponding expectation derived under the assumption of independence between the constituent elements of each combination. Unlike the method 'compare_enrichment', this method does NOT split the groups based on phenotypes. It simply treats the entire input as a single cohort and measures the magnitude of difference between the expected and observed frequencies of cooccurring events.




# Running *RareComb*
**Running *RareComb* involves the following considerations and steps,**

**Things to consider:**
> **1) Make sure the input variables in the input data are prefixed with 'Input_' and the output/outcome variables are prefixed with 'Output_'. If you prefer a different convention, please use the optional parameters 'input_format' and 'output_format' to specify the prefix of your choice**

> **2) Prior to invoking the function 'analyze_in_out_simultaneity', make sure the input file has more than one output/outcome variables**

> **3) Prior to invoking the functions 'compare_enrichment', 'compare_enrichment_depletion' and 'compare_enrichment_modifiers', make sure the input file has EXACTLY one output/outcome variable**


**Steps involved in running the functions in *RareComb*:**
> **Step 1) Generate an input Boolean dataframe with 'n' samples (rows) and 'p' variables (columns).**

> **Step 2) Invoke the function of interest using the input file along with the additional mandatory and optional parameters.**

> **Step 3) The final output returned by all the functions is a dataframe that the users can choose save to an output comma or tab separated file.**


# Usage & Examples


**1)** analyze_in_out_simultaneity(boolean_input_mult_df, combo_length, min_output_count,
                            max_output_count, min_indv_threshold, max_freq_threshold,
                            input_format, output_format, pval_filter_threshold,
                            adj_pval_type)
                 
**Total input parameters :** 10 (6 mandatory + 4 optional)


**`analyze_in_out_simultaneity(boolean_input_mult_df, 3, 1, 2, 5, 0.25,
                            input_format = 'Input_', output_format = 'Output_',
                            pval_filter_threshold = 0.05, adj_pval_type = 'BH')`**


**2)** compare_enrichment(boolean_input_df, combo_length, min_indv_threshold,                                            max_freq_threshold, input_format, output_format,
                         pval_filter_threshold, adj_pval_type, min_power_threshold,
                         sample_names_ind)
                 
**Total input parameters :** 10 (4 mandatory + 6 optional)

**`compare_enrichment(boolean_input_df, 3, 5, 0.25, input_format = 'Input_',
                     output_format = 'Output_', adj_pval_type = 'bonferroni',
                     sample_names_ind = 'N')`**
                     
**3)** compare_enrichment_depletion(boolean_input_df, combo_length, min_indv_threshold,
                              max_freq_threshold, input_format, output_format,
                              pval_filter_threshold, adj_pval_type, min_power_threshold,
                              sample_names_ind)

**Total input parameters :** 10 (4 mandatory + 6 optional)

**`compare_enrichment_depletion(boolean_input_df, 3, 5, 0.25, input_format = 'Input_',
                                output_format = 'Output_', adj_pval_type = 'bonferroni',
                                sample_names_ind = 'N')`**



**4)** compare_enrichment_modifiers(boolean_input_df, combo_length, min_indv_threshold,
                                  max_freq_threshold, primary_input_entities, input_format,
                                  output_format, pval_filter_threshold, adj_pval_type,
                                  min_power_threshold, sample_names_ind)
                                 
**Total input parameters :** 11 (5 mandatory + 6 optional)

**`compare_enrichment_modifiers(boolean_input_df, 2, 4, 0.25, input_format = 'Input_',
                                 output_format = 'Output_', primary_input_entities = input_list,
                                 adj_pval_type = 'bonferroni', sample_names_ind = 'N')`**                              
                                 
**5)** compare_expected_vs_observed(boolean_input_df, combo_length, min_indv_threshold,
                              max_freq_threshold, input_format,
                              pval_filter_threshold, adj_pval_type)

**Total input parameters :** 7 (4 mandatory + 3 optional)

**`compare_expected_vs_observed(boolean_input_df, 2, 10, 0.25, 0.05,
                                 input_format = 'Input_', adj_pval_type = 'BH')`**

**Further details on the list of parameters applicable to each function can be found in the documentation for the R package in the [CRAN website](https://cran.r-project.org/web/packages/available_packages_by_date.html).**


# Output files and columns
Each function returns a dataframe with the list of statistically significant combinations that meet the user-specified input criteria as the output. Since the definition of 'statistical significance' is defined differently for each function, the number and types of columns in the output file will vary for each function depending on the size of the requested combination, number of groups under analysis, if the supporting sample names are requested or not etc. For example, for functions that analyze the data based on a single binary outcome (Case/Control), the output file will contain frequency of individual and cooccurring events in each group separately, whereas the output file from analyzing multiple phenotypes together will only contain frequencies from the single group that is being analyzed. A list of output column names for each of the five functions along with their descriptions are provided below,

#### Output from *'analyze_in_out_simultaneity'*,
Column Names | Column Descriptions
-----------|-------------------
Item_1 | Name of the first item in the combination.
Item_2 | Name of the second item in the combination.
.. | Other items in the combination.
Item_***N*** | Name of the 'N'th item in the combination.
Obs_Count_Combo | Observed frequency of the cooccurring event within the combination.
Case_Obs_Count_I1 | Observed frequency of the individual item 'Item_1' in cases.
Case_Obs_Count_I2 | Observed frequency of the individual item 'Item_2' in cases.
.. |
Case_Obs_Count_I***N*** | Observed frequency of the individual item 'Item_N' in cases.
Output_Count | Number of output variables in the combination.
Exp_Prob_Combo | Expected probability of events to cooccur.
Obs_Prob_Combo | Observed probability of cooccurring events.
pvalue_more | p-values from the **one-tailed binomial test** to evaluate if the observed frequency is **greater** than the expected frequency of cooccurring events.
input_only_pvalue_more | p-values from the **one-tailed binomial test** considering only the input variables (genotype). This p-value can be used to evaluate if the genotypes in a combination by themselves are strongly associated with eachother.
Adj_Pval_bonf | p-values of the genotype-phenotype combination adjusted for multiple testing using the 'bonferroni' method.
Adj_Pval_BH | p-values of the genotype-phenotype combination adjusted for multiple testing using the 'Benjamini-Hochberg' method.


#### Output from *'compare_expected_vs_observed'*,
Column Names | Column Descriptions
-----------|-------------------
Item_1 | Name of the first item in the combination.
Item_2 | Name of the second item in the combination.
.. | Other items in the combination.
Item_***N*** | Name of the 'N'th item in the combination.
Obs_Count_Combo | Observed frequency of the cooccurring event within the combination.
Obs_Count_I1 | Observed frequency of the individual item 'Item_1' in cases.
Obs_Count_I2 | Observed frequency of the individual item 'Item_2' in cases.
.. |
Obs_Count_I***N*** | Observed frequency of the individual item 'Item_N' in cases.
Exp_Prob_Combo | Expected probability of events to cooccur.
Obs_Prob_Combo | Observed probability of cooccurring events.
pvalue_more | p-values from the **one-tailed binomial test** to evaluate if the observed frequency is **greater** than the expected frequency of cooccurring events.
Adj_Pval_bonf | p-values of the genotype-phenotype combination adjusted for multiple testing using the 'bonferroni' method.
Adj_Pval_BH | p-values of the genotype-phenotype combination adjusted for multiple testing using the 'Benjamini-Hochberg' method.


#### Output from *'compare_enrichment'*, *'compare_enrichment_depletion'* & *'compare_enrichment_modifiers'*,
Column Names | Column Descriptions
-----------|-------------------
Item_1 | Name of the first item in the combination.
Item_2 | Name of the second item in the combination.
.. | 
Item_***N*** | Name of the 'N'th item in the combination.
Case_Obs_Count_I1 | Observed frequency of the individual item 'Item_1' in cases.
Case_Obs_Count_I2 | Observed frequency of the individual item 'Item_2' in cases.
.. |
Case_Obs_Count_I***N*** | Observed frequency of the individual item 'Item_N' in cases.
Case_Exp_Prob_Combo | Expected probability of the cooccurring event within the combination in cases.
Case_Obs_Prob_Combo | Observed probability of the cooccurring event within the combination in cases.
Case_Exp_Count_Combo | Expected frequency of the cooccurring event within the combination in cases.
Case_Obs_Count_Combo | Observed frequency of the cooccurring event within the combination in cases.
Case_pvalue_more | p-values from the **one-tailed binomial test** to evaluate if the observed frequency is **greater** than the expected frequency of cooccurring events in cases.
Cont_Obs_Count_I1 | Observed frequency of the individual item 'Item_1' in controls.
Cont_Obs_Count_I2 | Observed frequency of the individual item 'Item_2' in controls.
.. |
Cont_Obs_Count_I***N*** | Observed frequency of the individual item 'Item_N' in controls.
Cont_Exp_Prob_Combo | Expected probability of the cooccurring event within the combination in controls.
Cont_Obs_Prob_Combo | Observed probability of the cooccurring event within the combination in controls.
Cont_Exp_Count_Combo | Expected frequency of the cooccurring event within the combination in controls.
Cont_Obs_Count_Combo | Observed frequency of the cooccurring event within the combination in controls.
Cont_pvalue_more | p-values from the **one-tailed binomial test** to evaluate if the observed frequency is **greater** than the expected frequency of cooccurring events in controls.
Control_pvalue_less (***applies only to 'compare_enrichment_depletion'***) | This output column replaces *Cont_pvalue_more* when the function *compare_enrichment_depletion* is invoked. This column provides the p-values from the **one-tailed binomial test** to evaluate if the observed frequency is **lesser** than the expected frequency of cooccurring events in controls.
Case_Adj_Pval_bonf | p-values of the combination in cases adjusted for multiple testing using the 'bonferroni' method.
Case_Adj_Pval_BH | p-values of the combination in cases adjusted for multiple testing using the 'Benjamini-Hochberg' method.
Effect_Size | Effect size measured as Cohen's d capturing the magnitude of difference in frequency of cooccurring events between cases and controls.
Power_One_Pct | Available statistical power for the 2-sample 2-proportion test to compare the frequencies of cooccurring events in cases and controls at ***1%*** significance threshold.
Power_Five_Pct | Available statistical power for the 2-sample 2-proportion test to compare the frequencies of cooccurring events in cases and controls at ***5%*** significance threshold.
Case_Samples | A list of sample names from ***cases*** that carry the significant combination identified by the function. This column is part of the output only when the function is invoked with the input parameter *'sample_names_ind'* set to *'Y'*.
Control_Samples | A list of sample names from ***controls*** that carry the significant combination identified by the function. This column is part of the output only when the function is invoked with the input parameter *'sample_names_ind'* set to *'Y'*.




# MIT License

Copyright (c) 2021 Vijay Kumar Pounraja

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


# Contact
For questions or comments, please contact Vijay Kumar Pounraja (vxm915@psu.edu) or Santhosh Girirajan (sxg47@psu.edu).
