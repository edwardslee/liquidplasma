## Introduction

This is the repository for the paper title “Thawed Plasma vs Liquid
Plasma: Tracking Coagulation Factor Activity Changes During Storage” by
Yurtsever et al. We have included the coagulation factor activity data
(fibrinogen, factor V, factor VII, factor VIII, protein C, and protein
S) measured in liquid plasma (on day 15, 26, and 27 of storage) and in
thawed plasma (on day 5).

## Data

Two data sets are included in the `figures` folder. The main data set
has factory activities for liquid plasma and thawed plasma units and is
in files named `liquid_plasma_data_updated.xlsx` and
`liquid_plasma.rds`. We provide the data in both .xlsx and .rds formats
for Excel and R users respectively. The data have the following columns:

-   `sample_id`: IDs of the plasma units used;
-   `product`: describes if the plasma unit is either “liquid plasma” or
    “thawed plasma”
-   `day`: the day of storage when the factor activity measurement was
    performed
-   `factor`: which factor was measured
-   `value`: the activity measurement; fibrinogen has the units mg/dL,
    and the others have the units %
-   `abo`: the ABO type of the donor for the unit
-   `comments`: any relevant comments regarding the measurements

The second data set is called `liquid_correlation.csv` and contains the
data for the analysis that compared measured factor activities in
segments vs main bags. The data have the following columns:

-   `sample`: IDs of the samples used. `_seg` represents segments of the
    unit, and `_bag` represents the main bag of the unit.
-   `bag`: ID of the plasma unit used.
-   `source`: indicates whether the sample is from a segment or a bag
-   `factor`: which factor was measured
-   `value`: the activity measurement; fibrinogen has the units mg/dL,
    and the others have the units %

## Code

The statistical models are coded in R and uses Stan (see
[here](https://mc-stan.org)) for sampling the posteriors of the Bayesian
models. We use both Richard Mcelreath’s rethinking package and directly
coded Stan models in this study. The R scripts are contained in the
`scripts` folder, and the Stan code are contained in the `stan` folder.
The following scripts implement the various parts of the analysis:

-   `lqp_coag_factors.R` runs the Student-t analysis to estimate the
    posterior distribution of coagulation factor activity for each
    product at each time point
-   `compare_lqp_tp5.R` runs the Bayesian two-sample t test analysis
-   `regression_all.R` runs the multilevel piecewise linear regression
    analysis
-   `student-t-analysis.R` contains the analysis of comparing a fixed nu
    vs a free nu when fitting the Student-t likelihoods
-   `segment_correlation.R` contains the analysis for the comparison of
    activities in segments and in their corresponding bags and
    quantifying the correlation between them.

Running the code requires installation of `tidyverse`, `rethinking`,
`cmdstanr`, and `readxl` packages, which can be installed by running the
following commands

    install.packages(c("tidyverse", "cmdstanr", "readxl"))
    devtools::install_github("rmcelreath/rethinking")
