NCC: R package for simulation and analysis of platform trials with non-concurrent controls
================

datasim\_cont()
===============

The function `datasim_cont()` enables data simulation of a platform trial with continuous endpoint and an arbitrary number of treatment arms entering at different time points.

Assumptions
-----------

-   equal sample sizes across all treatment arms
-   block randomization is used, a factor to multiply the number of active arms with in order to get the block size in each period can be specified as input argument (`period_blocks`, default=2)

Notation
--------

<table>
    <tr>
        <td><b>Paper</b></td>
        <td><b>Software</b></td>
    </tr>
    <tr>
        <td>$N$</td>
        <td>n_total</td>
    </tr>
    <tr>
        <td>$K$</td>
        <td>num_arms</td>
    </tr>
    <tr>
        <td>$\delta$</td>
        <td>t_arms</td>
    </tr>
    <tr>
        <td>$n$</td>
        <td>n_arm</td>
    </tr>
    <tr>
        <td>$\eta_0$</td>
        <td>mu0</td>
    </tr>
    <tr>
        <td>$\theta$</td>
        <td>delta</td>
    </tr>
    <tr>
        <td>$\lambda$</td>
        <td>lambda</td>
    </tr>
    <tr>
        <td>$\sigma$</td>
        <td>sigma</td>
    </tr>
    <tr>
        <td>$N_p$</td>
        <td>N_peak</td>
    </tr>

</table>
Usage options
-------------

Input
-----

The user has two options to define the sample sizes in the trial:

#### Option 1)

Specify total sample size in the trial (`n_total`), number of treatment arms (`num_arms`) and the timing of adding new arms in terms of number of patients allocated to the control arm (`t_arm`). The sample size per period is then computed using the `get_ss_matrix()` function. The sample sizes in treatment and control arms are varying, based on the 3 specified parameters and the assumptions of equal sample sizes for all treatment arms and 1:1:...:1 allocation ratio in all periods.

#### Option 2)

Specify the sample size per treatment arm (`n_arm`, assumed equal) and a matrix with allocation ratios (`alloc_ratios`) in each period. The matrix with sample sizes per period is then computed based on these input parameters.

#### Other input parameters

-   `period_blocks` - number to multiply the number of active arms with in order to get the block size per period (block size = `period_blocks` ⋅ \#active arms)
-   `mu0` - response in the control arm. Default=0
-   `delta` - vector with treatment effects for each arm
-   `lambda` - vector with strength of time trend in each arm
-   `sigma` - residual variance
-   `trend` - indicates the time trend pattern ("linear", "stepwise" or "inv\_u")
-   `N_peak` - point at which the inverted-u time trend switches direction in terms of overall sample size
-   `full` - boolean. Indicates whether the full dataset should be returned. Default=`FALSE`

### Output

The function outputs a dataframe with the simulated trial data. If the parameter `full` is set to `FALSE` (default), only variables needed for the analysis are included. If `full=TRUE`, the dataframe also contains additional information (lambdas, underlying responses).

### Examples

``` r
# Usage option 1

head(datasim_cont(n_total = 1000, num_arms = 3, t_arm = 120, 
                   n_arm = NULL, alloc_ratios = NULL, 
                   period_blocks = 2, delta = rep(0.25, 3), lambda = rep(0.15, 4), 
                   sigma = 1, trend = "linear"))
```

        response treatment period
    1 -0.5403472         1      1
    2 -0.4110947         0      1
    3 -0.6374202         0      1
    4 -0.1747773         1      1
    5  1.7030518         0      1
    6 -0.1522278         1      1

``` r
# Usage option 2

head(datasim_cont(n_total = NULL, num_arms = NULL, t_arm = NULL, 
                   n_arm = 200, alloc_ratios = matrix(c(1,1,1,
                                                        1,1,NA,
                                                        NA,1,1), ncol = 3, byrow = T),
                   period_blocks = 2, delta = rep(0.25, 3), lambda = rep(0.15, 4), 
                   sigma = 1, trend = "linear"))
```

        response treatment period
    1  2.3016223         0      1
    2  1.3360588         1      1
    3  0.4002106         0      1
    4  0.9501642         1      1
    5 -1.3864670         0      1
    6  2.0977060         0      1

``` r
# Usage option 1 - full dataset

head(datasim_cont(n_total = 1000, num_arms = 3, t_arm = 120, 
                   n_arm = NULL, alloc_ratios = NULL, 
                   period_blocks = 2, delta = rep(0.25, 3), lambda = rep(0.15, 4), 
                   sigma = 1, trend = "linear", full = T))
```

      j   response treatment period        means lambda0 lambda1 lambda2 lambda3
    1 1 -1.0765237         1      1 0.2500000000    0.15    0.15    0.15    0.15
    2 2 -0.6347640         0      1 0.0001501502    0.15    0.15    0.15    0.15
    3 3  1.1518523         0      1 0.0003003003    0.15    0.15    0.15    0.15
    4 4 -0.4457669         1      1 0.2504504505    0.15    0.15    0.15    0.15
    5 5  0.1921266         0      1 0.0006006006    0.15    0.15    0.15    0.15
    6 6  1.1135093         0      1 0.0007507508    0.15    0.15    0.15    0.15

------------------------------------------------------------------------

**Funding**

[EU-PEARL](https://eu-pearl.eu/) (EU Patient-cEntric clinicAl tRial pLatforms) project has received funding from the Innovative Medicines Initiative (IMI) 2 Joint Undertaking (JU) under grant agreement No 853966. This Joint Undertaking receives support from the European Union’s Horizon 2020 research and innovation programme and EFPIA andChildren’s Tumor Foundation, Global Alliance for TB Drug Development non-profit organisation, Spring works Therapeutics Inc. This publication reflects the authors’ views. Neither IMI nor the European Union, EFPIA, or any Associated Partners are responsible for any use that may be made of the information contained herein.

<!-- ## Example?? -->
<!-- This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated. -->
<!-- ... Including Plots -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo=FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->
