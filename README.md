Simulations of platform trials with non-concurrent controls
================

data\_sim\_cont()
-----------------

The function `data_sim_cont()` enables data simulation of a platform trial with an arbitrary number of treatment arms entering at different time points.

### Usage options

The user has two options to specify input arguments:

#### 1)

Specify total sample size in the trial (n\_total), number of treatment arms (num\_arms) and the timing of adding new arms in terms of number of patients allocated to the control arm (t\_arm). The sample size per period is then computed using the `get_ss_matrix()` function. The sample sizes in treatment and control arms are varying, based on the 3 specified parameters and the assumptions of equal sample sizes for all treatment arms and 1:1:...:1 allocation ratio in all periods.

``` r
# Example

head(data_sim_cont(n_total = 1000, num_arms = 3, t_arm = 120, 
                   n_arm = NULL, alloc_ratios = NULL, 
                   period_blocks = 2, delta = rep(0.25, 3), lambda = rep(0.15, 4), 
                   sigma = 1, trend = "linear"))
```

    ## Loading required package: rlang

    ## Warning: package 'rlang' was built under R version 3.6.3

    ##     response treatment period j        means lambda0 lambda1 lambda2 lambda3
    ## 1  0.2322430         0      1 1 0.0000000000    0.15    0.15    0.15    0.15
    ## 2 -0.8054627         1      1 2 0.2501501502    0.15    0.15    0.15    0.15
    ## 3  0.6682034         1      1 3 0.2503003003    0.15    0.15    0.15    0.15
    ## 4  0.6622246         0      1 4 0.0004504505    0.15    0.15    0.15    0.15
    ## 5 -1.1549538         1      1 5 0.2506006006    0.15    0.15    0.15    0.15
    ## 6 -1.8217935         0      1 6 0.0007507508    0.15    0.15    0.15    0.15

#### 2)

Specify the sample size per treatment arm (assumed equal) and a matrix with allocation ratios in each period. The matrix with sample sizes per period is then computed based on these input parameters.

``` r
# Example

head(data_sim_cont(n_total = NULL, num_arms = NULL, t_arm = NULL, 
                   n_arm = 200, alloc_ratios = , matrix(c(1,1,1,
                                                          1,1,NA,
                                                          NA,1,1), ncol = 3, byrow = T),
                   period_blocks = 2, delta = rep(0.25, 3), lambda = rep(0.15, 4), 
                   sigma = 1, trend = "linear"))
```

    ##      response treatment period j        means lambda0 lambda1 lambda2
    ## 1  0.17560861         1      1 1 0.2500000000    0.15    0.15    0.15
    ## 2  0.01861719         0      1 2 0.0002145923    0.15    0.15    0.15
    ## 3  1.05355776         0      1 3 0.0004291845    0.15    0.15    0.15
    ## 4  0.25512900         1      1 4 0.2506437768    0.15    0.15    0.15
    ## 5 -0.38940674         0      1 5 0.0008583691    0.15    0.15    0.15
    ## 6  0.75342553         1      1 6 0.2510729614    0.15    0.15    0.15

### Assumptions

-   sample sizes are equal across all treatment arms
-   block randomization is used, a factor to multiply the number of active arms with in order to get the block size in each period can be specified as input argument (period\_blocks, default=2)

get\_ss\_matrix()
-----------------

``` r
# Example 

get_ss_matrix(1000, 4, 100)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    ## [1,]  100   40   60   40   60   40  100
    ## [2,]  100   40   NA   NA   NA   NA   NA
    ## [3,]   NA   40   60   40   NA   NA   NA
    ## [4,]   NA   NA   NA   40   60   40   NA
    ## [5,]   NA   NA   NA   NA   NA   40  100

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
