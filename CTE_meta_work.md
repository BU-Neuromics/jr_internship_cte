Import data

    library(haven)
    library(plyr)
    library(DESeq2)

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     rename

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     desc

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::arrange()    masks plyr::arrange()
    ## ✖ dplyr::collapse()   masks IRanges::collapse()
    ## ✖ dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ purrr::compact()    masks plyr::compact()
    ## ✖ dplyr::count()      masks matrixStats::count(), plyr::count()
    ## ✖ dplyr::desc()       masks IRanges::desc(), plyr::desc()
    ## ✖ tidyr::expand()     masks S4Vectors::expand()
    ## ✖ dplyr::failwith()   masks plyr::failwith()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::first()      masks S4Vectors::first()
    ## ✖ dplyr::id()         masks plyr::id()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ dplyr::mutate()     masks plyr::mutate()
    ## ✖ ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()     masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()     masks S4Vectors::rename(), plyr::rename()
    ## ✖ dplyr::slice()      masks IRanges::slice()
    ## ✖ dplyr::summarise()  masks plyr::summarise()
    ## ✖ dplyr::summarize()  masks plyr::summarize()

    library(knitr)
    library(striprtf)
    spss_meta <- read_sav("/restricted/projectnb/cteseq/projects/CTE_data_processing/09.28.2023_dataset/CTE RNAseq phenotype data.sav")
    #spss_meta <- as.data.frame(spss_meta)
    csv_meta <- read.csv("/restricted/projectnb/cteseq/projects/CTE_data_processing/all_meta.csv")
    counts <- read.csv("/restricted/projectnb/cteseq/projects/CTE_data_processing/all_counts.csv")

    library(dplyr)

    spss_col <- spss_meta
    spss_col <- select(spss_col, -subjid, -locratio, -subratio, -CA1ratio, -CA23ratio, -CA4ratio, -AT8sulcus, -AT8crest)
    spss_names <- colnames(spss_col)

    val_tbl <- lapply(
      spss_names,
      function(col_name) {
        spss_col %>% count(.data[[col_name]])
      }
    )
    val_tbl

    ## [[1]]
    ## # A tibble: 7 × 2
    ##    apoe     n
    ##   <dbl> <int>
    ## 1    22     3
    ## 2    23    10
    ## 3    24     7
    ## 4    33   115
    ## 5    34    43
    ## 6    44     6
    ## 7    NA    37
    ## 
    ## [[2]]
    ## # A tibble: 30 × 2
    ##    footyrs     n
    ##      <dbl> <int>
    ##  1       1     1
    ##  2       2     1
    ##  3       3     7
    ##  4       4     5
    ##  5       5     1
    ##  6       6     7
    ##  7       7     5
    ##  8       8    16
    ##  9       9    14
    ## 10      10    11
    ## # … with 20 more rows
    ## 
    ## [[3]]
    ## # A tibble: 18 × 2
    ##      AFE     n
    ##    <dbl> <int>
    ##  1     0     1
    ##  2     4     1
    ##  3     5     7
    ##  4     6     8
    ##  5     7    12
    ##  6     8    13
    ##  7     9     8
    ##  8    10    28
    ##  9    11     9
    ## 10    12    26
    ## 11    13    23
    ## 12    14    37
    ## 13    15    14
    ## 14    16    10
    ## 15    18     3
    ## 16    20     1
    ## 17    22     1
    ## 18    NA    19
    ## 
    ## [[4]]
    ## # A tibble: 32 × 2
    ##    totyrs     n
    ##     <dbl> <int>
    ##  1      1     3
    ##  2      2     1
    ##  3      3     4
    ##  4      4     7
    ##  5      5     2
    ##  6      6     7
    ##  7      7     4
    ##  8      8    15
    ##  9      9    12
    ## 10     10    10
    ## # … with 22 more rows
    ## 
    ## [[5]]
    ## # A tibble: 6 × 2
    ##                                    race     n
    ##                               <dbl+lbl> <int>
    ## 1  1 [White]                              196
    ## 2  2 [Black/African American]              18
    ## 3  3 [American Indian/Alaska Native]        2
    ## 4  4 [Native Hawaiian/Pacific Islander]     1
    ## 5  8 [Other:]                               3
    ## 6 NA                                        1
    ## 
    ## [[6]]
    ## # A tibble: 9 × 2
    ##                           cod     n
    ##                     <dbl+lbl> <int>
    ## 1  1 [Suicide]                   13
    ## 2  2 [Accidental Overdose]        1
    ## 3  3 [Cardiovascular Disease]    42
    ## 4  4 [Neurodegenerative]        107
    ## 5  5 [Motor Neuron Disease]       1
    ## 6  6 [Cancer]                    21
    ## 7  7 [Other]                     25
    ## 8  8 [Injury]                     5
    ## 9 NA                              6
    ## 
    ## [[7]]
    ## # A tibble: 2 × 2
    ##   DementiaHx     n
    ##    <dbl+lbl> <int>
    ## 1    0 [No]     73
    ## 2    1 [Yes]   148
    ## 
    ## [[8]]
    ## # A tibble: 3 × 2
    ##                                         nphemo     n
    ##                                      <dbl+lbl> <int>
    ## 1  0 [0 No (SKIP TO QUESTION 12c)]               209
    ## 2  1 [1 Yes (COMPLETE QUESTION Q12b1 - Q12b3)]     2
    ## 3 NA                                              10
    ## 
    ## [[9]]
    ## # A tibble: 5 × 2
    ##            npavas     n
    ##         <dbl+lbl> <int>
    ## 1  0 [0 None]       107
    ## 2  1 [1 Mild]        52
    ## 3  2 [2 Moderate]    39
    ## 4  3 [3 Severe]      20
    ## 5 NA                  3
    ## 
    ## [[10]]
    ## # A tibble: 5 × 2
    ##             npwmr     n
    ##         <dbl+lbl> <int>
    ## 1  0 [0 None]        29
    ## 2  1 [1 Mild]        70
    ## 3  2 [2 Moderate]    91
    ## 4  3 [3 Severe]      28
    ## 5 NA                  3
    ## 
    ## [[11]]
    ## # A tibble: 9 × 2
    ##                                                             npbraak     n
    ##                                                           <dbl+lbl> <int>
    ## 1  0 [0 Stage 0: AD-type neurofibrillary degeneration not present]     41
    ## 2  1 [1 Stage I (B1)]                                                  11
    ## 3  2 [2 Stage II (B1)]                                                 28
    ## 4  3 [3 Stage III (B2)]                                                60
    ## 5  4 [4 Stage IV (B2)]                                                 28
    ## 6  5 [5 Stage V (B3)]                                                  20
    ## 7  6 [6 Stage VI (B3)]                                                 31
    ## 8  7 [7 The presence of a tauopathy (other than aging/AD) preclude]     1
    ## 9 NA                                                                    1
    ## 
    ## [[12]]
    ## # A tibble: 5 × 2
    ##                                  npneur     n
    ##                               <dbl+lbl> <int>
    ## 1  0 [0 No neuritic plaques (C0)]         115
    ## 2  1 [1 Sparse neuritic plaques (C1)]      60
    ## 3  2 [2 Moderate neuritic plaques (C2)]    34
    ## 4  3 [3 Frequent neuritic plaques (C3)]    11
    ## 5 NA                                        1
    ## 
    ## [[13]]
    ## # A tibble: 5 × 2
    ##                     npadnc     n
    ##                  <dbl+lbl> <int>
    ## 1  0 [0 Not AD]               83
    ## 2  1 [1 Low ADNC]             32
    ## 3  2 [2 Intermediate ADNC]    60
    ## 4  3 [3 High ADNC]            38
    ## 5 NA                           8
    ## 
    ## [[14]]
    ## # A tibble: 5 × 2
    ##                            npdiff     n
    ##                         <dbl+lbl> <int>
    ## 1  0 [0 No diffuse plaques]          72
    ## 2  1 [1 Sparse diffuse plaques]      46
    ## 3  2 [2 Moderate diffuse plaques]    41
    ## 4  3 [3 Frequent diffuse plaques]    59
    ## 5 NA                                  3
    ## 
    ## [[15]]
    ## # A tibble: 5 × 2
    ##             npamy     n
    ##         <dbl+lbl> <int>
    ## 1  0 [0 None]       119
    ## 2  1 [1 Mild]        44
    ## 3  2 [2 Moderate]    33
    ## 4  3 [3 Severe]      24
    ## 5 NA                  1
    ## 
    ## [[16]]
    ## # A tibble: 3 × 2
    ##                                        npold     n
    ##                                    <dbl+lbl> <int>
    ## 1  0 [0 No (SKIP TO QUESTION Q12d)]            164
    ## 2  1 [1 Yes (COMPLETE QUESTIONS Q12c1-12c4)]    56
    ## 3 NA                                             1
    ## 
    ## [[17]]
    ## # A tibble: 5 × 2
    ##           npold1     n
    ##        <dbl+lbl> <int>
    ## 1  0 [0]           187
    ## 2  1 [1]            15
    ## 3  2 [2]            11
    ## 4  3 [3 or more]     7
    ## 5 NA                 1
    ## 
    ## [[18]]
    ## # A tibble: 5 × 2
    ##           npold2     n
    ##        <dbl+lbl> <int>
    ## 1  0 [0]           200
    ## 2  1 [1]            11
    ## 3  2 [2]             2
    ## 4  3 [3 or more]     7
    ## 5 NA                 1
    ## 
    ## [[19]]
    ## # A tibble: 5 × 2
    ##           npold3     n
    ##        <dbl+lbl> <int>
    ## 1  0 [0]           201
    ## 2  1 [1]            12
    ## 3  2 [2]             5
    ## 4  3 [3 or more]     2
    ## 5 NA                 1
    ## 
    ## [[20]]
    ## # A tibble: 5 × 2
    ##           npold4     n
    ##        <dbl+lbl> <int>
    ## 1  0 [0]           202
    ## 2  1 [1]            13
    ## 3  2 [2]             3
    ## 4  3 [3 or more]     1
    ## 5 NA                 2
    ## 
    ## [[21]]
    ## # A tibble: 3 × 2
    ##                                         npoldd     n
    ##                                      <dbl+lbl> <int>
    ## 1  0 [0 No (SKIP TO QUESTION Q12e)]              216
    ## 2  1 [1 Yes (COMPLETE QUESTION Q12d1 - Q12d4)]     4
    ## 3 NA                                               1
    ## 
    ## [[22]]
    ## # A tibble: 4 × 2
    ##          npoldd1     n
    ##        <dbl+lbl> <int>
    ## 1  0 [0]           217
    ## 2  1 [1]             2
    ## 3  3 [3 or more]     1
    ## 4 NA                 1
    ## 
    ## [[23]]
    ## # A tibble: 2 × 2
    ##     npoldd2     n
    ##   <dbl+lbl> <int>
    ## 1     0 [0]   220
    ## 2    NA         1
    ## 
    ## [[24]]
    ## # A tibble: 3 × 2
    ##     npoldd3     n
    ##   <dbl+lbl> <int>
    ## 1     0 [0]   219
    ## 2     1 [1]     1
    ## 3    NA         1
    ## 
    ## [[25]]
    ## # A tibble: 2 × 2
    ##     npoldd4     n
    ##   <dbl+lbl> <int>
    ## 1     0 [0]   220
    ## 2    NA         1
    ## 
    ## [[26]]
    ## # A tibble: 5 × 2
    ##           nparter     n
    ##         <dbl+lbl> <int>
    ## 1  0 [0 None]        37
    ## 2  1 [1 Mild]        48
    ## 3  2 [2 Moderate]   101
    ## 4  3 [3 Severe]      33
    ## 5 NA                  2
    ## 
    ## [[27]]
    ## # A tibble: 3 × 2
    ##                                       nppath     n
    ##                                    <dbl+lbl> <int>
    ## 1  0 [0 No (SKIP TO QUESTION 13)]              183
    ## 2  1 [1 Yes (COMPLETE QUESTION 12g1- 12g12)]    30
    ## 3 NA                                             8
    ## 
    ## [[28]]
    ## # A tibble: 3 × 2
    ##      nppath6     n
    ##    <dbl+lbl> <int>
    ## 1  0 [0 No]    212
    ## 2  1 [1 Yes]     1
    ## 3 NA             8
    ## 
    ## [[29]]
    ## # A tibble: 7 × 2
    ##                         nplbod     n
    ##                      <dbl+lbl> <int>
    ## 1  0 [0 No]                      166
    ## 2  1 [1 Brainstem predominant]    18
    ## 3  2 [2 Limbic (transitional)]    10
    ## 4  3 [3 Neocortical (diffuse)]    14
    ## 5  4 [4 Amygdala predominant]      4
    ## 6  5 [5 Olfactory bulb]            8
    ## 7 NA                               1
    ## 
    ## [[30]]
    ## # A tibble: 3 × 2
    ##   npftdtau     n
    ##      <dbl> <int>
    ## 1        0    62
    ## 2        1   158
    ## 3       NA     1
    ## 
    ## [[31]]
    ## # A tibble: 2 × 2
    ##      nppick     n
    ##   <dbl+lbl> <int>
    ## 1  0 [0 No]   220
    ## 2 NA            1
    ## 
    ## [[32]]
    ## # A tibble: 2 × 2
    ##     npftdt2     n
    ##   <dbl+lbl> <int>
    ## 1  0 [0 No]   213
    ## 2 NA            8
    ## 
    ## [[33]]
    ## # A tibble: 3 × 2
    ##       npcort     n
    ##    <dbl+lbl> <int>
    ## 1  0 [0 No]    219
    ## 2  1 [1 Yes]     1
    ## 3 NA             1
    ## 
    ## [[34]]
    ## # A tibble: 3 × 2
    ##       npprog     n
    ##    <dbl+lbl> <int>
    ## 1  0 [0 No]    213
    ## 2  1 [1 Yes]     7
    ## 3 NA             1
    ## 
    ## [[35]]
    ## # A tibble: 3 × 2
    ##      npftdt5     n
    ##    <dbl+lbl> <int>
    ## 1  0 [0 No]    212
    ## 2  1 [1 Yes]     2
    ## 3 NA             7
    ## 
    ## [[36]]
    ## # A tibble: 2 × 2
    ##     npftdt8     n
    ##   <dbl+lbl> <int>
    ## 1  0 [0 No]   213
    ## 2 NA            8
    ## 
    ## [[37]]
    ## # A tibble: 3 × 2
    ##      npftdt9     n
    ##    <dbl+lbl> <int>
    ## 1  0 [0 No]    211
    ## 2  1 [1 Yes]     2
    ## 3 NA             8
    ## 
    ## [[38]]
    ## # A tibble: 2 × 2
    ##   npftdt10     n
    ##      <dbl> <int>
    ## 1        0   214
    ## 2       NA     7
    ## 
    ## [[39]]
    ## # A tibble: 3 × 2
    ##   npftdtdp     n
    ##      <dbl> <int>
    ## 1        0   211
    ## 2        1     9
    ## 3       NA     1
    ## 
    ## [[40]]
    ## # A tibble: 2 × 2
    ##                             npoftd     n
    ##                          <dbl+lbl> <int>
    ## 1  0 [0 No (SKIP TO QUESTION 18a)]   213
    ## 2 NA                                   8
    ## 
    ## [[41]]
    ## # A tibble: 61 × 2
    ##    agedeath     n
    ##       <dbl> <int>
    ##  1       18     1
    ##  2       21     1
    ##  3       22     1
    ##  4       23     1
    ##  5       25     2
    ##  6       27     2
    ##  7       29     1
    ##  8       30     1
    ##  9       32     1
    ## 10       34     1
    ## # … with 51 more rows
    ## 
    ## [[42]]
    ## # A tibble: 3 × 2
    ##      PathAD     n
    ##   <dbl+lbl> <int>
    ## 1   0 [No]    166
    ## 2   1 [Yes]    53
    ## 3  NA           2
    ## 
    ## [[43]]
    ## # A tibble: 4 × 2
    ##                          PathLBD     n
    ##                        <dbl+lbl> <int>
    ## 1  0 [Not present]                 178
    ## 2  1 [LBD Brainstem predominant]    18
    ## 3  2 [LBD limbic/neocortical]       24
    ## 4 NA                                 1
    ## 
    ## [[44]]
    ## # A tibble: 3 × 2
    ##        PathFTD     n
    ##      <dbl+lbl> <int>
    ## 1  0 [Absent]    195
    ## 2  1 [Present]    18
    ## 3 NA               8
    ## 
    ## [[45]]
    ## # A tibble: 3 × 2
    ##   PathMND     n
    ##     <dbl> <int>
    ## 1       0   212
    ## 2       1     1
    ## 3      NA     8
    ## 
    ## [[46]]
    ## # A tibble: 2 × 2
    ##   PathPrion     n
    ##   <dbl+lbl> <int>
    ## 1    0 [No]   213
    ## 2   NA          8
    ## 
    ## [[47]]
    ## # A tibble: 3 × 2
    ##         CTE     n
    ##   <dbl+lbl> <int>
    ## 1   0 [No]     68
    ## 2   1 [Yes]   151
    ## 3  NA           2
    ## 
    ## [[48]]
    ## # A tibble: 6 × 2
    ##         CTEStage     n
    ##        <dbl+lbl> <int>
    ## 1  0 [None]         62
    ## 2  1 [Stage I]      13
    ## 3  2 [Stage II]     19
    ## 4  3 [Stage III]    48
    ## 5  4 [Stage IV]     70
    ## 6 NA                 9
    ## 
    ## [[49]]
    ## # A tibble: 2 × 2
    ##   micdorfront     n
    ##         <dbl> <int>
    ## 1           0   213
    ## 2          NA     8
    ## 
    ## [[50]]
    ## # A tibble: 3 × 2
    ##   micsuptemp     n
    ##        <dbl> <int>
    ## 1          0   212
    ## 2          1     1
    ## 3         NA     8
    ## 
    ## [[51]]
    ## # A tibble: 2 × 2
    ##   micinfpar     n
    ##       <dbl> <int>
    ## 1         0   218
    ## 2        NA     3
    ## 
    ## [[52]]
    ## # A tibble: 2 × 2
    ##   micalc     n
    ##    <dbl> <int>
    ## 1      0   219
    ## 2     NA     2
    ## 
    ## [[53]]
    ## # A tibble: 60 × 2
    ##    agecogsx     n
    ##       <dbl> <int>
    ##  1        8     1
    ##  2       10     1
    ##  3       15     1
    ##  4       16     3
    ##  5       20     4
    ##  6       23     2
    ##  7       28     1
    ##  8       30     3
    ##  9       32     1
    ## 10       33     1
    ## # … with 50 more rows
    ## 
    ## [[54]]
    ## # A tibble: 61 × 2
    ##    cogdecage     n
    ##        <dbl> <int>
    ##  1        13     1
    ##  2        15     1
    ##  3        17     1
    ##  4        19     1
    ##  5        20     1
    ##  6        21     2
    ##  7        22     3
    ##  8        23     2
    ##  9        25     3
    ## 10        26     3
    ## # … with 51 more rows
    ## 
    ## [[55]]
    ## # A tibble: 67 × 2
    ##    mincogage     n
    ##        <dbl> <int>
    ##  1         8     1
    ##  2        10     1
    ##  3        13     1
    ##  4        15     1
    ##  5        16     2
    ##  6        17     1
    ##  7        19     1
    ##  8        20     5
    ##  9        21     2
    ## 10        22     2
    ## # … with 57 more rows
    ## 
    ## [[56]]
    ## # A tibble: 3 × 2
    ##     suicide     n
    ##   <dbl+lbl> <int>
    ## 1   0 [No]    202
    ## 2   1 [Yes]    13
    ## 3  NA           6
    ## 
    ## [[57]]
    ## # A tibble: 46 × 2
    ##    disdur     n
    ##     <dbl> <int>
    ##  1      0     2
    ##  2      1     5
    ##  3      2     3
    ##  4      3     5
    ##  5      4    12
    ##  6      5    13
    ##  7      6    10
    ##  8      7     6
    ##  9      8     6
    ## 10      9     9
    ## # … with 36 more rows
    ## 
    ## [[58]]
    ## # A tibble: 9 × 2
    ##                    sport     n
    ##                <dbl+lbl> <int>
    ## 1  1 [Football]            176
    ## 2  2 [Hockey]                4
    ## 3  3 [Boxing]                7
    ## 4  4 [Soccer]                6
    ## 5  5 [Rugby]                 2
    ## 6  6 [Pro Wrestling]         1
    ## 7  7 [Amateur Wrestling]     3
    ## 8  8 [Lacrosse]              2
    ## 9 NA                        20
    ## 
    ## [[59]]
    ## # A tibble: 5 × 2
    ##                                                     Group     n
    ##                                                 <dbl+lbl> <int>
    ## 1  1 [RHI without CTE]                                       54
    ## 2  2 [Mild CTE (Stages I or II)]                             31
    ## 3  3 [Severe CTE (stages III or IV)]                        118
    ## 4  4 [Non-contact sports w/out Neurodegenerative disease]    16
    ## 5 NA                                                          2
    ## 
    ## [[60]]
    ## # A tibble: 5 × 2
    ##   rs1990622     n
    ##   <chr>     <int>
    ## 1 ""          106
    ## 2 "C:C"        16
    ## 3 "C:T"        46
    ## 4 "T:C"        10
    ## 5 "T:T"        43
    ## 
    ## [[61]]
    ## # A tibble: 4 × 2
    ##   rs3173615     n
    ##   <chr>     <int>
    ## 1 ""          148
    ## 2 "C"          27
    ## 3 "CG"         35
    ## 4 "G"          11
