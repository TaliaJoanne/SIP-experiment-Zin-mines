---
title: "Zin restoration SIP"
subtitle: "05 Phylogentic analysis"
description: "V1.1"
author: "Roey Angel"
date: "2021-10-21"
bibliography: references.bib
link-citations: yes
csl: fems-microbiology-ecology.csl
output:
  rmarkdown::html_document:
    toc: true
    toc_float: true
    keep_md: true
    number_sections: false
    highlight: "pygments"
    theme: "flatly"
    dev: "png"
    df_print: "kable"
    fig_caption: true
    code_folding: "show"
---






[roey.angel@bc.cas.cz](mailto: roey.angel@bc.cas.cz)  

## Phylogenetic analysis
This analysis explores the phylogenetic distribution patters in the different samples, based on the DADA2-produced sequences. Large parts of this script are based on [this protocol](https://f1000research.com/articles/5-1492/v2) and the accompanying publication by Callahan and colleagues [-@callahan_bioconductor_2016]. This script should be run after [05_calc_tree_V2.0.sh](05_calc_tree_V2.0.sh).

### Setting general parameters:

```r
set.seed(1000)
min_lib_size <- 2500
data_path <- "./DADA2_pseudo/"
Proj_name <- "Zin_SIP"
Ps_file <- paste0(Proj_name, "_seq_prev_filt.Rds")
Tree_file <- "Tree/DADA2_reps_seq_prev_filt.filtered.align.treefile"
```

### Reading in raw data
Read abundance table, taxonomic classification and metadata into a phyloseq object.


```r
# read phyloseq object from data file
Ps_obj <- readRDS(paste0(data_path, Ps_file))

# Load phylogenetic tree
Tree <- read_tree(paste0(data_path, Tree_file))

# generate phyloseq object
Ps_obj <- merge_phyloseq(Ps_obj,
                        phy_tree(Tree)
                        )
```

### Exploring dataset features

### Plot phylogenetic trees
Now let's try to simplify the phylogenetic tree by agglomerating genera

```r
Ps_obj_subset <-
  subset_samples(Ps_obj, sample_sums(Ps_obj) > min_lib_size)
Ps_obj_subset <-
  filter_taxa(Ps_obj_subset, function(x)
    sum(x) > 0, TRUE)

# How many genera would be present after filtering?
length(get_taxa_unique(Ps_obj_subset, taxonomic.rank = "Genus"))
```

```
## [1] 209
```

```r
Ps_obj_glom <- tax_glom(Ps_obj_subset, 
                             "Genus", 
                             NArm = TRUE)

multiPlotTitleTextSize = 8
p2tree <- plot_tree(Ps_obj_subset, method = "treeonly",
                     ladderize = "left",
                     title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree <- plot_tree(Ps_obj_glom, method = "treeonly",
                     ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together
grid.arrange(nrow = 1, p2tree, p3tree)
```

![](05_Phylogeny_figures/tree-1.svg)<!-- -->

#### Save filtered phyloseq object

```r
saveRDS(Ps_obj, file = paste0(data_path, Proj_name, "_filt_wTree.Rds"))
saveRDS(Ps_obj_glom, file = paste0(data_path, Proj_name, "_filt_glom_wTree.Rds"))
```


```r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
  )
```

<details open>
<summary> <span title='Click to Expand'> Current session info </span> </summary>

```r

─ Session info ─────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.1.1 (2021-08-10)
 os       Ubuntu 18.04.6 LTS          
 system   x86_64, linux-gnu           
 ui       X11                         
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Prague               
 date     2021-10-21                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package          * version    date       lib source                           
 ade4               1.7-18     2021-09-16 [1] CRAN (R 4.1.1)                   
 ape                5.5        2021-04-25 [1] CRAN (R 4.0.3)                   
 assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.0.2)                   
 backports          1.2.1      2020-12-09 [1] CRAN (R 4.0.2)                   
 Biobase            2.52.0     2021-05-19 [1] Bioconductor                     
 BiocGenerics     * 0.38.0     2021-05-19 [1] Bioconductor                     
 biomformat         1.20.0     2021-05-19 [1] Bioconductor                     
 Biostrings       * 2.60.2     2021-08-05 [1] Bioconductor                     
 bitops             1.0-7      2021-04-24 [1] CRAN (R 4.0.3)                   
 broom              0.7.9      2021-07-27 [1] CRAN (R 4.1.0)                   
 bslib              0.3.1      2021-10-06 [1] CRAN (R 4.1.1)                   
 cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.0.2)                   
 cli                3.0.1      2021-07-17 [1] CRAN (R 4.1.0)                   
 clipr              0.7.1      2020-10-08 [1] CRAN (R 4.0.2)                   
 cluster            2.1.2      2021-04-17 [1] CRAN (R 4.0.3)                   
 codetools          0.2-18     2020-11-04 [1] CRAN (R 4.0.2)                   
 colorspace         2.0-2      2021-06-24 [1] CRAN (R 4.1.0)                   
 crayon             1.4.1      2021-02-08 [1] CRAN (R 4.0.3)                   
 data.table         1.14.2     2021-09-27 [1] CRAN (R 4.1.1)                   
 DBI                1.1.1      2021-01-15 [1] CRAN (R 4.0.3)                   
 dbplyr             2.1.1      2021-04-06 [1] CRAN (R 4.0.3)                   
 desc               1.4.0      2021-09-28 [1] CRAN (R 4.1.1)                   
 details            0.2.1      2020-01-12 [1] CRAN (R 4.0.2)                   
 digest             0.6.28     2021-09-23 [1] CRAN (R 4.1.1)                   
 dplyr            * 1.0.7      2021-06-18 [1] CRAN (R 4.1.0)                   
 ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.0.3)                   
 evaluate           0.14       2019-05-28 [1] CRAN (R 4.0.2)                   
 extrafont        * 0.17       2014-12-08 [1] CRAN (R 4.1.0)                   
 extrafontdb        1.0        2012-06-11 [1] CRAN (R 4.0.2)                   
 fansi              0.5.0      2021-05-25 [1] CRAN (R 4.0.3)                   
 farver             2.1.0      2021-02-28 [1] CRAN (R 4.0.3)                   
 fastmap            1.1.0      2021-01-25 [1] CRAN (R 4.0.3)                   
 forcats          * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)                   
 foreach            1.5.1      2020-10-15 [1] CRAN (R 4.0.2)                   
 fs                 1.5.0      2020-07-31 [1] CRAN (R 4.0.2)                   
 generics           0.1.0      2020-10-31 [1] CRAN (R 4.0.2)                   
 GenomeInfoDb     * 1.28.4     2021-09-05 [1] Bioconductor                     
 GenomeInfoDbData   1.2.6      2021-05-25 [1] Bioconductor                     
 ggplot2          * 3.3.5      2021-06-25 [1] CRAN (R 4.1.0)                   
 glue               1.4.2      2020-08-27 [1] CRAN (R 4.0.2)                   
 gridExtra        * 2.3        2017-09-09 [1] CRAN (R 4.0.2)                   
 gtable             0.3.0      2019-03-25 [1] CRAN (R 4.0.2)                   
 haven              2.4.3      2021-08-04 [1] CRAN (R 4.1.0)                   
 highr              0.9        2021-04-16 [1] CRAN (R 4.0.3)                   
 hms                1.1.1      2021-09-26 [1] CRAN (R 4.1.1)                   
 htmltools          0.5.2      2021-08-25 [1] CRAN (R 4.1.1)                   
 httr               1.4.2      2020-07-20 [1] CRAN (R 4.0.2)                   
 igraph             1.2.7      2021-10-15 [1] CRAN (R 4.1.1)                   
 IRanges          * 2.26.0     2021-05-19 [1] Bioconductor                     
 iterators          1.0.13     2020-10-15 [1] CRAN (R 4.0.2)                   
 jquerylib          0.1.4      2021-04-26 [1] CRAN (R 4.0.3)                   
 jsonlite           1.7.2      2020-12-09 [1] CRAN (R 4.0.2)                   
 knitr              1.36       2021-09-29 [1] CRAN (R 4.1.1)                   
 labeling           0.4.2      2020-10-20 [1] CRAN (R 4.0.2)                   
 lattice            0.20-45    2021-09-22 [1] CRAN (R 4.1.1)                   
 lifecycle          1.0.1      2021-09-24 [1] CRAN (R 4.1.1)                   
 lubridate          1.8.0      2021-10-07 [1] CRAN (R 4.1.1)                   
 magrittr         * 2.0.1      2020-11-17 [1] CRAN (R 4.0.2)                   
 MASS               7.3-54     2021-05-03 [1] CRAN (R 4.0.3)                   
 Matrix             1.3-4      2021-06-01 [1] CRAN (R 4.1.0)                   
 mgcv               1.8-38     2021-10-06 [1] CRAN (R 4.1.1)                   
 modelr             0.1.8      2020-05-19 [1] CRAN (R 4.0.2)                   
 multtest           2.48.0     2021-05-19 [1] Bioconductor                     
 munsell            0.5.0      2018-06-12 [1] CRAN (R 4.0.2)                   
 nlme               3.1-153    2021-09-07 [1] CRAN (R 4.1.1)                   
 permute            0.9-5      2019-03-12 [1] CRAN (R 4.0.2)                   
 phyloseq         * 1.36.0     2021-05-19 [1] Bioconductor                     
 pillar             1.6.4      2021-10-18 [1] CRAN (R 4.1.1)                   
 pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.0.2)                   
 plyr               1.8.6      2020-03-03 [1] CRAN (R 4.0.2)                   
 png                0.1-7      2013-12-03 [1] CRAN (R 4.0.2)                   
 purrr            * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)                   
 R6                 2.5.1      2021-08-19 [1] CRAN (R 4.1.1)                   
 ragg             * 1.1.3      2021-06-09 [1] CRAN (R 4.1.0)                   
 Rcpp               1.0.7      2021-07-07 [1] CRAN (R 4.1.0)                   
 RCurl              1.98-1.5   2021-09-17 [1] CRAN (R 4.1.1)                   
 readr            * 2.0.2      2021-09-27 [1] CRAN (R 4.1.1)                   
 readxl             1.3.1      2019-03-13 [1] CRAN (R 4.0.2)                   
 reprex             2.0.1      2021-08-05 [1] CRAN (R 4.1.0)                   
 reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.0.2)                   
 rhdf5              2.36.0     2021-05-19 [1] Bioconductor                     
 rhdf5filters       1.4.0      2021-05-19 [1] Bioconductor                     
 Rhdf5lib           1.14.2     2021-07-06 [1] Bioconductor                     
 rlang              0.4.12     2021-10-18 [1] CRAN (R 4.1.1)                   
 rmarkdown          2.11       2021-09-14 [1] CRAN (R 4.1.1)                   
 rprojroot          2.0.2      2020-11-15 [1] CRAN (R 4.0.2)                   
 rstudioapi         0.13       2020-11-12 [1] CRAN (R 4.0.2)                   
 Rttf2pt1           1.3.9      2021-07-22 [1] CRAN (R 4.1.0)                   
 rvest              1.0.2      2021-10-16 [1] CRAN (R 4.1.1)                   
 S4Vectors        * 0.30.2     2021-10-03 [1] Bioconductor                     
 sass               0.4.0      2021-05-12 [1] CRAN (R 4.0.3)                   
 scales           * 1.1.1      2020-05-11 [1] CRAN (R 4.0.2)                   
 sessioninfo        1.1.1      2018-11-05 [1] CRAN (R 4.0.2)                   
 speedyseq        * 0.5.3.9018 2021-08-11 [1] Github (mikemc/speedyseq@ceb941f)
 stringi            1.7.5      2021-10-04 [1] CRAN (R 4.1.1)                   
 stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.0.2)                   
 survival           3.2-13     2021-08-24 [1] CRAN (R 4.1.1)                   
 svglite          * 2.0.0      2021-02-20 [1] CRAN (R 4.1.0)                   
 systemfonts        1.0.3      2021-10-13 [1] CRAN (R 4.1.1)                   
 textshaping        0.3.6      2021-10-13 [1] CRAN (R 4.1.1)                   
 tibble           * 3.1.5      2021-09-30 [1] CRAN (R 4.1.1)                   
 tidyr            * 1.1.4      2021-09-27 [1] CRAN (R 4.1.1)                   
 tidyselect         1.1.1      2021-04-30 [1] CRAN (R 4.0.3)                   
 tidyverse        * 1.3.1      2021-04-15 [1] CRAN (R 4.0.3)                   
 tzdb               0.1.2      2021-07-20 [1] CRAN (R 4.1.0)                   
 utf8               1.2.2      2021-07-24 [1] CRAN (R 4.1.0)                   
 vctrs              0.3.8      2021-04-29 [1] CRAN (R 4.0.3)                   
 vegan              2.5-7      2020-11-28 [1] CRAN (R 4.0.3)                   
 withr              2.4.2      2021-04-18 [1] CRAN (R 4.0.3)                   
 xfun               0.27       2021-10-18 [1] CRAN (R 4.1.1)                   
 xml2               1.3.2      2020-04-23 [1] CRAN (R 4.0.2)                   
 XVector          * 0.32.0     2021-05-19 [1] Bioconductor                     
 yaml               2.2.1      2020-02-01 [1] CRAN (R 4.0.2)                   
 zlibbioc           1.38.0     2021-05-19 [1] Bioconductor                     

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library

```

</details>
<br>

```r
## References
```
