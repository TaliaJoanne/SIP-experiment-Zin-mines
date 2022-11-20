---
title: "Zin restoration SIP"
subtitle: "02 Decontaminate dataset"
description: "V1.2"
author: "Roey Angel"
email: "roey.angel@bc.cas.cz"
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






## Identify and remove contaminant ASVs
Decontamination of sequence library based on [Introduction to decontam](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) and Davis and colleagues [-@davis_simple_2017]. Decontamination is based on correlating sequence abundance frequencies to initial DNA concentrations used for PCR and also on examining sequence abundance prevalence patterns in the negative control samples.

### Setting general parameters:

```r
set.seed(1000)
data_path <- "./DADA2_pseudo/"
Metadata_table <- "./Metadata.csv"
Seq_table <- "DADA2_seqtab_nochim.tsv"
Tax_table <- "DADA2_taxa_silva.tsv"
Seq_file <- "DADA2_reps.fa"
Ps_file <- "DADA2_seqtab_nochim.RDS"
Proj_name <- "Zin_SIP"
Var1 = "Density (g ml-1)" # e.g sampling point / replicate
Var2 = "Treatment" # e.g. a treatment or a manipulation
Var3 = "Label (18O)" # e.g. a treatment/manipulation or an important covariant
Var4 = "Comparison_pair" # e.g. an important covariant
```

### Reading in raw data and generate phyloseq object

```r
# read OTU mat from data file
read_tsv(paste0(data_path, Seq_table),
         trim_ws = TRUE) %>% 
  rename_with(., ~ str_remove(.x, "_L001")) %>% 
  identity() ->
  abundance_mat # in tibble format
```

```
## Rows: 10275 Columns: 219
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr   (1): ASV
## dbl (218): 661Roey001-N1-2_S1_L001, 661Roey002-N1-3_S2_L001, 661Roey003-N1-4...
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
# get short names of samples
# abundance_mat %>% 
#   rownames() %>% 
#   str_remove("^Roey[0-9]{3,}-?") %>% 
#   str_split("_", simplify = T) %>% 
#   .[, 1] ->
#   short_names

# Read metadata file
read_csv(Metadata_table, 
                        trim_ws = TRUE) %>% 
  mutate(`16S copies` = replace(`16S copies`, which(`16S copies` == 0 | is.na(`16S copies`)), 1)) %>%  # add pseudo count
  filter(str_remove(Read1_file, "_L001_R1_001_noPrimers.fastq.gz") %in% colnames(abundance_mat)) %>% # remove metadata rows if the samples did not go through qual processing
  mutate(`Library size` = colSums(select(abundance_mat, -ASV))) %>% # Add lib size
  mutate(to_names = str_remove(Read1_file, "_L001_R1_001_noPrimers.fastq.gz"), .before = 1) %>% 
  mutate(across(c(!!sym(Var2), 
                  !!sym(Var3),
                  !!sym(Var4)), ~factor(.))) %>% 
  mutate(across(!!sym(Var3), ~fct_relevel(., "Labelled", after = Inf))) %>% 
  mutate(across(!!sym(Var1), ~as.numeric(.))) %>% 
  # mutate(Identifier = paste(Site, 
  #                           `Plant species`, 
  #                           `Preservation method`,
                            # Replicate,
  #                           sep = "_"))  %>%   # optional. For merging samples
  identity() ->
  Metadata
```

```
## Rows: 218 Columns: 26
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (10): merged_sample_name, sample_name, Read1_file, Sample, Treatment, La...
## dbl (14): input, filtered, denoised, merged, tabled, nonchim, Hours, TNA ext...
## lgl  (2): Density (old_calc), Control
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```
## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion
```

```r
# Order abundance_mat samples according to the metadata
sample_order <- match(Metadata$to_names, str_remove(colnames(select(abundance_mat, -ASV)), "_L001"))
abundance_mat %<>% select(c(1, sample_order + 1))


# read taxonomy from data file
Raw_tax_data <- read_tsv(paste0(data_path, Tax_table), 
                        trim_ws = TRUE, col_names = TRUE)
```

```
## Rows: 10275 Columns: 15
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (8): ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species
## dbl (7): Kingdom (BS), Phylum (BS), Class (BS), Order (BS), Family (BS), Gen...
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
Raw_tax_data %<>%
  replace(., is.na(.), "Unclassified") 

# Potentially store tax classification BS values
# Raw_tax_data %>%
#   dplyr::select(.,
#          `Kingdom (BS)`,
#          `Phylum (BS)`,
#          `Class (BS)`,
#          `Order (BS)`,
#          `Family (BS)`,
#          `Genus (BS)`) %>%
#   cbind(Name = colnames(abundance_mat),. ) ->
#   Taxonomy.bs

# merge it downstream with the PS object
# taxTraits <- tax_table(cbind(tax_table(Ps_obj), taxTraits))
# Ps_obj <- merge_phyloseq(Ps_obj, taxTraits)
# Taxonomy.bs %<>% 
#   filter(Taxonomy.bs$Name %in% row.names(Ps_obj_filt@tax_table))

Raw_tax_data %>%
  dplyr::select(.,
         Kingdom,
         Phylum,
         Class,
         Order,
         Family,
         Genus) %>% 
  # map_dfr(., as_factor) %>% 
  # map_dfr(fct_expand, "Rare")  %>%  
  mutate(to_names = abundance_mat$ASV, .before = 1)-> # must be a matrix or phyloseq drops row names and gives and error
  Taxonomy
# row.names(Taxonomy) <- colnames(abundance_mat)
# colnames(Taxonomy) <-
#   c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# read sequence data
ASV_seqs <- readDNAStringSet(
  file = paste0(data_path, Seq_file),
  format = "fasta", 
  nrec = -1L, 
  skip = 0L, 
  seek.first.rec = FALSE,
  use.names = TRUE)

# generate phyloseq object. **Note: only speedyseq allows constructing phyloseq from tibble!**
Ps_obj <- phyloseq(otu_table(abundance_mat, taxa_are_rows = TRUE),
                   sample_data(Metadata),
                   tax_table(Taxonomy),
                   refseq(ASV_seqs))
```

```
## Assuming first column, `ASV`, contains the taxa names
```

```
## Assuming first column, `to_names`, contains the sample names
```

```
## Assuming first column, `to_names`, contains the taxa names
```

```r
# rename_with_sample_data(Ps_obj, colnames(select(Metadata, -to_names)))

# merge samples in case the company re-sequenced some samples. **Note! Don't use it now to merge real DNA samples. Instead, do it downstream in one of the next scripts.**
Ps_obj %>% 
  phyloseq_merge_samples("Read1_file") %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) ->
  Ps_obj_merged 
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(fct_vars)` instead of `fct_vars` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(lgl_vars)` instead of `lgl_vars` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```

### Inspect Library Sizes

```r
Ps_obj_df <- as_tibble(sample_data(Ps_obj_merged)) # Put sample_data into a ggplot-friendly data.frame
Ps_obj_df <- Ps_obj_df[order(pull(Ps_obj_df, Library.size)), ]
Ps_obj_df$Index <- seq(nrow(Ps_obj_df))
ggplot(data = Ps_obj_df, 
       aes(x = Index, y = Library.size, color = Control)) + 
  geom_point() +
  scale_y_log10(breaks = c(
    min((pull(Ps_obj_df, Library.size))),
    10,
    100,
    1000,
    5000,
    10000,
    ceiling(max((pull(Ps_obj_df, Library.size))) / 10000) * 10000
    )) + 
  scale_color_brewer(type = 'qual', 
                     palette = 'Set1', 
                     direction = -1)
```

![](02_Decontamination_V1.2_files/figure-html/Library Sizes-1.png)<!-- -->

```r
summary(sample_sums(Ps_obj_merged))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     602   18987   23739   47311   42068  291129
```

```r
summary(taxa_sums(Ps_obj_merged))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       1      45     109    1004     326  777988
```


```r
Ps_obj_merged %<>%
  prune_samples(names(which(sample_sums(Ps_obj_merged) != 0)), .)
```

###  Identify contaminants - Frequency
Use the distribution of the frequency of each sequence feature as a function of the input DNA concentration to identify contaminants.


```r
contamdf_freq <-
  isContaminant(Ps_obj_merged, method = "frequency", conc = "X16S.copies")
# print(contamdf_freq)
# How many contaminants are found?
table(contamdf_freq$contaminant)
```

```
## 
## FALSE  TRUE 
## 10131   144
```

```r
# Which ones
which(contamdf_freq$contaminant)
```

```
##   [1]     6    45    85   152   214   242   280   285   346   387   392   413
##  [13]   434   463   506   514   523   536   541   547   566   567   620   670
##  [25]   691   725   735   783   848   917   929   935   994  1035  1043  1050
##  [37]  1054  1075  1077  1127  1161  1168  1171  1178  1180  1184  1308  1316
##  [49]  1358  1371  1394  1448  1475  1493  1511  1512  1581  1597  1598  1612
##  [61]  1621  1675  1680  1684  1688  1800  1802  1837  1875  1936  1970  1974
##  [73]  2004  2007  2027  2032  2068  2154  2198  2247  2266  2270  2271  2277
##  [85]  2428  2460  2563  2605  2617  2679  2720  2747  2848  2878  2951  2984
##  [97]  3040  3068  3069  3086  3141  3182  3225  3249  3297  3307  3346  3365
## [109]  3531  3655  3667  3752  3755  3902  3960  3991  4082  4112  4135  4281
## [121]  4556  4563  4575  4627  4746  4862  5026  5040  5372  5474  5526  5578
## [133]  5627  6283  6321  6432  6480  6794  7354  7869  8295  8547  9593 10028
```

Plot the frequency of sequence 1 and 3 (non-contaminants) against the DNA concentration, as an example.

```r
plot_frequency(Ps_obj_merged, taxa_names(Ps_obj_merged)[c(1, 3)], conc = "X16S.copies")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

![](02_Decontamination_V1.2_files/figure-html/plot frequency 1-1.png)<!-- -->

Plot the frequency of the contaminant sequences against the DNA concentration.

```r
plot_frequency(Ps_obj_merged, taxa_names(Ps_obj_merged)[which(contamdf_freq$contaminant)[1:20]], conc = "X16S.copies")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

![](02_Decontamination_V1.2_files/figure-html/plot frequency 2-1.png)<!-- -->

The frequency analysis detected $144$ sequences as contaminants.

###  Identify contaminants - Prevalence
Use the prevalence of sequences found in the control samples (no-template controls) to identify contaminants.

```r
contamdf_prev <- isContaminant(Ps_obj_merged, method = "prevalence", neg = "Control")
# How many contaminants are found?
table(contamdf_prev$contaminant)
```

```
## 
## FALSE  TRUE 
## 10242    33
```

```r
# Which ones
which(contamdf_prev$contaminant)
```

```
##  [1]   18   21   90  101  125  212  299  363  500  535  553  725  892  952  953
## [16] 1012 1160 1317 1330 1493 1522 1551 2647 2776 2824 2971 3245 3292 3779 3885
## [31] 3953 5466 9772
```

```r
# And using a more aggressive threshold
contamdf_prev05 <- isContaminant(Ps_obj_merged, method = "prevalence", neg = "Control", threshold = 0.5)
table(contamdf_prev05$contaminant)
```

```
## 
## FALSE  TRUE 
## 10147   128
```

```r
# Make phyloseq object of presence-absence in negative controls
Ps_obj_pa <-
  transform_sample_counts(Ps_obj, function(abund)
    1 * (abund > 0))
Ps_obj_pa_neg <-
  prune_samples(sample_data(Ps_obj_pa)$Control == "TRUE", Ps_obj_pa)
Ps_obj_pa_pos <-
  prune_samples(sample_data(Ps_obj_pa)$Control == "FALSE", Ps_obj_pa)
# Make data.frame of prevalence in positive and negative samples
df_pa <-
  data.frame(
    pa_pos = taxa_sums(Ps_obj_pa_pos),
    pa_neg = taxa_sums(Ps_obj_pa_neg),
    contaminant = contamdf_prev$contaminant
  )
ggplot(data = df_pa, aes(x = pa_neg, y = pa_pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```

![](02_Decontamination_V1.2_files/figure-html/prevalence-1.png)<!-- -->

The frequency analysis detected $33$ sequences as contaminants.
In total $177$ were detected as contaminants and will be removed.

### Save contaminant sequence names and decontaminated data

```r
c(taxa_names(Ps_obj)[which(contamdf_freq$contaminant)],
  taxa_names(Ps_obj)[which(contamdf_prev$contaminant)]) ->
  contaminant_seqs
  
write_csv(as_tibble(contaminant_seqs), 
            paste0(data_path, "decontam_contaminants.csv"), 
            col_names = FALSE)


good_seqs <- setdiff(taxa_names(Ps_obj), contaminant_seqs)
Ps_obj_clean <- prune_taxa(good_seqs, Ps_obj)

# save decontaminated seqtab
Ps_obj_clean %>% 
  # t() %>%
  get_taxa() %>% 
  as_tibble(rownames = "ASV") %>%
  write_tsv(., 
            paste0(data_path, str_remove(Seq_table, ".tsv"), "_decontam.tsv"), 
            col_names = TRUE)

Ps_obj_clean %>% 
  t() %>% 
  tax_table() %>% 
  as_tibble() %>%
  rename(.otu = "ASV") %>% 
  write_tsv(., 
            paste0(data_path, str_remove(Tax_table, ".tsv"), "_decontam.tsv"), 
            col_names = TRUE)

# save decontaminated metadata (just in case some samples were dropped)
Ps_obj_clean %>% 
  t() %>%
  sample_data() %>% 
  # setNames(., colnames(Metadata)) %>% 
  # as_tibble(rownames = "ASV") %>%
  write_csv(., 
            paste0("./", str_remove(Metadata_table, ".csv"), "_decontam.csv"), 
            col_names = TRUE)

# save decontaminated seqs
Ps_obj_clean %>% 
  refseq() %>% 
  writeXStringSet(., filepath = paste0(data_path, "DADA2_reps_decontam.fa"), format = "fasta", width = 1000)
 
# save R obj **saving taxa as columns!**
saveRDS(t(Ps_obj_clean), file = paste0(data_path, Proj_name, "_decontam.Rds"))
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

─ Session info ───────────────────────────────────────────────────────────────
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

─ Packages ───────────────────────────────────────────────────────────────────
 package          * version    date       lib source                           
 ade4               1.7-18     2021-09-16 [1] CRAN (R 4.1.1)                   
 ape                5.5        2021-04-25 [1] CRAN (R 4.0.3)                   
 assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.0.2)                   
 backports          1.2.1      2020-12-09 [1] CRAN (R 4.0.2)                   
 Biobase            2.52.0     2021-05-19 [1] Bioconductor                     
 BiocGenerics     * 0.38.0     2021-05-19 [1] Bioconductor                     
 biomformat         1.20.0     2021-05-19 [1] Bioconductor                     
 Biostrings       * 2.60.2     2021-08-05 [1] Bioconductor                     
 bit                4.0.4      2020-08-04 [1] CRAN (R 4.0.2)                   
 bit64              4.0.5      2020-08-30 [1] CRAN (R 4.0.2)                   
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
 decontam         * 1.12.0     2021-05-19 [1] Bioconductor                     
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
 RColorBrewer       1.1-2      2014-12-07 [1] CRAN (R 4.0.2)                   
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
 scales             1.1.1      2020-05-11 [1] CRAN (R 4.0.2)                   
 sessioninfo        1.1.1      2018-11-05 [1] CRAN (R 4.0.2)                   
 speedyseq        * 0.5.3.9018 2021-08-11 [1] Github (mikemc/speedyseq@ceb941f)
 stringi            1.7.5      2021-10-04 [1] CRAN (R 4.1.1)                   
 stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.0.2)                   
 survival           3.2-13     2021-08-24 [1] CRAN (R 4.1.1)                   
 svglite          * 2.0.0      2021-02-20 [1] CRAN (R 4.1.0)                   
 systemfonts        1.0.3      2021-10-13 [1] CRAN (R 4.1.1)                   
 tibble           * 3.1.5      2021-09-30 [1] CRAN (R 4.1.1)                   
 tidyr            * 1.1.4      2021-09-27 [1] CRAN (R 4.1.1)                   
 tidyselect         1.1.1      2021-04-30 [1] CRAN (R 4.0.3)                   
 tidyverse        * 1.3.1      2021-04-15 [1] CRAN (R 4.0.3)                   
 tzdb               0.1.2      2021-07-20 [1] CRAN (R 4.1.0)                   
 utf8               1.2.2      2021-07-24 [1] CRAN (R 4.1.0)                   
 vctrs              0.3.8      2021-04-29 [1] CRAN (R 4.0.3)                   
 vegan              2.5-7      2020-11-28 [1] CRAN (R 4.0.3)                   
 vroom              1.5.5      2021-09-14 [1] CRAN (R 4.1.1)                   
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

## References

