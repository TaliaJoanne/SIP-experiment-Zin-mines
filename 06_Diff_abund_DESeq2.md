Zin SIP
================
Roey Angel
2021-10-22

-   [Differential abundance modelling of SIP
    gradients](#differential-abundance-modelling-of-sip-gradients)
    -   [Setting general parameters:](#setting-general-parameters)
    -   [Load phyloseq object](#load-phyloseq-object)
    -   [Beta diversity analysis](#beta-diversity-analysis)
    -   [Subset the dataset](#subset-the-dataset)
    -   [Differential abundance models](#differential-abundance-models)
        -   [Inspect results](#inspect-results)
        -   [Plot differential abundance
            models](#plot-differential-abundance-models)
    -   [Plot labelled ASVs](#plot-labelled-asvs)
        -   [Plot phylogenetic trees with
            heatmaps](#plot-phylogenetic-trees-with-heatmaps)
-   [References](#references)

## Differential abundance modelling of SIP gradients

Here we attempt to detect ASVs that were labelled with <sup>18</sup>O in
the soil incubations using differential abundance modelling. Using
DESeq2 ([Love, Huber and Anders 2014](#ref-love_moderated_2014)) we
compare the relative abundance of each ASV in the fractions where
<sup>18</sup>O-labelled DNA is expected to be found (&gt;1.70 g
ml<sup>-1</sup>; AKA ‘heavy’ fractions) in the labelled gradients to the
heavy fractions in the unlabelled gradient. The method has been
previously described in Angel ([2019](#ref-angel_stable_2019)).

### Setting general parameters:

``` r
set.seed(2021)
alpha_thresh <- 0.1
LFC_thresh <- 0.2
samples_prep_path <- "./"
data_path <- "./DADA2_pseudo/"
Proj_name <- "Zin_SIP"
Ps_file <- paste0(Proj_name, "_filt_wTree.Rds")
Tree_file <- "./Tree/DADA2.Seqs_decontam_filtered.filtered.align.treefile"
```

### Load phyloseq object

This phyloseq object was created in
[05\_Taxonomical\_analysis.html](05_Taxonomical_analysis.html) by
including the iqtree-calculated tree. The phyloseq object excludes
contaminants, all sequences classified as eukaryota, chloroplast,
mitochondria or unknown, taxa with low prevalence.

``` r
Ps_obj <- readRDS(file = paste0(data_path, Ps_file))
# Ps_obj <- phyloseq_replace_zero(Ps_obj)
```

### Beta diversity analysis

Let us look first at the dissimilarity in community composition between
the different fractions. If the labelling was strong enough we should
see a deviation of (some of) the heavy fractions from the light ones.
However, a lack of a significant deviation does not mean unsuccessful
labelling because if only a small minority of the community was labelled
we might not see it here (but we will, hopefully, see it using DESeq2
modelling).

``` r
(mod1 <- adonis(vegdist(otu_table(Ps_obj), method = "horn") ~ Treatment + Library.size,
  data = as(sample_data(Ps_obj), "data.frame"),
  permutations = 999
))
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj), method = "horn") ~      Treatment + Library.size, data = as(sample_data(Ps_obj),      "data.frame"), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Treatment      1     4.195  4.1948  13.817 0.05794  0.001 ***
    ## Library.size   1     8.399  8.3989  27.665 0.11600  0.001 ***
    ## Residuals    197    59.808  0.3036         0.82606           
    ## Total        199    72.401                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot_lib_dist(Ps_obj)
```

![](06_Diff_abund_DESeq2_figures/beta%20div%20joint-1.png)<!-- -->

``` r
Ps_obj %>%
  scale_libraries(round = "round") ->
  Ps_obj_SIP_scaled
  
plot_lib_dist(Ps_obj_SIP_scaled)
```

![](06_Diff_abund_DESeq2_figures/beta%20div%20joint-2.png)<!-- -->

``` r
(mod2 <- adonis(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn") ~ Treatment + Library.size,
  data = as(sample_data(Ps_obj_SIP_scaled), "data.frame"),
  permutations = 999
))
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn") ~      Treatment + Library.size, data = as(sample_data(Ps_obj_SIP_scaled),      "data.frame"), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Treatment      1     4.184  4.1841   13.77 0.05778  0.001 ***
    ## Library.size   1     8.374  8.3741   27.56 0.11564  0.001 ***
    ## Residuals    197    59.859  0.3039         0.82658           
    ## Total        199    72.417                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(mod3 <- adonis(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn") ~ Treatment * Density.zone,
  data = as(sample_data(Ps_obj_SIP_scaled), "data.frame"),
  permutations = 999
))
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn") ~      Treatment * Density.zone, data = as(sample_data(Ps_obj_SIP_scaled),      "data.frame"), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Treatment                1     4.184  4.1841 12.9358 0.05778  0.001 ***
    ## Density.zone             1     3.397  3.3969 10.5018 0.04691  0.001 ***
    ## Treatment:Density.zone   1     1.439  1.4390  4.4488 0.01987  0.001 ***
    ## Residuals              196    63.397  0.3235         0.87544           
    ## Total                  199    72.417                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
Treatment_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Treatment"))
permutest(Treatment_disp)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
    ## Groups      1 0.08608 0.086082 7.3589    999  0.007 **
    ## Residuals 198 2.31616 0.011698                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(Treatment_disp)
```

![](06_Diff_abund_DESeq2_figures/beta%20div%20joint-3.png)<!-- -->

``` r
#Comparison_pair_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Comparison_pair"))
#permutest(Oxygen_disp)
#plot(Oxygen_disp)
#Hours_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Hours"))
#permutest(Hours_disp)
#plot(Hours_disp)
Density_disp <- betadisper(vegdist(otu_table(Ps_obj_SIP_scaled), method = "horn"), get_variable(Ps_obj_SIP_scaled, "Density.zone"))
permutest(Density_disp)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)  
    ## Groups      1 0.0581 0.058122 3.0648    999   0.09 .
    ## Residuals 198 3.7549 0.018964                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(Density_disp)
```

![](06_Diff_abund_DESeq2_figures/beta%20div%20joint-4.png)<!-- -->

``` r
Ord <- ordinate(Ps_obj_SIP_scaled, "CAP", "horn", 
                formula =  ~ Treatment * Density.zone)
explained <- as.numeric(format(round(eigenvals(Ord)/sum(eigenvals(Ord)) * 100, 1), nsmall = 1))
Ord_plt <- plot_ordination(Ps_obj, Ord, type = "samples", color = "Label..18O.", justDF = TRUE)

p_ord_joint <- ggplot(Ord_plt) +
  geom_point(aes(
               x = CAP1,
               y = CAP2,
               color = Label..18O.,
               size = Density..g.ml.1.,
               shape = Treatment
             ), alpha = 2 / 3) +
  guides(colour = guide_legend(title = "Labelling"), 
         size = guide_legend(title = "Density (g ml<sup>-1</sup>)"),
         shape = guide_legend(title = "Treatment")) +
  scale_colour_locuszoom() +
  # scale_colour_manual(values = Gradient.colours) +
  # scale_fill_manual(values = Gradient.colours, guide = "none") +
  labs(x = sprintf("CAP1 (%s%%)", explained[1]),
  y = sprintf("CAP2 (%s%%)", explained[2])) +
  coord_fixed(ratio = sqrt(explained[2] / explained[1])) +
   theme(legend.justification = "top",
         legend.title = element_markdown(size = 11)
         ) +
  scale_size_continuous(breaks = round(c(seq(min(Ord_plt$Density..g.ml.1.), 
                                       max(Ord_plt$Density..g.ml.1.), 
                                       length.out = 5), 
                                   1), 4),
                        range = c(0.1, 5)) +
  facet_grid("Treatment") +
  # ggtitle("Joint analysis") +
  NULL

save_figure(paste0(fig.path, "Oridnation"), 
            p_ord_joint, 
            pwidth = 10, 
            pheight = 8,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Oridnation", ".png"))
```

<img src="06_Diff_abund_DESeq2_figures/Oridnation.png" width="1920" />

### Subset the dataset

Because the DESeq2 models will be run on each gradient-pair separately,
we need to subset This is easily done using `HTSSIP::phyloseq_subset`
([Youngblut, Barnett and Buckley 2018](#ref-youngblut_htssip_2018))

``` r
# split, ignore replicates (for labelled ASV plots)
# test_expr_1 <- "(Site == '${Site}' & Oxygen == '${Oxygen}' & Label..13C. == 'Unlabelled') | (Site == '${Site}'  & Oxygen == '${Oxygen}' & Label..13C. == '${Label..13C.}')"
# params_1 <- get_treatment_params(Ps_obj_SIP, c("Site",
#                                    "Oxygen",
#                                    "Glucose",
#                                    "Label..13C."),
#                      "Label..13C. != 'Unlabelled'")

test_expr_1 <- "Treatment == '${Treatment}'"

params_1 <- get_treatment_params(Ps_obj, c("Treatment"))

Ps_obj %>% 
  subset_samples(Density.zone == "Heavy") %>% 
  phyloseq_subset(., params_1, test_expr_1) %>% 
  mclapply(.,
           function(x) {filter_taxa(x, function(y) sum(y) > 0, TRUE)},
           mc.cores = nrow(params_1)) -> # remove 0-summed ASVs
  Ps_obj_SIP_all_reps_l

names(Ps_obj_SIP_all_reps_l) %<>%
  map(., ~str_remove_all(.x, "\\s\\|\\s.*")) %>%
  map(., ~str_remove_all(.x, "\\(|\\)|Treatment == |'"))

#split, include time points (for DESeq2 modelling)
test_expr_2 <- "(Treatment == '${Treatment}' & Comparison_pair  == '${Comparison_pair}')"
params_2 <- get_treatment_params(Ps_obj, c("Treatment", 
                                   "Comparison_pair"))

# Generate a list of subsetted phyloseq objects
Ps_obj %>% 
  subset_samples(Density.zone == "Heavy") %>% 
  phyloseq_subset(., params_2, test_expr_2) %>% 
    mclapply(.,
           function(x) {filter_taxa(x, function(y) sum(y) > 0, TRUE)},
           mc.cores = nrow(params_2)) -> # remove 0-summed ASVs
  Ps_obj_SIP_each_rep_l

names(Ps_obj_SIP_each_rep_l) %<>%
  map(., ~str_remove_all(.x, "\\s\\|\\s.*")) %>%
  map(., ~str_remove_all(.x, "\\(|\\)|Treatment == |Comparison_pair  == |'"))
```

### Differential abundance models

Now run the differential abundance models using DESeq2. We then filter
the resutls to include only ASVs with Log\_2\_ fold change
&gt;`LFC_thresh` and significant at P&lt;`alpha_thresh`. Lastly, we run
‘LFC-shrinking’ based on Stephens
([**stephens\_fdr\_2016?**](#ref-stephens_fdr_2016)).

``` r
# generate a DESeq2 object
DESeq_obj_SIP_each_rep_l <- mclapply(Ps_obj_SIP_each_rep_l, 
                                   function(x) {phyloseq_to_deseq2_safe(x, 
                                                                        test_condition = "Label..18O.", 
                                                                        ref_level = "Unlabelled")}, 
                                   mc.cores = nrow(params_2))

# DESeq_obj_SIP_each_rep_zb_l <- mclapply(DESeq_obj_SIP_each_rep_l, 
#                                         function(x) {zinbwave(x,
#                                                               X="~ 1",
#                                                               epsilon = 1e10,
#                                                               verbose = TRUE,
#                                                               K = 0,
#                                                               observationalWeights = TRUE,
#                                                               BPPARAM = BiocParallel::SerialParam())},
#                                         mc.cores = nrow(params_2))

# DESeqDataSet(dds_zinbwave, design = ~ Label..18O.)


# run dds pipeline
DESeq_obj_SIP_each_rep_l %<>% mclapply(., 
                                     function(x) {DESeq(x, 
                                                        test = "Wald",
                                                        # test = "LRT",
                                                        # reduced = ~1, 
                                                        fitType = "local",
                                                        sfType = "poscounts")}, 
                                     mc.cores = nrow(params_2)) # run dds pipeline

map(seq(length(DESeq_obj_SIP_each_rep_l)), 
                        ~plotDispEsts(DESeq_obj_SIP_each_rep_l[[.x]]))
```

![](06_Diff_abund_DESeq2_figures/DESeq2%20models-1.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-2.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-3.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-4.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-5.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-6.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-7.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-8.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-9.png)<!-- -->![](06_Diff_abund_DESeq2_figures/DESeq2%20models-10.png)<!-- -->

    ## [[1]]
    ## [[1]]$rect
    ## [[1]]$rect$w
    ## [1] 1.025522
    ## 
    ## [[1]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[1]]$rect$left
    ## [1] 3.147725
    ## 
    ## [[1]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[1]]$text
    ## [[1]]$text$x
    ## [1] 3.446905 3.446905 3.446905
    ## 
    ## [[1]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[2]]
    ## [[2]]$rect
    ## [[2]]$rect$w
    ## [1] 0.971229
    ## 
    ## [[2]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[2]]$rect$left
    ## [1] 2.846002
    ## 
    ## [[2]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[2]]$text
    ## [[2]]$text$x
    ## [1] 3.129343 3.129343 3.129343
    ## 
    ## [[2]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[3]]
    ## [[3]]$rect
    ## [[3]]$rect$w
    ## [1] 0.9998511
    ## 
    ## [[3]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[3]]$rect$left
    ## [1] 3.254806
    ## 
    ## [[3]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[3]]$text
    ## [[3]]$text$x
    ## [1] 3.546497 3.546497 3.546497
    ## 
    ## [[3]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[4]]
    ## [[4]]$rect
    ## [[4]]$rect$w
    ## [1] 0.950096
    ## 
    ## [[4]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[4]]$rect$left
    ## [1] 2.763033
    ## 
    ## [[4]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[4]]$text
    ## [[4]]$text$x
    ## [1] 3.040209 3.040209 3.040209
    ## 
    ## [[4]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[5]]
    ## [[5]]$rect
    ## [[5]]$rect$w
    ## [1] 0.8653113
    ## 
    ## [[5]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[5]]$rect$left
    ## [1] 2.691079
    ## 
    ## [[5]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[5]]$text
    ## [[5]]$text$x
    ## [1] 2.94352 2.94352 2.94352
    ## 
    ## [[5]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[6]]
    ## [[6]]$rect
    ## [[6]]$rect$w
    ## [1] 0.8841422
    ## 
    ## [[6]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[6]]$rect$left
    ## [1] 2.504911
    ## 
    ## [[6]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[6]]$text
    ## [[6]]$text$x
    ## [1] 2.762845 2.762845 2.762845
    ## 
    ## [[6]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[7]]
    ## [[7]]$rect
    ## [[7]]$rect$w
    ## [1] 0.9312939
    ## 
    ## [[7]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[7]]$rect$left
    ## [1] 2.747865
    ## 
    ## [[7]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[7]]$text
    ## [[7]]$text$x
    ## [1] 3.019555 3.019555 3.019555
    ## 
    ## [[7]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[8]]
    ## [[8]]$rect
    ## [[8]]$rect$w
    ## [1] 1.019995
    ## 
    ## [[8]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[8]]$rect$left
    ## [1] 2.884701
    ## 
    ## [[8]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[8]]$text
    ## [[8]]$text$x
    ## [1] 3.182269 3.182269 3.182269
    ## 
    ## [[8]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[9]]
    ## [[9]]$rect
    ## [[9]]$rect$w
    ## [1] 0.9839436
    ## 
    ## [[9]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[9]]$rect$left
    ## [1] 2.906529
    ## 
    ## [[9]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[9]]$text
    ## [[9]]$text$x
    ## [1] 3.193579 3.193579 3.193579
    ## 
    ## [[9]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481
    ## 
    ## 
    ## 
    ## [[10]]
    ## [[10]]$rect
    ## [[10]]$rect$w
    ## [1] 0.935808
    ## 
    ## [[10]]$rect$h
    ## [1] 2.460759
    ## 
    ## [[10]]$rect$left
    ## [1] 2.523998
    ## 
    ## [[10]]$rect$top
    ## [1] -5.899241
    ## 
    ## 
    ## [[10]]$text
    ## [[10]]$text$x
    ## [1] 2.797005 2.797005 2.797005
    ## 
    ## [[10]]$text$y
    ## [1] -6.51443 -7.12962 -7.74481

``` r
# extract results from a DESeq analysis
# DESeq_res_SIP_each_rep_l  <- mclapply(DESeq_obj_SIP_each_rep_l, 
#                              function(x) {
#                                results(x, 
#                                        altHypothesis = "greater",
#                                        alpha = alpha_thresh, 
#                                        # filterFun = ihw,
#                                        contrast = c("Label..18O.", "Labelled", "Unlabelled"))}, # redundant if phyloseq_to_deseq2_safe() was used but doesn't hurt
#                              mc.cores = nrow(params_2)) 

DESeq_res_SIP_each_rep_LFC_l <- mclapply(DESeq_obj_SIP_each_rep_l, 
                                     function(x) {
                                       results(x,
                                               lfcThreshold = LFC_thresh,
                                               altHypothesis = "greater",
                                               alpha = alpha_thresh,
                                               # filterFun = ihw, # optional alternative to BH (package IHW)
                                               contrast = c("Label..18O.", "Labelled", "Unlabelled"))}, # redundant if phyloseq_to_deseq2_safe() was used but doesn't hurt
                                     mc.cores = nrow(params_2)) # Extract results from a DESeq analysis


DESeq_res_SIP_each_rep_LFC_shrink_l <- map(seq(length(DESeq_obj_SIP_each_rep_l)), 
                                             ~lfcShrink(DESeq_obj_SIP_each_rep_l[[.x]],
                                                         res = DESeq_res_SIP_each_rep_LFC_l[[.x]],
                                                         coef = "Label..18O._Labelled_vs_Unlabelled",
                                                         type = "ashr"))
names(DESeq_res_SIP_each_rep_LFC_shrink_l) <- names(DESeq_res_SIP_each_rep_LFC_l)

# Compare
# plotMA(DESeq_res_SIP_each_rep_l[[2]], ylim = c(-2,2))
plotMA(DESeq_res_SIP_each_rep_LFC_l[[1]])
```

![](06_Diff_abund_DESeq2_figures/DESeq2%20models-11.png)<!-- -->

``` r
plotMA(DESeq_res_SIP_each_rep_LFC_shrink_l[[1]])
```

![](06_Diff_abund_DESeq2_figures/DESeq2%20models-12.png)<!-- -->

``` r
# summarise results (lfcShrink doesn't change the values)
for (i in seq(1, length(DESeq_res_SIP_each_rep_LFC_l))) { # didn't manage with map
  print(names(DESeq_res_SIP_each_rep_LFC_l[i]))
  summary(DESeq_res_SIP_each_rep_LFC_l[[i]])
}
```

    ## [1] "Natural & 1"
    ## 
    ## out of 1093 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 414, 38%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 361, 33%
    ## low counts [2]     : 102, 9.3%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Natural & 2"
    ## 
    ## out of 560 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 5, 0.89%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 238, 42%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Natural & 3"
    ## 
    ## out of 208 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 26, 12%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 97, 47%
    ## low counts [2]     : 53, 25%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Natural & 4"
    ## 
    ## out of 108 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 6, 5.6%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 58, 54%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Natural & 5"
    ## 
    ## out of 99 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 1, 1%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 62, 63%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Restored & 1"
    ## 
    ## out of 169 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 19, 11%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 97, 57%
    ## low counts [2]     : 17, 10%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Restored & 2"
    ## 
    ## out of 125 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 5, 4%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 69, 55%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Restored & 3"
    ## 
    ## out of 181 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 5, 2.8%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 94, 52%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Restored & 4"
    ## 
    ## out of 186 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 4, 2.2%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 98, 53%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## [1] "Restored & 5"
    ## 
    ## out of 213 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0.20 (up)    : 35, 16%
    ## LFC < -0.20 (down) : 0, 0%
    ## outliers [1]       : 98, 46%
    ## low counts [2]     : 57, 27%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
DESeq_res_SIP_each_rep_LFC_shrink_l %>% 
  map(., ~subset(.x, padj < alpha_thresh & log2FoldChange > LFC_thresh)) %>% 
  map(., ~as.data.frame(.x)) %>% 
  map(., ~rownames_to_column(.x, "ASV")) %>% 
  bind_rows(., .id = "Comparison") %>% 
  arrange(Comparison, desc(baseMean))  %>% 
  separate(., "Comparison" ,c("Treatment","Comparison_pair"), sep = " & ") %T>% 
  write_csv(., file = paste0("DESeq2_each_rep_a-", alpha_thresh, "_LFC0-", LFC_thresh, ".txt")) ->
  DESeq_res_SIP_each_rep_LFC_sig_df
```

#### Inspect results

``` r
DESeq_res_SIP_each_rep_LFC_sig_df %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_value()
```

![](06_Diff_abund_DESeq2_figures/vis%20DES%20res-1.png)<!-- -->

``` r
DESeq_res_SIP_each_rep_LFC_sig_df %>% 
  get_variable() %>% 
  select_if(is.numeric) %>% 
  vis_cor()
```

![](06_Diff_abund_DESeq2_figures/vis%20DES%20res-2.png)<!-- -->

#### Plot differential abundance models

``` r
# ps_obj <- Ps_obj
# DESeq_results <- DESeq_res_SIP_byTime_LFC0.322_l[9]
# plot_DESeq(DESeq_results, ps_obj, plot_title = names(DESeq_results))

DESeq_plots <- map(seq(length(DESeq_res_SIP_each_rep_LFC_shrink_l)), 
                        ~plot_DESeq(DESeq_res_SIP_each_rep_LFC_shrink_l[.x],  
                                                Ps_obj, plot_title = names(DESeq_res_SIP_each_rep_LFC_shrink_l[.x])))

Natural_DESeq2 <- ((DESeq_plots[[1]] + 
                     theme(legend.position = "none") +
                     theme(axis.text.x = element_blank())) +
                    (DESeq_plots[[2]] + 
                       theme(legend.position = "none", 
                             axis.text.x = element_blank(), 
                             axis.title.y = element_blank())) +
                    (DESeq_plots[[3]] + 
                       theme(legend.position = "none",
                             axis.text.x = element_blank())) +
                    (DESeq_plots[[4]] + 
                       theme(legend.position = "none", 
                             axis.text.x = element_blank(), 
                             axis.title.y = element_blank())) +
                    (DESeq_plots[[5]]) + 
                    plot_layout(ncol = 2, guides = "collect") & 
                    theme(legend.position = 'bottom'))

save_figure(paste0(fig.path, "Natural_DESeq2"), 
            Natural_DESeq2, 
            pwidth = 14, 
            pheight = 12,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Natural_DESeq2", ".png"))
```

<img src="06_Diff_abund_DESeq2_figures/Natural_DESeq2.png" width="2688" />

``` r
Restored_DESeq2 <- ((DESeq_plots[[6]] + 
                     theme(legend.position = "none") +
                     theme(axis.text.x = element_blank())) +
                    (DESeq_plots[[7]] + 
                       theme(legend.position = "none", 
                             axis.text.x = element_blank(), 
                             axis.title.y = element_blank())) +
                    (DESeq_plots[[8]] + 
                       theme(legend.position = "none",
                             axis.text.x = element_blank())) +
                    (DESeq_plots[[9]] + 
                       theme(legend.position = "none", 
                             axis.text.x = element_blank(), 
                             axis.title.y = element_blank())) +
                    (DESeq_plots[[10]] + 
                       theme(legend.position = "none", 
                             axis.title.y = element_blank())) + 
                    plot_layout(ncol = 2, guides = "collect") & 
                    theme(legend.position = 'bottom'))

save_figure(paste0(fig.path, "Restored_DESeq2"), 
            Restored_DESeq2, 
            pwidth = 14, 
            pheight = 12,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Restored_DESeq2", ".png"))
```

<img src="06_Diff_abund_DESeq2_figures/Restored_DESeq2.png" width="2688" />

### Plot labelled ASVs

``` r
plot_combintions <- crossing(Treatment = c("Natural", "Restored"))

Labelled_ASVs <- map(seq(length(Ps_obj_SIP_all_reps_l)), ~plot_otus_by_density(Ps_obj_SIP_all_reps_l[[.x]], 
                     ASV2plot = filter(DESeq_res_SIP_each_rep_LFC_sig_df, Treatment == plot_combintions$Treatment[.x])))

map(seq(length(Ps_obj_SIP_all_reps_l)), 
    ~save_figure(paste0(fig.path, "Labelled_ASVs_", paste(plot_combintions[.x, ], collapse = "_")), 
                 Labelled_ASVs[[.x]], 
                 pwidth = 16, 
                 pheight = 12,
                 dpi = 600))
```

    ## [[1]]
    ## [1] "06_Diff_abund_DESeq2_figures/Labelled_ASVs_Natural.svgz"
    ## 
    ## [[2]]
    ## [1] "06_Diff_abund_DESeq2_figures/Labelled_ASVs_Restored.svgz"

``` r
plots2display <- list.files(path = paste0(fig.path), 
                    pattern = "^Labelled_ASVs_(.*).png$",
                    full.names = TRUE)

knitr::include_graphics(plots2display)
```

<img src="06_Diff_abund_DESeq2_figures//Labelled_ASVs_Natural.png" width="3072" /><img src="06_Diff_abund_DESeq2_figures//Labelled_ASVs_Restored.png" width="3072" />

#### Plot phylogenetic trees with heatmaps

``` r
c("Natural", "Restored")  ->
  col_order

DESeq_res_SIP_each_rep_LFC_shrink_l %>% 
  map(., ~as.data.frame(.x)) %>% 
  map(., ~rownames_to_column(.x, "ASV")) %>% 
  bind_rows(., .id = "Comparison") %>% 
  # filter(str_detect(Comparison, "Labelled")) %>% # remove unlabelled samples [c(-5, -10, -15, -20)]
  mutate(Labelled = ifelse(padj < alpha_thresh & log2FoldChange > LFC_thresh, "Labelled", "Unlabelled")) %>% 
  # arrange(Comparison, desc(baseMean)) %>% 
  separate(., "Comparison" ,c("Treatment","Comparison_pair"), sep = " & ") %>% 
  mutate(Treatment_pair = paste(Treatment, Comparison_pair)) %>% 
  mutate(across(Treatment_pair, ~as.factor(.))) %>% 
  # mutate(Site_Oxygen = factor(paste0(Site, "-", Oxygen),
  #                             levels = c("Plesne-Oxic", "Plesne-Anoxic", "Certovo-Oxic", "Certovo-Anoxic"),
  #                             labels = c("Pl-Ox", "Pl-Anox", "Ct-Ox", "Ct-Anox"))) %>%
  # mutate(across(c("Hours"), ~factor(., 
  #                                   levels = c("12 h", "24 h", "48 h", "72 h", "216 h"),
  #                                   labels = c("12", "24", "48", "72", "216")))) 
  identity() ->
  # mutate(Site_oxygen = paste(Site, Oxygen)) ->
  DESeq_res_SIP_each_rep_df

# Summarise number of labelled and unlabelled ASVs
DESeq_res_SIP_each_rep_df %>% 
  group_by(Labelled) %>% 
  summarise(n = n()) 
```

<div class="kable-table">

<table>
<thead>
<tr>
<th style="text-align:left;">
Labelled
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Labelled
</td>
<td style="text-align:right;">
520
</td>
</tr>
<tr>
<td style="text-align:left;">
Unlabelled
</td>
<td style="text-align:right;">
1500
</td>
</tr>
<tr>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
922
</td>
</tr>
</tbody>
</table>

</div>

``` r
# detect taxa with NA from DESeq analysis
DESeq_res_SIP_each_rep_df %<>% 
  filter(!is.na(Labelled)) #%>% 
  # pull(Labelled) -> 
  # bad_seqs

# remove NA taxa from PS obj
Ps_obj %>% 
  # prune_taxa(setdiff(taxa_names(Ps_obj), "Seq_2375"), .) %>% # outlier
  prune_taxa(DESeq_res_SIP_each_rep_df$ASV, .) ->
  Ps_obj_SIP4tree_plot


# Remove long name
tax_table(Ps_obj_SIP4tree_plot)[, "Order"] %<>%  str_replace_all(., "Gammaproteobacteria Incertae Sedis", "Incertae Sedis")


taxa2plot <- tibble(rank = c(rep("Class", 3), rep("Phylum", 4)), 
                    subrank = c(rep("Order", 3), rep("Class", 4)), 
                    Taxa2plot = c("Actinobacteria", 
                                  "Alphaproteobacteria", 
                                  "Gammaproteobacteria", 
                                  "Gemmatimonadota",
                                  "Cyanobacteria",
                                  "Bacteroidota",
                                  "Firmicutes"),
                    # lab_rows = c(4, 5, 6, 3, 3, 3, 3),
                    # pwidth = c(5, 6, 8, 3, 3, 3, 3), 
                    # pheight = c(rep(10, 7)),)
                    lab_rows = c(rep(4, 7)),
                    pwidth = c(rep(5, 7)), 
                    pheight = c(rep(8, 7)),)

tree_p_l <- map(seq(nrow(taxa2plot)), 
                ~wrap_ggtree_heatmap(ps_obj = Ps_obj_SIP4tree_plot,
                                     DESeq_res_df = DESeq_res_SIP_each_rep_df,
                                     rank = taxa2plot$rank[.x],
                                     subrank = taxa2plot$subrank[.x],
                                     Taxa2plot = taxa2plot$Taxa2plot[.x],
                                     lab_rows = 4,
                                     pwidth = 5,
                                     pheight = 8))

trees2display <- list.files(path = paste0(fig.path), 
                    pattern = "^Tree_HM_(.*).png$",
                    full.names = TRUE)

knitr::include_graphics(trees2display)
```

<img src="06_Diff_abund_DESeq2_figures//Tree_HM_Actinobacteria.png" width="960" /><img src="06_Diff_abund_DESeq2_figures//Tree_HM_Alphaproteobacteria.png" width="960" /><img src="06_Diff_abund_DESeq2_figures//Tree_HM_Bacteroidota.png" width="960" /><img src="06_Diff_abund_DESeq2_figures//Tree_HM_Cyanobacteria.png" width="960" /><img src="06_Diff_abund_DESeq2_figures//Tree_HM_Firmicutes.png" width="960" /><img src="06_Diff_abund_DESeq2_figures//Tree_HM_Gammaproteobacteria.png" width="960" /><img src="06_Diff_abund_DESeq2_figures//Tree_HM_Gemmatimonadota.png" width="960" />

``` r
all_trees <- ((tree_p_l[[1]] | tree_p_l[[2]] + guides(fill = FALSE) | tree_p_l[[3]] + guides(fill = FALSE) | tree_p_l[[6]] + guides(fill = FALSE)) / (tree_p_l[[4]] + guides(fill = FALSE) | tree_p_l[[5]] + guides(fill = FALSE) | tree_p_l[[7]] + guides(fill = FALSE) | plot_spacer())) + plot_layout(heights = c(2, 1))

save_figure(paste0(fig.path, "all_trees"), 
            all_trees, 
            pwidth = 16, 
            pheight = 18,
            dpi = 900)
```

``` r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
 )
```

<details open>
<summary>
<span title="Click to Expand"> Current session info </span>
</summary>

``` r
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
 date     2021-10-22                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package              * version    date       lib
 ade4                   1.7-18     2021-09-16 [1]
 annotate               1.70.0     2021-05-19 [1]
 AnnotationDbi          1.54.1     2021-06-08 [1]
 ape                    5.5        2021-04-25 [1]
 aplot                  0.1.1      2021-09-22 [1]
 ashr                   2.2-47     2020-02-20 [1]
 assertthat             0.2.1      2019-03-21 [1]
 backports              1.2.1      2020-12-09 [1]
 Biobase              * 2.52.0     2021-05-19 [1]
 BiocGenerics         * 0.38.0     2021-05-19 [1]
 BiocParallel           1.26.2     2021-08-22 [1]
 biomformat             1.20.0     2021-05-19 [1]
 Biostrings           * 2.60.2     2021-08-05 [1]
 bit                    4.0.4      2020-08-04 [1]
 bit64                  4.0.5      2020-08-30 [1]
 bitops                 1.0-7      2021-04-24 [1]
 blob                   1.2.2      2021-07-23 [1]
 broom                  0.7.9      2021-07-27 [1]
 cachem                 1.0.6      2021-08-19 [1]
 cellranger             1.1.0      2016-07-27 [1]
 cli                    3.0.1      2021-07-17 [1]
 clipr                  0.7.1      2020-10-08 [1]
 cluster                2.1.2      2021-04-17 [1]
 codetools              0.2-18     2020-11-04 [1]
 colorspace             2.0-2      2021-06-24 [1]
 crayon                 1.4.1      2021-02-08 [1]
 data.table             1.14.2     2021-09-27 [1]
 DBI                    1.1.1      2021-01-15 [1]
 dbplyr                 2.1.1      2021-04-06 [1]
 DelayedArray           0.18.0     2021-05-19 [1]
 desc                   1.4.0      2021-09-28 [1]
 DESeq2               * 1.32.0     2021-05-19 [1]
 details                0.2.1      2020-01-12 [1]
 digest                 0.6.28     2021-09-23 [1]
 dplyr                * 1.0.7      2021-06-18 [1]
 ellipsis               0.3.2      2021-04-29 [1]
 evaluate               0.14       2019-05-28 [1]
 extrafont            * 0.17       2014-12-08 [1]
 extrafontdb            1.0        2012-06-11 [1]
 fansi                  0.5.0      2021-05-25 [1]
 farver                 2.1.0      2021-02-28 [1]
 fastmap                1.1.0      2021-01-25 [1]
 forcats              * 0.5.1      2021-01-27 [1]
 foreach                1.5.1      2020-10-15 [1]
 fs                     1.5.0      2020-07-31 [1]
 genefilter             1.74.1     2021-10-12 [1]
 geneplotter            1.70.0     2021-05-19 [1]
 generics               0.1.0      2020-10-31 [1]
 GenomeInfoDb         * 1.28.4     2021-09-05 [1]
 GenomeInfoDbData       1.2.6      2021-05-25 [1]
 GenomicRanges        * 1.44.0     2021-05-19 [1]
 ggfun                  0.0.4      2021-09-17 [1]
 ggplot2              * 3.3.5      2021-06-25 [1]
 ggplotify              0.1.0      2021-09-02 [1]
 ggpomological        * 0.1.2      2020-08-13 [1]
 ggrepel              * 0.9.1      2021-01-15 [1]
 ggsci                * 2.9        2018-05-14 [1]
 ggtext               * 0.1.1      2020-12-17 [1]
 ggtree               * 3.0.4      2021-08-22 [1]
 glue                 * 1.4.2      2020-08-27 [1]
 gridExtra              2.3        2017-09-09 [1]
 gridGraphics           0.5-1      2020-12-13 [1]
 gridtext               0.1.4      2020-12-10 [1]
 gtable                 0.3.0      2019-03-25 [1]
 haven                  2.4.3      2021-08-04 [1]
 highr                  0.9        2021-04-16 [1]
 hms                    1.1.1      2021-09-26 [1]
 htmltools              0.5.2      2021-08-25 [1]
 HTSSIP               * 1.4.1      2021-01-15 [1]
 httr                   1.4.2      2020-07-20 [1]
 igraph                 1.2.7      2021-10-15 [1]
 invgamma               1.1        2017-05-07 [1]
 IRanges              * 2.26.0     2021-05-19 [1]
 irlba                  2.3.3      2019-02-05 [1]
 iterators              1.0.13     2020-10-15 [1]
 jsonlite               1.7.2      2020-12-09 [1]
 kableExtra           * 1.3.4      2021-02-20 [1]
 KEGGREST               1.32.0     2021-05-19 [1]
 knitr                  1.36       2021-09-29 [1]
 labeling               0.4.2      2020-10-20 [1]
 lattice              * 0.20-45    2021-09-22 [1]
 lazyeval               0.2.2      2019-03-15 [1]
 lifecycle              1.0.1      2021-09-24 [1]
 locfit                 1.5-9.4    2020-03-25 [1]
 lubridate              1.8.0      2021-10-07 [1]
 magrittr             * 2.0.1      2020-11-17 [1]
 markdown               1.1        2019-08-07 [1]
 MASS                   7.3-54     2021-05-03 [1]
 Matrix                 1.3-4      2021-06-01 [1]
 MatrixGenerics       * 1.4.3      2021-08-26 [1]
 matrixStats          * 0.61.0     2021-09-17 [1]
 memoise                2.0.0      2021-01-26 [1]
 mgcv                   1.8-38     2021-10-06 [1]
 mixsqp                 0.3-43     2020-05-14 [1]
 modelr                 0.1.8      2020-05-19 [1]
 multtest               2.48.0     2021-05-19 [1]
 munsell                0.5.0      2018-06-12 [1]
 nlme                   3.1-153    2021-09-07 [1]
 patchwork            * 1.1.1      2020-12-17 [1]
 permute              * 0.9-5      2019-03-12 [1]
 phyloseq             * 1.36.0     2021-05-19 [1]
 pillar                 1.6.4      2021-10-18 [1]
 pkgconfig              2.0.3      2019-09-22 [1]
 plyr                   1.8.6      2020-03-03 [1]
 png                    0.1-7      2013-12-03 [1]
 purrr                * 0.3.4      2020-04-17 [1]
 R6                     2.5.1      2021-08-19 [1]
 ragg                 * 1.1.3      2021-06-09 [1]
 RColorBrewer         * 1.1-2      2014-12-07 [1]
 Rcpp                   1.0.7      2021-07-07 [1]
 RCurl                  1.98-1.5   2021-09-17 [1]
 readr                * 2.0.2      2021-09-27 [1]
 readxl                 1.3.1      2019-03-13 [1]
 reprex                 2.0.1      2021-08-05 [1]
 reshape2               1.4.4      2020-04-09 [1]
 rhdf5                  2.36.0     2021-05-19 [1]
 rhdf5filters           1.4.0      2021-05-19 [1]
 Rhdf5lib               1.14.2     2021-07-06 [1]
 rlang                  0.4.12     2021-10-18 [1]
 rmarkdown              2.11       2021-09-14 [1]
 rprojroot              2.0.2      2020-11-15 [1]
 RSQLite                2.2.8      2021-08-21 [1]
 rstudioapi             0.13       2020-11-12 [1]
 Rttf2pt1               1.3.9      2021-07-22 [1]
 rvest                  1.0.2      2021-10-16 [1]
 S4Vectors            * 0.30.2     2021-10-03 [1]
 scales               * 1.1.1      2020-05-11 [1]
 sessioninfo            1.1.1      2018-11-05 [1]
 speedyseq            * 0.5.3.9018 2021-08-11 [1]
 SQUAREM                2021.1     2021-01-13 [1]
 stringi                1.7.5      2021-10-04 [1]
 stringr              * 1.4.0      2019-02-10 [1]
 SummarizedExperiment * 1.22.0     2021-05-19 [1]
 survival               3.2-13     2021-08-24 [1]
 svglite              * 2.0.0      2021-02-20 [1]
 systemfonts            1.0.3      2021-10-13 [1]
 textshaping            0.3.6      2021-10-13 [1]
 tibble               * 3.1.5      2021-09-30 [1]
 tidyr                * 1.1.4      2021-09-27 [1]
 tidyselect             1.1.1      2021-04-30 [1]
 tidytree               0.3.5      2021-09-08 [1]
 tidyverse            * 1.3.1      2021-04-15 [1]
 treeio                 1.16.2     2021-08-17 [1]
 truncnorm              1.0-8      2018-02-27 [1]
 tzdb                   0.1.2      2021-07-20 [1]
 utf8                   1.2.2      2021-07-24 [1]
 vctrs                  0.3.8      2021-04-29 [1]
 vegan                * 2.5-7      2020-11-28 [1]
 viridis              * 0.6.2      2021-10-13 [1]
 viridisLite          * 0.4.0      2021-04-13 [1]
 visdat               * 0.6.0.9000 2021-06-07 [1]
 vroom                  1.5.5      2021-09-14 [1]
 webshot                0.5.2      2019-11-22 [1]
 withr                  2.4.2      2021-04-18 [1]
 xfun                   0.27       2021-10-18 [1]
 XML                    3.99-0.8   2021-09-17 [1]
 xml2                   1.3.2      2020-04-23 [1]
 xtable                 1.8-4      2019-04-21 [1]
 XVector              * 0.32.0     2021-05-19 [1]
 yaml                   2.2.1      2020-02-01 [1]
 yulab.utils            0.0.4      2021-10-09 [1]
 zlibbioc               1.38.0     2021-05-19 [1]
 source                                  
 CRAN (R 4.1.1)                          
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 Bioconductor                            
 CRAN (R 4.1.1)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.1)                          
 Github (gadenbuie/ggpomological@69f3815)
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 Github (buckleylab/HTSSIP@29ec56b)      
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 Bioconductor                            
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.0)                          
 Bioconductor                            
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 Bioconductor                            
 Bioconductor                            
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.1)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Github (mikemc/speedyseq@ceb941f)       
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.3)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.1.0)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.3)                          
 Github (ropensci/visdat@8121dfe)        
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.3)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.1.1)                          
 CRAN (R 4.0.2)                          
 CRAN (R 4.0.2)                          
 Bioconductor                            
 CRAN (R 4.0.2)                          
 CRAN (R 4.1.1)                          
 Bioconductor                            

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```

</details>

<br>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-angel_stable_2019" class="csl-entry">

Angel R. Stable Isotope Probing Techniques and Methodological
Considerations Using <sup>15</sup>N. In: Dumont MG, Hernández García M
(eds.). *Stable Isotope Probing: Methods and Protocols*. New York, NY:
Springer New York, 2019, 175–87.

</div>

<div id="ref-love_moderated_2014" class="csl-entry">

Love MI, Huber W, Anders S. Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2. *Genome Biol* 2014;**15**:550.

</div>

<div id="ref-youngblut_htssip_2018" class="csl-entry">

Youngblut ND, Barnett SE, Buckley DH. HTSSIP: An R package for analysis
of high throughput sequencing data from nucleic acid stable isotope
probing (SIP) experiments. *PLOS ONE* 2018;**13**:e0189616.

</div>

</div>
