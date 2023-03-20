Zin SIP
================
Talia Gabay and Roey Angel
2023-03-20

- <a href="#set-parameters-and-load-phyloseq-object"
  id="toc-set-parameters-and-load-phyloseq-object">Set parameters and load
  phyloseq object</a>
- <a href="#subset-phyloseq-object" id="toc-subset-phyloseq-object">subset
  phyloseq object</a>
- <a href="#subsetting-asvs-to-include-only-active-ones"
  id="toc-subsetting-asvs-to-include-only-active-ones">Subsetting ASVs to
  include only active ones</a>
- <a href="#make-a-new-asv-table-of-only-active-sequences"
  id="toc-make-a-new-asv-table-of-only-active-sequences">Make a new ASV
  table of only active sequences</a>
- <a href="#create-a-venn-diagram" id="toc-create-a-venn-diagram">Create a
  venn diagram</a>
- <a href="#create-a-new-phyloseq-object-of-only-active-asvs"
  id="toc-create-a-new-phyloseq-object-of-only-active-asvs">Create a new
  phyloseq object of only active ASVs</a>
- <a href="#beta-diversity-analysis" id="toc-beta-diversity-analysis">Beta
  diversity analysis</a>

### Set parameters and load phyloseq object

``` r
data_path <- "C:/Users/Talia Gabay/OneDrive - post.bgu.ac.il/PhD/SIP experiment/Analyses/Analysis_With_Roey/DADA2_pseudo/"
Proj_name <- "Zin_SIP"
Tree_file <- paste0(data_path, "Tree/DADA2_reps_seq_prev_filt.filtered.align.treefile")
#Load phyloseq object
Ps_file <- paste0(Proj_name, "_filt_wTree.Rds")
Ps_obj <- readRDS(file = paste0(data_path, Ps_file))
```

### subset phyloseq object

First subset phyloseq object created in \_06_Diff_abund_DESeq2 to remove
negative controls, light fractions and unlabeled samples

``` r
Ps_obj %>% 
  subset_samples(Label..18O.== "Labelled" & Density.zone == "Heavy") -> ps_obj_only_labelled 
```

### Subsetting ASVs to include only active ones

In this section, we will choose active sequences that appear in the
table ‘DESeq_res_SIP_each_rep_LFC_sig_df’ We do this for each sample
separately, since in each sample, different sequences were identified as
active.

``` r
#For each samples - look at DESeq table and see which sequences are active in a 
#specific sample 
labelled_seqs_Na_2 <- c("Seq_57", "Seq_235", "Seq_460", "Seq_127", "Seq_623") #select active sequence from table for a spesific sample, in this case "Na_2"
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Natural" & Comparison_pair == "2") %>% # select relevant samples from the ps object based on the treatment and comparison pair
  otu_table() %>% # takes the OTU table from the phyloseq object
  as(., "matrix") %>%
  t() -> mat_Na_2 #creates a matrix
mat_Na_2[!(rownames(mat_Na_2) %in% labelled_seqs_Na_2), ] = 0 # gives zero value for sequences that are not in the labelled list
Na2_summed <- as.matrix(rowSums(mat_Na_2)) #Sum fractions
colnames(Na2_summed) <- c("Na_2") #Rename variable to fit sample names in metadata

#Do the same for the rest of the samples
labelled_seqs_Na_3 <- c("Seq_248",
                        "Seq_29",
                        "Seq_250",
                        "Seq_100",
                        "Seq_1052",
                        "Seq_10",
                        "Seq_8030",
                        "Seq_3",
                        "Seq_19",
                        "Seq_215",
                        "Seq_8",
                        "Seq_31",
                        "Seq_254",
                        "Seq_46",
                        "Seq_1812",
                        "Seq_133",
                        "Seq_39",
                        "Seq_734",
                        "Seq_91",
                        "Seq_49",
                        "Seq_1365",
                        "Seq_175",
                        "Seq_979",
                        "Seq_32",
                        "Seq_108",
                        "Seq_95")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Natural" & Comparison_pair == "3") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Na_3
mat_Na_3[!(rownames(mat_Na_3) %in% labelled_seqs_Na_3), ] = 0
Na3_summed <- as.matrix(rowSums(mat_Na_3)) 
colnames(Na3_summed) <- c("Na_3") 


labelled_seqs_Na_4 <- c("Seq_1",
                        "Seq_9",
                        "Seq_3",
                        "Seq_11",
                        "Seq_5",
                        "Seq_139")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Natural" & Comparison_pair == "4") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Na_4
mat_Na_4[!(rownames(mat_Na_4) %in% labelled_seqs_Na_4), ] = 0
Na4_summed <- as.matrix(rowSums(mat_Na_4)) 
colnames(Na4_summed) <- c("Na_4")

labelled_seqs_Na_5 <- c("Seq_116")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Natural" & Comparison_pair == "5") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Na_5
mat_Na_5[!(rownames(mat_Na_5) %in% labelled_seqs_Na_5), ] = 0
Na5_summed <- as.matrix(rowSums(mat_Na_5)) 
colnames(Na5_summed) <- c("Na_7")

labelled_seqs_Res_1 <- c("Seq_97",
                         "Seq_4",
                         "Seq_42",
                         "Seq_19",
                         "Seq_20",
                         "Seq_1",
                         "Seq_269",
                         "Seq_148",
                         "Seq_9",
                         "Seq_124",
                         "Seq_30",
                         "Seq_187",
                         "Seq_661",
                         "Seq_411",
                         "Seq_562",
                         "Seq_659",
                         "Seq_1183",
                         "Seq_1808",
                         "Seq_7632")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Restored" & Comparison_pair == "1") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Res_1
mat_Res_1[!(rownames(mat_Res_1) %in% labelled_seqs_Res_1), ] = 0
Res1_summed <- as.matrix(rowSums(mat_Res_1)) 
colnames(Res1_summed) <- c("Res_3")


labelled_seqs_Res_2 <- c("Seq_122",
                         "Seq_63",
                         "Seq_4",
                         "Seq_275",
                         "Seq_68")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Restored" & Comparison_pair == "2") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Res_2
mat_Res_2[!(rownames(mat_Res_2) %in% labelled_seqs_Res_2), ] = 0
Res2_summed <- as.matrix(rowSums(mat_Res_2)) 
colnames(Res2_summed) <- c("Res_4")


labelled_seqs_Res_3 <- c("Seq_239",
                         "Seq_504",
                         "Seq_11",
                         "Seq_57",
                         "Seq_5")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Restored" & Comparison_pair == "3") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Res_3
mat_Res_3[!(rownames(mat_Res_3) %in% labelled_seqs_Res_3), ] = 0
Res3_summed <- as.matrix(rowSums(mat_Res_3)) 
colnames(Res3_summed) <- c("Res_5")

labelled_seqs_Res_4 <- c("Seq_24",
                         "Seq_3",
                         "Seq_67",
                         "Seq_73")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Restored" & Comparison_pair == "4") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Res_4
mat_Res_4[!(rownames(mat_Res_4) %in% labelled_seqs_Res_4), ] = 0
Res4_summed <- as.matrix(rowSums(mat_Res_4)) 
colnames(Res4_summed) <- c("Res_7")


labelled_seqs_Res_5 <- c("Seq_15",
                         "Seq_24",
                         "Seq_11",
                         "Seq_3",
                         "Seq_124",
                         "Seq_171",
                         "Seq_138",
                         "Seq_347",
                         "Seq_612",
                         "Seq_10",
                         "Seq_283",
                         "Seq_97",
                         "Seq_562",
                         "Seq_1221",
                         "Seq_9",
                         "Seq_4",
                         "Seq_34",
                         "Seq_16",
                         "Seq_1065",
                         "Seq_62",
                         "Seq_543",
                         "Seq_1616",
                         "Seq_52",
                         "Seq_1002",
                         "Seq_59",
                         "Seq_2396",
                         "Seq_37",
                         "Seq_1149",
                         "Seq_217",
                         "Seq_161",
                         "Seq_67",
                         "Seq_55",
                         "Seq_47",
                         "Seq_462",
                         "Seq_2260")
ps_obj_only_labelled %>%
  subset_samples(Treatment == "Restored" & Comparison_pair == "5") %>%
  otu_table() %>%
  as(., "matrix") %>%
  t() -> mat_Res_5
mat_Res_5[!(rownames(mat_Res_5) %in% labelled_seqs_Res_5), ] = 0
Res5_summed <- as.matrix(rowSums(mat_Res_5)) 
colnames(Res5_summed) <- c("Res_8")
```

### Make a new ASV table of only active sequences

``` r
#Bind all of the tables for each sample together
Only_active_ASV_table <- cbind(Na2_summed, Na3_summed, Na4_summed, Na5_summed, Res1_summed,
                               Res2_summed, Res3_summed, Res4_summed, Res5_summed)
```

### Create a venn diagram

For visualizing which sequences appear in each treatment and which
overlap

``` r
#First, make df of only natural samples
natural_only <- cbind(Na2_summed, Na3_summed, Na4_summed, Na5_summed)
#sum abundances of all samples and change column name
natural_summed <- as.data.frame(rowSums(natural_only))
colnames(natural_summed) <- c("Natural")
#filter out zero values
natural_no_zeros <- filter(natural_summed, Natural > 0)
#Make a column of sequence names
natural_no_zeros <- tibble::rownames_to_column(natural_no_zeros, "Natural_seqs")
#Do the same for restored samples
restored_only <- cbind(Res1_summed, Res2_summed, Res3_summed, Res4_summed, Res5_summed)
restored_summed <- as.data.frame(rowSums(restored_only))
colnames(restored_summed) <- c("Restored")
restored_no_zeros <- filter(restored_summed, Restored > 0)
restored_no_zeros <- tibble::rownames_to_column(restored_no_zeros, "Restored_seqs")
# Make a list of Natural and Restored active seqs to compare fo venn diagram
Venn <- list(Natural_seqs = natural_no_zeros$Natural_seqs, Restored_seqs = restored_no_zeros$Restored_seqs)
#Make venn diagram
Venn_actice <- ggvenn(
  Venn, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5
)
save_figure(paste0(fig.path, "Venn_active"), 
            Venn_actice, 
            pwidth = 10, 
            pheight = 8,
            dpi = 600)
```

    ## Loading required package: ragg

``` r
knitr::include_graphics(paste0(fig.path, "Venn_active", ".png"))
```

<img src="07_Active_community_analysis_figures/Venn_active.png" width="6000" />
\#Make a consolidated metadata table create a metadata table where all
fractions are consolidated to fit the subsetted ASV table

``` r
ps_metadata <- as.data.frame(sample_data(ps_obj_only_labelled))
metadata_active = subset(ps_metadata, select = c(Sample, Treatment, Label..18O., Density.zone, Fraction.no.) )
metadata_active <- subset(metadata_active, Fraction.no. == '2')
#Remove sample number 1 (is an outlier)
metadata_active <- metadata_active[-c(1), ] 
metadata_active <- as.data.frame(metadata_active)
rownames(metadata_active) <- metadata_active$Sample
metadata_active <- metadata_active[ ,-c(1)]
```

### Create a new phyloseq object of only active ASVs

``` r
#Build phyloseq object of only active ASVs
asvtable <- otu_table(Only_active_ASV_table,taxa_are_rows = TRUE)
taxtable <-tax_table(tax_table(Ps_obj))
metadata <- sample_data(metadata_active, sample)
#create phyloseq file
ps_active_only <- phyloseq(asvtable,taxtable,metadata)
#Add phylogenetic tree to phyloseq object
tree_active <- read_tree(Tree_file)
ps_active_only <- merge_phyloseq(ps_active_only,
                                 phy_tree(tree_active))
#Filter out ASVs that don't appear in any sample
ps_active_only_filtered = filter_taxa(ps_active_only, function(x) sum(x) > 1, TRUE)
```

\##Relative abuncances of active community

``` r
#Abundances at phylum level
active.phylum <- ps_active_only_filtered %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at chosen taxonomic level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()                                          # Melt to long format
active.agg <- aggregate(active.phylum$Abundance, FUN=mean, by=list(Treatment=active.phylum$Treatment,Phylum=active.phylum$Phylum))
active.phylum.filtered <- filter(active.agg,active.agg$x > 0.01)
write.csv(active.agg, file= "active phyla table.csv")

#Comparing relative abundance - Natural and Post-mining (Mann-Whitney)
actino <- filter(active.phylum, Phylum == "Actinobacteriota")
with(actino, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.67328, p-value = 0.005418

``` r
with(actino, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.94931, p-value = 0.7322

``` r
bartlett.test(Abundance ~ Treatment, data=actino)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 2.0704, df = 1, p-value = 0.1502

``` r
wilcox.test(Abundance ~ Treatment, data=actino)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 5, p-value = 0.2683
    ## alternative hypothesis: true location shift is not equal to 0

``` r
bacteroi <- filter(active.phylum, Phylum == "Bacteroidota")
with(bacteroi, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.64387, p-value = 0.002074

``` r
with(bacteroi, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.55523, p-value = 0.0001456

``` r
bartlett.test(Abundance ~ Treatment, data=bacteroi)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 1.1814, df = 1, p-value = 0.2771

``` r
wilcox.test(Abundance ~ Treatment, data=bacteroi)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 11, p-value = 0.8933
    ## alternative hypothesis: true location shift is not equal to 0

``` r
cyano <- filter(active.phylum, Phylum == "Cyanobacteria")
with(cyano, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.78213, p-value = 0.07384

``` r
with(cyano, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.90734, p-value = 0.4518

``` r
bartlett.test(Abundance ~ Treatment, data=cyano)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 0.002009, df = 1, p-value = 0.9642

``` r
wilcox.test(Abundance ~ Treatment, data=cyano)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 5, p-value = 0.2622
    ## alternative hypothesis: true location shift is not equal to 0

``` r
firmi <- filter(active.phylum, Phylum == "Firmicutes")
with(firmi, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.75236, p-value = 0.04063

``` r
with(firmi, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.86403, p-value = 0.2431

``` r
bartlett.test(Abundance ~ Treatment, data=firmi)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 5.0879, df = 1, p-value = 0.02409

``` r
wilcox.test(Abundance ~ Treatment, data=firmi)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 10, p-value = 1
    ## alternative hypothesis: true location shift is not equal to 0

``` r
gammati <- filter(active.phylum, Phylum == "Gemmatimonadota")
with(gammati, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.62978, p-value = 0.001241

``` r
with(gammati, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.59375, p-value = 0.0005163

``` r
bartlett.test(Abundance ~ Treatment, data=gammati)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 4.5729, df = 1, p-value = 0.03248

``` r
wilcox.test(Abundance ~ Treatment, data=gammati)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 9.5, p-value = 1
    ## alternative hypothesis: true location shift is not equal to 0

``` r
myxo <- filter(active.phylum, Phylum == "Myxococcota")
with(myxo, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.62978, p-value = 0.001241

``` r
with(myxo, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.76379, p-value = 0.03974

``` r
bartlett.test(Abundance ~ Treatment, data=myxo)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 21.225, df = 1, p-value = 4.084e-06

``` r
wilcox.test(Abundance ~ Treatment, data=myxo)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 7.5, p-value = 0.5605
    ## alternative hypothesis: true location shift is not equal to 0

``` r
patesci <- filter(active.phylum, Phylum == "Patescibacteria")
with(patesci, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.62978, p-value = 0.001241

``` r
with(patesci, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.8929, p-value = 0.3719

``` r
bartlett.test(Abundance ~ Treatment, data=patesci)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 9.7543, df = 1, p-value = 0.001789

``` r
wilcox.test(Abundance ~ Treatment, data=patesci)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 5, p-value = 0.2272
    ## alternative hypothesis: true location shift is not equal to 0

``` r
proteo <- filter(active.phylum, Phylum == "Proteobacteria")
with(proteo, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.82955, p-value = 0.1666

``` r
with(proteo, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.88452, p-value = 0.3303

``` r
bartlett.test(Abundance ~ Treatment, data=proteo)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 1.2899, df = 1, p-value = 0.2561

``` r
wilcox.test(Abundance ~ Treatment, data=proteo)
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  Abundance by Treatment
    ## W = 11, p-value = 0.9048
    ## alternative hypothesis: true location shift is not equal to 0

``` r
unclass <- filter(active.phylum, Phylum == "Unclassified")
with(unclass, shapiro.test(Abundance[Treatment == "Natural"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Natural"]
    ## W = 0.62978, p-value = 0.001241

``` r
with(unclass, shapiro.test(Abundance[Treatment == "Restored"]))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  Abundance[Treatment == "Restored"]
    ## W = 0.55218, p-value = 0.000131

``` r
bartlett.test(Abundance ~ Treatment, data=unclass)
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Abundance by Treatment
    ## Bartlett's K-squared = 9.9357, df = 1, p-value = 0.001621

``` r
wilcox.test(Abundance ~ Treatment, data=unclass)
```

    ## Warning in wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...): cannot
    ## compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Abundance by Treatment
    ## W = 11, p-value = 0.8668
    ## alternative hypothesis: true location shift is not equal to 0

``` r
#Bar plot - all phyla
cl <- colors(distinct = TRUE)
set.seed(15887)
colours <- sample(cl,18)
bp.active <- ggplot(active.agg, aes(x = Treatment, y = x, fill = Phylum)) + 
  geom_bar(stat = "identity", width = 0.4) + theme_light(base_size = 18) + 
  theme_classic() +
  theme(legend.position="bottom", legend.text=element_text(size=20), 
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),) +
  labs(x="", y="Relative abundance (%)") + 
  scale_fill_manual(values=c("#6A3D9A","#A6CEE3", "#1F78B4", 
                        "#FFFF99", "#B2DF8A", "#33A02C", 
                       "#FB9A99", "#E31A1C", "#FDBF6F", 
                       "#FF7F00", "#CAB2D6")) 
                       

save_figure(paste0(fig.path, "Phyla_abundances_active_community"), 
            bp.active, 
            pwidth = 13, 
            pheight = 9,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "Phyla_abundances_active_community", ".png"))
```

<img src="07_Active_community_analysis_figures/Phyla_abundances_active_community.png" width="7800" />

### Beta diversity analysis

we will examine the active community using weighed Unifrac

``` r
#weighed Unifrac
active_unifrac <- UniFrac(ps_active_only_filtered, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE) #Calculate Unifrac
ordu = ordinate(ps_active_only_filtered, "PCoA", "unifrac", weighted=TRUE) #Ordinate
ord_unifrac <- plot_ordination(ps_active_only_filtered, ordu, color="Treatment") #Plot ordination

#Adonis
adonis2(active_unifrac ~ Treatment, method = "unifrac", data = as(sample_data(ps_active_only_filtered), "data.frame"), permutations = 999)
```

<div class="kable-table">

|           |  Df |  SumOfSqs |        R2 |       F | Pr(\>F) |
|:----------|----:|----------:|----------:|--------:|--------:|
| Treatment |   1 | 0.2826436 | 0.1502988 | 1.23819 |   0.238 |
| Residual  |   7 | 1.5979013 | 0.8497012 |      NA |      NA |
| Total     |   8 | 1.8805449 | 1.0000000 |      NA |      NA |

</div>

``` r
#PCoA plot with ggplot
Unifrac_ggp = ord_unifrac + 
  geom_point(size = 4, aes(colour = Treatment)) + 
  stat_ellipse(level=0.95, aes(color = Treatment, group = Treatment)) +
  scale_colour_manual(values = c("#00AFBB", "#FC4E07")) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) 
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.

``` r
Unifrac_ggp
```

![](07_Active_community_analysis_V1.1_files/figure-gfm/ordination%20active-1.png)<!-- -->

``` r
save_figure(paste0(fig.path, "unifrac_PCoA"), 
            Unifrac_ggp, 
            pwidth = 10, 
            pheight = 8,
            dpi = 600)

knitr::include_graphics(paste0(fig.path, "unifrac_PCoA", ".png"))
```

<img src="07_Active_community_analysis_figures/unifrac_PCoA.png" width="6000" />

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

─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.2 (2022-10-31 ucrt)
 os       Windows 10 x64 (build 22621)
 system   x86_64, mingw32
 ui       RTerm
 language (EN)
 collate  English_Israel.utf8
 ctype    English_Israel.utf8
 tz       Asia/Jerusalem
 date     2023-03-20
 pandoc   2.19.2 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────
 ! package          * version    date (UTC) lib source
   ade4               1.7-20     2022-11-01 [1] CRAN (R 4.2.2)
   ape                5.6-2      2022-03-02 [1] CRAN (R 4.2.1)
   assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.2.1)
   backports          1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
   Biobase            2.56.0     2022-04-26 [1] Bioconductor
   BiocGenerics     * 0.42.0     2022-04-26 [1] Bioconductor
   biomformat         1.24.0     2022-04-26 [1] Bioconductor
   Biostrings       * 2.64.1     2022-08-25 [1] Bioconductor
   bitops             1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
   broom              1.0.2      2022-12-15 [1] CRAN (R 4.2.2)
   cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.2.1)
   cli                3.5.0      2022-12-20 [1] CRAN (R 4.2.2)
   clipr              0.8.0      2022-02-22 [1] CRAN (R 4.2.1)
   cluster            2.1.4      2022-08-22 [2] CRAN (R 4.2.2)
   codetools          0.2-18     2020-11-04 [2] CRAN (R 4.2.2)
   colorspace         2.0-3      2022-02-21 [1] CRAN (R 4.2.1)
   crayon             1.5.2      2022-09-29 [1] CRAN (R 4.2.2)
   data.table         1.14.6     2022-11-16 [1] CRAN (R 4.2.2)
   DBI                1.1.3      2022-06-18 [1] CRAN (R 4.2.1)
   dbplyr             2.2.1      2022-06-27 [1] CRAN (R 4.2.1)
   desc               1.4.2      2022-09-08 [1] CRAN (R 4.2.2)
   details            0.3.0      2022-03-27 [1] CRAN (R 4.2.1)
   digest             0.6.31     2022-12-11 [1] CRAN (R 4.2.2)
   dplyr            * 1.0.10     2022-09-01 [1] CRAN (R 4.2.1)
   ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.2.1)
   evaluate           0.19       2022-12-13 [1] CRAN (R 4.2.2)
   extrafont        * 0.18       2022-04-12 [1] CRAN (R 4.2.1)
   extrafontdb        1.0        2012-06-11 [1] CRAN (R 4.2.0)
   fansi              1.0.3      2022-03-24 [1] CRAN (R 4.2.1)
   farver             2.1.1      2022-07-06 [1] CRAN (R 4.2.1)
   fastmap            1.1.0      2021-01-25 [1] CRAN (R 4.2.1)
   forcats          * 0.5.2      2022-08-19 [1] CRAN (R 4.2.1)
   foreach            1.5.2      2022-02-02 [1] CRAN (R 4.2.1)
   fs                 1.5.2      2021-12-08 [1] CRAN (R 4.2.1)
   gargle             1.2.1      2022-09-08 [1] CRAN (R 4.2.1)
   generics           0.1.3      2022-07-05 [1] CRAN (R 4.2.1)
   GenomeInfoDb     * 1.32.3     2022-08-11 [1] Bioconductor
   GenomeInfoDbData   1.2.8      2022-09-06 [1] Bioconductor
   ggplot2          * 3.4.0      2022-11-04 [1] CRAN (R 4.2.2)
   ggsci            * 2.9        2018-05-14 [1] CRAN (R 4.2.1)
   ggtext           * 0.1.2      2022-09-16 [1] CRAN (R 4.2.2)
   ggvenn           * 0.1.9      2022-09-16 [1] Github (yanlinlin82/ggvenn@b7ff54b)
   glue             * 1.6.2      2022-02-24 [1] CRAN (R 4.2.1)
   googledrive        2.0.0      2021-07-08 [1] CRAN (R 4.2.1)
   googlesheets4      1.0.1      2022-08-13 [1] CRAN (R 4.2.1)
   gridExtra          2.3        2017-09-09 [1] CRAN (R 4.2.1)
   gridtext           0.1.5      2022-09-16 [1] CRAN (R 4.2.2)
   gtable             0.3.1      2022-09-01 [1] CRAN (R 4.2.1)
   haven              2.5.1      2022-08-22 [1] CRAN (R 4.2.1)
   highr              0.10       2022-12-22 [1] CRAN (R 4.2.2)
   hms                1.1.2      2022-08-19 [1] CRAN (R 4.2.1)
   htmltools          0.5.4      2022-12-07 [1] CRAN (R 4.2.2)
   httr               1.4.4      2022-08-17 [1] CRAN (R 4.2.1)
   igraph             1.3.5      2022-09-22 [1] CRAN (R 4.2.2)
   IRanges          * 2.30.1     2022-08-25 [1] Bioconductor
   iterators          1.0.14     2022-02-05 [1] CRAN (R 4.2.1)
   jsonlite           1.8.4      2022-12-06 [1] CRAN (R 4.2.2)
   knitr              1.41       2022-11-18 [1] RSPM (R 4.2.0)
   labeling           0.4.2      2020-10-20 [1] CRAN (R 4.2.0)
   lattice          * 0.20-45    2021-09-22 [2] CRAN (R 4.2.2)
   lifecycle          1.0.3      2022-10-07 [1] CRAN (R 4.2.2)
   lubridate          1.9.0      2022-11-06 [1] CRAN (R 4.2.2)
   magrittr         * 2.0.3      2022-03-30 [1] CRAN (R 4.2.1)
   MASS               7.3-58.1   2022-08-03 [2] CRAN (R 4.2.2)
   Matrix             1.5-3      2022-11-11 [1] CRAN (R 4.2.2)
   mgcv               1.8-41     2022-10-21 [2] CRAN (R 4.2.2)
   modelr             0.1.10     2022-11-11 [1] CRAN (R 4.2.2)
   multtest           2.52.0     2022-04-26 [1] Bioconductor
   munsell            0.5.0      2018-06-12 [1] CRAN (R 4.2.1)
   nlme               3.1-160    2022-10-10 [2] CRAN (R 4.2.2)
   patchwork        * 1.1.2      2022-08-19 [1] CRAN (R 4.2.1)
   permute          * 0.9-7      2022-01-27 [1] CRAN (R 4.2.1)
   phyloseq         * 1.40.0     2022-04-26 [1] Bioconductor
   pillar             1.8.1      2022-08-19 [1] CRAN (R 4.2.1)
   pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.2.1)
   plyr               1.8.8      2022-11-11 [1] CRAN (R 4.2.2)
   png                0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
   purrr            * 1.0.0      2022-12-20 [1] CRAN (R 4.2.2)
   R6                 2.5.1      2021-08-19 [1] CRAN (R 4.2.1)
   ragg             * 1.2.4      2022-10-24 [1] CRAN (R 4.2.2)
   RColorBrewer     * 1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
   Rcpp               1.0.9      2022-07-08 [1] CRAN (R 4.2.1)
   RCurl              1.98-1.9   2022-10-03 [1] CRAN (R 4.2.2)
   readr            * 2.1.3      2022-10-01 [1] CRAN (R 4.2.2)
   readxl             1.4.1      2022-08-17 [1] CRAN (R 4.2.1)
   reprex             2.0.2      2022-08-17 [1] CRAN (R 4.2.1)
   reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.2.1)
   rhdf5              2.40.0     2022-04-26 [1] Bioconductor
 D rhdf5filters       1.8.0      2022-04-26 [1] Bioconductor
   Rhdf5lib           1.18.2     2022-05-15 [1] Bioconductor
   rlang              1.0.6      2022-09-24 [1] CRAN (R 4.2.2)
   rmarkdown          2.19       2022-12-15 [1] CRAN (R 4.2.2)
   rprojroot          2.0.3      2022-04-02 [1] CRAN (R 4.2.1)
   rstudioapi         0.14       2022-08-22 [1] CRAN (R 4.2.1)
   Rttf2pt1           1.3.11     2022-10-08 [1] CRAN (R 4.2.1)
   rvest              1.0.3      2022-08-19 [1] CRAN (R 4.2.1)
   S4Vectors        * 0.34.0     2022-04-26 [1] Bioconductor
   scales             1.2.1      2022-08-20 [1] CRAN (R 4.2.2)
   sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.2.1)
   speedyseq        * 0.5.3.9018 2022-09-09 [1] Github (mikemc/speedyseq@ceb941f)
   stringi            1.7.8      2022-07-11 [1] CRAN (R 4.2.1)
   stringr          * 1.5.0      2022-12-02 [1] CRAN (R 4.2.2)
   survival           3.4-0      2022-08-09 [2] CRAN (R 4.2.2)
   svglite          * 2.1.0      2022-02-03 [1] CRAN (R 4.2.1)
   systemfonts        1.0.4      2022-02-11 [1] CRAN (R 4.2.1)
   textshaping        0.3.6      2021-10-13 [1] CRAN (R 4.2.1)
   tibble           * 3.1.8      2022-07-22 [1] CRAN (R 4.2.1)
   tidyr            * 1.2.1      2022-09-08 [1] CRAN (R 4.2.2)
   tidyselect         1.2.0      2022-10-10 [1] CRAN (R 4.2.2)
   tidyverse        * 1.3.2      2022-07-18 [1] CRAN (R 4.2.2)
   timechange         0.1.1      2022-11-04 [1] CRAN (R 4.2.2)
   tzdb               0.3.0      2022-03-28 [1] CRAN (R 4.2.1)
   utf8               1.2.2      2021-07-24 [1] CRAN (R 4.2.1)
   vctrs              0.5.1      2022-11-16 [1] CRAN (R 4.2.2)
   vegan            * 2.6-4      2022-10-11 [1] CRAN (R 4.2.2)
   viridis          * 0.6.2      2021-10-13 [1] CRAN (R 4.2.1)
   viridisLite      * 0.4.1      2022-08-22 [1] CRAN (R 4.2.1)
   visdat           * 0.5.3      2019-02-15 [1] CRAN (R 4.2.1)
   withr              2.5.0      2022-03-03 [1] CRAN (R 4.2.1)
   xfun               0.36       2022-12-21 [1] CRAN (R 4.2.2)
   xml2               1.3.3      2021-11-30 [1] CRAN (R 4.2.1)
   XVector          * 0.36.0     2022-04-26 [1] Bioconductor
   yaml               2.3.6      2022-10-18 [1] CRAN (R 4.2.2)
   zlibbioc           1.42.0     2022-04-26 [1] Bioconductor

 [1] C:/Users/Talia Gabay/AppData/Local/R/win-library/4.2
 [2] C:/Program Files/R/R-4.2.2/library

 D ── DLL MD5 mismatch, broken installation.

──────────────────────────────────────────────────────────────────────────────
```

</details>

<br>
