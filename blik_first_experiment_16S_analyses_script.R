---
title: "Land use drives prokaryotic community composition of directly adjacent grasslands"
author: "Martinez-Cuesta, R."
date: "05.09.2024"
---
###DADA2 processing of 16S rRNA gene amplicon sequences:
library(dada2)
library(plotly)

fnFs = sort(list.files("blik_raw_samples", pattern="_R1_", full.names = TRUE))
fnRs = sort(list.files("blik_raw_samples", pattern="_R2_", full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#plot quality of forward reads
forwplot = ggplotly(plotQualityProfile(fnFs[1:length(fnFs)], aggregate = TRUE) + geom_hline(yintercept=c(15,25,35), color=c("red","blue","green"), size=0.5), width =750)
forwplot

#plot quality of reverse reads
revqplot = ggplotly(plotQualityProfile(fnRs[1:length(fnRs)], aggregate=TRUE) + geom_hline(yintercept=c(15,25,35), color=c("red","blue","green"), size=0.5), width =800)
revqplot

#create folder for filtered reads and determine sample names
filtFs = file.path("blik_raw_samples/filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs = file.path("blik_raw_samples/filtered", paste0(sample.names, "_R2_filt.fastq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

#filtering of reads based on quality plots
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), truncLen = c(272,205),
                     maxN=0, maxEE=c(2,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 

#learn error rates
errF = learnErrors(filtFs, multithread=TRUE)
errR = learnErrors(filtRs, multithread=TRUE)

#plot the error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication step
derepFs = derepFastq(filtFs[exists], verbose=TRUE)
derepRs = derepFastq(filtRs[exists2], verbose=TRUE)
names(derepFs) = sample.names
names(derepRs) = sample.names

#sample inference
dadaFs = dada(filtFs, err=errF, multithread=TRUE)
dadaRs = dada(filtRs, err=errF, multithread=TRUE)

print("Forward Reads")
dadaFs[[1]]
print("Reverse Reads")
dadaRs[[1]]

#merge forward and reverse reads
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#construct sequence table and inspect distribution of sequence lengths
seqtab = makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))#

#filter for expected sequenced length and inspect the output
seqtab = seqtab[,nchar(colnames(seqtab)) %in% 401:430]
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")

#remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)*100

#inspect output of DADA2 pipeline along the previous steps
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample.names

#export output to QIIME2
library(openssl)
seqtab.nochimmd5 = seqtab.nochim
sequences = colnames(seqtab.nochimmd5)
sequencesmd5 = md5(sequences)
colnames(seqtab.nochimmd5) = sequencesmd5
dir.create("qiime2")
write.table(t(seqtab.nochimmd5), "qiime2/seqtab-nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, fout='qiime2/rep-seqs.fna', ids=sequencesmd5)
write.table(t(track), "qiime2/stats.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
````

```{bash}
  
#Taxonomic assignation with QIIME2 in linux terminal
qiime tools import \
--input-path rep-seqs.fna \
--type "FeatureData[Sequence]" \
--output-path rep-seqs.qza

echo -n "#OTU Table" | cat - seqtab-nochim.txt > biom-table.txt

biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5

qiime tools import \
--input-path table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path table.qza

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file metadata.txt

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

#phylogeny part
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#after creating the classifier
qiime feature-classifier classify-sklearn --i-reads rep-seqs.qza \
--i-classifier classifier.qza \
--p-n-jobs 6 \
--output-dir taxa
```

```{r}
  
###export results into a phyloseq object and processing
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(ggvenn)
library(microbiome)
library(microbiomeMarker)
library(readr)
library(decontam)

metadata = read_tsv("metadata.tsv")
tree = read_qza("taxa/rooted-tree.qza")
asvs = read_qza(ASVs)$data %>% as.data.frame()
taxtable =  read_qza("taxa/classification.qza")$data %>%
  as.tibble() %>% separate(Taxon,
                           sep = ";",
                           c("Domain",
                             "Phylum",
                             "Class",
                             "Order",
                             "Family",
                             "Genus",
                             "Species"),
                           fill = "right",
                           extra = "drop") %>%
  as.matrix() %>%
  gsub(x = .,
       pattern = "[A-Z]_[0-9]__",
       replacement = "") %>% 
  as.data.frame() %>%
  map_df(~ gsub(
    pattern = "metagenome|uncultured|unidentified|Unknown",
    replacement = NA,
    .x)) %>%
  mutate_if(is_character, str_trim) %>%
  column_to_rownames("Feature.ID") %>%
  mutate(Domain = ifelse(is.na(Domain),
                         "U. Domain",
                         Domain),
         Phylum = coalesce(Phylum,
                           ifelse(grepl("^U.", Domain),
                                  Domain,
                                  paste("U.", Domain))),
         Class = coalesce(Class,
                          ifelse(grepl("^U.", Phylum),
                                 Phylum,
                                 paste("U.", Phylum))),
         Order = coalesce(Order,
                          ifelse(grepl("^U.", Class),
                                 Class,
                                 paste("U.", Class))),
         Family = coalesce(Family,
                           ifelse(grepl("^U.", Order),
                                  Order,
                                  paste("U.", Order))),
         Genus = coalesce(Genus,
                          ifelse(grepl("^U.", Family),
                                 Family,
                                 paste("U.", Family))),
         Species = coalesce(Species,
                            ifelse(grepl("^U.", Genus),
                                   Genus,
                                   paste("U.", Genus))))
#create the PHYLOSEQ object
physeq_blik = phyloseq(otu_table(asvs, taxa_are_rows = T),
                           phy_tree(tree$data),
                           tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% as.matrix()),
                           sample_data(metadata %>% as.data.frame() %>% column_to_rownames("#SampleID")))

##remove ASVs detected in the negative controls plus sequences not belonging to bacteria nor archaea
sample_data(physeq_blik)$is.neg <- sample_data(physeq_blik)$Field_type == "negative"

#identify the contaminants
contamdf.prev05 = isContaminant(physeq_blik, method="prevalence", neg="is.neg", threshold=0.05)
table(contamdf.prev05$contaminant)

#remove the contaminants
physeq_blik_filt = prune_taxa(!contamdf.prev05$contaminant, physeq_blik)
physeq_blik_filt = physeq_blik_filt %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) )
physeq_blik_filt = subset_samples(physeq_blik_filt, is.neg !=TRUE)
physeq_blik_2_filt = subset_samples(physeq_blik_filt, site !="Undetermined")

sample_data(physeq_blik_filt)$is.neg = NULL

physeq_blik_filt = prune_taxa(taxa_sums(physeq_blik_filt) > 0, physeq_blik_filt)

save(physeq_blik_filt,
     file = paste("phy_blik_1st_exp_filt",
                  ".RData",
                  sep = "")) #save the object

##obtain rarefaction curves
tab = otu_table(physeq_blik_filt) %>% t()
rare = rarecurve(tab, step=10000, lwd=1, ylab="ASVs",  label=F)

physeq_blik_filt = normalize(physeq_blik_filt, method = "CSS")

##obtain the number of observed ASVs per sample
observed_ASVs = alpha(physeq_blik_filt, index = "observed") %>% as.data.frame()

##bind with metadata information
observed_asv_merged = cbind(sample_data(physeq_blik_filt)$Sampling_time, observed_ASVs)
observed_asv_merged = cbind(sample_data(physeq_blik_filt)$Field_type, observed_asv_merged)
observed_asv_merged = cbind(sample_data(physeq_blik_filt)$Soil, observed_asv_merged)

##plotting the results

##add comparisons for Wilcoxon's rank sum tests

my_comparisons = list(c("G1", "C1"),
                      c("G2", "C2"),
                      c("G3", "C3"),
                      c("FP-G", "FP-C"))


my_colors = c("FP-C" = "#FFCC33", "C" = "forestgreen", "FP-G" = "purple", "G" = "royalblue1")

custom_order = c("G1", "G2", "G3", "FP-G", "C1", "C2", "C3", "FP-C") #order for the x axis

st.order = c("March", "June", "November") #order for the facet

observed_asv_merged$Sampling_time = factor(observed_asv_merged$Sampling_time, levels = st.order) #order of the facet

observed.asvs.plot = ggplot(observed_asv_merged, aes(x = factor(Soil, levels = custom_order), y = observed, color = Field_type)) +
  geom_boxplot() +
  geom_jitter(aes(color = Field_type), size = 2, alpha = 0.5) +
  facet_grid("Sampling_time" ,drop = FALSE) +
  theme(ylab(label = "Total ASVs number")) + guides(color = FALSE)

observed.asvs.plot + theme_pubclean() + 
  scale_color_manual(values = my_colors) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", p.adjust.method = "fdr",
                     comparisons = my_comparisons, size = 4, 
                     vjust = 1.5, hide.ns = T) #add wilcox's statistics to the plot

##Plot results of 16S rRNA gene qPCRs
biomass.blik1 = read.csv("16SqPCR_first_exp_0523.csv", header = T) #load the .csv with the results

facet_order = c("March", "June", "November")
plot.biomass.blik1 = ggplot(biomass.blik1, aes(x = factor(Soil, levels = custom_order), y = copies.soil, fill = Field_type)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  facet_grid( ~ factor(Sampling_time, levels = facet_order), scales = "free_y") +
  scale_fill_manual(values = my_colors) +  # Set custom colors
  ggpubr::theme_pubclean() +
ggpubr::stat_compare_means(comparisons = list(c("G1", "C1"), c("G2", "C2"), c("G3", "C3"), c("FP-G", "FP-C")), label = "p.adjust", hide.ns = T, method = "wilcox.test") +
labs(y = "16S rRNA gene copies per g-1 dry soil", x = "Sampling points") + guides(fill = guide_legend(title = "Site"))

##Ordination with PCoAs
#subsetting by sampling time for the stats
blik.sp = subset_samples(physeq_blik_filt, Sampling_time == "March")
blik.su = subset_samples(physeq_blik_filt, Sampling_time == "June")
blik.au = subset_samples(physeq_blik_filt, Sampling_time == "November")

##PERMOVA tests for differences between field types in microbial composition
bray = distance(blik.<sampling_time>, method = "bray")
sample.df = data.frame(sample_data(blik.<sampling_time))
adonis2(bray ~ Field_type, data = sample.df)

##Plotting the PCoAs
ord_plot = plot_ordination(physeq_blik_filt, ordination, type="samples", color="Field_type",
                         shape = "catena_level") +
  geom_point(size=4) + scale_color_manual(values = my_colors) + facet_grid(~ factor(Sampling_time, levels = facet_order)) +
  stat_ellipse(type = "t") + ggpubr::theme_pubr()

ord_plot + 
  annotate("text", x = 0.6, y = 0.7, label = "P-value < 0.001 ***", 
           data = data.frame(Time = "Time1"), size = 2, color = "black") +
  annotate("text", x = 0.6, y = 0.7, label = "P-value < 0.001 ***", 
           data = data.frame(Time = "Time2"), size = 2, color = "black") +
  annotate("text", x = 0.6, y = 0.7, label = "P-value < 0.001 ***", 
           data = data.frame(Time = "Time3"), size = 2, color = "black") #to include the results of the PERMANOVA tests

##Core microbiome pre-processing to filter for ASVs that only occur in 60% of each sample collection (3 out of 5 replicates)
melted = physeq_blik_filt %>%
  psmelt()

filtered =
  melted %>%
  group_by(Soil, Field_type, Sampling_time, OTU) %>%
  summarise(prev = sum(Abundance > 0)/n()) %>%
  filter(prev >= 0.6)

shared_asv = filtered$OTU

#Remove the ASVs that are did not pass the filtering step
phy60 = physeq_blik_filt %>% prune_taxa(shared_asv, .)

#subset for sampling time
blik_sp = subset_samples(phy60, Sampling_time == "March")
blik_su = subset_samples(phy60, Sampling_time == "June")
blik_au = subset_samples(phy60, Sampling_time == "November")

#determine parameters for the core
detection_core = "0.01"
prevalence_core = "0.6"
venn_list = list()

Field_type = c("G", "C", "FP-G", "FP-C")

for (x in Field_type){
  core = core_members(subset_samples(blik_<sampling_time>, Field_type == x),
                      detection = detection_core,
                      prevalence = prevalence_core)
  venn_list[[x]] = core
}

venn_sp = ggvenn(venn_list, fill_color = my_colors, stroke_linetype = 0.8) + ggplot2::ggtitle(label = "March")
venn_su = ggvenn(venn_list, fill_color = my_colors, stroke_linetype = 0.8) + ggplot2::ggtitle(label = "June")
venn_au = ggvenn(venn_list, fill_color = my_colors, stroke_linetype = 0.8) + ggplot2::ggtitle(label = "November")

library(patchwork)
library(gridExtra)

grid.arrange(venn_sp, venn_su, venn_au, ncol = 3) #combine the Venn diagrams

##LDA analysis with edgeR comparing the two floodplains and the grassland with the cropland using the three seasonal subsets

sub.fps = subset_samples(physeq_blik_filt, Field_type %in% c("FP-C","FP-G")) #subset for only the floodplain samples
sub.gc = subset_samples(physeq_blik_filt, Field_type %in% c("C","G"))

#sub-setting for sampling time
gc.sp = subset_samples(sub.gc, Sampling_time == "March")
gc.su = amp_subset_samples(sub.gc, Sampling_time == "June")
gc.au = amp_subset_samples(sub.gc, Sampling_time == "November")

fps.sp = subset_samples(sub.fps, Sampling_time == "March")
fps.su = amp_subset_samples(sub.fps, Sampling_time == "June")
fps.au = amp_subset_samples(sub.fps, Sampling_time == "November")

#filter out the taxa with all 0s after the subset
gc.sp = prune_taxa(taxa_sums(gc.sp) > 0, gc.sp) 
gc.su = prune_taxa(taxa_sums(gc.su) > 0, gc.su) 
gc.au = prune_taxa(taxa_sums(gc.au) > 0, gc.au)

fps.sp = prune_taxa(taxa_sums(fps.sp) > 0, fps.sp) 
fps.su = prune_taxa(taxa_sums(fps.su) > 0, fps.su) 
fps.au = prune_taxa(taxa_sums(fps.au) > 0, fps.au)

edgeR_gc_<sampling_time> = run_edger(
  gc.<st>,
  group = "Field_type",
  pvalue_cutoff = 0.001,
  p_adjust = "fdr"
  norm = "none")

edgeR_fp_<sampling_time> = run_edger(
  fp.<st>,
  group = "Field_type",
  pvalue_cutoff = 0.001,
  p_adjust = "fdr"
  norm = "none")

fc_cutoff  = 2

#remove duplicates
edgeR_gc_<sampling_time> = edgeR_gc_<sampling_time>[!duplicated(edgeR_gc_<sampling_time>$padj, fromLast = TRUE), ]
edgeR_fp_<sampling_time> = edgeR_fp_<sampling_time>[!duplicated(edgeR_fp_<sampling_time>$padj, fromLast = TRUE), ]

#set the threshold for log2FoldChange
edge.fp.<sampling_time> = subset(edgeR_fp_<sampling_time>, abs(ef_logFC) > fc_cutoff)
edge.gc.<sampling_time> = subset(edgeR_gc_<sampling_time>, abs(ef_logFC) > fc_cutoff)

#order the markers for log2FoldChange
edge.gc.sp$feature = fct_reorder(edge.gc.sp$feature, edge.gc.sp$ef_logFC, .desc = TRUE)
edge.gc.su$feature = fct_reorder(edge.gc.su$feature, edge.gc.su$ef_logFC, .desc = TRUE)
edge.gc.au$feature = fct_reorder(edge.gc.au$feature, edge.gc.au$ef_logFC, .desc = TRUE)

colors2 = c("forestgreen", "royalblue1")

#plotting the results of the cropland and grassland
bp.gc.sp = ggplot(edge.gc.sp, aes(x = ef_logFC, y = feature)) + 
  geom_point(aes(size = pvalue, fill = factor(enrich_group)), shape = 21) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        panel.background = element_blank()) +
  ylab("Taxa") +
  xlab("log2 Fold Change") +
  scale_fill_manual(values = colors2) +
  scale_size_continuous(name = "P.value", breaks = c(1e-05, 2e-05, 3e-05, 4e-05)) +
  theme_bw() +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Field type")) + 
  labs(title = "March")

bp.gc.su = ggplot(edge.gc.su, aes(x = ef_logFC, y = feature)) + 
  geom_point(aes(size = pvalue, fill = factor(enrich_group)), shape = 21) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        panel.background = element_blank()) +
  ylab("Taxa") +
  xlab("log2 Fold Change") +
  scale_fill_manual(values = colors2) +
  scale_size_continuous(name = "P.value", breaks = c(1e-05, 2e-05, 3e-05, 4e-05)) +
  theme_bw() +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Field type")) + 
  labs(title = "June")

bp.gc.au = ggplot(edge.gc.au, aes(x = ef_logFC, y = feature)) + 
  geom_point(aes(size = pvalue, fill = factor(enrich_group)), shape = 21) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        panel.background = element_blank()) +
  ylab("Taxa") +
  xlab("log2 Fold Change") +
  scale_fill_manual(values = colors2) +
  scale_size_continuous(name = "P.value", breaks = c(1e-05, 2e-05, 3e-05, 4e-05)) +
  theme_bw() +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Field type")) + 
  labs(title = "November")

#combining the plots
all.gc = ggpubr::ggarrange(bp.gc.sp, bp.gc.su, bp.gc.au, ncol = 3, align = "hv", common.legend = T)

#same for the floodplain samples
edge.fp.sp$feature = fct_reorder(edge.fp.sp$feature, edge.fp.sp$ef_logFC, .desc = TRUE)
edge.fp.su$feature = fct_reorder(edge.fp.su$feature, edge.fp.su$ef_logFC, .desc = TRUE)
edge.fp.au$feature = fct_reorder(edge.fp.au$feature, edge.fp.au$ef_logFC, .desc = TRUE)

colors = c("#FFCC33", "purple")

bp.fp.sp = ggplot(edge.fp.sp, aes(x = ef_logFC, y = feature)) + 
  geom_point(aes(size = pvalue, fill = factor(enrich_group)), shape = 21) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        panel.background = element_blank()) +
  ylab("Taxa") +
  xlab("log2 Fold Change") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(name = "P.value", breaks = c(2e-04, 4e-04, 6e-04, 8e-04)) +
  theme_bw() +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Field type")) + 
  labs(title = "March")

bp.fp.su = ggplot(edge.fp.su, aes(x = ef_logFC, y = feature)) + 
  geom_point(aes(size = pvalue, fill = factor(enrich_group)), shape = 21) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        panel.background = element_blank()) +
  ylab("Taxa") +
  xlab("log2 Fold Change") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(name = "P.value", breaks = c(2e-04, 4e-04, 6e-04, 8e-04)) +
  theme_bw() +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Field type")) + 
  labs(title = "June")

bp.fp.au = ggplot(edge.fp.au, aes(x = ef_logFC, y = feature)) + 
  geom_point(aes(size = pvalue, fill = factor(enrich_group)), shape = 21) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        panel.background = element_blank()) +
  ylab("Taxa") +
  xlab("log2 Fold Change") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(name = "P.value", breaks = c(2e-04, 4e-04, 6e-04, 8e-04)) +
  theme_bw() +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "Field type")) + 
  labs(title = "November")

#combining the three plots of the floodplains
all.fp = ggpubr::ggarrange(bp.fp.sp, bp.fp.su, bp.fp.au, ncol = 3, align = "hv", common.legend = T)
#combining all the plots
all.gc.fp = ggpubr::ggarrange(all.gc, all.fp, nrow = 2, align = "hv", common.legend = F, labels = c("A", "B"))

##Microbial Network Analysis
library(NetCoMi)

fpg = phyloseq::subset_samples(physeq_blik_filt_0, Field_type == "FP-G")
fpc = phyloseq::subset_samples(physeq_blik_filt_0, Field_type == "FP-C")

meta.fpc = sample_data(fpc)
meta.fpg = sample_data(fpg)

fpc = prune_taxa(taxa_sums(fpc) > 0, fpc)
fpg = prune_taxa(taxa_sums(fpg) > 0, fpg)

fpc.gen = tax_glom(fpc, taxrank = "Genus")
fpg.gen = tax_glom(fpg, taxrank = "Genus")

net.floodplains = netConstruct(fpc.gen, fpg.gen,
                       taxRank = "Genus",
                       measure = "spearman",
                       normMethod = "clr",
                       zeroMethod = "multRepl",
                       sparsMethod = "t-test",
                       adjust = "lfdr",
                       alpha = 0.001)

netprops = netAnalyze(net.floodplains, clustMethod = "cluster_fast_greedy")
netprops

##plot the two networks
plot(netprops, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     #featVecCol = phyla,
     labelScale = FALSE,
     groupNames = c("FP-C", "FP-G"),
     cexLabels = 0.3,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "black")

summary(net.floodplains, groupNames = c("FP-C", "FP-G"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 5) #obtain table with the main network properties
