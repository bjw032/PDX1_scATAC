library(monocle3)
library(Signac)
library(Seurat)
library(cicero)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(rtracklayer)
library(harmony)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)


setwd("~/Hislet_aggr/")
counts <- Read10X_h5(filename = "~/Hislet_aggr/outs/filtered_peak_bc_matrix.h5")
# metadata <- read.csv(
#   file = "~/Hislet_aggr/outs/singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )
# 
# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c(":", "-"),
#   fragments = '~/Hislet_aggr/outs/fragments.tsv.gz',
#   genome = "hg38",
#   min.cells = 10,
#   min.features = 200
# )
# 
# #
# # # extract gene annotations from EnsDb
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# 
# # change to UCSC style since the data was mapped to hg38
# seqlevelsStyle(annotation) <- "UCSC"
# 
# genome(annotation) <- "hg38"
# 
# 
# saveRDS(annotation,"~/Hislet_aggr/outs/annotations.rds")
# # add the gene information to the object
# 
# scATAC <- CreateSeuratObject(
#   counts = chrom_assay,
#   assay = "peaks",
#   meta.data = metadata
# )
# 
# Annotation(scATAC) <- annotation
# 
# saveRDS(scATAC,"~/Hislet_aggr/outs/scATAC.rds")
# 
# # scATAC <- readRDS("~/Hislet_aggr/outs/scATAC.rds")
# 
# # compute nucleosome signal score per cell
# scATAC <- NucleosomeSignal(object = scATAC)
# 
# # compute TSS enrichment score per cell
# scATAC <- TSSEnrichment(object = scATAC, fast = FALSE)
# 
# # add blacklist ratio and fraction of reads in peaks
# scATAC$pct_reads_in_peaks <- scATAC$peak_region_fragments / scATAC$passed_filters * 100
# scATAC$blacklist_ratio <- scATAC$blacklist_region_fragments / scATAC$peak_region_fragments
# 
# scATAC$nucleosome_group <- ifelse(scATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# FragmentHistogram(object = scATAC, group.by = 'nucleosome_group')
# 
# scATAC.filt <- subset(
#   x = scATAC,
#   subset = peak_region_fragments > 3000 &
#     peak_region_fragments < 20000 &
#     pct_reads_in_peaks > 15 &
#     #blacklist_ratio < 0.05 &
#     nucleosome_signal < 4 &
#     TSS.enrichment > 2
# )
# 
# scATAC.filt <- RunTFIDF(scATAC.filt)
# scATAC.filt <- FindTopFeatures(scATAC.filt, min.cutoff = 'q0')
# scATAC.filt <- RunSVD(scATAC.filt)
# scATAC.filt <- RunUMAP(object = scATAC.filt, reduction = 'lsi', dims = 2:30)
# scATAC.filt <- FindNeighbors(object = scATAC.filt, reduction = 'lsi', dims = 2:30)
# scATAC.filt <- FindClusters(object = scATAC.filt, verbose = FALSE, algorithm = 3)
# scATAC.filt$sample <- gsub(".*-","",colnames(scATAC.filt))
# DimPlot(object = scATAC.filt, label = TRUE,split.by = "sample")# + NoLegend()
# table(scATAC.filt$sample)
# 
# 
# ### Gene labels and stuff
# # scATAC.filt <- readRDS("~/Hislet_aggr/outs/scATAC.filt.rds")
# gene.activities <- GeneActivity(scATAC.filt)
# 
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# scATAC.filt[['RNA']] <- CreateAssayObject(counts = gene.activities)
# scATAC.filt <- NormalizeData(
#   object = scATAC.filt,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(scATAC.filt$nCount_RNA)
# )
# 
# DefaultAssay(scATAC.filt) <- 'peaks'
# 
# # saveRDS(scATAC.filt,"~/Hislet_aggr/outs/scATAC.filt.rds")
# 
# library(harmony)
# scATAC_hm <- RunHarmony(
#   object = scATAC.filt,
#   group.by.vars = 'sample',
#   reduction = 'lsi',
#   assay.use = 'peaks',
#   project.dim = FALSE
# )
# scATAC_hm <- scATAC_hm %>%
#   RunUMAP(reduction = "harmony", dims = 1:10) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:10) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# 
# DimPlot(scATAC_hm, reduction = "umap", pt.size = 0.1,label = T) + ggplot2::ggtitle("Harmony integration")
# DimPlot(scATAC_hm, reduction = "umap", pt.size = 0.1,label = T,split.by = "sample") + ggplot2::ggtitle("Harmony integration")
# 
# saveRDS(scATAC_hm,"~/Hislet_aggr/outs/scATAC_hm.rds")
# 
# FeaturePlot(
#   object = scATAC_hm,
#   features = c('INS-IGF2',"GCG","SST"),
#   pt.size = 0.1,
#   max.cutoff = 'q90',
#   ncol = 3
# )
# markers <- FindAllMarkers(scATAC_hm,assay = "RNA",min.pct = .3)
# for (f in unique(markers$cluster)) {
#   print(head(subset(markers,cluster==f)))
# }
# head(markers)
# saveRDS(markers,"~/Hislet_aggr/outs/markers.rds")
# readRDS(markers,"~/Hislet_aggr/outs/markers.rds")
# 
# top_genes <- vector()
# for (f in unique(markers$cluster)) {
#   print(rownames(head(subset(markers,cluster==f),n=10)))
# }
# 
# f=4
# top_genes_select <- c("GCG","INS-IGF2","CFTR","BCAT1","ZMAT4",
#                       "TMEM178B","SLCO3A1","COL1A2","SST","KCNK5",
#                       "PAM","CCL5","PTGIR","RGS1")
# 
# 
# ggsave(filename = "/home/bjw032/Hislet_aggr/DimPlot.pdf",plot = DimPlot(scATAC_hm),height = 8,width = 8)
# 
# FeaturePlot(
#   object = scATAC_hm,
#   features = c('INS-IGF2',"GCG","SST"),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol = 3
# )
# DefaultAssay(scATAC_hm) <- "RNA"
# 
# FeaturePlot(
#   object = scATAC_hm,
#   features = str_split_fixed(top_genes_select,pattern = "[.]",n = 2)[,1] ,
#   pt.size = 0.1,
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   ncol = 5
# )
# 
# FeaturePlot(
#   object = scATAC_hm,
#   features = str_split_fixed(top_genes_select,pattern = "[.]",n = 2)[,1] ,
#   pt.size = 0.1,
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   ncol = 5
# )
# 
# Gene_groups <- FeaturePlot(
#   object = scATAC_hm,
#   features = top_genes_select,
#   pt.size = 0.1,
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   ncol = 5)
#   
# ggsave(filename = "/home/bjw032/Hislet_aggr/cell_markers.pdf",plot = Gene_groups,height = 8,width = 16)
# 


# scATAC.filt <- subset(scATAC.filt, idents = 13, invert = TRUE)
# scATAC.filt.renamed <- RenameIdents(
#   object = scATAC.filt,
#   '0' = 'Beta-1',
#   '1' = 'Beta-2',
#   '2' = 'Alpha-1',
#   '3' = 'Stellate',
#   '4' = 'Alpha-2',
#   '5' = 'Alpha-3',
#   '6' = 'Beta-3',
#   '7' = 'Ductal-1',
#   '8' = 'Delta',
#   '9' = 'Unknown-1',
#   '10' = 'Unknown-2',
#   '11' = 'Unknown-3',
#   '12' = 'Acinar')

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = c(9606),all_versions = FALSE) # human
)

pwm.mm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = c(10090),all_versions = FALSE) # mouse
)

pwm.all <- c(pwm,pwm.mm)# all


# pwm_combined <- c(pwm,pwm.mm) # human and mouse

# add motif information
# 
# main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
# keep.peaks <- which(as.character(seqnames(granges(scATAC_hm))) %in% main.chroms)
# scATAC_hm[["peaks"]] <- subset(scATAC_hm[["peaks"]], features = rownames(scATAC_hm[["peaks"]])[keep.peaks])
# 
# scATAC_hm <- AddMotifs(scATAC_hm, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm.all)
# 
# scATAC_hm <- RunChromVAR(
#   object = scATAC_hm,
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )
# 
# # saveRDS(scATAC_hm,"~/Hislet_aggr/outs/scATAC_hm.rds")
# scATAC_hm <- readDS("~/Hislet_aggr/outs/scATAC_hm.rds")
# 
# umap = cbind("Barcode" = rownames(Embeddings(object = scATAC_hm, reduction = "umap")), Embeddings(object = scATAC_hm, reduction = "umap"))
# scATAC_hm$umap_1 <- as.numeric(umap[,2])
# scATAC_hm$umap_2 <- as.numeric(umap[,3])
# 
# 
# scATAC_hm_filt <- subset(scATAC_hm,umap_1>-10 & umap_1<10 & umap_2 <10)
# scATAC_hm_filt <- subset(scATAC_hm_filt,umap_1<0 | umap_2 < -1)
# scATAC_hm_filt <- subset(scATAC_hm_filt,seurat_clusters %in% c(0,1,4,5,6,7,9,11))
# # saveRDS(scATAC_hm_filt,"~/Hislet_aggr/outs/scATAC_hm_filt.rds")
# scATAC_hm_filt <- readRDS("~/Hislet_aggr/outs/scATAC_hm_filt.rds")
# 
# DimPlot(scATAC_hm_filt)
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "MA0132.2", #PDX1
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
# ) + ggplot2::ggtitle("PDX1 Motif")
# 
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "MA0132.2", #PDX1
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
#   ) + ggplot2::ggtitle("PDX1 Motif")
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "MA0107.1", #RELA
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
# ) + ggplot2::ggtitle("RELA Motif")
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "INS-IGF2", 
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("grey","red","darkred")
# ) + ggplot2::ggtitle("INS-IGF2 Locus")
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "GCG", 
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("grey","red","darkred")
# ) + ggplot2::ggtitle("GCG Locus")
# 
# p_RELA_Endo <- VlnPlot(scATAC_hm_filt,features = "MA0107.1",idents = c(1,4,6,11)) + ggtitle("RELA")
# p_PDX1_Endo <- VlnPlot(scATAC_hm_filt,features = "MA0132.2",idents = c(1,4,6,11)) + ggtitle("PDX1")
# ggsave("~/Hislet_aggr/p_RELA_Endo.pdf",p_RELA_Endo,
#        height=4,width=3)
# ggsave("~/Hislet_aggr/p_PDX1_Endo.pdf",p_PDX1_Endo,
#        height=4,width=3)
# 
# VlnPlot(scATAC_hm_filt,features = "MA0107.1",idents = c(1,4,6,11),pt.size = 0,split.by = "sample") + ggtitle("RELA")
# VlnPlot(scATAC_hm_filt,features = "MA0132.2",idents = c(1,4,6,11),pt.size = 0,split.by = "sample") + ggtitle("PDX1")
# 
# DimPlot(scATAC_hm_filt,label = T,split.by = "sample")
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/whole_islet_PDX1.pdf",
#        FeaturePlot(
#          object = scATAC.filt,
#          features = "MA0132.2", #PDX1
#          #min.cutoff = 'q01',
#          max.cutoff = 'q90',
#          pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")))
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/10X_Pdx1_BetaCells.pdf",
#        FeaturePlot(
#          object = subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10),
#          features = "MA0132.2", #PDX1
#          #min.cutoff = 'q01',
#          max.cutoff = 'q90',
#          pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")))
# DimPlot(subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10))
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/10X_NFKB1_BetaCells.pdf",
#        FeaturePlot(
#          object = subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10),
#          features = "MA0105.4", #NFKB1
#          #min.cutoff = 'q01',
#          max.cutoff = 'q90',
#          pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
#        ) + ggplot2::ggtitle("NFKB1 Motif"))
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/10X_RELA_BetaCells.pdf",
#        FeaturePlot(
#          object = subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10),
#          features = "MA0107.1", #RELA
#          #min.cutoff = 'q01',
#          max.cutoff = 'q90',
#          pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
#        ) + ggplot2::ggtitle("RELA Motif"))
# 
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/10X_clustering.pdf",DimPlot(scATAC.filt,label = T))
# DimPlot(scATAC.filt,label = T)
# 
# p_PDX1_10X <- VlnPlot(subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10 & seurat_clusters %in% c(0,1,6,9)),features = "MA0132.2",sort="increasing",pt.size = 0)
# #VlnPlot(subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10 & seurat_clusters %in% c(0,1,6,9)),features = "MA0105.4",sort="decreasing")
# p_RELA_10X <- VlnPlot(subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10 & seurat_clusters %in% c(0,1,6,9)),features = "MA0107.1",sort="decreasing",pt.size = 0)
# 
# CoveragePlot(subset(scATAC.filt,umap2>0 & umap2<10 & umap1 < 0 & umap1 > -10 & seurat_clusters %in% c(0,1,6,9)),region = "RELB")
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/10X_PDX1.pdf",p_PDX1_10X)
# ggsave("~/Documents/3-pdx1A4/Paper_figures/10X_RELA.pdf",p_RELA_10X)
# 
# 
# 
# INS_p <- FeaturePlot(
#   object = scATAC.filt,
#   features = "INS-IGF2", #RELA
#   min.cutoff = 'q05',
#   max.cutoff = 'q95',
#   pt.size = 0.5,cols = c("grey","blue","darkblue")
# ) + ggplot2::ggtitle("INS-IGF2")
# 
# GCG_p <- FeaturePlot(
#   object = scATAC.filt,
#   features = "GCG", #RELA
#   min.cutoff = 'q05',
#   max.cutoff = 'q95',
#   pt.size = 0.5,cols = c("grey","blue","darkblue")
# ) + ggplot2::ggtitle("GCG")
# 
# SST_p <- FeaturePlot(
#   object = scATAC.filt,
#   features = "SST", #RELA
#   #min.cutoff = 'q05',
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("grey","blue","darkblue")
# ) + ggplot2::ggtitle("SST")
# 
# 
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/cell_markers.pdf",INS_p + GCG_p + SST_p)
# 
# ggsave("~/Documents/3-pdx1A4/Paper_figures/other_cell_markers.pdf",
#        FeaturePlot(
#          object = scATAC.filt,
#          features = "MAFA", #RELA
#          min.cutoff = 'q05',
#          max.cutoff = 'q90',
#          pt.size = 0.5,cols = c("grey","blue","darkblue")
#        ) + ggplot2::ggtitle("MAFA") +
#          FeaturePlot(
#            object = scATAC.filt,
#            features = "ARX", #RELA
#            min.cutoff = 'q05',
#            max.cutoff = 'q90',
#            pt.size = 0.5,cols = c("grey","blue","darkblue")
#          ) + ggplot2::ggtitle("ARX") +
#          FeaturePlot(
#            object = scATAC.filt,
#            features = "LEPR", #RELA
#            min.cutoff = 'q05',
#            max.cutoff = 'q99',
#            pt.size = 0.5,cols = c("grey","blue","darkblue")
#          ) + ggplot2::ggtitle("HHEX"))
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "MA0107.1", #RELA
#   #min.cutoff = 'q05',
#   #max.cutoff = 'q99',
#   pt.size = 0.5,cols = c("grey","blue","darkblue")
# ) + ggplot2::ggtitle("MA0107.1")
# 
# FeaturePlot(
#   object = scATAC_hm_filt,
#   features = "MA0101.1", #REL
#   #min.cutoff = 'q05',
#   #max.cutoff = 'q99',
#   pt.size = 0.5,cols = c("grey","blue","darkblue")
# ) + ggplot2::ggtitle("MA0107.1")
# 
# MA0101.1
# 
# 
# FeaturePlot(
#   object = subset(scATAC_hm),
#   features = "MA0603.1", #BMAL1
#   min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   pt.size = 0.1,cols = c("blue","darkgrey","red")
# ) + ggplot2::ggtitle("BMAL1 Motif")
# 
# ggsave("/home/bjw032/Hislet_aggr/Hislet_endo_dimplot.pdf",plot = DimPlot(scATAC_hm_filt,label = T),height = 4,width = 5)
# DimPlot(scATAC_hm_filt,label = T)
# DefaultAssay(scATAC_hm) <- "peaks"
# CoveragePlot(subset(scATAC_hm,seurat_clusters %in% c(1,4,6,11)),region = "MAFA",
#              peaks = TRUE,
#              tile = TRUE,extend.upstream = 10000,extend.downstream = 2000)
# 
# CoveragePlot(subset(scATAC_hm,seurat_clusters %in% c(1,4,6,11)),region = "TNFAIP3",
#              peaks = TRUE,
#              tile = TRUE,extend.upstream = 5000,extend.downstream = 2000)
# CoveragePlot(subset(scATAC_hm,seurat_clusters %in% c(1,4,6,11)),region = "TNFRSF11B",
#              peaks = TRUE,
#              tile = TRUE,extend.upstream = 5000,extend.downstream = 2000)
# 
# CoveragePlot(subset(scATAC_hm,seurat_clusters %in% c(1,4,6,11)),region = "INS-IGF2",
#              peaks = TRUE,
#              tile = TRUE,extend.upstream = 5000,extend.downstream = 2000)
# 
# 
# CoveragePlot(subset(scATAC_hm,seurat_clusters %in% c(1,0,7,12)),region = "INS-IGF2",
#              peaks = TRUE,
#              tile = TRUE,extend.upstream = 5000,extend.downstream = 2000)
# 
# 
# saveRDS(scATAC.filt,"~/Hislet_aggr/outs/scATAC.filt.rds")
# 
# 
# VlnPlot(object = scATAC.filt, features = "MA0107.1",idents = c(6,9,0,1))# + NoLegend()
# VlnPlot(object = scATAC.filt, features = "INS-IGF2",idents = c(0,1,6,9))# + NoLegend()
# VlnPlot(object = scATAC.filt, features = "MA0132.2",idents = c(6,9,0,1))# + NoLegend()
# 
# 
# DimPlot(object = scATAC.filt, label = TRUE)# + NoLegend()
# 
# 
# 
# da_genes_alpha.v.beta <- FindMarkers(
#   object = scATAC.filt.renamed,
#   ident.1 = c("Alpha-1","Alpha-2","Alpha-3"),
#   ident.2 = c("Beta-1","Beta-2","Beta-3"),
#   min.pct = 0.2,
#   test.use = 'wilcox'
# )
# 
# DefaultAssay(scATAC.filt.renamed) <- 'peaks'
# 
# da_peaks_alpha.v.beta <- FindMarkers(
#   object = scATAC.filt.renamed,
#   ident.1 = c("Alpha-1","Alpha-2","Alpha-3"),
#   ident.2 = c("Beta-1","Beta-2","Beta-3"),
#   min.pct = 0.2,
#   test.use = 'LR',
#   latent.vars = 'peak_region_fragments'
# )
# 
# levels(scATAC.filt.renamed) <-   c('Beta-1','Beta-2','Beta-3',
#                                    'Alpha-1','Alpha-2','Alpha-3',
#                                    'Delta', 'Acinar','Stellate',
#                                    'Ductal-1', 'Unknown-1', 'Unknown-2','Unknown-3')
# 
# 
# CoveragePlot(
#   object = scATAC.filt.renamed,
#   region = rownames(da_peaks_alpha.v.beta)[1],
#   extend.upstream = 40000,
#   extend.downstream = 20000
# )
# 
# CoveragePlot(
#   object = scATAC.filt.renamed,
#   region = "LINGO2",
#   extend.upstream = 40000,
#   extend.downstream = 40000
# )
# 
# 
# #### Export UMAP and cell-identity for Loupe ####
# 
# umap = cbind("Barcode" = rownames(Embeddings(object = scATAC.filt.renamed, reduction = "umap")), Embeddings(object = scATAC.filt.renamed, reduction = "umap"))
# write.table(umap, file="./umap.csv", sep = ",", quote = F, row.names = F, col.names = T)
# 
# cell.idents <- data.frame(barcode=names(Idents(scATAC.filt.renamed)),
#                           id=Idents(scATAC.filt.renamed))
# write.csv(cell.idents,file = "~/Documents/8-HumanIslets/scATAC/cell.idents.csv",quote = F,row.names = F)


# #### Cicero for cis-acting regulatory regions ####
# DefaultAssay(scATAC.filt.renamed) <- "peaks"
# 
# scATAC.cds <- as.cell_data_set(scATAC.filt.footprint)
# scATAC.cicero <- make_cicero_cds(scATAC.cds, reduced_coordinates = reducedDims(scATAC.cds)$UMAP)
# 
# # get the chromosome sizes from the Seurat object
# genome <- seqlengths(scATAC.filt.renamed)
# 
# # use chromosome 1 to save some time
# # omit this step to run on the whole genome
# 
# # convert chromosome sizes to a dataframe
# genome.df <- data.frame("chr" = names(genome), "length" = genome)
# 
# # run cicero
# conns <- run_cicero(scATAC.cicero, genomic_coords = genome.df, sample_num = 100)
# ccans <- generate_ccans(conns)
# links <- ConnectionsToLinks(conns = conns, ccans = ccans)
# scATAC.filt.footprint_ccans <- scATAC.filt.footprint
# Links(scATAC.filt.footprint_ccans) <- links
# 
# #saveRDS(scATAC.filt.footprint_ccans, file = "~/Documents/3-pdx1A4/scATAC.filt.footprint_ccans.rds")
# #saveRDS(scATAC.filt.renamed, file = "~/Documents/8-HumanIslets/scATAC/scATAC_filtered.rds")
# 
# scATAC.filt.footprint_ccans <- readRDS(file = "~/Documents/3-pdx1A4/scATAC.filt.footprint_ccans.rds")
# 
# 
# CoveragePlot(scATAC.filt.renamed, region = "GLP1R")
# 
# phase1_interactions <- read.delim2("~/Downloads/1-h1-H3K27Ac-phase1.interaction.txt",header = F)
# phase1_interactions_end <- paste(str_split_fixed(str_split_fixed(phase1_interactions$V4,pattern = ",",n = 2)[,1],":",2)[,1],
#                                  str_split_fixed(str_split_fixed(str_split_fixed(phase1_interactions$V4,pattern = ",",n = 2)[,1],":",2)[,2],"-",2)[,1],
#                                  str_split_fixed(str_split_fixed(str_split_fixed(phase1_interactions$V4,pattern = ",",n = 2)[,1],":",2)[,2],"-",2)[,2],
#                                  sep="-")
# phase1_interactions_begin <- paste(phase1_interactions[,1],phase1_interactions[,2],phase1_interactions[,3],sep = "-")
# 
# phase1_loops <- read.table("~/Downloads/1-h1-H3K27Ac-phase1.allValidPairs.hic")
# 
# 
# HiChIP_conns <- data.frame(Peak1=phase1_interactions_begin,Peak2=phase1_interactions_end,coaccess=1)
# 
# 
# conns$in_HiChIP <- compare_connections(conns, HiChIP_conns)
# conns$in_HiChIP_500 <- compare_connections(conns, HiChIP_conns, maxgap=500)
# refined.conns <- conns[which(conns$in_HiChIP_500),]
# 
# HiChIP_ccans <- generate_ccans(HiChIP_conns,coaccess_cutoff_override = 0)
# ccans <- generate_ccans(conns)
# refined.ccans <- generate_ccans(refined.conns)
# 
# links <- ConnectionsToLinks(conns = conns, ccans = ccans)
# HiChIP_links <- ConnectionsToLinks(conns = HiChIP_conns, ccans = HiChIP_ccans)
# refined.links <- ConnectionsToLinks(conns = refined.conns, ccans = refined.ccans)
# 
# Links(scATAC.filt.renamed) <- links
# Links(scATAC.filt.renamed) <- HiChIP_links
# Links(scATAC.filt.renamed) <- refined.links
# 
# CoveragePlot(scATAC.filt.renamed, region = "GLP1R",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# CoveragePlot(scATAC.filt.renamed, region = "INS",extend.upstream = 20000,extend.downstream = 20000)
# CoveragePlot(scATAC.filt.renamed,region = "PER1",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "MAFA",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "NFKBIA",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "TNFAIP3",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "ADARB1",extend.upstream = 200000,extend.downstream = 200000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "NR1D1",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "ADRB1",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# library(motifmatchr)
# library(JASPAR2020)
# library(TFBSTools)
# library(BSgenome.Hsapiens.UCSC.hg38)
# 
# # gather the footprinting information for sets of motifs
# scATAC.filt.footprint <- Footprint(
#   object = scATAC.filt.renamed,
#   motif.name = c("Arntl","Atf1","BHLHE41","CLOCK","CREB1","CREB3",
#                  "CREB3L1","CREM","CTCF","DBP","ESRRA","ETS1","ETS2",
#                  "FOS","FOS::JUNB","FOSL1","FOXA1","FOXA2","HES6","HES7",
#                  "HNF1A","HNF4A","ISL2","JUNB","MAFA","Mafb","MEF2A",
#                  "MEF2C","MLXIPL","NEUROD1","NEUROG1","NFATC2","NFKB1",
#                  "NR1D1","NR1D2","NKX6-1","NR4A2","PAX4","PAX6","PDX1",
#                  "PPARA::RXRA","PPARG","RELA","RFX1","RORA","RORB",
#                  "RORC","Rxra","RUNX1","SREBF1","STAT1","STAT3","TEAD1",
#                  "TEF","THRB","XBP1"),
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )
# 
# scATAC.filt.footprint <- Footprint(
#   object = scATAC.filt.renamed,
#   motif.name = c("Arntl","Atf1","BHLHE41","CLOCK","CREB1","CREB3",
#                  "CREB3L1","CREM","DBP","ESRRA","ETS1","ETS2",
#                  "FOS","FOS::JUNB","FOSL1","FOXA1","FOXA2","HES6","HES7",
#                  "HNF1A","HNF4A","ISL2","JUNB","MAFA","Mafb","MEF2A",
#                  "MEF2C","MLXIPL","NEUROD1","NEUROG1","NFATC2","NFKB1",
#                  "NR1D1","NR1D2","NKX6-1","NR4A2","PAX4","PAX6","PDX1",
#                  "PPARA::RXRA","PPARG","RELA","RORA","RORB",
#                  "TEF","THRB"),
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )
# 
# #CTCF is bad
# #RFX1 is bad
# #XBP1 is bad
# 
# #scATAC.filt.renamed <- readRDS("~/Documents/8-HumanIslets/scATAC/scATAC_filtered.rds")
# #saveRDS(scATAC.filt.renamed, file = "~/Documents/8-HumanIslets/scATAC/scATAC_filtered.rds")
# #saveRDS(scATAC.filt.footprint, file = "~/Documents/8-HumanIslets/scATAC/scATAC_filtered_footprint.rds")
# #scATAC.filt.footprint <- readRDS("~/Documents/8-HumanIslets/scATAC/scATAC_filtered_footprint.rds")
# 
# scATAC.filt.footprint_ab <- RenameIdents(
#   object = scATAC.filt.footprint,
#   'Beta-1' = 'Beta',
#   'Beta-2' = 'Beta',
#   'Alpha-1' = 'Alpha',
#   'Stellate' = 'Stellate',
#   'Alpha-2' = 'Alpha',
#   'Alpha-3' = 'Alpha',
#   'Beta-3' = 'Beta',
#   'Ductal-1' = 'Ductal',
#   'Delta' = 'Delta',
#   'Unknown-1' = 'Unknown-1',
#   'Unknown-2' = 'Unknown-2',
#   'Unknown-3' = 'Unknown-3',
#   'Acinar' = 'Acinar')
# 
# # plot the footprint data for each group of cells
# p2 <- PlotFootprint(scATAC.filt.footprint_ab, features = c("Atf1"),
#                     idents = c("Beta","Alpha"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint_ab, features = c("RELA"),
#                     idents = c("Beta","Alpha"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint, features = c("RELA"),
#                     idents = c("Beta-1","Beta-3"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint, features = c("PDX1"),
#                     idents = c("Beta-1","Beta-3"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint, features = c("MAFA"),
#                     idents = c("Beta-1","Beta-3"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint, features = c("PDX1"),
#                     idents = c("Beta-1","Beta-3"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint, features = c("CLOCK"),
#                     idents = c("Beta-1","Beta-3"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# p2 <- PlotFootprint(scATAC.filt.footprint, features = c("FOXA2"),
#                     idents = c("Beta-1","Beta-3"))
# p2 + patchwork::plot_layout(ncol = 1)
# 
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "ADRB1",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# CoveragePlot(scATAC.filt.footprint_ccans,region = "ADRA1B",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "P2RY1",extend.upstream = 1000000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "TNFAIP3",extend.upstream = 500,extend.downstream = 30000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ab,region = "P2RY1",extend.upstream = 2000000,extend.downstream = 200000,
#              idents = c("Beta","Alpha"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "RORC",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "RORC",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ab,region = "RORC",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta","Alpha"))
# 
# CoveragePlot(scATAC.filt.footprint_ab,region = "RORC",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta","Alpha"))
# 
# CoveragePlot(scATAC.filt.footprint_ab,region = "IL1R1",extend.upstream = 200000,extend.downstream = 20000,
#              idents = c("Beta","Alpha"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "IL1R1",extend.upstream = 20000,extend.downstream = 20000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "PER2",extend.upstream = 20000,extend.downstream = 200000,
#              idents = c("Beta-1","Beta-2","Beta-3","Alpha-1","Alpha-2","Alpha-3"))
# 
# CoveragePlot(scATAC.filt.footprint_ccans,region = "FMO3",extend.upstream = 10000,extend.downstream = 0,
#              idents = c('Beta-1','Beta-2','Beta-3',
#                         'Alpha-1','Alpha-2','Alpha-3',
#                         'Delta', 'Stellate',
#                         'Ductal-1'))

# Islet1_scATAC.filt <- readRDS(file = "~/Hislet_aggr/Islet1.rds")
# Islet2_scATAC.filt <- readRDS(file = "~/Hislet_aggr/Islet2.rds")
# Islet3_scATAC.filt <- readRDS(file = "~/Hislet_aggr/Islet3.rds")
# 
# 
# frags <- Fragments(Islet1_scATAC.filt)  # get list of fragment objects
# Fragments(Islet1_scATAC.filt) <- NULL  # remove fragment information from assay
# new.paths <- '~/Hislet_aggr/gaulton/summary/S13.fragments.tsv.gz'
# for (i in seq_along(frags)) {
#   frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
# }
# Fragments(Islet1_scATAC.filt) <- frags # assign updated list back to the object
# 
# frags <- Fragments(Islet2_scATAC.filt)  # get list of fragment objects
# Fragments(Islet2_scATAC.filt) <- NULL  # remove fragment information from assay
# new.paths <- '~/Hislet_aggr/gaulton/summary/S14.fragments.tsv.gz'
# for (i in seq_along(frags)) {
#   frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
# }
# Fragments(Islet2_scATAC.filt) <- frags # assign updated list back to the object
# 
# frags <- Fragments(Islet3_scATAC.filt)  # get list of fragment objects
# Fragments(Islet3_scATAC.filt) <- NULL  # remove fragment information from assay
# new.paths <- '~/Hislet_aggr/gaulton/summary/S15.fragments.tsv.gz'
# for (i in seq_along(frags)) {
#   frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
# }
# Fragments(Islet3_scATAC.filt) <- frags # assign updated list back to the object
# 
# scATAC_hm <- readRDS("~/Hislet_aggr/outs/scATAC_hm.rds")
# 
# ## This keeps timing out...
# combined_10X <- merge(
#   x = scATAC_hm,
#   y = list(Islet1_scATAC.filt,Islet2_scATAC.filt, Islet3_scATAC.filt),
#   add.cell.ids = c("10X","Islet1", "Islet2", "Islet3")
# )
# UpdatePath(Fragments(combined_10X)[[2]], )
# saveRDS(combined_10X,"~/Hislet_aggr/combined_10X.rds")
# combined_10X <- readRDS("~/Hislet_aggr/combined_10X.rds")
# 
# combined_10X <- RunTFIDF(combined_10X)
# combined_10X <- FindTopFeatures(combined_10X, min.cutoff = 'q0')
# combined_10X <- RunSVD(combined_10X)
# combined_10X <- RunUMAP(object = combined_10X, reduction = 'lsi', dims = 2:30)
# combined_10X <- FindNeighbors(object = combined_10X, reduction = 'lsi', dims = 2:30)
# combined_10X <- FindClusters(object = combined_10X, verbose = FALSE, algorithm = 3)
# library(stringr)
# combined_10X$sample <- str_split_fixed(colnames(combined_10X),pattern="_",n=2)[,1]
# DimPlot(object = combined_10X, label = TRUE,split.by = "sample")# + NoLegend()
# 
# DefaultAssay(combined_10X) <- 'peaks'
# 
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"
# genome(annotation) <- "hg38"
# Annotation(combined_10X) <- annotation
# gene.activities <- GeneActivity(combined_10X)
# 
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# combined_10X[['RNA']] <- CreateAssayObject(counts = gene.activities)
# combined_10X <- NormalizeData(
#   object = combined_10X,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(combined_10X$nCount_RNA)
# )
# 
# library(harmony)
# combined_10X_hm <- RunHarmony(
#   object = combined_10X,
#   group.by.vars = 'sample',
#   reduction = 'lsi',
#   assay.use = 'peaks',
#   project.dim = FALSE
# )
# 
# combined_10X_hm <- combined_10X_hm %>%
#   RunUMAP(reduction = "harmony", dims = 1:10) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:10) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# 
# DimPlot(combined_10X_hm,split.by = "sample",reduction = "umap")
# 
# 
# FeaturePlot(combined_10X_hm,"GCG",max.cutoff = "q90")
# 
# 
# 
# 
# DimPlot(combined_10X_hm)
# # saveRDS(combined_10X,file = "~/Documents/3-pdx1A4/gaulton/Islet_10X_combined.rds")
# 
# 
# #
# # # extract gene annotations from EnsDb
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- "UCSC"
# genome(annotation) <- "hg38"
# Annotation(combined_10X_hm) <- annotation
# gene.activities <- GeneActivity(combined_10X_hm)
# 
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# combined_10X_hm[['RNA']] <- CreateAssayObject(counts = gene.activities)
# combined_10X_hm <- NormalizeData(
#   object = combined_10X_hm,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(combined_10X_hm$nCount_RNA)
# )
# DefaultAssay(combined_10X_hm) <- 'peaks'
# 
# FeaturePlot(combined_10X_hm,"INS-IGF2",max.cutoff = "q90", split.by="sample")
# 
# 
# pwm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(species = c(9606),all_versions = FALSE) # human
# )
# 
# pwm.mm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(species = c(10090),all_versions = FALSE) # mouse
# )
# 
# pwm.all <- c(pwm,pwm.mm)# all
# 
# 
# # pwm_combined <- c(pwm,pwm.mm) # human and mouse
# 
# # add motif information
# 
# main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
# keep.peaks <- which(as.character(seqnames(granges(combined_10X_hm))) %in% main.chroms)
# combined_10X_hm[["peaks"]] <- subset(combined_10X_hm[["peaks"]], features = rownames(combined_10X_hm[["peaks"]])[keep.peaks])
# combined_10X_hm <- AddMotifs(combined_10X_hm, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm.all)
# combined_10X_hm <- RunChromVAR(
#   object = combined_10X_hm,
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )
# 
# 
# FeaturePlot(
#   object = combined_10X_hm,
#   features = "MA0132.2", #PDX1
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',split.by="sample",
#   pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
# ) + ggplot2::ggtitle("PDX1 Motif")
# 
# FeaturePlot(
#   object = combined_10X_hm,
#   features = "MA0107.1", #RELA
#   #min.cutoff = 'q01',
#   max.cutoff = 'q90',
#   pt.size = 0.5,cols = c("darkblue","blue","grey","red","darkred")
# ) + ggplot2::ggtitle("RELA Motif")

#### Heatmaps ####

# 
# scATAC_hm <- readRDS("~/Hislet_aggr/outs/scATAC_hm.rds")
# scATAC_hm_filt <- readRDS("~/Hislet_aggr/outs/scATAC_hm_filt.rds")
# DimPlot(scATAC_hm)

# DefaultAssay(scATAC_hm) <- "chromvar"
# all_motif_markers <- FindAllMarkers(scATAC_hm,assay = "chromvar",min.pct = .05)
# saveRDS(all_motif_markers,"~/Hislet_aggr/outs/all_motif_markers_PC.rds")
# 
# DefaultAssay(scATAC_hm) <- "peaks"
# all_peak_markers <- FindAllMarkers(scATAC_hm,assay = "peaks",min.pct = .05)
# saveRDS(all_peak_markers,"~/Hislet_aggr/outs/all_peak_markers_PC.rds")
# 
# DefaultAssay(scATAC_hm) <- "RNA"
# all_RNA_markers <- FindAllMarkers(scATAC_hm,assay = "RNA",min.pct = .05)
# saveRDS(all_RNA_markers,"~/Hislet_aggr/outs/all_RNA_markers_PC.rds")

# DefaultAssay(scATAC_hm_filt) <- "chromvar"
# endo_motif_markers <- FindAllMarkers(scATAC_hm_filt,assay = "chromvar",min.pct = .05)
# saveRDS(endo_motif_markers,"~/Hislet_aggr/outs/endo_motif_markers_PC.rds")
# 
# DefaultAssay(scATAC_hm_filt) <- "peaks"
# endo_peak_markers <- FindAllMarkers(scATAC_hm_filt,assay = "peaks",min.pct = .05)
# saveRDS(endo_peak_markers,"~/Hislet_aggr/outs/endo_peak_markers_PC.rds")
# 
# DefaultAssay(scATAC_hm_filt) <- "RNA"
# endo_RNA_markers <- FindAllMarkers(scATAC_hm_filt,assay = "RNA",min.pct = .05)
# saveRDS(endo_RNA_markers,"~/Hislet_aggr/outs/endo_RNA_markers_PC.rds")

all_motif_markers <- readRDS("~/Hislet_aggr/outs/all_motif_markers_PC.rds")
all_peak_markers <- readRDS("~/Hislet_aggr/outs/all_peak_markers_PC.rds")
all_RNA_markers <- readRDS("~/Hislet_aggr/outs/all_RNA_markers_PC.rds")

endo_motif_markers <- readRDS("~/Hislet_aggr/outs/endo_motif_markers_PC.rds")
endo_peak_markers <- readRDS("~/Hislet_aggr/outs/endo_peak_markers_PC.rds")
endo_RNA_markers <- readRDS("~/Hislet_aggr/outs/endo_RNA_markers_PC.rds")

# 
# DefaultAssay(scATAC_hm) <- "peaks"
# all_peak_markers <- FindAllMarkers(scATAC_hm,assay = "peaks",min.pct = .05)
# saveRDS(all_peak_markers,"~/Hislet_aggr/outs/all_peak_markers_PC.rds")
# 
# DefaultAssay(scATAC_hm) <- "RNA"
# all_RNA_markers <- FindAllMarkers(scATAC_hm,assay = "RNA",min.pct = .05)
# saveRDS(all_RNA_markers,"~/Hislet_aggr/outs/all_RNA_markers_PC.rds")
# 
# DefaultAssay(scATAC_hm) <- 'RNA'
# scATAC_hm <- ScaleData(scATAC_hm)
# top10_RNA <- all_RNA_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# scATAC_hm$celltype <- scATAC_hm$seurat_clusters
# Idents(scATAC_hm) <- "celltype"
# # levels(scATAC_hm) <- c(0,1,2,3,5,4,6)
# p_RNA_all <- DoHeatmap(scATAC_hm,features=top10_RNA$gene,disp.min = -2,disp.max = 2) + scale_fill_gradientn(colors = c("cornflowerblue", "black", "yellow"))
# 
# DefaultAssay(scATAC_hm) <- 'chromvar'
# scATAC_hm <- ScaleData(scATAC_hm)
# top10_motif <- all_motif_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# scATAC_hm$celltype <- scATAC_hm$seurat_clusters
# Idents(scATAC_hm) <- "celltype"
# # levels(scATAC_hm) <- c(0,1,2,3,5,4,6)
# p_motif_all <- DoHeatmap(scATAC_hm,features=top10_motif$gene,disp.min = -2,disp.max = 2) + scale_fill_gradientn(colors = c("cornflowerblue", "black", "yellow"))
# top10_motif$gene
# 
# name(pfm_all)
# top10_motif_df <- data.frame(unique(top10_motif$gene))
# top10_motif_df$name <- NA
# for(f in 1:nrow(top10_motif_df)){
#   top10_motif_df$name[f] <- name(pfm_all[top10_motif$gene[f],])
# }
# top10_motif_df$col1 <- runif(nrow(top10_motif_df))
# top10_motif_df$col2 <- runif(nrow(top10_motif_df))
# pheatmap(top10_motif_df[,3:4],cluster_cols = F,cluster_rows = F,
#          labels_row =top10_motif_df$name)
# 
# 
# 
# 
# 
# DefaultAssay(scATAC_hm_filt) <- 'chromvar'
# scATAC_hm_filt <- ScaleData(scATAC_hm_filt)
# top10_motif_endo <- endo_motif_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
# scATAC_hm_filt$celltype <- scATAC_hm_filt$seurat_clusters
# Idents(scATAC_hm_filt) <- "celltype"
# # levels(scATAC_hm) <- c(0,1,2,3,5,4,6)
# p_motif_endo <- DoHeatmap(scATAC_hm_filt,features=top10_motif_endo$gene,disp.min = -2,disp.max = 2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
# 
# DefaultAssay(scATAC_hm_filt) <- 'RNA'
# scATAC_hm_filt <- ScaleData(scATAC_hm_filt)
# top10_RNA_endo <- endo_RNA_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# scATAC_hm_filt$celltype <- scATAC_hm_filt$seurat_clusters
# Idents(scATAC_hm_filt) <- "celltype"
# # levels(scATAC_hm) <- c(0,1,2,3,5,4,6)
# p_RNA_endo <- DoHeatmap(scATAC_hm_filt,features=top10_RNA_endo$gene,disp.min = -1.5,disp.max = 1.5) + scale_fill_gradientn(colors = c("cornflowerblue", "black", "yellow"))
# 
# DefaultAssay(scATAC_hm_filt) <- 'chromvar'
# motif_1.v.4 <- FindMarkers(scATAC_hm_filt,ident.1 = 1,ident.2 = 4,assay = "chromvar")
# #motif_1.v.4 <- subset(motif_1.v.4,p_val_adj < 0.05)
# motif_1.v.4$names <- NA
# for(f in 1:nrow(motif_1.v.4)){
#   motif_1.v.4$names[f] <- name(pfm_all[substr(rownames(motif_1.v.4)[f],start = 0,stop = 8),])
# }
# 
# motif_1.v.4$diff <- motif_1.v.4$pct.2-motif_1.v.4$pct.1
# sub_motif_1.v.4 <- subset(motif_1.v.4,abs(avg_log2FC) > 1)
# motif_1.v.4_up <- subset(sub_motif_1.v.4,diff>0)
# motif_1.v.4_up <- motif_1.v.4_up[order(abs(motif_1.v.4_up$p_val_adj),decreasing = F),]
# motif_1.v.4_down <- subset(sub_motif_1.v.4,diff < 0)
# motif_1.v.4_down <- motif_1.v.4_down[order(abs(motif_1.v.4_down$p_val_adj),decreasing = F),]
# 
# plot(motif_1.v.4$diff,y = -log10(motif_1.v.4$p_val_adj))
# volcanoplot_markers(motif_1.v.4)
# volcanoplot_markers <- function (res, lfcthresh=.25, sigthresh=0.05, main="Volcano Plot", labelsig=TRUE, textcx=1, ...) {
#   with(res, plot(diff, -log10(p_val_adj), pch=20, main=main, ...))
#   with(subset(res, p_val_adj<sigthresh ), points(diff, -log10(p_val_adj), pch=20, col="red", ...))
#   with(subset(res, abs(diff)>lfcthresh), points(diff, -log10(p_val_adj), pch=20, col="orange", ...))
#   with(subset(res, p_val_adj<sigthresh & abs(diff)>lfcthresh), points(diff, -log10(p_val_adj), pch=20, col="green", ...))
#   if (labelsig) {
#     require(calibrate)
#     with(subset(res, p_val_adj<sigthresh & abs(diff)>lfcthresh), textxy(diff, -log10(p_val_adj), labs=names, cex=textcx, ...))
#   }
# }
# 
# motif.name_1.v.4 <- name(pfm_all[rownames(motif_1.v.4)[which(rownames(motif_1.v.4) %in% names(pfm_all))]])
# motif.name_1.v.4 <- data.frame(motif.name_1.v.4)
# motif.name_1.v.4[c(rownames(motif_1.v.4_down)),]
# motif.name_1.v.4[c(rownames(motif_1.v.4_up)),]
# 
# motif.rename_1.v.4 <- motif.name_1.v.4
# motif.rename_1.v.4$JASPAR_name <- rownames(motif.rename_1.v.4)
# rownames(motif.rename_1.v.4) <- motif.rename_1.v.4$motif.name_1.v.4
# TFs_of_interest_up <- c("PDX1","NEUROD1","POU6F2","Ptf1a","NKX6-1",
#                         "TCF3","Ahr::Arnt","ATOH1(var.2)","BHLHE22(var.2)","PAX6","RORC","TEF","DBP","HLF",
#                         "RELA","NFKB1","FOS::JUNB","ELF3","KLF14","Smad2::Smad3","HNF4A","THRB","IKZF1","RXRG","IRF1","NR4A2")
# df_annotate <- data.frame(row.names = TFs_of_interest_up,col1 = runif(length(TFs_of_interest_up)), col2 = runif(length(TFs_of_interest_up)))
# 
# 
# select_motifs <- DoHeatmap(subset(scATAC_hm_filt,seurat_clusters %in% c(1,4)),
#           features = motif.rename_1.v.4[TFs_of_interest_up,"JASPAR_name"])
# select_motif_labels <- pheatmap(df_annotate,cluster_rows = F,cluster_cols = F)
# ggsave("~/Hislet_aggr/select_motifs_1.v.4.pdf",plot = select_motifs,device = "pdf",
#       width=6,height = 6)
# ggsave("~/Hislet_aggr/select_motifs_labels_1.v.4.pdf",plot = select_motif_labels,device = "pdf",
#        width=6,height = 6)
# 
# 
# DoHeatmap(subset(scATAC_hm_filt,seurat_clusters %in% c(1,4)),
#           features = c(rownames(head(motif_1.v.4_down,n=50)),
#                        rownames(head(motif_1.v.4_up,n=50))),)
# 
# 
# DefaultAssay(scATAC_hm_filt) <- 'RNA'
# RNA_1.v.4 <- FindMarkers(scATAC_hm_filt,ident.1 = 1,ident.2 = 4)
# 
# ggsave("/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_RNA_heatmap_PC.pdf",plot = p_endo_RNA)
# 
# motif.name <- name(pfm_all[subset(top10_motif_endo,gene %in% names(pfm_all))$gene])

