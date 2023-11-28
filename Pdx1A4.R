library(Signac)
library(Seurat)
library(cicero)
library(monocle)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(rtracklayer)
library(harmony)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(dplyr)


# BiocManager::install(c("Signac","Seurat","cicero","monocle","GenomeInfoDb",
#                        "EnsDb.Mmusculus.v79",
#                        "ggplot2",
#                        "patchwork",
#                        "hdf5r",
#                        "rtracklayer",
#                        "harmony",
#                        "JASPAR2020",
#                        "TFBSTools",
#                        "BSgenome.Mmusculus.UCSC.mm10",
#                        "chromVAR",
#                        "dplyr"), force=T
#                      )


set.seed(1234)

setwd("/home/bjw032/Documents/scATAC/01.10.2022_analysis/")

aggr_counts <- Read10X_h5(filename = "/media/bjw032/seagate/scATAC/Pdx1A4/Pdx1A4_aggr_2/outs/filtered_peak_bc_matrix.h5")
aggr_metadata <- read.csv(
  file = "/media/bjw032/seagate/scATAC/Pdx1A4/Pdx1A4_aggr_2/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

aggr_chrom_assay <- CreateChromatinAssay(
  counts = aggr_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = '/media/bjw032/seagate/scATAC/Pdx1A4/Pdx1A4_aggr_2/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

aggr_scATAC <- CreateSeuratObject(
  counts = aggr_chrom_assay,
  assay = "peaks",
  meta.data = aggr_metadata
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

#DOESN'T WORK... UCSC IS BROKEN AS OF 7/10/2021#
# change to UCSC style since the data was mapped to mm10 ...
seqlevelsStyle(annotations) <- 'UCSC'
#DOESN'T WORK... UCSC IS BROKEN AS OF 7/10/2021#

# add the gene information to the object
Annotation(aggr_scATAC) <- annotations

# compute nucleosome signal score per cell
aggr_scATAC <- NucleosomeSignal(object = aggr_scATAC)

# compute TSS enrichment score per cell
aggr_scATAC <- TSSEnrichment(object = aggr_scATAC, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
aggr_scATAC$pct_reads_in_peaks <- aggr_scATAC$peak_region_fragments / aggr_scATAC$passed_filters * 100
aggr_scATAC$blacklist_ratio <- aggr_scATAC$blacklist_region_fragments / aggr_scATAC$peak_region_fragments


# scATAC$high.tss <- ifelse(scATAC$TSS.enrichment > 2, 'High', 'Low')
# TSSPlot(scATAC, group.by = 'high.tss') + NoLegend()

aggr_scATAC$nucleosome_group <- ifelse(aggr_scATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# Con_scATAC$high.tss <- ifelse(Con_scATAC$TSS.enrichment > 2, 'High', 'Low')
# TSSPlot(Con_scATAC, group.by = 'high.tss') + NoLegend()

aggr_scATAC_filt <- subset(
  x = aggr_scATAC,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)



# aggr_scATAC_filt <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/aggr_scATAC_filt.rds")
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- which(as.character(seqnames(granges(aggr_scATAC_filt))) %in% main.chroms)
aggr_scATAC_filt[["peaks"]] <- subset(aggr_scATAC_filt[["peaks"]], features = rownames(aggr_scATAC_filt[["peaks"]])[keep.peaks])

aggr_scATAC_filt <- RunTFIDF(aggr_scATAC_filt)
aggr_scATAC_filt <- FindTopFeatures(aggr_scATAC_filt, min.cutoff = 'q0')
aggr_scATAC_filt <- RunSVD(aggr_scATAC_filt)
aggr_scATAC_filt <- RunUMAP(object = aggr_scATAC_filt, reduction = 'lsi', dims = 2:30)
aggr_scATAC_filt <- FindNeighbors(object = aggr_scATAC_filt, reduction = 'lsi', dims = 2:30)
aggr_scATAC_filt <- FindClusters(object = aggr_scATAC_filt, verbose = FALSE, algorithm = 3)
aggr_scATAC_filt$dataset <- gsub(pattern = "..*.-","",colnames(aggr_scATAC_filt))
aggr_scATAC_filt$genotype <- ifelse(aggr_scATAC_filt$dataset < 3,yes = "Con",no="Mut")
aggr_scATAC_filt$day <- ifelse(aggr_scATAC_filt$dataset == 1 | aggr_scATAC_filt$dataset == 3,yes = 1,no=2)

umap = cbind("Barcode" = rownames(Embeddings(object = aggr_scATAC_filt, reduction = "umap")), Embeddings(object = aggr_scATAC_filt, reduction = "umap"))
write.table(umap, file="/media/bjw032/seagate/scATAC/Pdx1A4/Pdx1A4_aggr_2/umap.csv", sep = ",", quote = F, row.names = F, col.names = T)
cell.idents <- data.frame(barcode=umap[,1],
                          id=Idents(aggr_scATAC_filt))
write.csv(cell.idents,file = "/home/bjw032/Documents/scATAC/01.10.2022_analysis/aggr_scATAC_filt.idents.csv",quote = F,row.names = F)

DimPlot(object = aggr_scATAC_filt, split.by = "day",label = TRUE)# + NoLegend()
DimPlot(object = aggr_scATAC_filt, split.by = "dataset",label = TRUE)# + NoLegend()
DimPlot(object = aggr_scATAC_filt, split.by = "genotype",label = TRUE)# + NoLegend()

aggr_scATAC_filt_gene.actvities <- GeneActivity(aggr_scATAC_filt)
aggr_scATAC_filt[['RNA']] <- CreateAssayObject(counts = aggr_scATAC_filt_gene.actvities)
aggr_scATAC_filt <- NormalizeData(
  object = aggr_scATAC_filt,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(aggr_scATAC_filt$nCount_RNA)
)

FeaturePlot(aggr_scATAC_filt,"Ins1",split.by="genotype",max.cutoff = "q95",order =T,pt.size = .5)
DotPlot(aggr_scATAC_filt,assay = "RNA","Ins1",split.by = "genotype",cols = c("grey","red"), dot.scale = 7) + RotatedAxis()
p_Con <- DotPlot(subset(aggr_scATAC_filt,genotype=="Con"),assay = "RNA","Mafa",cols = c("grey","red"), dot.scale = 7,scale.max = 30,scale.min = 20,col.min = -1.5,col.max = 2) + RotatedAxis()
p_Mut <- DotPlot(subset(aggr_scATAC_filt,genotype=="Mut"),assay = "RNA","Mafa",cols = c("grey","red"), dot.scale = 7,scale.max = 30,scale.min = 20,col.min = -1.5,col.max = 2) + RotatedAxis()
plot_grid(p_Con, p_Mut)
#saveRDS(aggr_scATAC_filt,file = "/home/bjw032/Documents/scATAC/01.10.2022_analysis/aggr_scATAC_filt.rds")
#aggr_scATAC_filt <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/aggr_scATAC_filt.rds")

hm_aggr <- RunHarmony(
  object = aggr_scATAC_filt,
  group.by.vars = 'day',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
hm_aggr <- hm_aggr %>%
  RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

DimPlot(hm_aggr, reduction = "umap", pt.size = 0.1,label = T) + ggplot2::ggtitle("Harmony integration")

umap <- cbind("Barcode" = rownames(Embeddings(object = hm_aggr, reduction = "umap")), Embeddings(object = hm_aggr, reduction = "umap"))
write.table(umap, file="/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_umap.csv", sep = ",", quote = F, row.names = F, col.names = T)
cell.idents <- data.frame(barcode=umap[,1],
                          id=Idents(hm_aggr))
write.csv(cell.idents,file = "/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_cell.idents.csv",quote = F,row.names = F)

# FeaturePlot(hm_aggr,"Ins1",max.cutoff = "q95",order = T,split.by="genotype",pt.size = .5)
# FeaturePlot(hm_aggr,"Gcg",max.cutoff = "q95",order = T,split.by="genotype",pt.size = .5,)
# 
# FeaturePlot(hm_aggr,"Pdx1",order = T,split.by="genotype",pt.size = .5,cols = c("Light Grey","Black"),keep.scale = "all")
# 
# FeaturePlot(hm_aggr,"Pdx1",order = T,split.by="genotype",pt.size = .5,
#             cols = c("Light Grey","Black"), keep.scale = "all")
# FeaturePlot(hm_aggr,"Pdx1",order = T,split.by="genotype",max.cutoff = "q95",pt.size = .5,cols = c("Light Grey","Black"),)
# FeaturePlot(hm_aggr,"Ins2",order = T,split.by="genotype",max.cutoff = "q95",pt.size = .5,cols = c("Light Grey","Black"),)
# 
# 
# FeaturePlot(subset(hm_aggr,genotype=="Con"),"Pdx1",max.cutoff = "q95",order = T,pt.size = .5,cols = c("Light Grey","Black"),keep.scale = "all")
# FeaturePlot(subset(hm_aggr,genotype=="Mut"),"Pdx1",max.cutoff = "q95",order = T,pt.size = .5,cols = c("Light Grey","Black"),keep.scale = "all")
# DimPlot(hm_aggr, reduction = "umap", pt.size = 0.1) + ggplot2::ggtitle("Harmony integration",label = T)
# DimPlot(hm_aggr, reduction = "umap", pt.size = 0.1,label = T) + ggplot2::ggtitle("Harmony integration")

#### Motif Analysis ####

#hm_aggr_motifs <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_motifs.rds")
pfm_human <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = F)
)

pfm_mouse <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = F)
)

# pfm_pdx1 <- PFMatrixList(getMatrixByID(JASPAR2020, ID="MA0132.1"))

pfm_all <- c(pfm_human,pfm_mouse)

DefaultAssay(hm_aggr) <- 'peaks'
hm_aggr_motifs <- AddMotifs(
  object = hm_aggr,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm_all
)
#### Start
# saveRDS(hm_aggr_motifs,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_motifs.rds")
# hm_aggr_motifs <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_motifs.rds")
DimPlot(hm_aggr_motifs)

hm_aggr_motifs <- RunChromVAR(
  object = hm_aggr_motifs,
  genome = BSgenome.Mmusculus.UCSC.mm10
)


umap = cbind("Barcode" = rownames(Embeddings(object = hm_aggr_motifs, reduction = "umap")), Embeddings(object = hm_aggr_motifs, reduction = "umap"))
hm_aggr_motifs$umap_1 <- as.numeric(umap[,2])
hm_aggr_motifs$umap_2 <- as.numeric(umap[,3])

#hm_aggr@meta.data <- hm_aggr@meta.data[colnames(hm_aggr),]
#DefaultAssay(hm_aggr) <- "peaks"
hm_aggr_cds <- as.CellDataSet(x = hm_aggr)
hm_aggr_cicero <- make_cicero_cds(hm_aggr_cds, reduced_coordinates = cbind(hm_aggr_motifs$umap_1,hm_aggr_motifs$umap_2))
genome <- seqlengths(hm_aggr)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns <- run_cicero(hm_aggr_cicero, genome.df, sample_num = 100)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns,ccans = ccans)
Links(hm_aggr_motifs) <- links

DimPlot(hm_aggr_motifs)

DimPlot(hm_aggr_motifs, reduction = "umap", pt.size = 0.1,label = T) + ggplot2::ggtitle("Harmony integration")
DimPlot(hm_endocrine, reduction = "umap", pt.size = 0.1,label = T) + ggplot2::ggtitle("Harmony integration")

# 
# DefaultAssay(hm_endocrine) <- "chromvar"
# 
# 
# FeaturePlot(hm_aggr_motifs,features = "MA0132.2",split.by="genotype",max.cutoff = "q95",min.cutoff = "q05",order =T,pt.size = .5)
# FeaturePlot(hm_aggr_motifs,features = "MA0107.1",split.by="genotype",max.cutoff = "q95",min.cutoff = "q05",order =T,pt.size = .5)
# 
# 
# FeaturePlot(hm_endocrine,features = "MA0132.2",split.by="genotype",max.cutoff = "q95",min.cutoff = "q05",order =T,pt.size = .5) #pdx1
# FeaturePlot(hm_endocrine,features = "MA0107.1",split.by="genotype",max.cutoff = "q95",min.cutoff = "q05",order =T,pt.size = .5) #rela
# 
# VlnPlot(hm_endocrine,features = "MA0132.2",split.by="genotype")
# VlnPlot(hm_endocrine,features = "MA0107.1",split.by="genotype")
# VlnPlot(hm_endocrine,features = "MA0874.1",split.by="genotype")
# 
# saveRDS(hm_aggr_motifs,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_motifs_PC.rds")
# saveRDS(hm_endocrine,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_endocrine_motifs_PC.rds")
# 
# 
# hm_endocrine$new_idents <- ifelse(hm_endocrine$seurat_clusters %in% c(0,1,2,4,5,7,8),"beta","alpha")
# 
# VlnPlot(hm_endocrine,features = "MA0107.1",split.by="genotype",group.by ="new_idents" )
# VlnPlot(hm_endocrine,features = "MA0132.2",split.by="genotype",group.by ="new_idents" )
# 
# 
# DimPlot(hm_aggr, reduction = "umap", pt.size = 0.1,label = T) + ggplot2::ggtitle("Harmony integration")
# 
# Mut_b_cells <- colnames(subset(hm_endocrine,genotype=="Mut" & new_idents=="beta"))
# Con_b_cells <- colnames(subset(hm_endocrine,genotype=="Con" & new_idents=="beta"))
# 
# 
# DefaultAssay(hm_endocrine) <- "chromvar"
# 
# for(f in 1:8){
#   assign(paste0("i",0,"_","i",f),
#          FindMarkers(
#            object = hm_endocrine,
#            ident.1 = '0',
#            ident.2 = f,
#            min.pct = 0.05,
#            test.use = 'LR',
#            latent.vars = 'peak_region_fragments'
#          ))
# }
# 
# hm_endocrine$celltype <- hm_endocrine$seurat_clusters
# hm_endocrine$celltype.dx <- paste(hm_endocrine$celltype, hm_endocrine$genotype, sep = "_")
# Idents(hm_endocrine) <- "celltype.dx"
# DefaultAssay(hm_endocrine) <- "peaks"
# for(f in unique(hm_endocrine$seurat_clusters)){
#   assign(paste0("endo_",f,"_Mut_v_Con"),
#          FindMarkers(
#            object = hm_endocrine,
#            ident.1 = paste0(f,"_","Mut"),
#            ident.2 = paste0(f,"_","Con"),
#            assay="peaks",
#            min.pct = 0.05,
#            test.use = 'LR',logfc.threshold = 0.1,
#            latent.vars = 'peak_region_fragments'
#          ))
# }
# 
# DefaultAssay(hm_endocrine) <- "chromvar"
# Idents(hm_endocrine)
# for(f in unique(hm_endocrine$seurat_clusters)){
#   assign(paste0("endo_",f,"_Mut_v_Con_chromvar"),
#          FindMarkers(
#            object = hm_endocrine,
#            ident.1 = paste0(f,"_","Mut"),
#            ident.2 = paste0(f,"_","Con"),
#            assay="chromvar",
#            min.pct = 0.05,
#            test.use = 'LR',logfc.threshold = 0.1,
#            latent.vars = 'peak_region_fragments'
#          ))
# }
# 
# DefaultAssay(hm_endocrine) <- "RNA"
# Idents(hm_endocrine)
# for(f in unique(hm_endocrine$seurat_clusters)){
#   assign(paste0("endo_",f,"_Mut_v_Con_RNA"),
#          FindMarkers(
#            object = hm_endocrine,
#            ident.1 = paste0(f,"_","Mut"),
#            ident.2 = paste0(f,"_","Con"),
#            assay="RNA",
#            min.pct = 0.05,
#            test.use = 'LR',logfc.threshold = 0.1,
#            latent.vars = 'peak_region_fragments'
#          ))
# }
# 
# 
# Idents(hm_endocrine) <- "celltype"
# DimPlot(hm_endocrine)
# head(endo_0_Mut_v_Con)
# head(endo_1_Mut_v_Con)
# head(endo_2_Mut_v_Con)
# head(endo_3_Mut_v_Con)
# head(endo_4_Mut_v_Con)
# head(endo_5_Mut_v_Con)
# head(endo_6_Mut_v_Con)
# head(endo_7_Mut_v_Con)
# head(endo_8_Mut_v_Con)
# 
# features_Mut_v_Con_chromvar <- unique(c(
#   rownames(head(endo_0_Mut_v_Con_chromvar)),
#   rownames(head(endo_1_Mut_v_Con_chromvar)),
#   rownames(head(endo_2_Mut_v_Con_chromvar)),
#   rownames(head(endo_3_Mut_v_Con_chromvar)),
#   rownames(head(endo_4_Mut_v_Con_chromvar)),
#   rownames(head(endo_5_Mut_v_Con_chromvar)),
#   rownames(head(endo_6_Mut_v_Con_chromvar)),
#   rownames(head(endo_7_Mut_v_Con_chromvar)),
#   rownames(head(endo_8_Mut_v_Con_chromvar))))
# 
# 
# 
# features_Mut_v_Con_RNA <- unique(c(
#   rownames(head(endo_0_Mut_v_Con_RNA)),
#   rownames(head(endo_1_Mut_v_Con_RNA)),
#   rownames(head(endo_2_Mut_v_Con_RNA)),
#   rownames(head(endo_3_Mut_v_Con_RNA)),
#   rownames(head(endo_4_Mut_v_Con_RNA)),
#   rownames(head(endo_5_Mut_v_Con_RNA)),
#   rownames(head(endo_6_Mut_v_Con_RNA)),
#   rownames(head(endo_7_Mut_v_Con_RNA)),
#   rownames(head(endo_8_Mut_v_Con_RNA))))
# 
# DimPlot(hm_endocrine)
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_0_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_1_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_2_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_3_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_4_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_5_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_6_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_7_Mut_v_Con_chromvar)))
# MotifPlot(object = hm_endocrine,assay = 'peaks',
#           motifs = head(rownames(endo_8_Mut_v_Con_chromvar)))
# 


# FeaturePlot(hm_endocrine,"MA1642.1",split.by="genotype",
#             max.cutoff = "q95",min.cutoff = "q05",order =T,pt.size = .5)
# 
# FeaturePlot(hm_endocrine,"MA1642.1",split.by="genotype",
#             max.cutoff = "q95",min.cutoff = "q05", order =T,pt.size = .5)
# 
# FeaturePlot(hm_endocrine,"MA0105.4",split.by="genotype",
#             max.cutoff = "q95",min.cutoff = "q05", order =T,pt.size = .5)
# 
# FeaturePlot(hm_aggr_motifs,"MA0105.4",split.by="genotype",
#             max.cutoff = "q95",min.cutoff = "q05", order =T,pt.size = .5)
# 
# FeaturePlot(hm_aggr_motifs,"MA0132.2",split.by="genotype",
#             max.cutoff = "q90",min.cutoff = "q10", order =T,pt.size = .5)
# 
# 
# 
saveRDS(hm_endocrine,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_endocrine.rds")
# saveRDS(hm_aggr_motifs,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_motifs.rds")
# 
# 
# DimPlot(object = hm_endocrine, reduction = 'umap',label = TRUE)# + NoLegend()
# FeaturePlot(hm_endocrine,"Pdx1",order = T,split.by="genotype",pt.size = .5,cols = c("Light Grey","Black"),)
# FeaturePlot(hm_endocrine,"Ins2",order = T,split.by="genotype",pt.size = .5,cols = c("Light Grey","Black"),keep.scale = "all")
# 
# 
# hm_endocrine <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_endocrine.rds")


AddMotifs(hm_aggr_motifs)

hm_aggr_motifs <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/hm_aggr_motifs.rds")
DimPlot(hm_endocrine,split.by = "celltype")
DimPlot(hm_aggr_motifs)


DefaultAssay(hm_aggr_motifs) <- "chromvar"
all_motif_markers <- FindAllMarkers(hm_aggr_motifs,assay = "chromvar",min.pct = .05)
saveRDS(all_motif_markers,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/all_motif_markers.rds")

DefaultAssay(hm_aggr_motifs) <- "peaks"
all_peak_markers <- FindAllMarkers(hm_aggr_motifs,assay = "peaks",min.pct = .05)
saveRDS(all_peak_markers,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/all_peak_markers.rds")

DefaultAssay(hm_aggr_motifs) <- "RNA"
all_RNA_markers <- FindAllMarkers(hm_aggr_motifs,assay = "RNA",min.pct = .05)
saveRDS(all_RNA_markers,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/all_RNA_markers.rds")

DefaultAssay(hm_endocrine) <- "chromvar"
Idents(hm_endocrine) <- "seurat_clusters"
endo_motif_markers <- FindAllMarkers(hm_endocrine,assay = "chromvar",min.pct = .05)
saveRDS(endo_motif_markers,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_motif_markers.rds")

DefaultAssay(hm_endocrine) <- "peaks"
endo_peak_markers <- FindAllMarkers(hm_endocrine,assay = "peaks",min.pct = .05)
saveRDS(endo_peak_markers,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_peak_markers.rds")

DefaultAssay(hm_endocrine) <- "RNA"
endo_RNA_markers <- FindAllMarkers(hm_endocrine,assay = "RNA",min.pct = .05)
saveRDS(endo_RNA_markers,"/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_RNA_markers.rds")



# all_motif_markers <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/all_motif_markers.rds")
# all_peak_markers <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/all_peak_markers.rds")
# all_RNA_markers <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/all_RNA_markers.rds")
# endo_motif_markers <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_motif_markers.rds")
# endo_peak_markers <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_peak_markers.rds")
# endo_RNA_markers <- readRDS("/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_RNA_markers.rds")

DefaultAssay(hm_endocrine) <- 'RNA'
hm_endocrine <- ScaleData(hm_endocrine)
top10_RNA <- endo_RNA_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
hm_endocrine$celltype <- hm_endocrine$seurat_clusters
hm_endocrine$celltype.dx <- paste(hm_endocrine$celltype, hm_endocrine$genotype, sep = "_")
hm_endocrine$celltype.dx <- as.factor(hm_endocrine$celltype.dx)
# Idents(hm_endocrine) <- "celltype.dx"
# levels(hm_endocrine) <- paste0(rep(c(0,1,2,4,5,3,6),each=2),"_",rep(c("Con","Mut"),5))
Idents(hm_endocrine) <- "celltype"
levels(hm_endocrine) <- c(0,1,2,4,5,3,6)
p_endo_RNA <- DoHeatmap(hm_endocrine,features=top10_RNA$gene,disp.min = -2,disp.max = 2) + scale_fill_gradientn(colors = c("cornflowerblue", "black", "yellow"))
ggsave("/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_RNA_heatmap.pdf",plot = p_endo_RNA)


DefaultAssay(hm_endocrine) <- 'chromvar'
hm_endocrine <- ScaleData(hm_endocrine)
top10_chromvar <- endo_motif_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Idents(hm_endocrine) <- "celltype"
levels(hm_endocrine) <- c(0,1,2,4,5,3,6)
p_endo_chromvar <- DoHeatmap(hm_endocrine, features = top10_chromvar$gene,disp.min = -2,disp.max = 2)
ggsave("/home/bjw032/Documents/scATAC/01.10.2022_analysis/endo_chromvar_heatmap.pdf",plot = p_endo_chromvar)
FeaturePlot(hm_endocrine,features = "Rela")
VlnPlot(hm_endocrine,"MA0107.1")
VlnPlot(hm_endocrine,"MA0132.2")


motif.name <- name(pfm_all[features_Mut_v_Con_chromvar])


DefaultAssay(hm_endocrine) <- 'peaks'
hm_endocrine <- ScaleData(hm_endocrine)
top10 <- endo_peak_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(hm_endocrine, features = top10$gene) + scale_fill_gradientn(colors = c("white", "red"))

DoHeatmap(hm_endocrine, assay="peaks",features = top10$gene,group.by = "celltype.dx") + NoLegend()
DoHeatmap(hm_endocrine, features = top10$gene,group.by = "celltype") + NoLegend()
rownames(hm_endocrine)

signficance_df <- data.frame(row.names = features_Mut_v_Con_chromvar,
                             endo_0_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_1_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_2_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_3_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_4_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_5_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_6_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_7_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"],
                             endo_8_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"p_val_adj"])


signficance_log2FC <- data.frame(row.names = features_Mut_v_Con_chromvar,
                                 endo_0_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_1_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_2_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_3_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_4_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_5_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_6_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_7_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"],
                                 endo_8_Mut_v_Con_chromvar[features_Mut_v_Con_chromvar,"avg_log2FC"])


paletteLength <- 50
myColor <- colorRampPalette(c("white", "red", "darkred"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(0, -log10(.05), length.out=ceiling(paletteLength/max(-log10(signficance_df),na.rm = T)) + 1),
              seq(max(-log10(signficance_df),na.rm = T)/paletteLength, max(-log10(signficance_df),na.rm = T), length.out=floor(paletteLength)))
pheatmap(-log10(signficance_df),cluster_cols = F,cluster_rows = F,show_colnames = F,
         color = myColor,breaks = myBreaks)

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/4 + 1)),
              seq(max(signficance_log2FC,na.rm = T)/paletteLength, max(signficance_log2FC,na.rm = T), length.out=floor(paletteLength)))
myBreaks <- seq(-4,4, length.out=floor(paletteLength))

pheatmap(signficance_log2FC,cluster_cols = F,cluster_rows = F,show_colnames = F,
         color = myColor,breaks = myBreaks)


# frags <- Fragments(hm_endocrine)  # get list of fragment objects
# Fragments(hm_endocrine) <- NULL  # remove fragment information from assay
# new.paths <- '/media/bjw032/seagate/scATAC/Pdx1A4/Pdx1A4_aggr_2/outs/fragments.tsv.gz'
# for (i in seq_along(frags)) {
#   frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
# }
# Fragments(hm_endocrine) <- frags # assign updated list back to the object
# 

VlnPlot(hm_endocrine,"MA0105.4")
VlnPlot(hm_endocrine,"MA0132.2")
VlnPlot(hm_endocrine,"MA1618.1")

Idents()
hm_endocrine$geno_celltype <- paste0(hm_endocrine$genotype,"-",hm_endocrine$new_idents)
VlnPlot(hm_endocrine,"Dusp5",group.by = "geno_celltype")
VlnPlot(hm_endocrine,"Dusp5",group.by = "celltype.dx")
DefaultAssay(hm_endocrine) <- "peaks"
CoveragePlot(hm_endocrine,"Dusp5",idents = c("1_Con","1_Mut"),extend.upstream = 1000,extend.downstream = 50000)

VlnPlot(hm_endocrine,"Fosb",group.by = "genotype")

DotPlot(hm_endocrine,feature="Ins1")

hm_aggr_motifs$celltype <- hm_aggr_motifs$seurat_clusters
hm_aggr_motifs$celltype.dx <- paste(hm_aggr_motifs$celltype, hm_aggr_motifs$genotype, sep = "_")
Idents(hm_aggr_motifs) <- "celltype.dx"
Idents(hm_aggr_motifs) <- "celltype"
DefaultAssay(hm_aggr_motifs) <- 'RNA'
hm_aggr_motifs <- ScaleData(hm_aggr_motifs)
top10 <- all_RNA_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(hm_aggr_motifs,features=top10$gene, group.by = "celltype") + NoLegend()

Idents(hm_aggr_motifs) <- "celltype"
DefaultAssay(hm_aggr_motifs) <- 'chromvar'
hm_aggr_motifs <- ScaleData(hm_aggr_motifs)
top10 <- all_motif_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(hm_aggr_motifs,features=top10$gene, group.by = "celltype") + NoLegend()


Idents(hm_endocrine) <- "celltype"
ScaleData(hm_endocrine)
hm_endocrine <- ScaleData(hm_endocrine)
DefaultAssay(hm_endocrine) <- "peaks"
DoHeatmap(hm_endocrine,assay = "peaks")

VlnPlot(hm_endocrine,"MA1642.1")
VlnPlot(hm_endocrine,"MA0603.1")
VlnPlot(hm_endocrine,"MA0819.1",split.by = "genotype")
VlnPlot(hm_endocrine,"MA0136.2",ident=0,split.by = "genotype")
VlnPlot(hm_endocrine,"MA0819.1",ident=0,split.by = "genotype")

DotPlot(hm_endocrine,feature="Ins1")
GetMo

MA0603.1

MA0136.2
MA1642.1



for(f in unique(hm_aggr_motifs$seurat_clusters)){
  assign(paste0("agrr_i",0,"_","i",f),
         FindMarkers(
           object = hm_aggr_motifs,
           ident.1 = '0',
           ident.2 = f,
           min.pct = 0.05,
           test.use = 'LR',
           latent.vars = 'peak_region_fragments'
         ))
}



all_chromvar <- FindAllMarkers(hm_endocrine,test.use = 'LR',latent.vars = "peak_region_fragments")

shared_motifs <- intersect(rownames(i0_i1),
                           intersect(rownames(i0_i2),
                                     intersect(rownames(i0_i3),
                                               intersect(rownames(i0_i4),
                                                         intersect(rownames(i0_i5),
                                                                   intersect(rownames(i0_i6),
                                                                             intersect(rownames(i0_i7),rownames(i0_i8))))))))

cross_compare_chromvar <- data.frame(row.names = shared_motifs,
                                     i0_i1[shared_motifs,"pct.1"],
                                     i0_i1[shared_motifs,"pct.2"],
                                     i0_i2[shared_motifs,"pct.2"],
                                     i0_i3[shared_motifs,"pct.2"],
                                     i0_i4[shared_motifs,"pct.2"],
                                     i0_i5[shared_motifs,"pct.2"],
                                     i0_i6[shared_motifs,"pct.2"],
                                     i0_i7[shared_motifs,"pct.2"],
                                     i0_i8[shared_motifs,"pct.2"]
)
library(pheatmap)
pheatmap(cross_compare_chromvar,scale="row",cluster_rows = F)  



assign("test",12)
i1_i2 <- FindMarkers(
  object = hm_endocrine,
  ident.1 = '2',
  ident.2 = '3',
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)


remove(hm_endocrine)


Idents(hm_endocrine)
# add motif information

setwd("~")


hm_aggr_motifs$celltype.dx <- paste(Idents(hm_aggr_motifs), hm_aggr_motifs$genotype, sep = "_")
hm_aggr_motifs$celltype <- Idents(hm_aggr_motifs)
for(f in unique(hm_aggr_motifs$seurat_clusters)){
  assign(paste0("aggr_",f,"_Mut_v_Con"),
         FindMarkers(
           object = hm_aggr_motifs,
           ident.1 = paste0(f,"_","Mut"),
           ident.2 = paste0(f,"_","Con"),
           min.pct = 0.05,
           test.use = 'LR',logfc.threshold = 0.1,
           latent.vars = 'peak_region_fragments'
         ))
}
for(f in unique(hm_aggr_motifs$seurat_clusters)){
  assign(paste0("aggr_",f,"_Mut_v_Con_motifs"),
         FindMarkers(
           object = hm_aggr_motifs,
           assay = "chromvar",
           ident.1 = paste0(f,"_","Mut"),
           ident.2 = paste0(f,"_","Con"),
           min.pct = 0.05,
           test.use = 'LR',logfc.threshold = 0.0,
           latent.vars = 'peak_region_fragments'
         ))
}
hm_endocrine$geno_celltype <- paste0(hm_endocrine$genotype,"-",hm_endocrine$new_idents)

CoveragePlot(hm_aggr_motifs,idents = c("0_Mut","0_Con"),region = "chr19-53836409-53837210 ",
             extend.downstream = 1000,extend.upstream = 1000)

Idents(hm_endocrine) <- "geno_celltype"

CoveragePlot(hm_endocrine,"Ins1",
             extend.downstream = 1000,extend.upstream = 1000)

CoveragePlot(hm_endocrine,"Pdx1",
             extend.downstream = 1000,extend.upstream = 10000)


for(f in 1:10){
  ggsave(paste0("/media/bjw032/seagate/scATAC/Pdx1A4/Pdx1A4_aggr_2/pdfs/",rownames(da_peaks)[f],".jpeg"),
         CoveragePlot(
           object = hm_aggr,
           region = rownames(da_peaks)[f],
           extend.upstream = 100000,
           extend.downstream = 200000,
           idents = c("0","3")),
         device="jpeg")
}

DefaultAssay(hm_aggr_motifs) <- "peaks"
hm_aggr_cds <- as.CellDataSet(x = hm_aggr_motifs)
hm_aggr_cicero <- make_cicero_cds(hm_aggr_cds, reduced_coordinates = cbind(hm_aggr_motifs$umap_1,hm_aggr_motifs$umap_2))
genome <- seqlengths(hm_aggr_motifs)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns <- run_cicero(hm_aggr_cicero, genome.df, sample_num = 100) 
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns,ccans = ccans)
Links(hm_aggr_motifs) <- links


CoveragePlot(hm_aggr_motifs,"Ins1",
             extend.downstream = 1000,extend.upstream = 1000)

CoveragePlot(hm_endocrine,"Pdx1",
             extend.downstream = 1000,extend.upstream = 10000)




