### Single cell atac UVA

library(Seurat)
gtfv32 <-'~/Desktop/scatac/gencode.v32.annotation.gtf'
Chromosomes <-c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

##### scA #####

scA = "~/Desktop/scatac/scA/filtered_peak_bc_matrix.h5"
scA <-Read10X_h5(scA)
scA.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scA, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scA.atac <- CreateSeuratObject(counts = scA, assay = "ATAC", project = "scATAC")
scA.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scA.activity.matrix)
scA.meta <- read.table("~/Desktop/scatac/scA/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scA.meta <- scA.meta[colnames(scA.atac), ]
scA.atac <- AddMetaData(scA.atac, metadata = scA.meta)
scA.atac <- subset(scA.atac, subset = nCount_ATAC > 5000)
scA.atac$tech <- "atac"

#### scC
scC = "~/Desktop/scatac/scC/filtered_peak_bc_matrix.h5"
scC <-Read10X_h5(scC)
scC.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scC, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scC.atac <- CreateSeuratObject(counts = scC, assay = "ATAC", project = "scCTAC")
scC.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scC.activity.matrix)
scC.meta <- read.table("~/Desktop/scatac/scC/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scC.meta <- scC.meta[colnames(scC.atac), ]
scC.atac <- AddMetaData(scC.atac, metadata = scC.meta)
scC.atac <- subset(scC.atac, subset = nCount_ATAC > 5000)
scC.atac$tech <- "atac"

#### scD
scD = "~/Desktop/scatac/scD/filtered_peak_bc_matrix.h5"
scD <-Read10X_h5(scD)
scD.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scD, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scD.atac <- CreateSeuratObject(counts = scD, assay = "ATAC", project = "scDTAC")
scD.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scD.activity.matrix)
scD.meta <- read.table("~/Desktop/scatac/scD/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scD.meta <- scD.meta[colnames(scD.atac), ]
scD.atac <- AddMetaData(scD.atac, metadata = scD.meta)
scD.atac <- subset(scD.atac, subset = nCount_ATAC > 5000)
scD.atac$tech <- "atac"

#### scE
scE = "~/Desktop/scatac/scE/filtered_peak_bc_matrix.h5"
scE <-Read10X_h5(scE)
scE.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scE, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scE.atac <- CreateSeuratObject(counts = scE, assay = "ATAC", project = "scETAC")
scE.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scE.activity.matrix)
scE.meta <- read.table("~/Desktop/scatac/scE/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scE.meta <- scE.meta[colnames(scE.atac), ]
scE.atac <- AddMetaData(scE.atac, metadata = scE.meta)
scE.atac <- subset(scE.atac, subset = nCount_ATAC > 5000)
scE.atac$tech <- "atac"

#### scF

scF = "~/Desktop/scatac/scF/filtered_peak_bc_matrix.h5"
scF <-Read10X_h5(scF)
scF.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scF, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scF.atac <- CreateSeuratObject(counts = scF, assay = "ATAC", project = "scFTAC")
scF.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scF.activity.matrix)
scF.meta <- read.table("~/Desktop/scatac/scF/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scF.meta <- scF.meta[colnames(scF.atac), ]
scF.atac <- AddMetaData(scF.atac, metadata = scF.meta)
scF.atac <- subset(scF.atac, subset = nCount_ATAC > 5000)
scF.atac$tech <- "atac"

####scG

scG = "~/Desktop/scatac/scG/filtered_peak_bc_matrix.h5"
scG <-Read10X_h5(scG)
scG.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scG, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scG.atac <- CreateSeuratObject(counts = scG, assay = "ATAC", project = "scGTAC")
scG.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scG.activity.matrix)
scG.meta <- read.table("~/Desktop/scatac/scG/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scG.meta <- scG.meta[colnames(scG.atac), ]
scG.atac <- AddMetaData(scG.atac, metadata = scG.meta)
scG.atac <- subset(scG.atac, subset = nCount_ATAC > 5000)
scG.atac$tech <- "atac"

##### scH #####

scH = "~/Desktop/scatac/scH/filtered_peak_bc_matrix.h5"
scH <-Read10X_h5(scH)
scH.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scH, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scH.atac <- CreateSeuratObject(counts = scH, assay = "ATAC", project = "scHTAC")
scH.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scH.activity.matrix)
scH.meta <- read.table("~/Desktop/scatac/scH/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                       stringsAsFactors = FALSE)
scH.meta <- scH.meta[colnames(scH.atac), ]
scH.atac <- AddMetaData(scH.atac, metadata = scH.meta)
scH.atac <- subset(scH.atac, subset = nCount_ATAC > 5000)
scH.atac$tech <- "atac"

##### scI_A11 #####

scI_A11 = "~/Desktop/scatac/scI_A11/filtered_peak_bc_matrix.h5"
scI_A11 <-Read10X_h5(scI_A11)
scI_A11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scI_A11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scI_A11.atac <- CreateSeuratObject(counts = scI_A11, assay = "ATAC", project = "scI_A11TAC")
scI_A11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scI_A11.activity.matrix)
scI_A11.meta <- read.table("~/Desktop/scatac/scI_A11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scI_A11.meta <- scI_A11.meta[colnames(scI_A11.atac), ]
scI_A11.atac <- AddMetaData(scI_A11.atac, metadata = scI_A11.meta)
scI_A11.atac <- subset(scI_A11.atac, subset = nCount_ATAC > 5000)
scI_A11.atac$tech <- "atac"

##### scJ_B11 #####

scJ_B11 = "~/Desktop/scatac/scJ_B11/filtered_peak_bc_matrix.h5"
scJ_B11 <-Read10X_h5(scJ_B11)
scJ_B11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scJ_B11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scJ_B11.atac <- CreateSeuratObject(counts = scJ_B11, assay = "ATAC", project = "scJ_B11TAC")
scJ_B11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scJ_B11.activity.matrix)
scJ_B11.meta <- read.table("~/Desktop/scatac/scJ_B11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scJ_B11.meta <- scJ_B11.meta[colnames(scJ_B11.atac), ]
scJ_B11.atac <- AddMetaData(scJ_B11.atac, metadata = scJ_B11.meta)
scJ_B11.atac <- subset(scJ_B11.atac, subset = nCount_ATAC > 5000)
scJ_B11.atac$tech <- "atac"

##### scK_C11 #####

scK_C11 = "~/Desktop/scatac/scK_C11/filtered_peak_bc_matrix.h5"
scK_C11 <-Read10X_h5(scK_C11)
scK_C11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scK_C11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scK_C11.atac <- CreateSeuratObject(counts = scK_C11, assay = "ATAC", project = "scK_C11TAC")
scK_C11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scK_C11.activity.matrix)
scK_C11.meta <- read.table("~/Desktop/scatac/scK_C11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scK_C11.meta <- scK_C11.meta[colnames(scK_C11.atac), ]
scK_C11.atac <- AddMetaData(scK_C11.atac, metadata = scK_C11.meta)
scK_C11.atac <- subset(scK_C11.atac, subset = nCount_ATAC > 5000)
scK_C11.atac$tech <- "atac"

##### scL_D11 #####

scL_D11 = "~/Desktop/scatac/scL_D11/filtered_peak_bc_matrix.h5"
scL_D11 <-Read10X_h5(scL_D11)
scL_D11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scL_D11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scL_D11.atac <- CreateSeuratObject(counts = scL_D11, assay = "ATAC", project = "scL_D11TAC")
scL_D11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scL_D11.activity.matrix)
scL_D11.meta <- read.table("~/Desktop/scatac/scL_D11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scL_D11.meta <- scL_D11.meta[colnames(scL_D11.atac), ]
scL_D11.atac <- AddMetaData(scL_D11.atac, metadata = scL_D11.meta)
scL_D11.atac <- subset(scL_D11.atac, subset = nCount_ATAC > 1000)
scL_D11.atac$tech <- "atac"

##### scM_E11 #####

scM_E11 = "~/Desktop/scatac/scM_E11/filtered_peak_bc_matrix.h5"
scM_E11 <-Read10X_h5(scM_E11)
scM_E11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scM_E11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scM_E11.atac <- CreateSeuratObject(counts = scM_E11, assay = "ATAC", project = "scM_E11TAC")
scM_E11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scM_E11.activity.matrix)
scM_E11.meta <- read.table("~/Desktop/scatac/scM_E11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scM_E11.meta <- scM_E11.meta[colnames(scM_E11.atac), ]
scM_E11.atac <- AddMetaData(scM_E11.atac, metadata = scM_E11.meta)
scM_E11.atac <- subset(scM_E11.atac, subset = nCount_ATAC > 5000)
scM_E11.atac$tech <- "atac"

##### scO_G11 #####

scO_G11 = "~/Desktop/scatac/scO_G11/filtered_peak_bc_matrix.h5"
scO_G11 <-Read10X_h5(scO_G11)
scO_G11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scO_G11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scO_G11.atac <- CreateSeuratObject(counts = scO_G11, assay = "ATAC", project = "scO_G11TAC")
scO_G11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scO_G11.activity.matrix)
scO_G11.meta <- read.table("~/Desktop/scatac/scO_G11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scO_G11.meta <- scO_G11.meta[colnames(scO_G11.atac), ]
scO_G11.atac <- AddMetaData(scO_G11.atac, metadata = scO_G11.meta)
scO_G11.atac <- subset(scO_G11.atac, subset = nCount_ATAC > 5000)
scO_G11.atac$tech <- "atac"

##### scP_H11 #####

scP_H11 = "~/Desktop/scatac/scP_H11/filtered_peak_bc_matrix.h5"
scP_H11 <-Read10X_h5(scP_H11)
scP_H11.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scP_H11, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scP_H11.atac <- CreateSeuratObject(counts = scP_H11, assay = "ATAC", project = "scP_H11TAC")
scP_H11.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scP_H11.activity.matrix)
scP_H11.meta <- read.table("~/Desktop/scatac/scP_H11/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scP_H11.meta <- scP_H11.meta[colnames(scP_H11.atac), ]
scP_H11.atac <- AddMetaData(scP_H11.atac, metadata = scP_H11.meta)
scP_H11.atac <- subset(scP_H11.atac, subset = nCount_ATAC > 5000)
scP_H11.atac$tech <- "atac"

##### scQ_A12 #####

scQ_A12 = "~/Desktop/scatac/scQ_A12/filtered_peak_bc_matrix.h5"
scQ_A12 <-Read10X_h5(scQ_A12)
scQ_A12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scQ_A12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scQ_A12.atac <- CreateSeuratObject(counts = scQ_A12, assay = "ATAC", project = "scQ_A12TAC")
scQ_A12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scQ_A12.activity.matrix)
scQ_A12.meta <- read.table("~/Desktop/scatac/scQ_A12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scQ_A12.meta <- scQ_A12.meta[colnames(scQ_A12.atac), ]
scQ_A12.atac <- AddMetaData(scQ_A12.atac, metadata = scQ_A12.meta)
scQ_A12.atac <- subset(scQ_A12.atac, subset = nCount_ATAC > 5000)
scQ_A12.atac$tech <- "atac"

##### scR_B12 #####

scR_B12 = "~/Desktop/scatac/scR_B12/filtered_peak_bc_matrix.h5"
scR_B12 <-Read10X_h5(scR_B12)
scR_B12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scR_B12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scR_B12.atac <- CreateSeuratObject(counts = scR_B12, assay = "ATAC", project = "scR_B12TAC")
scR_B12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scR_B12.activity.matrix)
scR_B12.meta <- read.table("~/Desktop/scatac/scR_B12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scR_B12.meta <- scR_B12.meta[colnames(scR_B12.atac), ]
scR_B12.atac <- AddMetaData(scR_B12.atac, metadata = scR_B12.meta)
scR_B12.atac <- subset(scR_B12.atac, subset = nCount_ATAC > 5000)
scR_B12.atac$tech <- "atac"

##### scS_C12 #####

scS_C12 = "~/Desktop/scatac/scS_C12/filtered_peak_bc_matrix.h5"
scS_C12 <-Read10X_h5(scS_C12)
scS_C12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scS_C12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scS_C12.atac <- CreateSeuratObject(counts = scS_C12, assay = "ATAC", project = "scS_C12TAC")
scS_C12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scS_C12.activity.matrix)
scS_C12.meta <- read.table("~/Desktop/scatac/scS_C12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scS_C12.meta <- scS_C12.meta[colnames(scS_C12.atac), ]
scS_C12.atac <- AddMetaData(scS_C12.atac, metadata = scS_C12.meta)
scS_C12.atac <- subset(scS_C12.atac, subset = nCount_ATAC > 5000)
scS_C12.atac$tech <- "atac"

##### scT_D12 #####

scT_D12 = "~/Desktop/scatac/scT_D12/filtered_peak_bc_matrix.h5"
scT_D12 <-Read10X_h5(scT_D12)
scT_D12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scT_D12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scT_D12.atac <- CreateSeuratObject(counts = scT_D12, assay = "ATAC", project = "scT_D12TAC")
scT_D12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scT_D12.activity.matrix)
scT_D12.meta <- read.table("~/Desktop/scatac/scT_D12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scT_D12.meta <- scT_D12.meta[colnames(scT_D12.atac), ]
scT_D12.atac <- AddMetaData(scT_D12.atac, metadata = scT_D12.meta)
scT_D12.atac <- subset(scT_D12.atac, subset = nCount_ATAC > 5000)
scT_D12.atac$tech <- "atac"

##### scU_E12 #####

scU_E12 = "~/Desktop/scatac/scU_E12/filtered_peak_bc_matrix.h5"
scU_E12 <-Read10X_h5(scU_E12)
scU_E12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scU_E12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scU_E12.atac <- CreateSeuratObject(counts = scU_E12, assay = "ATAC", project = "scU_E12TAC")
scU_E12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scU_E12.activity.matrix)
scU_E12.meta <- read.table("~/Desktop/scatac/scU_E12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scU_E12.meta <- scU_E12.meta[colnames(scU_E12.atac), ]
scU_E12.atac <- AddMetaData(scU_E12.atac, metadata = scU_E12.meta)
scU_E12.atac <- subset(scU_E12.atac, subset = nCount_ATAC > 5000)
scU_E12.atac$tech <- "atac"

##### scV_F12 #####

scV_F12 = "~/Desktop/scatac/scV_F12/filtered_peak_bc_matrix.h5"
scV_F12 <-Read10X_h5(scV_F12)
scV_F12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scV_F12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scV_F12.atac <- CreateSeuratObject(counts = scV_F12, assay = "ATAC", project = "scV_F12TAC")
scV_F12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scV_F12.activity.matrix)
scV_F12.meta <- read.table("~/Desktop/scatac/scV_F12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scV_F12.meta <- scV_F12.meta[colnames(scV_F12.atac), ]
scV_F12.atac <- AddMetaData(scV_F12.atac, metadata = scV_F12.meta)
scV_F12.atac <- subset(scV_F12.atac, subset = nCount_ATAC > 5000)
scV_F12.atac$tech <- "atac"

##### scW_G12 #####

scW_G12 = "~/Desktop/scatac/scW_G12/filtered_peak_bc_matrix.h5"
scW_G12 <-Read10X_h5(scW_G12)
scW_G12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scW_G12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scW_G12.atac <- CreateSeuratObject(counts = scW_G12, assay = "ATAC", project = "scW_G12TAC")
scW_G12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scW_G12.activity.matrix)
scW_G12.meta <- read.table("~/Desktop/scatac/scW_G12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scW_G12.meta <- scW_G12.meta[colnames(scW_G12.atac), ]
scW_G12.atac <- AddMetaData(scW_G12.atac, metadata = scW_G12.meta)
scW_G12.atac <- subset(scW_G12.atac, subset = nCount_ATAC > 5000)
scW_G12.atac$tech <- "atac"

##### scX_H12 #####

scX_H12 = "~/Desktop/scatac/scX_H12/filtered_peak_bc_matrix.h5"
scX_H12 <-Read10X_h5(scX_H12)
scX_H12.activity.matrix <- CreateGeneActivityMatrix(peak.matrix = scX_H12, annotation.file = gtfv32 , seq.levels =Chromosomes, upstream = 2000, verbose = TRUE)

scX_H12.atac <- CreateSeuratObject(counts = scX_H12, assay = "ATAC", project = "scX_H12TAC")
scX_H12.atac[["ACTIVITY"]] <- CreateAssayObject(counts = scX_H12.activity.matrix)
scX_H12.meta <- read.table("~/Desktop/scatac/scX_H12/singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                           stringsAsFactors = FALSE)
scX_H12.meta <- scX_H12.meta[colnames(scX_H12.atac), ]
scX_H12.atac <- AddMetaData(scX_H12.atac, metadata = scX_H12.meta)
scX_H12.atac <- subset(scX_H12.atac, subset = nCount_ATAC > 5000)
scX_H12.atac$tech <- "atac"


### combine all atac runs into one seurat object
atac.combined <- merge(scA.atac, y = c(scC.atac, scD.atac, scE.atac, scF.atac, scG.atac, scH.atac, scI_A11.atac,
                                  scJ_B11.atac, scK_C11.atac,scL_D11.atac, scM_E11.atac,  scO_G11.atac, scP_H11.atac,
                                  scQ_A12.atac,scR_B12.atac,  scS_C12.atac,scT_D12.atac, scU_E12.atac, scV_F12.atac, scW_G12.atac, scX_H12.atac ),
                                  add.cell.ids = c("scA", "scC", "scD", "scE", "scF", "scG", "scH", "scI_A11",
                                  "scJ_B11", "scK_C11","scL_D11", "scM_E11",  "scO_G11", "scP_H11",
                                  "scQ_A12","scR_B12",  "scS_C12","scT_D12", "scU_E12", "scV_F12", "scW_G12", "scX_H12" ), project = "UVA-ATAC")

### Save combined seurat object for later use 
saveRDS(atac.combined, file="atac.combined.rds")

