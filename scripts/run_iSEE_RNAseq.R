library(iSEE)
## BiocManager::install("Gviz") # if not already installed

sce <- readRDS("/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/outputR/shiny_sce.rds")
sce <- sce$sce_gene
# rownames(sce) <- paste0(rowData(sce)$gene_id, "__", rowData(sce)$symbol)
# rowData(sce) <- tidyr::unnest(as.data.frame(rowData(sce)))
# rowData(sce)$edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.mlog10PValue <- 
#   -log10(rowData(sce)$edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.PValue)
# 
# reddim <- redDimPlotDefaults(sce, 5)
# reddim$PointSize <- 5
# reddim$ColorBy <- "Column data"
# reddim$ColorByColData <- "condition"
# 
# rowdata <- rowDataPlotDefaults(sce, 5)
# rowdata$YAxis <- "edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.mlog10PValue"
# rowdata$XAxis <- "Row data"
# rowdata$XAxisRowData <- "edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.logFC"
# 
# rowstat <- rowStatTableDefaults(sce, 5)
# rowstat[['SelectByPlot']] <- c("Row data plot 1", "---", "---", "---", "---")
# rowstat[['Selected']] <- c(1L, 2238L, 1L, 1L, 1L)
# 
# featassay <- featAssayPlotDefaults(sce, 5)
# featassay$XAxis <- "Column data"
# featassay$XAxisColData <- "names"
# featassay$Assay <- 4L
# featassay$PointSize <- 4
# featassay[['YAxisRowTable']] <- c("Row statistics table 2", "---", "---", "---", "---")

source(here("scripts/custom_iSEE_panels.R"))

gtf <- prepareGtf("/Volumes/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf")
saveRDS(gtf, file = "/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/reference/Mus_musculus.GRCm38.90.gtf.rds")

cdp <- customDataPlotDefaults(sce, 2)
cdp$Function <- c("customGviz")

cdp$Arguments <- c("bigwig_files 1_mo_KI_H1_lib350686_6534_Aligned.sortedByCoord.out.bw,1_mo_KI_H2_lib350687_6534_Aligned.sortedByCoord.out.bw,1_mo_KI_H3_lib350688_6534_Aligned.sortedByCoord.out.bw,1_mo_KI_H4_lib350689_6534_Aligned.sortedByCoord.out.bw,1_mo_KI_H5_lib350690_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_H6_lib350691_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_SNS1_lib350692_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_SNS2_lib350693_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_SNS3_lib350694_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_SNS4_lib350695_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_SNS5_lib350696_6537_Aligned.sortedByCoord.out.bw,1_mo_KI_SNS6_lib350697_6537_Aligned.sortedByCoord.out.bw,1_mo_WT_H1_lib350670_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_H2_lib350671_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_H3_lib350672_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_H4_lib350673_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_H5_lib350674_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_H6_lib350675_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_N1_lib350682_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_N2_lib350683_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_N3_lib350684_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_N4_lib350685_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_SNS1_lib350676_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_SNS2_lib350677_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_SNS3_lib350678_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_SNS4_lib350679_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_SNS5_lib350680_6534_Aligned.sortedByCoord.out.bw,1_mo_WT_SNS6_lib350681_6534_Aligned.sortedByCoord.out.bw,6_mo_KI_H1_lib350710_6537_Aligned.sortedByCoord.out.bw,6_mo_KI_H2_lib350711_6537_Aligned.sortedByCoord.out.bw,6_mo_KI_H4_lib350713_6537_Aligned.sortedByCoord.out.bw,6_mo_KI_H5_lib350714_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_H7_lib350712_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_H8_lib350715_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_SNS1_lib350716_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_SNS2_lib350717_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_SNS4_lib350718_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_SNS5_lib350719_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_SNS7_lib350720_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_KI_SNS8_lib350721_6534_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_H1_lib350698_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_H2_lib350699_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_H4_lib350701_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_H5_lib350702_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_H6_lib350703_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_H7_lib351457_6534_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_SNS1_lib350704_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_SNS2_lib350705_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_SNS4_lib350706_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_SNS5_lib350707_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_SNS6_lib350708_6537_Aligned.sortedByCoord.out.bw,/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/STARbigwig/6_mo_WT_SNS7_lib350709_6537_Aligned.sortedByCoord.out.bw\nbigwig_names 1_mo_KI_H1,1_mo_KI_H2,1_mo_KI_H3,1_mo_KI_H4,1_mo_KI_H5,1_mo_KI_H6,1_mo_KI_SNS1,1_mo_KI_SNS2,1_mo_KI_SNS3,1_mo_KI_SNS4,1_mo_KI_SNS5,1_mo_KI_SNS6,1_mo_WT_H1,1_mo_WT_H2,1_mo_WT_H3,1_mo_WT_H4,1_mo_WT_H5,1_mo_WT_H6,1_mo_WT_N1,1_mo_WT_N2,1_mo_WT_N3,1_mo_WT_N4,1_mo_WT_SNS1,1_mo_WT_SNS2,1_mo_WT_SNS3,1_mo_WT_SNS4,1_mo_WT_SNS5,1_mo_WT_SNS6,6_mo_KI_H1,6_mo_KI_H2,6_mo_KI_H4,6_mo_KI_H5,6_mo_KI_H7,6_mo_KI_H8,6_mo_KI_SNS1,6_mo_KI_SNS2,6_mo_KI_SNS4,6_mo_KI_SNS5,6_mo_KI_SNS7,6_mo_KI_SNS8,6_mo_WT_H1,6_mo_WT_H2,6_mo_WT_H4,6_mo_WT_H5,6_mo_WT_H6,6_mo_WT_H7,6_mo_WT_SNS1,6_mo_WT_SNS2,6_mo_WT_SNS4,6_mo_WT_SNS5,6_mo_WT_SNS6,6_mo_WT_SNS7\ngranges /Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/reference/Mus_musculus.GRCm38.90.gtf.rds\nchr 1\nstart 6.1e6\nend 6.2e6\nshowgene FUS")
  


app <- iSEE(sce, 
            # redDimArgs = reddim,
            # rowDataArgs = rowdata,
            # rowStatArgs = rowstat, 
            # featAssayArgs = featassay, 
            customDataArgs = cdp, 
            customDataFun = list(customGviz = customGviz),
            initialPanels = DataFrame(
              Name = c("Reduced dimension plot 1", "Custom data plot 1",
                       "Row data plot 1", "Row statistics table 1",
                       "Feature assay plot 1", "Row statistics table 2"),
              Width = c(4, 8, 3, 3, 3, 3) ))
shiny::runApp(app)
