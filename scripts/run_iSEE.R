library(iSEE)
library(here)
################################################################################
# Settings for reduced dimension plots
################################################################################

redDimPlotArgs <- new('DataFrame', nrows=5L, rownames=sprintf('redDimPlot%i', seq_len(5)))
redDimPlotArgs[['Type']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['XAxis']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['YAxis']] <- c(2L, 2L, 2L, 2L, 2L)
redDimPlotArgs[['DataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['VisualBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SelectBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
redDimPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
redDimPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
redDimPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
redDimPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- list()
redDimPlotArgs[['MultiSelectHistory']] <- tmp
redDimPlotArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
redDimPlotArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

tmp <- vector('list', 5)
tmp[[1]] <- c("Color", "Shape", "Points", "Other")
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
redDimPlotArgs[['VisualChoices']] <- tmp
redDimPlotArgs[['PointSize']] <- c(4, 1, 1, 1, 1)
redDimPlotArgs[['PointAlpha']] <- c(0.9, 1, 1, 1, 1)
redDimPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
redDimPlotArgs[['FontSize']] <- c(2, 1, 1, 1, 1)
redDimPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
redDimPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
redDimPlotArgs[['LassoData']] <- tmp
redDimPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['ContourColor']] <- c("#0000FF", "blue", "blue", "blue", "blue")
redDimPlotArgs[['ColorBy']] <- c("Column data", "None", "None", "None", "None")
redDimPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
redDimPlotArgs[['ColorByColData']] <- c("group", "ID", "ID", "ID", "ID")
redDimPlotArgs[['ShapeBy']] <- c("Column data", "None", "None", "None", "None")
redDimPlotArgs[['ShapeByColData']] <- c("genotype", "type", "type", "type", "type")
redDimPlotArgs[['SizeBy']] <- c("None", "None", "None", "None", "None")
redDimPlotArgs[['SizeByColData']] <- c("ID", "ID", "ID", "ID", "ID")
redDimPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
redDimPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['ColorByFeatNameAssay']] <- c(1L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
redDimPlotArgs[['ColorBySampName']] <- c(53L, 1L, 1L, 1L, 1L)
redDimPlotArgs[['ColorBySampNameColor']] <- c("#FF0000", "red", "red", "red", "red")
redDimPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
redDimPlotArgs[['RowFacetByColData']] <- c("type", "type", "type", "type", "type")
redDimPlotArgs[['ColumnFacetByColData']] <- c("type", "type", "type", "type", "type")

################################################################################
# Settings for column data plots
################################################################################

colDataPlotArgs <- new('DataFrame', nrows=5L, rownames=sprintf('colDataPlot%i', seq_len(5)))
colDataPlotArgs[['YAxis']] <- c("age", "ID", "ID", "ID", "ID")
colDataPlotArgs[['XAxis']] <- c("Column data", "None", "None", "None", "None")
colDataPlotArgs[['XAxisColData']] <- c("group", "sample", "sample", "sample", "sample")
colDataPlotArgs[['DataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['VisualBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
colDataPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
colDataPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
colDataPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
colDataPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- list()
colDataPlotArgs[['MultiSelectHistory']] <- tmp
colDataPlotArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
colDataPlotArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

tmp <- vector('list', 5)
tmp[[1]] <- c("Color", "Shape", "Points")
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
colDataPlotArgs[['VisualChoices']] <- tmp
colDataPlotArgs[['PointSize']] <- c(2, 1, 1, 1, 1)
colDataPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
colDataPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
colDataPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
colDataPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
colDataPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
colDataPlotArgs[['LassoData']] <- tmp
colDataPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['ContourColor']] <- c("#0000FF", "blue", "blue", "blue", "blue")
colDataPlotArgs[['ColorBy']] <- c("Column data", "None", "None", "None", "None")
colDataPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
colDataPlotArgs[['ColorByColData']] <- c("fraction", "ID", "ID", "ID", "ID")
colDataPlotArgs[['ShapeBy']] <- c("None", "None", "None", "None", "None")
colDataPlotArgs[['ShapeByColData']] <- c("type", "type", "type", "type", "type")
colDataPlotArgs[['SizeBy']] <- c("None", "None", "None", "None", "None")
colDataPlotArgs[['SizeByColData']] <- c("ID", "ID", "ID", "ID", "ID")
colDataPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
colDataPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
colDataPlotArgs[['ColorByFeatNameAssay']] <- c(1L, 1L, 1L, 1L, 1L)
colDataPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
colDataPlotArgs[['ColorBySampName']] <- c(1L, 1L, 1L, 1L, 1L)
colDataPlotArgs[['ColorBySampNameColor']] <- c("#FF0000", "red", "red", "red", "red")
colDataPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
colDataPlotArgs[['RowFacetByColData']] <- c("type", "type", "type", "type", "type")
colDataPlotArgs[['ColumnFacetByColData']] <- c("type", "type", "type", "type", "type")

################################################################################
# Settings for feature assay plots
################################################################################

featAssayPlotArgs <- new('DataFrame', nrows=5L, rownames=sprintf('featAssayPlot%i', seq_len(5)))
featAssayPlotArgs[['Assay']] <- c(4L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['XAxis']] <- c("Column data", "None", "None", "None", "None")
featAssayPlotArgs[['XAxisColData']] <- c("group", "ID", "ID", "ID", "ID")
featAssayPlotArgs[['XAxisFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['XAxisRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['YAxisFeatName']] <- c(8144L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['YAxisRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['DataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['VisualBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SelectByPlot']] <- c("Reduced dimension plot 1", "---", "---", "---", "---")
featAssayPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
featAssayPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
featAssayPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
featAssayPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- list()
featAssayPlotArgs[['MultiSelectHistory']] <- tmp
featAssayPlotArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
featAssayPlotArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

tmp <- vector('list', 5)
tmp[[1]] <- c("Color", "Shape", "Points", "Other")
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
featAssayPlotArgs[['VisualChoices']] <- tmp
featAssayPlotArgs[['PointSize']] <- c(2, 1, 1, 1, 1)
featAssayPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
featAssayPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
featAssayPlotArgs[['FontSize']] <- c(1.5, 1, 1, 1, 1)
featAssayPlotArgs[['LegendPosition']] <- c("Right", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
featAssayPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
featAssayPlotArgs[['LassoData']] <- tmp
featAssayPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['ContourColor']] <- c("#0000FF", "blue", "blue", "blue", "blue")
featAssayPlotArgs[['ColorBy']] <- c("Column data", "None", "None", "None", "None")
featAssayPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
featAssayPlotArgs[['ColorByColData']] <- c("age", "ID", "ID", "ID", "ID")
featAssayPlotArgs[['ShapeBy']] <- c("Column data", "None", "None", "None", "None")
featAssayPlotArgs[['ShapeByColData']] <- c("genotype", "type", "type", "type", "type")
featAssayPlotArgs[['SizeBy']] <- c("None", "None", "None", "None", "None")
featAssayPlotArgs[['SizeByColData']] <- c("ID", "ID", "ID", "ID", "ID")
featAssayPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['ColorByFeatNameAssay']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
featAssayPlotArgs[['ColorBySampName']] <- c(1L, 1L, 1L, 1L, 1L)
featAssayPlotArgs[['ColorBySampNameColor']] <- c("#FF0000", "red", "red", "red", "red")
featAssayPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
featAssayPlotArgs[['RowFacetByColData']] <- c("type", "type", "type", "type", "type")
featAssayPlotArgs[['ColumnFacetByColData']] <- c("type", "type", "type", "type", "type")

################################################################################
# Settings for row statistics tables
################################################################################

rowStatTableArgs <- new('DataFrame', nrows=5L, rownames=sprintf('rowStatTable%i', seq_len(5)))
rowStatTableArgs[['Selected']] <- c(7989L, 1L, 1L, 1L, 1L)
rowStatTableArgs[['Search']] <- c("", "", "", "", "")

tmp <- vector('list', 5)
tmp[[1]] <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "[\"true\"]")
tmp[[2]] <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "")
tmp[[3]] <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "")
tmp[[4]] <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "")
tmp[[5]] <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
              "", "", "", "", "", "", "", "", "", "")
rowStatTableArgs[['SearchColumns']] <- tmp
rowStatTableArgs[['SelectBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
rowStatTableArgs[['SelectByPlot']] <- c("Row data plot 1", "---", "---", "---", "---")
rowStatTableArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
rowStatTableArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

################################################################################
# Settings for row data plots
################################################################################

rowDataPlotArgs <- new('DataFrame', nrows=5L, rownames=sprintf('rowDataPlot%i', seq_len(5)))
rowDataPlotArgs[['YAxis']] <- c("edgeR:groupSNS.WT.1_mo-groupH.WT.1_mo:logFC", "gene_id", "gene_id", "gene_id", 
                                "gene_id")
rowDataPlotArgs[['XAxis']] <- c("Row data", "None", "None", "None", "None")
rowDataPlotArgs[['XAxisRowData']] <- c("edgeR:groupSNS.WT.1_mo-groupH.WT.1_mo:logCPM", "gene_name", "gene_name", "gene_name", 
                                       "gene_name")
rowDataPlotArgs[['DataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['VisualBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['SelectBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
rowDataPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
rowDataPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
rowDataPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
tmp[[1]] <- list(xmin = 13.29670997673, xmax = 16.830792574485, ymin = 2.6875450324357, ymax = 4.0226033975128, 
                 coords_css = list(xmin = 496.03125, xmax = 587.03125, ymin = 107.734375, ymax = 131.734375), 
                 coords_img = list(xmin = 496.03125, xmax = 587.03125, ymin = 107.734375, ymax = 131.734375), 
                 img_css_ratio = list(x = 1L, y = 1L), mapping = list(x = "X", y = "Y", colour = "ColorBy"), 
                 domain = list(left = -4.4635, right = 17.3935, bottom = -12.984, top = 8.664), 
                 range = list(left = 38.7190443065069, right = 601.520547945205, bottom = 413.457699325771, 
                              top = 24.2971850062453), log = list(x = NULL, y = NULL), direction = "xy", 
                 brushId = "rowDataPlot1_Brush", outputId = "rowDataPlot1")
rowDataPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- list()
rowDataPlotArgs[['MultiSelectHistory']] <- tmp
rowDataPlotArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
rowDataPlotArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

tmp <- vector('list', 5)
tmp[[1]] <- "Color"
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
rowDataPlotArgs[['VisualChoices']] <- tmp
rowDataPlotArgs[['PointSize']] <- c(1, 1, 1, 1, 1)
rowDataPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
rowDataPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
rowDataPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
rowDataPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
rowDataPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
rowDataPlotArgs[['LassoData']] <- tmp
rowDataPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['ContourColor']] <- c("blue", "blue", "blue", "blue", "blue")
rowDataPlotArgs[['ColorBy']] <- c("Row data", "None", "None", "None", "None")
rowDataPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
rowDataPlotArgs[['ColorByRowData']] <- c("edgeR:groupSNS.WT.1_mo-groupH.WT.1_mo:mlog10PValue", "gene_id", "gene_id", "gene_id", 
                                         "gene_id")
rowDataPlotArgs[['ShapeBy']] <- c("None", "None", "None", "None", "None")
rowDataPlotArgs[['ShapeByRowData']] <- c("gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype")
rowDataPlotArgs[['SizeBy']] <- c("None", "None", "None", "None", "None")
rowDataPlotArgs[['SizeByRowData']] <- c("entrezid", "entrezid", "entrezid", "entrezid", "entrezid")
rowDataPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
rowDataPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
rowDataPlotArgs[['ColorByFeatNameColor']] <- c("#FF0000", "red", "red", "red", "red")
rowDataPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
rowDataPlotArgs[['ColorBySampName']] <- c(1L, 1L, 1L, 1L, 1L)
rowDataPlotArgs[['ColorBySampNameAssay']] <- c(1L, 1L, 1L, 1L, 1L)
rowDataPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
rowDataPlotArgs[['RowFacetByRowData']] <- c("gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype")
rowDataPlotArgs[['ColumnFacetByRowData']] <- c("gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype")

################################################################################
# Settings for sample assay plots
################################################################################

sampAssayPlotArgs <- new('DataFrame', nrows=5L, rownames=sprintf('sampAssayPlot%i', seq_len(5)))
sampAssayPlotArgs[['YAxisSampName']] <- c(7L, 1L, 1L, 1L, 1L)
sampAssayPlotArgs[['YAxisColTable']] <- c("---", "---", "---", "---", "---")
sampAssayPlotArgs[['Assay']] <- c(1L, 1L, 1L, 1L, 1L)
sampAssayPlotArgs[['XAxis']] <- c("Row data", "None", "None", "None", "None")
sampAssayPlotArgs[['XAxisRowData']] <- c("gene_biotype", "gene_id", "gene_id", "gene_id", "gene_id")
sampAssayPlotArgs[['XAxisSampName']] <- c(2L, 2L, 2L, 2L, 2L)
sampAssayPlotArgs[['XAxisColTable']] <- c("---", "---", "---", "---", "---")
sampAssayPlotArgs[['DataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['VisualBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['SelectBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['SelectByPlot']] <- c("Row data plot 1", "---", "---", "---", "---")
sampAssayPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
sampAssayPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
sampAssayPlotArgs[['SelectColor']] <- c("#FF0000", "red", "red", "red", "red")

tmp <- vector('list', 5)
sampAssayPlotArgs[['BrushData']] <- tmp

tmp <- vector('list', 5)
tmp[[1]] <- list()
sampAssayPlotArgs[['MultiSelectHistory']] <- tmp
sampAssayPlotArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
sampAssayPlotArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

tmp <- vector('list', 5)
tmp[[1]] <- "Color"
tmp[[2]] <- "Color"
tmp[[3]] <- "Color"
tmp[[4]] <- "Color"
tmp[[5]] <- "Color"
sampAssayPlotArgs[['VisualChoices']] <- tmp
sampAssayPlotArgs[['PointSize']] <- c(1, 1, 1, 1, 1)
sampAssayPlotArgs[['PointAlpha']] <- c(1, 1, 1, 1, 1)
sampAssayPlotArgs[['Downsample']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['SampleRes']] <- c(200, 200, 200, 200, 200)
sampAssayPlotArgs[['FontSize']] <- c(1, 1, 1, 1, 1)
sampAssayPlotArgs[['LegendPosition']] <- c("Bottom", "Bottom", "Bottom", "Bottom", "Bottom")

tmp <- vector('list', 5)
sampAssayPlotArgs[['ZoomData']] <- tmp

tmp <- vector('list', 5)
sampAssayPlotArgs[['LassoData']] <- tmp
sampAssayPlotArgs[['ContourAdd']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['ContourColor']] <- c("blue", "blue", "blue", "blue", "blue")
sampAssayPlotArgs[['ColorBy']] <- c("None", "None", "None", "None", "None")
sampAssayPlotArgs[['ColorByDefaultColor']] <- c("#000000", "black", "black", "black", "black")
sampAssayPlotArgs[['ColorByRowData']] <- c("gene_id", "gene_id", "gene_id", "gene_id", "gene_id")
sampAssayPlotArgs[['ShapeBy']] <- c("None", "None", "None", "None", "None")
sampAssayPlotArgs[['ShapeByRowData']] <- c("gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype")
sampAssayPlotArgs[['SizeBy']] <- c("None", "None", "None", "None", "None")
sampAssayPlotArgs[['SizeByRowData']] <- c("entrezid", "entrezid", "entrezid", "entrezid", "entrezid")
sampAssayPlotArgs[['ColorByRowTable']] <- c("---", "---", "---", "---", "---")
sampAssayPlotArgs[['ColorByFeatName']] <- c(1L, 1L, 1L, 1L, 1L)
sampAssayPlotArgs[['ColorByFeatNameColor']] <- c("#FF0000", "red", "red", "red", "red")
sampAssayPlotArgs[['ColorByColTable']] <- c("---", "---", "---", "---", "---")
sampAssayPlotArgs[['ColorBySampName']] <- c(1L, 1L, 1L, 1L, 1L)
sampAssayPlotArgs[['ColorBySampNameAssay']] <- c(1L, 1L, 1L, 1L, 1L)
sampAssayPlotArgs[['FacetByRow']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['FacetByColumn']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
sampAssayPlotArgs[['RowFacetByRowData']] <- c("gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype")
sampAssayPlotArgs[['ColumnFacetByRowData']] <- c("gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype", "gene_biotype")

################################################################################
# Settings for column statistics tables
################################################################################

colStatTableArgs <- new('DataFrame', nrows=5L, rownames=sprintf('colStatTable%i', seq_len(5)))
colStatTableArgs[['Selected']] <- c(1L, 1L, 1L, 1L, 1L)
colStatTableArgs[['Search']] <- c("", "", "", "", "")

tmp <- vector('list', 5)
tmp[[1]] <- c("", "", "", "", "", "", "", "", "", "")
tmp[[2]] <- c("", "", "", "", "", "", "", "", "", "")
tmp[[3]] <- c("", "", "", "", "", "", "", "", "", "")
tmp[[4]] <- c("", "", "", "", "", "", "", "", "", "")
tmp[[5]] <- c("", "", "", "", "", "", "", "", "", "")
colStatTableArgs[['SearchColumns']] <- tmp
colStatTableArgs[['SelectBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
colStatTableArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
colStatTableArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
colStatTableArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)

################################################################################
# Settings for custom data plots
################################################################################

customDataPlotArgs <- new('DataFrame', nrows=0L, rownames=sprintf('customDataPlot%i', seq_len(0)))
customDataPlotArgs[['Function']] <- character(0)
customDataPlotArgs[['Arguments']] <- character(0)
customDataPlotArgs[['VisibleArgs']] <- logical(0)
customDataPlotArgs[['ColumnSource']] <- character(0)
customDataPlotArgs[['RowSource']] <- character(0)
customDataPlotArgs[['DataBoxOpen']] <- logical(0)
customDataPlotArgs[['SelectBoxOpen']] <- logical(0)

################################################################################
# Settings for custom statistics tables
################################################################################

customStatTableArgs <- new('DataFrame', nrows=0L, rownames=sprintf('customStatTable%i', seq_len(0)))
customStatTableArgs[['Function']] <- character(0)
customStatTableArgs[['Arguments']] <- character(0)
customStatTableArgs[['VisibleArgs']] <- logical(0)
customStatTableArgs[['ColumnSource']] <- character(0)
customStatTableArgs[['RowSource']] <- character(0)
customStatTableArgs[['DataBoxOpen']] <- logical(0)
customStatTableArgs[['SelectBoxOpen']] <- logical(0)
customStatTableArgs[['Search']] <- character(0)

################################################################################
# Settings for heat maps
################################################################################

heatMapPlotArgs <- new('DataFrame', nrows=5L, rownames=sprintf('heatMapPlot%i', seq_len(5)))
heatMapPlotArgs[['Assay']] <- c(1L, 1L, 1L, 1L, 1L)
heatMapPlotArgs[['FeatNameBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)

tmp <- vector('list', 5)
tmp[[1]] <- c(17600L, 27123L, 17602L, 17605L, 28502L, 17611L, 28835L, 17601L, 17604L, 17609L, 
              17603L, 17606L, 17599L, 17598L, 17607L)
tmp[[2]] <- 1L
tmp[[3]] <- 1L
tmp[[4]] <- 1L
tmp[[5]] <- 1L
heatMapPlotArgs[['FeatName']] <- tmp
heatMapPlotArgs[['ColDataBoxOpen']] <- c(TRUE, FALSE, FALSE, FALSE, FALSE)

tmp <- vector('list', 5)
tmp[[1]] <- c("fraction", "genotype", "age")
tmp[[2]] <- "ID"
tmp[[3]] <- "ID"
tmp[[4]] <- "ID"
tmp[[5]] <- "ID"
heatMapPlotArgs[['ColData']] <- tmp
heatMapPlotArgs[['FeatNameSource']] <- c("Row data plot 1", "---", "---", "---", "---")

tmp <- vector('list', 5)
tmp[[1]] <- "Centered"
tmp[[2]] <- "Centered"
tmp[[3]] <- "Centered"
tmp[[4]] <- "Centered"
tmp[[5]] <- "Centered"
heatMapPlotArgs[['CenterScale']] <- tmp
heatMapPlotArgs[['Lower']] <- c(NA, -Inf, -Inf, -Inf, -Inf)
heatMapPlotArgs[['Upper']] <- c(NA, Inf, Inf, Inf, Inf)
heatMapPlotArgs[['ColorScale']] <- c("blue-white-orange", "purple-black-yellow", "purple-black-yellow", "purple-black-yellow", 
                                     "purple-black-yellow")

tmp <- vector('list', 5)
heatMapPlotArgs[['ZoomData']] <- tmp
heatMapPlotArgs[['SelectBoxOpen']] <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
heatMapPlotArgs[['SelectByPlot']] <- c("---", "---", "---", "---", "---")
heatMapPlotArgs[['SelectEffect']] <- c("Transparent", "Transparent", "Transparent", "Transparent", "Transparent")
heatMapPlotArgs[['SelectAlpha']] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
heatMapPlotArgs[['SelectColor']] <- c("red", "red", "red", "red", "red")
heatMapPlotArgs[['SelectMultiType']] <- c("Active", "Active", "Active", "Active", "Active")
heatMapPlotArgs[['SelectMultiSaved']] <- c(0L, 0L, 0L, 0L, 0L)


################################################################################
# Initial panel settings
################################################################################

initialPanels <- DataFrame(
  Name=c("Reduced dimension plot 1", "Column data plot 1", "Feature assay plot 1", "Row statistics table 1", 
         "Row data plot 1", "Sample assay plot 1", "Column statistics table 1", "Heat map 1"
  ),
  Width=c(6L, 4L, 4L, 4L, 4L, 4L, 4L, 4L),
  Height=c(600L, 500L, 500L, 500L, 500L, 500L, 500L, 500L)
)

################################################################################
# Start the app
################################################################################

# sce <- readRDS("/Volumes/Shared-1/data/seq/sonu_RNAseq/polyA_Aug2019/output/outputR/shiny_sce.rds")
sce <- readRDS("/Volumes/seq/sonu_RNAseq/polyA_Aug2019/output/outputR/shiny_sce.rds")



# set_here(path = "/Volumes/Shared/data/seq/sonu_RNAseq/polyA_Aug2019/")
# sce <- readRDS(here("output/outputR/shiny_sce.rds" ))


# sce on gene level
sce_g <- sce$sce_gene

app <- iSEE(sce_g, redDimArgs = redDimPlotArgs, colDataArgs = colDataPlotArgs, 
            featAssayArgs = featAssayPlotArgs, rowStatArgs = rowStatTableArgs, 
            rowDataArgs = rowDataPlotArgs, sampAssayArgs = sampAssayPlotArgs, 
            colStatArgs = colStatTableArgs, customDataArgs = customDataPlotArgs, 
            customStatArgs = customStatTableArgs, heatMapArgs = heatMapPlotArgs, 
            initialPanels = initialPanels)
            
shiny::runApp(app)
