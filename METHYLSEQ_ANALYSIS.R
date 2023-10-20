#!/usr/bin/env Rscript
################################################################################
### R script to compare two different conditions with count file (generated using corset and salmon) and DESeq2 package
### Aditya Narayan Sarangi
### Designed to be executed with bulkRNASeqPIPE
################################################################################

rm(list=ls())                                        # remove all the objects from the R session

suppressMessages(library(ggplot2))
suppressMessages(library(interacCircos))
suppressMessages(library(pagedown))
suppressMessages(library(stringr))
suppressMessages(library(RCircos))
suppressMessages(library(methylKit))
suppressMessages(library(htmlwidgets))
suppressMessages(library(genomation))
suppressMessages(library(xlsx))
suppressMessages(library(annotate))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(Rsubread))
suppressMessages(library(ggrepel))
suppressMessages(library(optparse))
suppressMessages(library(pandoc))
suppressMessages(library(viridis))
suppressMessages(library(viridisLite))

# options list with associated default value.
option_list <- list( 
  
  make_option(c("-p", "--project"),
              dest="projectName",
              default="REPORT",
              help="Project Name [default: %default]."),
  
  make_option(c("-m", "--metadata"),
              dest="metadata",
              help="path to the metadata file [default: %default]."),
  
  make_option(c("-i", "--input_dir"),
              dest="input_dir",
              help="path to the Methylkit Imput Directory [default: %default]."),
  
  make_option(c("-b", "--bed_file"),
              dest="bed_file",
              help="path to the tree.nwk file [default: %default].")
)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list, 
                       description="Genetate Images using microeco.",
                       epilogue="For comments, bug reports etc... please contact Darshan Dhakan <darshan.dhakan@basesolve.com>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

#Check mandetory inputs 
if ( is.null(opt$metadata) ) {
  stop("--metadata tsv file in QIIME2 format must be provided. See script usage (--help)")
}

if ( is.null(opt$input_dir) ) {
  stop("--path to the input Methylkit Directory needs to be provided must be provided. See script usage (--help)")
}

if ( is.null(opt$bed_file) ) {
  stop("--taxonmy path to silva taxonomy .qza file must be provided. See script usage (--help)")
}

# get options and arguments
workDir <- getwd()
metadata <- opt$metadata  
input_dir <- opt$input_dir  
bed_dir <- opt$bed_file
projectName <- opt$projectName

################################################################################
###                             running script                               ###
################################################################################

MethylSeq_Analysis <- function(Path, metadata, projectName, bed_file)
{
  library(ggplot2)
  library(interacCircos)
  library(pagedown)
  library(stringr)
  library(RCircos)
  library(methylKit)
  library(htmlwidgets)
  library(genomation)
  library(xlsx)
  library(annotate)
  library(org.Hs.eg.db)
  library(Rsubread)
  library(ggrepel)
  
  
  # - titleStyle : style object to use for title
  xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
    rows <-createRow(sheet,rowIndex=rowIndex)
    sheetTitle <-createCell(rows, colIndex=1)
    setCellValue(sheetTitle[[1,1]], title)
    setCellStyle(sheetTitle[[1,1]], titleStyle)
  }
  
  
  FILES <- list.files(path = Path, pattern = ".gz", all.files = TRUE, full.names = TRUE, recursive = FALSE)
  METADATA_TABLE <- read.csv(file = metadata, sep = "\t", header = TRUE, row.names = 1)
  
  ## Create the Directory ##
  dir_p <- function(dir)
    {
      if (dir.exists(dir))
      {
        print ("Directory already exists")
      }
      else {
        dir.create(dir, recursive = TRUE)
      }
    } 
  
  dir_p(projectName)
  
  ## Perform the Groupwise Comparison
  
  GP1 <- unique(METADATA_TABLE$Group)
  print (GP1)
  for (i in seq_along(GP1[-length(GP1)]))
  {
    print (GP1[i])
    GP2 <- GP1[-(1:i)]
    for (j in seq_along(GP2))
    {
      print (GP2[j])
      if (GP1[i] == GP2[j])
      {
        next
      }
      else {
        
        ## Now add tor the excel sheet
        wb<-createWorkbook(type="xlsx")
        
        # Title and sub title styles
        TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16,
                                           color="dodgerblue4", isBold=TRUE, underline=1)
        SUB_TITLE_STYLE <- CellStyle(wb) +
          Font(wb,  heightInPoints=14,
               isItalic=TRUE, isBold=FALSE)
        # Styles for the data table row/column names
        TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
        TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
          Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
          Border(color="black", position=c("TOP", "BOTTOM"),
                 pen=c("BORDER_THIN", "BORDER_THICK"))
        
        print (paste0("Group1 is ", GP1[i]))
        print (paste0("Group2 is ", GP2[j]))
        
        print (paste0("### Starting the Analysis for ", GP1[i], " and ", GP2[j]))
        Comp <- paste0(GP1[i], "_vs_", GP2[j])
        DIR <- file.path(projectName, "GROUPWISE", Comp)
        dir_p(DIR)
        plot_types <- c("Methylation_Stats", "Coverage_Stats", "Correlation", "Circos_Plot", "Diff_Meth_Plots", "Gene_Set_Plots")
        for (plot in plot_types){
          DIR <- file.path(projectName,"GROUPWISE",Comp,"figures",plot)
          dir_p(DIR)
        }
        DIR <- file.path(projectName, "GROUPWISE", Comp)
        TABLE_OUTPUT <- file.path(DIR, "tables")
        print (TABLE_OUTPUT)
        dir_p(TABLE_OUTPUT)
        
        TABLE_OUTPUT <- file.path(projectName, "GROUPWISE", Comp, "tables")
        META_TABLE <- METADATA_TABLE[METADATA_TABLE$Group %in% c(GP1[i], GP2[j]),]
        FILES2 <- FILES
        FILES2 <- gsub(x = FILES, pattern = ".*/",replacement = "")
        FILES2 <- gsub(x = FILES2, pattern = "_mbias.tsv_EXTRACTED.*", replacement = "")
        FILES2 <- FILES2[match(rownames(META_TABLE), FILES2)]
        NEW_FILES <- FILES[match(rownames(META_TABLE), FILES2)]
        NEW_FILES <- as.list(NEW_FILES)
        GROUPS <- META_TABLE$Group
        GROUPS <- factor(META_TABLE$Group, labels = c(0, 1))
        GROUPS <- as.numeric(as.character(GROUPS))
        context <- gsub(FILES, pattern = ".*_mbias.tsv_EXTRACTED_", replacement = "")
        context <- unique(gsub(context, pattern = ".methylKit.gz", replacement = ""))
        
        myobj <- methRead(NEW_FILES, sample.id = as.list(FILES2), assembly = "hg38", treatment = GROUPS, context = context, mincov = 5)
        assign("myobj", value = myobj, envir = .GlobalEnv)
        for (m in seq_along(myobj))
        {
          DIR <- file.path(projectName,"GROUPWISE",Comp,"figures","Methylation_Stats")
          png(paste0(DIR,"/",FILES2[m], "_Methylation_Coverage_plot.png"), height = 12, width = 16, units = 'in', res = 300)
          getMethylationStats(myobj[[m]],plot=TRUE,both.strands=FALSE)
          dev.off()
          
          DIR <- file.path(projectName,"GROUPWISE",Comp,"figures","Coverage_Stats")
          png(paste0(DIR,"/",FILES2[m], "_Coverage_plot.png"), height = 12, width = 16, units = 'in', res = 300)
          getCoverageStats(myobj[[m]],plot=TRUE,both.strands=FALSE)
          dev.off()
        }
        ## Merge the Samples for further analysis ##
        ## this will merge based on whether the regions are found even in a single sample also.
        meth.min=unite(myobj,min.per.group=1L)
        assign("meth.min", value = meth.min, envir = .GlobalEnv)
        
        # Get the Correlation between samples:
        DIR <- file.path(projectName,"GROUPWISE",Comp,"figures","Correlation")
        png(paste0(DIR,"/",FILES2[m], "_Correlation_plot.png"), height = 12, width = 16, units = 'in', res = 300)
        getCorrelation(object = meth.min, method = "pearson", plot = TRUE)
        dev.off()
        
        
        ## Get the Methylation levels at different tiles across the entire Chromosome ##
        
        tiles=tileMethylCounts(myobj,win.size=10000000,step.size=10000000)
        TABLE1 <- data.frame(tiles[[1]]@.Data)
        TABLE2 <- data.frame(tiles[[2]]@.Data)
        CIRCOS_TABLE1 <- data.frame(TABLE1[,c(1,2,3,6)],rep("TP53",nrow(TABLE1)), rep("NA", nrow(TABLE1)))
        colnames(CIRCOS_TABLE1) <- c("chr", "start", "end", "value", "Name", "link")
        CIRCOS_TABLE1 <- CIRCOS_TABLE1[CIRCOS_TABLE1$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),]
        CIRCOS_TABLE1 <- CIRCOS_TABLE1[!CIRCOS_TABLE1$chr %in% c("chr7.1", "chr8.1"),]
        CIRCOS_TABLE2 <- data.frame(TABLE2[,c(1,2,3,6)],rep("TP53",nrow(TABLE2)), rep("NA", nrow(TABLE2)))
        colnames(CIRCOS_TABLE2) <- c("chr", "start", "end", "value", "Name", "link")
        CIRCOS_TABLE2 <- CIRCOS_TABLE2[CIRCOS_TABLE2$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),]
        CIRCOS_TABLE2 <- CIRCOS_TABLE2[!CIRCOS_TABLE2$chr %in% c("chr7.1", "chr8.1"),]
        
        PLOT <- Circos(moduleList=CircosHistogram('HISTOGRAM01', data = CIRCOS_TABLE1, fillColor= "#ff7f0e",maxRadius = 200,minRadius = 160,animationDisplay = FALSE) + CircosHistogram('HISTOGRAM02', data = CIRCOS_TABLE2, fillColor= "forestgreen",maxRadius = 160,minRadius = 120,animationDisplay = FALSE), genome=list("chr1"=250000000,"chr2"=250000000,"chr3"=200000000,"chr4"=200000000,"chr5"=190000000,"chr6"=180000000,"chr7"=160000000,"chr8"=150000000,"chr7"=45000000,"chr8"=43000000,"chr9"=140000000, "chr10"=140000000,"chr11"=140000000,"chr12"=140000000,"chr13"=120000000,"chr14"=110000000,"chr15"=110000000,"chr16"=100000000,"chr17"=90000000,"chr18"=90000000,"chr19"=60000000,"chr20"=70000000,"chr21"=50000000,"chr22"=60000000,"chrX"=160000000,"chrY"=60000000), outerRadius = 240, width = 900,height = 750, innerRadius = 200, svgClassName = "PLOT.svg")
        DIR <- file.path(projectName,"GROUPWISE",Comp,"figures","Circos_Plot")
        saveWidget(PLOT, paste0(DIR, "/","Combined_Circos.html"), selfcontained = TRUE)
        chrome_print(input = paste0(DIR, "/", "Combined_Circos.html"), output = paste0(DIR, "/","Combined_Circos.pdf"))
        
        ## Calculate the Differential Methylation 
        
        DIR <- file.path(projectName,"GROUPWISE",Comp,"figures","Diff_Meth_Plots")
        myDiff=calculateDiffMeth(meth.min,num.cores=2)
        myDiff <- myDiff[myDiff$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),]
        png(paste0(DIR, "/", "DMR_DISTRIBUTION.png"), height = 12, width = 18, units = 'in', res = 300)
        diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)
        dev.off()
        
        ## Annotate the Differentially methylated bases ##
        gene.obj=readTranscriptFeatures(bed_file)
        myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
        myDiff25p <- myDiff25p[myDiff25p$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),]
        SUMMARY <- annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
        
        TABLE <- data.frame(SUMMARY@annotation)
        colnames(TABLE) <- "Percentage"
        COLORS <- c("dodgerblue", "indianred", "green", "purple")
        FTS <- c("promoter", "exon", "intron", "intergenic")
        names(COLORS) <- FTS
        png(paste0(DIR,"/", "REGIONWISE_DMR_DISTRIBUTION.png"), height = 12, width = 16, units = 'in', res = 300)
        pie(TABLE$Percentage, labels = paste0(rownames(TABLE), " - ", round(TABLE$Percentage, 2), "%"), col = COLORS[as.character(rownames(TABLE))], main = "Differentially Methylated Regions", cex = 2)
        legend("topright", legend = as.character(rownames(TABLE)), fill = COLORS[as.character(rownames(TABLE))], cex = 2)
        dev.off()
        
        ## Calculate the Distribution of Methylation near Promoters ##
        
        ## Add the Table:
        
        ## Get the Association of Differential Methylated Regions with TSS ##
        
        diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
        TSS_TABLE <- getAssociationWithTSS(diffAnn)
        Gene_names <- TSS_TABLE$feature.name
        Gene_names <- gsub(x = Gene_names, pattern = "\\..*", replacement = "")
        Genes <- select(org.Hs.eg.db, keys = Gene_names,columns = c('SYMBOL', 'ENTREZID'), keytype = "REFSEQ")
        TSS_TABLE <- cbind(TSS_TABLE, Genes)
        sheet <- createSheet(wb, sheetName = paste0("DMR_TSS_DISTANCE_TABLE"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("DISTANCES_TO_TSS"), titleStyle = TITLE_STYLE)
        addDataFrame(TSS_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
        
        library(tidyverse)
        library(dplyr)
        
        ## Plot the Differentially methylated regions across the Chromosome:
        
        DIFF_TABLE <- data.frame(myDiff25p@.Data)
        colnames(DIFF_TABLE) <- c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff")
        DIFF_TABLE <- DIFF_TABLE[DIFF_TABLE$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),]
        DIFF_TABLE$chr <- factor(DIFF_TABLE$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
        DIFF_TABLE$Methylation <- ifelse(DIFF_TABLE$meth.diff > 0, "HYPER", "HYPO")
        DIFF_TABLE <- DIFF_TABLE %>% group_by(chr) %>% summarise(chr_len=max(start)) %>% mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>% select(-chr_len) %>% left_join(DIFF_TABLE, ., by=c("chr"="chr")) %>% arrange(chr, start) %>% mutate( BPcum=start+tot)
        axisdf = DIFF_TABLE %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        PLOT <- ggplot(DIFF_TABLE, aes(x=BPcum, y=meth.diff, col = chr)) + scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + geom_point(stat = "identity") + scale_color_manual(values = rep(c("orange", "dodgerblue"), 22)) + xlab(label = "Chromosomes")  + ylab("Difference in Methylation (%)") + theme(axis.text.x = element_text(size = 14, face = "bold", angle = 60, hjust = 1), axis.text.y = element_text(size = 14, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.position = "none")
        png(filename = paste0(DIR, "/", "DMR_25P_CHROMOSOMEWISE_PLOT.png"), height = 12, width = 16, units = 'in', res = 300)
        plot(PLOT)
        dev.off()
        
        
        DIFF_TABLE$Gene <- TSS_TABLE$SYMBOL
        Top_Genes <- DIFF_TABLE[order(DIFF_TABLE$qvalue)[1:30],"BPcum"]
        DIFF_TABLE$is_annotate <- ifelse(DIFF_TABLE$BPcum %in% Top_Genes, "yes", "no")
        PLOT <- ggplot(DIFF_TABLE, aes(x=BPcum, y=meth.diff, col = chr)) + scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + geom_point(stat = "identity") + scale_color_manual(values = rep(c("orange", "dodgerblue"), 22)) + geom_point(data = subset(DIFF_TABLE, is_annotate=="yes"), color = "purple", size = 4) +  geom_label_repel(data = subset(DIFF_TABLE, is_annotate=="yes"), aes(label = Gene), size = 6) + xlab(label = "Chromosomes") + ylab(label = "Difference in Methylation (%)") + theme(axis.text.x = element_text(size = 14, face = "bold", angle = 60, hjust = 1), axis.text.y = element_text(size = 14, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.position = "none")
        png(filename = paste0(DIR, "/", "DMR_25P_CHROMOSOMEWISE_ANNOTATED.png"), height = 12, width = 16, units = 'in', res = 300)
        plot(PLOT)
        dev.off()
        
        detach("package:tidyverse")
        detach("package:dplyr")
        detach("package:tidyr")
        
        ## Regional Analysis of Differential Methylated Regions: ##
        
        ## Perform Differential Methylation Analysis for Promoter Regions
        
        promoters_count=regionCounts(myobj,gene.obj$promoters)
        promoters_count <- unite(promoters_count, min.per.group = 1L)
        myDiff_Promoters <- calculateDiffMeth(promoters_count, num.cores = 2)
        # Filter the Differential Methylated Regions
        myDiff25p_Promoters=getMethylDiff(myDiff_Promoters,difference=25,qvalue=0.05)
        # Get the DMR data
        PROMOTER_DMR_TABLE <- data.frame(myDiff25p_Promoters@.Data)
        colnames(PROMOTER_DMR_TABLE) <- c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff")
        ## Now Annotate the Promoter regions to find the Genes associated:
        PROMOTERS_ANNOTATED <- annotateWithGeneParts(as(myDiff25p_Promoters,"GRanges"),gene.obj)
        Gene_names <- PROMOTERS_ANNOTATED@dist.to.TSS$feature.name
        Gene_names <- gsub(x = Gene_names, pattern = "\\..*", replacement = "")
        Genes <- select(org.Hs.eg.db, keys = Gene_names,columns = c('SYMBOL', 'ENTREZID'), keytype = "REFSEQ")
        PROMOTER_DMR_TABLE <- cbind(PROMOTER_DMR_TABLE, Genes)
        sheet <- createSheet(wb, sheetName = paste0("DMR_PROMOTERS"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0("DMR_PROMOTERS"), titleStyle = TITLE_STYLE)
        addDataFrame(PROMOTER_DMR_TABLE, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)

        ## Plot the Genes with their Promoters as Bar plot ##
        
        PROMOTER_DMR_TABLE <- PROMOTER_DMR_TABLE[!is.na(PROMOTER_DMR_TABLE$SYMBOL),]
        PROMOTER_DMR_TABLE <- PROMOTER_DMR_TABLE[rev(order(abs(PROMOTER_DMR_TABLE$meth.diff))), ]
        PROMOTER_DMR_TABLE <- PROMOTER_DMR_TABLE[1:20,]
        PROMOTER_DMR_TABLE$Methylation <- ifelse(PROMOTER_DMR_TABLE$meth.diff > 0, "hyper", "hypo")
        
        PLOT <- ggplot(PROMOTER_DMR_TABLE, aes(x=SYMBOL, y = meth.diff, fill = Methylation)) + geom_bar(stat = "identity") + coord_flip() + ylab(label = "Differential Methylation (%)") + xlab(label = "Genes") + scale_fill_viridis_d() + theme(axis.text = element_text(size = 14, face = "bold"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 14, face = "bold"))
        png(paste0(DIR, "/", "DMR_PROMOTER_GENES_PLOT.png"), height = 15, width = 15, units = 'in', res = 300)
        plot(PLOT)
        dev.off()        
        
        ## Generate the TSS Distance Heatmap ##
        
        promoter <- promoterRegions(annotation = "hg38",upstream = 3000L,downstream = 3000L)
        promoter <- as(promoter, "GRanges")
        
        Sample1_Obj <- as(myobj[[1]], "GRanges")
        Sample2_Obj <- as(myobj[[2]], "GRanges")
        targets <- list(Control = Sample1_Obj, Treated = Sample2_Obj)
        sml = ScoreMatrixList(targets = targets, windows = promoter, bin.num = 50, strand.aware = TRUE)
        png(paste0(DIR, "/", "DISTRIBUTION_METHYLATION_NEAR_PROMOTERS.png"), height = 12, width = 16, units = 'in', res = 300)
        multiHeatMatrix(sml, xcoords = c(-10000, 10000), clustfun = function(x) kmeans(x, centers = 2)$cluster, grid = TRUE, col = viridis(5,direction = -1,option = "A"), group.col = c("skyblue", "indianred"))
        dev.off()
        
        ## Generate the GO and KEGG Enrichment Plots
        
        library(clusterProfiler)
        entrezids <- TSS_TABLE$ENTREZID %>% as.character() %>% unique()
        ego <- enrichGO(gene = entrezids, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, readable = TRUE)
        cluster_summary <- data.frame(ego)
        DIR <- file.path(projectName, "GROUPWISE", Comp, "figures","Gene_Set_Plots")
        wb2<-createWorkbook(type="xlsx")
        TITLE_STYLE2 <- CellStyle(wb2)+ Font(wb2,  heightInPoints=16,
                                             color="dodgerblue4", isBold=TRUE, underline=1)
        SUB_TITLE_STYLE2 <- CellStyle(wb2) +
          Font(wb2,  heightInPoints=14,
               isItalic=TRUE, isBold=FALSE)
        # Styles for the data table row/column names
        TABLE_ROWNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE)
        TABLE_COLNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE) +
          Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
          Border(color="black", position=c("TOP", "BOTTOM"),
                 pen=c("BORDER_THIN", "BORDER_THICK"))
        if (nrow(cluster_summary) > 0)
        {
          png(paste0(DIR, "/", "GENE_ONTOLOGY_PLOTS.png"), height = 16, width = 16, units = 'in', res = 300)
          plot(dotplot(ego, showCategory=30))
          dev.off()
          sheet <- createSheet(wb2, sheetName = paste0("GO_ENRICHMENT_TABLE"))
          xlsx.addTitle(sheet, rowIndex=1, title=paste0("GO Enrichment of the DMR regions Identified"), titleStyle = TITLE_STYLE2)
          addDataFrame(cluster_summary, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
        }
        # Plot the KEGG Enrichment:
        ekegg <- enrichKEGG(gene = entrezids, organism = 'hsa', pvalueCutoff = 1, keyType="kegg")
        if (nrow(ekegg@result) > 0)
        {
          png(paste0(DIR, "/", "KEGG_ENRICHMENT_PLOT.png"), height = 12, width = 14, units = 'in', res = 300)
          plot(dotplot(ekegg))
          dev.off()
          sheet <- createSheet(wb2, sheetName = paste0("KEGG_ENRICHMENT"))
          xlsx.addTitle(sheet, rowIndex=1, title=paste0("KEGG Enrichment of the DMR Regions Identified"), titleStyle = TITLE_STYLE2)
          addDataFrame(ekegg@result, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
        }
        detach("package:clusterProfiler", unload = TRUE)
        saveWorkbook(wb2, paste0(TABLE_OUTPUT, "/",  "FUNCTIONAL_ENRICHMENT_TABLE.xlsx"))
        saveWorkbook(wb, paste0(TABLE_OUTPUT, "/",  "DMR_TABLES.xlsx"))
      }
    }
  }
}

MethylSeq_Analysis(Path = input_dir, metadata = metadata, projectName = projectName, bed_file = bed_dir)