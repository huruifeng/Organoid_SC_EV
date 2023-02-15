###########################################
# R script for running pathway analysis
# Usage: Rscript Enrich.R DEresult.all_file Sample_mate_info_file output_folder
##########################################

suppressPackageStartupMessages(library('tidyverse',logical.return=T)) || {install.packages('tidyverse'); suppressPackageStartupMessages(library('tidyverse',logical.return=T))}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
suppressPackageStartupMessages(library('SPIA',logical.return=T)) || {BiocManager::install('SPIA',ask=F); suppressPackageStartupMessages(library('SPIA',logical.return=T))}
suppressPackageStartupMessages(library('org.Hs.eg.db',logical.return=T)) || {BiocManager::install("org.Hs.eg.db",ask=F); suppressPackageStartupMessages(library('org.Hs.eg.db',logical.return=T))}
suppressPackageStartupMessages(library('biomaRt',logical.return=T)) || {BiocManager::install('biomaRt',ask=F); suppressPackageStartupMessages(library('biomaRt',logical.return=T))}
suppressPackageStartupMessages(library("topGO",logical.return=T)) || {BiocManager::install("topGO",ask=F); suppressPackageStartupMessages(library("topGO",logical.return=T))}
suppressPackageStartupMessages(library('fgsea',logical.return=T)) || {BiocManager::install('fgsea',ask=F); suppressPackageStartupMessages(library('fgsea',logical.return=T))}
suppressPackageStartupMessages(library('pathview',logical.return=T)) || {BiocManager::install('pathview',ask=F); suppressPackageStartupMessages(library('pathview',logical.return=T))}

require('BiocParallel',quietly=T, warn.conflicts=F) || BiocManager::install('BiocParallel');


# args<-commandArgs(TRUE)
# 
# prefix=args[1]
# # prefix="PPMI_CaseAll_CtrlNonCarrier_allgene_DEanalysis"

output_dir = paste0("../results/Enrich/EV_20")

# Create folder if the directory doesn't exist
file.exists(output_dir) || dir.create(output_dir, recursive = T)


input_DE_result=  "/Volumes/bioinformatics/projects/fei2021/results/DE/DEresult.EV.Days.all.xls.gz"

topN=20
Q_CUTOFF=0.05
index="hg19"
filters="medium"
program="all"

if(is.null(program)) program="all"

message(paste("-i", input_DE_result,"-o",output_dir,"-N",topN,"-q",Q_CUTOFF,"-g",index,"-F",filters))

# check input
if(!file.exists(input_DE_result)) {stop(paste(input_DE_result, "doesn't exist. Exit!"), call.=FALSE);}

set.seed(12)
##======================================
message("[INFO] 1. Loading data...")
##=====================================

DE_result=read.delim(input_DE_result, stringsAsFactors=F, row.names = 1, header=T, check.names =F)

if(nrow(DE_result)<=5) {stop(paste(input_DE_result, "has too few (<=5) records. Exit!"), call.=FALSE);}

# get human ortholog if the DE result is from non-human species
DE_result$human_orth_symbol=DE_result$symbol

#head(DE_result)
dim(DE_result); DE_result=DE_result[!is.na(DE_result$human_orth_symbol),]; dim(DE_result)

# to save all genes for GSEA
DE_result0=DE_result;

pc_gene = DE_result0[DE_result0$geneType=="protein_coding",]
DE_result = DE_result0[DE_result0$geneType=="protein_coding",]

# filter
if(filters == 'no') DE_result=DE_result # all i.e. p<=1
if(filters == 'loose') DE_result=subset(DE_result, pvalue<=.05) # pvalue<0.05
if(filters == 'medium') DE_result=subset(DE_result, padj<=.05)
if(filters == 'strigent') DE_result=subset(DE_result, padj<=0.05 & lfcSE<2)

if(nrow(DE_result)<=5) {stop(paste(input_DE_result, "has too few (<=5) records. Exit!"), call.=FALSE);}

source("lib.R")

##==================================================
message("2. running GO terms enrichment anlaysis")
##==================================================
if(("symbol" %in% colnames(DE_result)) && (program=="all" || grepl("GO",program))) {
  ## require: only gene symbol

  df <- read.table("/Volumes/NEUROGEN/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F)
  gene_names = df$V7[df$V8=="protein_coding"]

  # all DE genes
  topGOenrichment(DE_result$human_orth_symbol, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all', output=file.path(output_dir,"all_gene"))
  if("log2FoldChange" %in% colnames(DE_result)){
    # UP regulated DE genes
    topGOenrichment(DE_result$human_orth_symbol[DE_result$log2FoldChange>0],
                    allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='up', output=file.path(output_dir,"all_gene"))
    # DOWN regulated DE genes
    topGOenrichment(DE_result$human_orth_symbol[DE_result$log2FoldChange<0],
                    allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='down', output=file.path(output_dir,"all_gene"))
  }

  ## simply over-representation analysis (ORA)
  ORA(DE_result$human_orth_symbol, allGenes=gene_names, topN=topN, nCutoff=3, pCutoff=0.05, output=file.path(output_dir,"all_gene"))
  if("log2FoldChange" %in% colnames(DE_result)){
    # UP regulated DE genes
    ORA(DE_result$human_orth_symbol[DE_result$log2FoldChange>0 ],
        allGenes=gene_names,topN=topN, output=file.path(output_dir,"all_gene.up"))
    # DOWN regulated DE genes
    ORA(DE_result$human_orth_symbol[DE_result$log2FoldChange<0],
        allGenes=gene_names,topN=topN, output=file.path(output_dir,"all_gene.down"))
  }
}

##==================================================
message("3. GSEA pathway analysis")
##==================================================
if(all(c('symbol', 'log2FoldChange') %in% colnames(DE_result)) && (program=="all" || grepl("GSEA",program))){
  ## require: gene symbol and log2FoldChange
  ## require no cutoff for the input, e.g. use all genes

  # # Extract the symbol and stat (as metric)
  # ranks = DE_result0 %>% dplyr::select(symbol=human_orth_symbol, stat=stat) %>%
  #   group_by(symbol) %>% top_n(n=1, wt=abs(stat)) %>% ungroup()

  # OR use -log10(p-value)*sign(2logFC) as metric  (08312018)
  # see discussion on the metric: https://www.biostars.org/p/159029/
  # ranks <- data.frame(symbol =  DE_result0$human_orth_symbol, stat = -log10(DE_resul0t$pval) * sign(DE_result0$log2FoldChange))
  ranks <- data.frame(symbol =  pc_gene$human_orth_symbol, stat = pc_gene$log2FoldChange)
  ranks=ranks[!is.na(ranks$stat),]

  ranks = setNames(ranks$stat, ranks$symbol)
  ranks <- ranks[!duplicated(names(ranks))]
  # genesTables <- genesTables %>% mutate(rank = rank(log2FoldChange,  ties.method = "random")) 
  
  gt=data.frame();
  pdf(file.path(output_dir, 'fGSEA_LFC.pdf'))
  for(i in c("c2.cgp","c2.cp.kegg", "c2.cp.wikipathways",
             "c3.mir.mirdb", "c3.mir.mir_legacy","c3.tft.gtrd","c3.tft.tft_legacy",
             "c5.go.bp","c5.go.cc","c5.go.mf","c5.hpo","c8.all")){
    message(i);

    gmt.file = paste0("/Volumes/bioinformatics/referenceGenome/Homo_sapiens/msigdb_v7.4_files/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/",i,".v7.4.symbols.gmt")

    pathways = gmtPathways(gmt.file)
    fgseaRes = fgseaMultilevel(pathways, stats=ranks, minSize=10, maxSize=500)
    fgseaRes = fgseaRes[order(fgseaRes$padj),]  # default is BH-adjusted p-value, same as FDR

    # get FDR < 0.05 pathways
    fgseaRes = subset(fgseaRes, padj < Q_CUTOFF)

    if(nrow(fgseaRes)>0) {
      # plot
      for(x in fgseaRes$pathway) {p=plotEnrichment(pathways[[x]], stats=ranks) + labs(title=x, x = "Rank based on log2FC(desc)"); print(p);}

      # flatten the list of leadingEdge
      fgseaRes$leadingEdge = vapply(fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))

      gt=rbind(gt, data.frame(gene_set=i, fgseaRes))
    }
  }
  dev.off()
  write.table(gt, file=file.path(output_dir, 'fGSEA_LFC.xls'), sep="\t", quote=F, row.names=F, col.names=T)

}

# ## require: only gene symbol and log2FoldChange
# if(all(c('symbol', 'log2FoldChange') %in% colnames(DE_result)) && (program=="all" || grepl("SPIA",program))){
#   ###################################################
#   message("#Step4: SPIA analysis")
#   ###################################################
#   ## Prerequisition
#   # download all KEGG pathway and install in SPIA
#   # cd /Library/Frameworks/R.framework/Versions/4.0/Resources/library/SPIA
#   # curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1,a,":"); print "curl http://rest.kegg.jp/get/"a[2]"/kgml -o extdata/keggxml/hsa/"a[2]".xml"}' | bash
#   # curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1,a,":"); print "curl http://rest.kegg.jp/get/"a[2]"/image -o extdata/keggxml/hsa/"a[2]".png"}' | bash
#   #library(SPIA)
#   #makeSPIAdata(kgml.path=system.file("extdata/keggxml/hsa",package="SPIA"),organism="hsa",out.path="./extdata")
# 
#   # add ENTREZ ID using the human_orth_symbol (as Pathway data might be more enriched in human data)
#   # Note: ortholog genes (e.g. PAX6 in human and Pax6 in mouse) have different ENTREZ ID
#   human_orth_Entrez = mapIds(org.Hs.eg.db,
#                              keys = DE_result$human_orth_symbol,
#                              column = "ENTREZID",
#                              keytype = "SYMBOL",
#                              multiVals='first')
#   DE_result$human_orth_Entrez = as.numeric(human_orth_Entrez) # NA for those genes without Entrez ID
#   degGIs <- as.vector(DE_result$log2FoldChange)
#   names(degGIs) <- DE_result$human_orth_Entrez
#   # remove NA
#   degGIs <- degGIs[!is.na(names(degGIs))]
#   # remove duplicate
#   degGIs=degGIs[!duplicated(names(degGIs))]
# 
#   res <- spia(de=degGIs,
#               all=mappedkeys(org.Hs.egGENENAME),
#               organism='hsa',
#               nB=2000, plots=FALSE,
#               beta=NULL,
#               combine="fisher")
# 
#   ## A data frame containing the ranked pathways and various statistics:
#   # pSize is the number of genes on the pathway;
#   # NDE is the number of DE genes per pathway;
#   # tA is the observed total preturbation accumulation in the pathway;
#   # pNDE is the probability to observe at least NDE genes on the pathway using a hypergeometric model;
#   # pPERT is the probability to observe a total accumulation more extreme than tA only by chance;
#   # pG is the p-value obtained by combining pNDE and pPERT;
#   # pGFdr and pGFWER are the False Discovery Rate and respectively Bonferroni adjusted global p-values;
#   # Status gives the direction in which the pathway is perturbed (activated or inhibited).
#   # KEGGLINK gives a web link to the KEGG website that displays the pathway image with the differentially expressed genes highlighted in red.
# 
#   # add color code to the KEGGLINK based on fold-change of each DE gene
#   # Potential issue: some pathways are not in the Pathview web server and it will return an error (see https://support.bioconductor.org/p/90506/)
# 
#   library(pathview) # BiocManager::install("pathview",ask=F);
#   sapply(filter(res, pGFdr<Q_CUTOFF, !(ID %in% c("04723","04320","04215","05206")), file.exists(paste0("~/Downloads/hsa/hsa",ID,".xml")),
#                 !file.exists(paste0("hsa",ID,".","AMPPD",".png"))) %>% dplyr::select(ID),
#          function(pid) pathview(gene.data  = degGIs,
#                                 pathway.id = pid,
#                                 species = "hsa",
#                                 kegg.dir = "~/Downloads/hsa/",
#                                 out.suffix = "AMPPD",
#                                 limit  = list(gene=max(abs(degGIs)), cpd=1))
#          )
# 
#   # make text file and plot
#   write.table(subset(res, pGFdr<Q_CUTOFF), file.path(output_dir, 'SPIA.xls'), sep="\t", quote=F, row.names=F, col.names=T)
# 
#   pdf(file.path(output_dir, 'SPIA.pdf'))
#   plotP(res, threshold=Q_CUTOFF)
#   res$ID=res$Name;
#   plotP(res, threshold=Q_CUTOFF)
#   dev.off()
# }


