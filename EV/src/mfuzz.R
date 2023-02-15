library(tidyverse)
################################################
## intersection of consecutive DE genes for 
################################################

args<-commandArgs(TRUE)

clusterN=as.integer(args[1])

# clusterN = 3

output_dir = paste0("../results/mufuzz_consecutive/cluster")
# Create folder if the directory doesn't exist
file.exists(output_dir) || dir.create(output_dir, recursive = T)

#### for EV
RNA_source = 'EV'
Timepoints_sorted_levels = c(0,seq(1,11,2),seq(12,30,2))
consecutiveDEgenes.EV=data.frame()

for(i in 2:length(Timepoints_sorted_levels)){
  time_i = Timepoints_sorted_levels[i];
  time_i_m_1 = Timepoints_sorted_levels[i-1]
  message(paste("# Running pair-wise comparison between consecutive Day",time_i,"vs",time_i_m_1, "..."))
  res= read.table(file=file.path("/data/bioinformatics/projects/fei2021/results/DE", paste0("DEresult.", RNA_source, ".consecutive_Day",time_i,"vs",time_i_m_1,".padj05.xls")), 
                  sep="\t", header = T)
  consecutiveDEgenes.EV = rbind(mutate(res, comparison = paste0("D",time_i,"vs",time_i_m_1)), consecutiveDEgenes.EV)
}

dim(consecutiveDEgenes.EV); 
# head(consecutiveDEgenes.EV)
consecutiveDEgenes.EV = consecutiveDEgenes.EV %>% mutate(value=1, geneID=X) %>% pivot_wider(id_cols=geneID, names_from=comparison, values_from=value, values_fill = 0)

#### for cell
RNA_source = 'cell'
Timepoints_sorted_levels = c(0,11,30)
consecutiveDEgenes.cell=data.frame()
for(i in 2:length(Timepoints_sorted_levels)){
  time_i = Timepoints_sorted_levels[i];
  time_i_m_1 = Timepoints_sorted_levels[i-1]
  message(paste("# Running pair-wise comparison between consecutive Day",time_i,"vs",time_i_m_1, "..."))
  res= read.table(file=file.path("/data/bioinformatics/projects/fei2021/results/DE", paste0("DEresult.", RNA_source, ".consecutive_Day",time_i,"vs",time_i_m_1,".padj05.xls")), 
                  sep="\t", header = T)
  consecutiveDEgenes.cell = rbind(mutate(res, comparison = paste0("D",time_i,"vs",time_i_m_1)), consecutiveDEgenes.cell)
}

dim(consecutiveDEgenes.cell); 
# head(consecutiveDEgenes.cell)
consecutiveDEgenes.cell = consecutiveDEgenes.cell %>% mutate(value=1, geneID=X) %>% pivot_wider(id_cols=geneID, names_from=comparison, values_from=value, values_fill = 0) 


################################################
## read fpkm for clustering
################################################
# read data
fpkm.all = read.table("/data/bioinformatics/projects/fei2021/results/merged/genes.fpkm.cufflinks.allSamples.xls", header = T, row.names = 1)
# head(fpkm.all)

# group rep into mean
fpkm.long = rownames_to_column(fpkm.all) %>% pivot_longer(cols = contains("_Rep"), names_to = "sampleName", values_to = "fpkm") %>% 
  separate(sampleName, c("RNAsource",NA,"sampleID","Day","Rep"), sep = "_", remove = F) %>% 
  mutate(Day = factor(Day, levels = paste0("D",c(0,seq(1,11,2),seq(12,30,2))))) %>% 
  group_by(rowname, RNAsource, Day) %>%
  summarise(fpkm=mean(fpkm)) 

fpkm.cell = filter(fpkm.long, RNAsource=="cell") %>% 
  pivot_wider(id_cols = 'rowname', names_from='Day', values_from='fpkm') %>% column_to_rownames()

fpkm.EV = filter(fpkm.long, RNAsource=="EV") %>% 
  pivot_wider(id_cols = 'rowname', names_from='Day', values_from='fpkm') %>% column_to_rownames()

# fpkm = fpkm.cell; consecutiveDEgenes=consecutiveDEgenes.cell
fpkm = fpkm.EV; consecutiveDEgenes=consecutiveDEgenes.EV

genes_annotation = read.table("/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed+2", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneSymbol","geneType"));
gene_names = genes_annotation$geneSymbol[genes_annotation$geneType=="protein_coding"]
source("lib.R")
################################################
## clustering time-series data with Mfuzz
## Ref: https://mp.weixin.qq.com/s/ubfTKFyUeBR585EfR9nfig
################################################
pdf(file.path(output_dir, paste0("MFuzz_",clusterN,".pdf")), width = 8,height = 8)
library(Mfuzz) # BiocManager::install("Mfuzz")
# only the consecutive DE genes
expr = fpkm[consecutiveDEgenes$geneID,]
expr = log10(expr + 0.01)
eset <- new("ExpressionSet",exprs = as.matrix(expr))
dim(eset)
# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset, min.std=0)
## 2. 标准化：聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
eset <- standardise(eset)

## 3. 聚类：Mfuzz中的聚类算法需要提供两个参数，第一个参数为希望最终得到的聚类的个数，这个参数由我们直接指定；
# 第二个参数称之为fuzzifier值，用小写字母m表示，可以通过函数评估一个最佳取值

cat("#Cluster N = ", clusterN,"\n")
c <- clusterN # 聚类个数
m <- mestimate(eset) #  评估出最佳的m值
cl <- mfuzz(eset, c = c, m = m) # 聚类

# cl$size # gene size of cluster
# cl$membership # membership value

file_name = paste0("mfuzz_clusterN_",clusterN,".txt")
write.table(cl$membership, paste0(output_dir,"/",file_name),sep="\t")

# num_ls = c()
# sum = 0
# for(j in 1:clusterN){
#   cat("#Cluster N = ", clusterN, " -->>ORA:", j,"\n")
#
#   # cl$cluster[cl$cluster == j] # member of cluster 1
#
#   ## only keep genes with cluster value > 0.5
#   gene_id_ls = rownames(cl$membership[(cl$membership[,i] > 0.5),])
#
#   # gene_symbol_ls = genes_annotation$geneSymbol[match(names(cl$cluster[cl$cluster == j]), genes_annotation$geneID)]
#   gene_symbol_ls = genes_annotation$geneSymbol[match(gene_id_ls, genes_annotation$geneID)]
#   res_ORA = ORA(gene_symbol_ls, allGenes=gene_names, topN=20, nCutoff=3, pCutoff=0.05,output=paste0(output_dir,"/mfuzz_cluster_",j))
#   sum = sum + nrow(res_ORA)
#   num_ls = c(num_ls,nrow(res_ORA));
# }
# cat("Sum:", sum, "\n")
# num_ls = c(num_ls, sum)
# write.table(num_ls, paste0(output_dir,"/go_number.txt"),sep="\t")

## 4. visualization
library(RColorBrewer)
mfuzz.plot2(eset,cl,min.mem=0, time.labels=gsub("D","",colnames(fpkm)),
            xlab = "Differentiation time (Day)", ylab = "Normalized expression values (log10 fpkm)",
            colo = "fancy",centre=T, x11=F, centre.lwd=3)
# mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
dev.off()




