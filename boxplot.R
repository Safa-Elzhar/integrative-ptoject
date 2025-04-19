install.packages("BiocManager")
BiocManager::install("DESeq2", force = TRUE)

library(BiocManager)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(NMF)


data = as.matrix(read.table("exp", row.names=1, header=T))
View(data)
pheno  <- read.table("clinical", sep = "\t",row.names = 1)
View(pheno)
dim(data)
dim(pheno)
rownames(pheno) <- gsub("-", ".", rownames(pheno))

# Check dimensions to confirm the sizes
ncol(data)
nrow(pheno)  

Filter the 'pheno' data frame to match the number of columns in 'data'
pheno_filtered <- pheno[1:ncol(data), ]
nrow(pheno_filtered)

# Ensure row names of pheno and column names of data are set correctly

# Identify common identifiers

common_ids <- intersect(rownames(pheno), colnames(data))

# Filter pheno and data to include only common identifiers

pheno_filtered <- pheno[common_ids, ]
data_filtered <- data[, common_ids]

# Check dimensions to ensure they match

nrow(pheno_filtered)  
ncol(data_filtered)   



pheno2=pheno_filtered[colnames(data),]
View(pheno2)



Remove rows with NA in the design variable (V22)

pheno_filtered <- pheno_filtered[!is.na(pheno_filtered$V22), ]


dds <- DESeqDataSetFromMatrix(countData = data_filtered, colData = pheno_filtered, design = ~ V22)









common_ids <- intersect(colnames(data), rownames(pheno2))
if(length(common_ids) == 0) {
  stop("No matching identifiers found between colnames(data) and rownames(pheno2).")
}



pheno_filtered <- pheno2[common_ids, , drop=FALSE]
data_filtered <- data[, common_ids, drop=FALSE]

data_filtered <- data_filtered[, rownames(pheno_filtered), drop=FALSE]

if (!all(rownames(pheno_filtered) == colnames(data_filtered))) {
  stop("Row names of pheno_filtered do not match column names of data_filtered.")
}
pheno_filtered <- pheno_filtered[!is.na(pheno_filtered$V22), ]
dds <- DESeqDataSetFromMatrix(countData = data_filtered, colData = pheno_filtered, design = ~ V29)




















# Check for common identifiers

common_ids <- intersect(colnames(data_filtered), rownames(pheno_filtered))

# Filter pheno_filtered and data_filtered to include only common identifiers

pheno_filtered <- pheno_filtered[common_ids, , drop=FALSE]
data_filtered <- data_filtered[, common_ids, drop=FALSE]

# Check for NAs in pheno_filtered$V29 and remove corresponding rows in pheno_filtered and columns in data_filtered

na_indices <- which(is.na(pheno_filtered$V29))

if (length(na_indices) > 0) {
  pheno_filtered <- pheno_filtered[-na_indices, , drop=FALSE]
  data_filtered <- data_filtered[, -na_indices, drop=FALSE]
}

# Check dimensions after filtering

cat("Dimensions of pheno_filtered:", dim(pheno_filtered), "\n")
cat("Dimensions of data_filtered:", dim(data_filtered), "\n")

# Verify that row names of pheno_filtered match column names of data_filtered

if (!all(rownames(pheno_filtered) == colnames(data_filtered))) {
  stop("Row names of pheno_filtered do not match column names of data_filtered.")
}

# Create DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = data_filtered, colData = pheno_filtered, design = ~ V29)




# Ensure 'V29' in pheno_filtered is a factor

pheno_filtered$V29 <- as.factor(pheno_filtered$V29)

# Create DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = data_filtered, colData = pheno_filtered, design = ~ V29)

View(dds)
dds.run = DESeq(dds)



cond1="DECEASED" 
cond2="LIVING"
res=results(dds.run, contrast = c("V29",cond1 ,cond2))
res=as.data.frame(res[complete.cases(res), ])
View(res)
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>1.2,]
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="LIVING  vs DECEASED DEGs"))

with(res, plot(log2FoldChange, -log10(padj), pch=20, main="LIVING  vs DECEASED DEGs"))
with(subset(res, padj<.05 & (log2FoldChange)>1.1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(res, padj<.05 & (log2FoldChange)< -1.1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

aheatmap(log2(exp.degs+1), annCol =pheno_filtered$V29, col = rev(brewer.pal(9,"RdBu")), main="mRNA Control vs infection")



aheatmap(log2(exp.degs+1),col = rev(brewer.pal(9,"RdBu")), main="mRNA LIVING  vs DECEASED DEGs", annCol = pheno_filtered$V29)

