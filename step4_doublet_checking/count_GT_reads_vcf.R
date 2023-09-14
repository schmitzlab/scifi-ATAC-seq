## counts SNPs ##
args <- commandArgs(trailingOnly=T)
if(length(args) != 2){stop("Rscript count_GT_reads.R [prefix] [vcf]")}

# load libraries
library(Matrix)

# functions
refVCF <- function(x){
    
    # split genotypes
    g0 <- as.data.frame(do.call(rbind, strsplit(as.character(x$V10),":")))
    g1 <- as.data.frame(do.call(rbind, strsplit(as.character(x$V11),":")))
    #df <- data.frame(g0=g0$V1,g1=g1$V1,row.names=paste(x$V1,x$V2,sep="_"))
    df <- data.frame(g0=g0$V1,g1=g1$V1)
    return(df)
    
}

# prefix
prefix <- as.character(args[1])
vcf <- as.character(args[2])
setwd(paste0("./",prefix))


# load ids
bcs <- read.table(paste0(prefix,".BARCODES.txt"))
#snps <- read.table("common_variants_covered.vcf")
snps <- read.table(vcf)

# reformat vcf
snps <- refVCF(snps)

# load data
ref <- readMM('ref.mtx')
alt <- readMM('alt.mtx')

# annotate
rownames(ref) <- rownames(snps)
rownames(alt) <- rownames(snps)
colnames(ref) <- as.character(bcs$V1)
colnames(alt) <- as.character(bcs$V1)

# count GT reads
snps.g0.ref <- snps[snps$g0=="0/0" | snps$g0=="0/1",]
snps.g0.alt <- snps[snps$g0=="1/1" | snps$g0=="0/1",]
snps.g1.ref <- snps[snps$g1=="0/0" | snps$g1=="0/1",]
snps.g1.alt <- snps[snps$g1=="1/1" | snps$g1=="0/1",]

# read counts
B73.ref.cnts <- Matrix::colSums(ref[rownames(snps.g0.ref),])
B73.alt.cnts <- Matrix::colSums(alt[rownames(snps.g0.alt),])
Mo17.ref.cnts <- Matrix::colSums(ref[rownames(snps.g1.ref),])
Mo17.alt.cnts <- Matrix::colSums(alt[rownames(snps.g1.alt),])
B73.cnts <- B73.ref.cnts + B73.alt.cnts
Mo17.cnts <- Mo17.ref.cnts + Mo17.alt.cnts

# make output
df <- data.frame(B73_reads=B73.cnts, Mo17_reads=Mo17.cnts, row.names=names(B73.cnts))
write.table(df, file=paste0(prefix,".genotype_read_counts.txt"), quote=F, row.names=T, col.names=T, sep="\t")
