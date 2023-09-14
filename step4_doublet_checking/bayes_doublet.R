###################################################################################################
## description
###################################################################################################
#module load R/4.1.0-foss-2019b
# load arguments
library("dplyr")
library("vioplot")
library("scales")


args <- commandArgs(trailingOnly=T)
if(length(args) < 2){stop("Rscript bayes_doublet-xz.R [reads_count] [name] [p.call]")}
#args
counts <- as.character(args[1])
name <- as.character(args[2])
p.call <- as.numeric(args[3])


#counts <- "B73Mo17_min1k_min1_max2500000.genotype_read_counts.txt"
#geno <- "clusters.clean.tsv"
#name <- "scifiM2_merge_clean"
#name <- "B73Mo17_min1k_min1_max2500000"

counts.m <- read.table(counts)
#counts.m$barcode <- rownames(counts.m)

#------assign cell status with bayes method based on Alex's method in Cell paper-----------
bayes <- function(data, p0 = 0.475, p1 = 0.475, pn = 0.05, p.call = 0.9) {
        data <- as.matrix(data)
        res <- apply(data, 1, function(x) {
                hap0 <- x[1]
                hap1 <- x[2]
                total <- hap0 + hap1
                n_call <- 0.99
                p0_bc.loge <- dbinom(hap0, size = total, prob = p.call, log = T)
                p1_bc.loge <- dbinom(hap1, size = total, prob = p.call, log = T)
                p2_bc.loge <- dbinom(hap1, size = total, prob = 0.5, log = T)
                
                p0_bc <- dbinom(hap0, size = total, prob = p.call)
                p1_bc <- dbinom(hap1, size = total, prob = p.call)
                p2_bc <- dbinom(hap1, size = total, prob = 0.5)
                
                p0_top <- p0_bc * p0
                p1_top <- p1_bc * p1
                p2_top <- p2_bc * pn
                total_prob <- p0_top + p1_top + p2_top
                prob_0 <- p0_top / total_prob
                prob_1 <- p1_top / total_prob
                prob_2 <- p2_top / total_prob
                return(list(B73.binom.loge = p0_bc.loge, Mo17.binom.loge = p1_bc.loge, Doublet.binom.loge = p2_bc.loge, B73.p = prob_0, Mo17.p = prob_1, Doublet.p = prob_2))
        })
        res <- as.data.frame(do.call(rbind,res))
        return(res)
}

res <- bayes(counts.m, p.call = p.call)
#res <- bayes(counts.m)

#---assign doublet type based on binamial test---
res$binom.type <- apply(res[,c(1:3)],1,function(x){names(x)[which.max(as.matrix(x))]})
res$binom.type <- gsub(".binom.loge","",res$binom.type)

#---assign doublet type based on bayes test---
res$bayes.type <- apply(res[,c(4:6)],1,function(x){as.character(names(x)[which.max(as.matrix(x))])})
res$bayes.type <- apply(res[,c(7:8)],1,function(x){
        x <- as.array(x)
        type <- ifelse(length(x[[2]]) == 0, x[[1]], x[[2]])
        return(type)
})
res$bayes.type <- gsub(".p","",res$bayes.type)
res <- apply(res,2,function(x){unlist(x, use.names = F)})

#merge
res <- cbind(counts.m,res)
write.table(res, paste0(name,"_pcall",p.call,"_bayes.doulblet.txt"),sep = "\t", quote = F)

#---------------Plot-------
pdf(paste0(name, "_pcall",p.call,".doubletPlot.bayes.pdf"), width=16, height=4)
layout(matrix(c(1:4), nrow=1))
#---plot violin plot for B73, Mo17 and doublet from binormial prediction---
res$total <- res$B73_reads + res$Mo17_reads + 1
ave.b73 <- round(mean(res[res$binom.type == "B73",]$total), digits = 2)
ave.mo17 <- round(mean(res[res$binom.type == "Mo17",]$total), digits = 2)
ave.dou <- round(mean(res[res$binom.type == "Doublet",]$total), digits = 2)

vioplot(log10(res$total) ~ res$binom.type ,col = c("#1874CD","#BFBFBF","#CD2626"),ylab = "Reads number(log10)",main = "Reads number distribution")

legend("topright",
       legend = c("B73", "Doublet", "Mo17"),
       pch = 15,
       col = c("#1874CD","#BFBFBF","#CD2626"),
       )
mtext(paste0("B73 mean=",ave.b73,";Mo17 mean=",ave.mo17,";Doulet mean=",ave.dou), side=3, cex = 0.5)

#---plot doublet dot plot with binormal result----
res$idxH.rate <- pmin(res$B73_reads, res$Mo17_reads)/(res$total)
b73.n <- nrow(res[res$binom.type == "B73",])
b73.r <- round(b73.n/nrow(res), digits = 2)
mo17.n <- nrow(res[res$binom.type == "Mo17",])
mo17.r <- round(mo17.n/nrow(res),digits = 2)
doublet.n <- nrow(res[res$binom.type == "Doublet",])
doublet.r <- round(doublet.n/nrow(res),digits = 2)
idxH.soup.mean <- round(mean(res[res$binom.type == "B73" | res$binom.type == "Mo17",]$idxH.rate),digits = 2)
idxH.soup.median <- round(median(res[res$binom.type == "B73" | res$binom.type == "Mo17",]$idxH.rate), digits = 2)
#unassign.n <- nrow(res[res$binom.type == "NA",])
#unassign.r <- round(unassign.n/nrow(res), digits = 2)

res1 <- res[res$binom.type != "NA",]
colors <- c("#1874CD","#BFBFBF","#CD2626")
myColors <- ifelse(res1$binom.type == "B73","#1874CD",ifelse(res1$binom.type == "Doublet", "#BFBFBF","#CD2626"))
#head(myColors)
plot(x = res1$B73_reads/1000, 
     y = res1$Mo17_reads/1000,
     #col = colors[factor(res1$binom.type)],
     col=alpha(myColors, 0.3),
     cex = 0.8,
     pch=16,
     xlim = c(0,20),
     ylim = c(0,20),
     xlab = "B73 reads number(x1000)",
     ylab = "Mo17 reads number(x1000)",
     main = "Doublet plot for binomial test"
     )
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)
legend("topright",
       legend = c("B73", "Doublet", "Mo17"),
       pch = 16,
       col = colors)
mtext(paste0("B73: ",b73.n,", ",b73.r,";Mo17: ",mo17.n,", ",mo17.r,";Dou: ",doublet.n,", ",doublet.r,";\nindexH mean: ",idxH.soup.mean,", median:",idxH.soup.median), side=3, cex = 0.5)

#---plot violin plot for B73, Mo17 and doublet from bayes prediction---
res$total <- res$B73_reads + res$Mo17_reads +1
ave.b73 <- round(mean(res[res$bayes.type == "B73",]$total), digits = 2)
ave.mo17 <- round(mean(res[res$bayes.type == "Mo17",]$total), digits = 2)
ave.dou <- round(mean(res[res$bayes.type == "Doublet",]$total), digits = 2)

vioplot(log10(res$total) ~ res$bayes.type ,col = c("#1874CD","#BFBFBF","#CD2626"),ylab = "Reads number(log10)",main = "Reads number distribution")

legend("topright",
       legend = c("B73", "Doublet", "Mo17"),
       pch = 15,
       col = c("#1874CD","#BFBFBF","#CD2626"),
)
mtext(paste0("B73 mean=",ave.b73,";Mo17 mean=",ave.mo17,";Doulet mean=",ave.dou), side=3, cex = 0.5)

#---plot doublet dot plot with binormal result----
res$idxH.rate <- pmin(res$B73_reads, res$Mo17_reads)/(res$total)
b73.n <- nrow(res[res$bayes.type == "B73",])
b73.r <- round(b73.n/nrow(res), digits = 2)
mo17.n <- nrow(res[res$bayes.type == "Mo17",])
mo17.r <- round(mo17.n/nrow(res),digits = 2)
doublet.n <- nrow(res[res$bayes.type == "Doublet",])
doublet.r <- round(doublet.n/nrow(res),digits = 2)
idxH.soup.mean <- round(mean(res[res$bayes.type == "B73" | res$bayes.type == "Mo17",]$idxH.rate),digits = 2)
idxH.soup.median <- round(median(res[res$bayes.type == "B73" | res$bayes.type == "Mo17",]$idxH.rate), digits = 2)
#unassign.n <- nrow(res[res$bayes.type == "NA",])
#unassign.r <- round(unassign.n/nrow(res), digits = 2)

res1 <- res[res$bayes.type != "NA",]
colors <- c("#1874CD","#BFBFBF","#CD2626")
myColors <- ifelse(res1$bayes.type == "B73","#1874CD",ifelse(res1$bayes.type == "Doublet", "#BFBFBF","#CD2626"))
#head(myColors)
plot(x = res1$B73_reads/1000, 
     y = res1$Mo17_reads/1000,
     #col = colors[factor(res1$bayes.type)],
     col=alpha(myColors, 0.3),
     cex = 0.8,
     pch=16,
     xlim = c(0,20),
     ylim = c(0,20),
     xlab = "B73 reads number(x1000)",
     ylab = "Mo17 reads number(x1000)",
     main = "Doublet plot for bayes test"
)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 1)
legend("topright",
       legend = c("B73", "Doublet", "Mo17"),
       pch = 16,
       col = colors)
mtext(paste0("B73: ",b73.n,", ",b73.r,";Mo17: ",mo17.n,", ",mo17.r,";Dou: ",doublet.n,", ",doublet.r,";\nindexH mean: ",idxH.soup.mean,", median:",idxH.soup.median), side=3, cex = 0.5)

dev.off()

