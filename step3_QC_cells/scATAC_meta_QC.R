###############################################################################
## analyze barcode statistics
###############################################################################
#-----------------------History--------------
#20230702 do not use knee as cutoff, using min tn5=1000.
library(MASS)
library(viridis)

# args
args <- commandArgs(T)
if(length(args)!=2){stop("Rscript scATAC_meta_QC-xz.R [input] [name]")}

# load args
input <- as.character(args[1])
#input <- "NCS3_Gm_Middle_maturation_stage_seeds-rep2.updated_metadata.txt"
name <- as.character(args[2])

# setwd
#setwd(getwd())

# functions
plotDist   <- function(x, main=""){
    x <- x[order(x$unique, decreasing=T),]
    rank <- log10(seq(1:nrow(x)))
    depth <- log10(x$unique+1)
    df <- data.frame(rank=rank, depth=depth)
    fit <- smooth.spline(rank, depth, spar=0.25)
    
    # find local minima slope
    #X <- data.frame(t=seq(min(rank),max(rank),length=nrow(df))) 
    #Y <- predict(fit, newdata=X, deriv=1) 
    #xvals <- Y$x
    #yvals <- Y$y
    #df.n <- data.frame(xvals=xvals,yvals=yvals,log10depth=depth,rank=rank)
    #df.n <- df.n[3000:100000,]
    #df.n <- subset(df.n, df.n$log10depth >= 3 & df.n$log10depth <= 5)
    #knee <- df.n$xvals[which.min(df.n$yvals)]
    #cells <- as.numeric(rownames(df.n[which.min(df.n$yvals),]))
    #reads <- as.integer(10^(as.numeric(df.n[which.min(df.n$yvals),]$log10depth)))
    
    df.n <- data.frame(log10depth=depth,rank=rank)
    df.n <- subset(df.n, df.n$log10depth >= 3)
    cells <- as.numeric(rownames(df.n[nrow(df.n),]))
    knee <- rank[cells]
    reads <- as.integer(10^(as.numeric(df.n[nrow(df.n),]$log10depth)))

    # plot
    plot(rank[(cells+1):length(rank)], depth[(cells+1):length(rank)], 
         type="l", lwd=2, col="grey75", main=main, 
         xlim=range(rank),
         ylim=range(depth),
         xlab="Barcode rank (log10)",
         ylab="Unique Tn5 insertions (log10)")
    lines(rank[1:cells], depth[1:cells], lwd=2, col="darkorchid4")
    grid()
    abline(v=knee, col="red", lty=2, lwd=1)
    abline(h=log10(reads), col="red", lty=2, lwd=1)
    text(2,1, labels=paste("# cells=",cells,", # reads=",reads,sep=""))
    return(head(x, n=cells))
}

# load data
# <- sub(".updated_metadata.txt","", input)
message(" - loading meta data for ", name)
a <- read.table(input)
a$unique <- a$total

pdf(paste(name,".QC_FIGURES.pdf",sep = ""), width=12, height=3)
layout(matrix(c(1:4), nrow=1))
###### plot potential knee and filter cells#######
message(" - calling cells ... for ", name)
#pdf(file=paste(name,"cell_calling.pdf", sep = "_"), width=4, height=4)
meta.v1 <- plotDist(a, main=name)
message(" - number of cells after cell calling v1 = ",nrow(meta.v1))
write.table(meta.v1, file=paste(name,".updated_metadata_v1.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)

###### plot TSS and filter cells ###########
#    meta.v1$prop <- meta.v1$tss/meta.v1$unique
#    threshs <- c()
    meta.v1$tss_z <- scale(meta.v1$pTSS)
    meta.v1 <- meta.v1[order(meta.v1$tss_z, decreasing=T),]
    below <- subset(meta.v1, meta.v1$tss_z < -2)
    if(!is.null(below)){
	    thresh <- max(0.2, max(below$pTSS))
    }else{
	    thresh <- 0.2
    }
    meta.v1 <- meta.v1[order(meta.v1$unique, decreasing=T),]
    meta.v2 <- subset(meta.v1, meta.v1$pTSS >= thresh)
#    threshs <- c(threshs, thresh)
    den <- kde2d(log10(meta.v1$unique), meta.v1$pTSS, 
                 n=300, h=c(0.2, 0.05),
                 lims=c(c(2.5, 6.5), c(0,1.0)))
    image(den, useRaster=T, col=c("white", rev(magma(100))),
          xlab="Unique Tn5 insertions (log10)", ylab="Fraction reads in TSS",
          main=name)
    grid(lty=1, lwd=0.5, col="grey90")
    abline(h=thresh, col="red", lty=2, lwd=1)
    legend("topright", legend=paste("# cells = ", nrow(meta.v2), sep=""), fill=NA, col=NA, border=NA)
    box()

	
message(" - number of cells after pTSS filter v2 = ",nrow(meta.v2))
write.table(meta.v2, file=paste(name,".updated_metadata_v2.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)

###### plot FRiP and filter cells ###########
#    meta.v2$prop <- meta.v2$tss/meta.v2$unique
    meta.v2 <- meta.v2[order(meta.v2$acr_z, decreasing=T),]
    below <- subset(meta.v2, meta.v2$acr_z < -2)
    if(!is.null(below)){
	    thresh <- max(0.2, max(below$FRiP))
    }else{
	    thresh <- 0.2
    }
    meta.v2 <- meta.v2[order(meta.v2$unique, decreasing=T),]
    meta.v3 <- subset(meta.v2, meta.v2$FRiP >= thresh)
#    threshs <- c(threshs, thresh)
    den <- kde2d(log10(meta.v2$unique), meta.v2$FRiP, 
                 n=300, h=c(0.2, 0.05),
                 lims=c(c(2.5, 6.5), c(0,1.0)))
    image(den, useRaster=T, col=c("white", rev(magma(100))),
          xlab="Unique Tn5 insertions (log10)", ylab="Fraction Tn5 insertions in ACRs",
          main=name)
    grid(lty=1, lwd=0.5, col="grey90")
    abline(h=thresh, col="red", lty=2, lwd=1)
    legend("topright", legend=paste("# cells = ", nrow(meta.v3), sep=""), fill=NA, col=NA, border=NA)
    box()

	
message(" - number of cells after FRiP filter v3 = ",nrow(meta.v3))
write.table(meta.v3, file=paste(name,".updated_metadata_v3.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)


###### plot pOrg and filter cells ###########
#    meta.v2$prop <- meta.v2$tss/meta.v2$unique
    meta.v3$pOrg_z <- as.numeric(scale(meta.v3$pOrg))
    meta.v3 <- meta.v3[order(meta.v3$pOrg_z, decreasing=T),]
    #below <- subset(meta.v3, meta.v3$pOrg_z > 2)
    #if(!is.null(below)){
	#    thresh <- min(0.3, min(below$pOrg))
    #}else{
	    thresh <- 0.3
    #}
    meta.v3 <- meta.v3[order(meta.v3$unique, decreasing=T),]
    meta.v4 <- subset(meta.v3, meta.v3$pOrg <= thresh)
#    threshs <- c(threshs, thresh)
    den <- kde2d(log10(meta.v3$unique), meta.v3$pOrg, 
                 n=300, h=c(0.2, 0.05),
                 lims=c(c(2.5, 6.5), c(0,0.5)))
    image(den, useRaster=T, col=c("white", rev(magma(100))),
          xlab="Unique Tn5 insertions (log10)", ylab="Fraction Tn5 insertions in Organelle",
          main=name)
    grid(lty=1, lwd=0.5, col="grey90")
    abline(h=thresh, col="red", lty=2, lwd=1)
    legend("topright", legend=paste("# cells = ", nrow(meta.v4), sep=""), fill=NA, col=NA, border=NA)
    box()

	
message(" - number of cells after pOrg filter v4 = ",nrow(meta.v4))
write.table(meta.v4, file=paste(name,".updated_metadata_v4.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)
meta.v5 <- subset(meta.v4, meta.v4$call == 1)
write.table(meta.v5, file=paste(name,".updated_metadata_v5.txt", sep = ""), sep="\t", quote=F, row.names=T, col.names=T)
dev.off()




# layout(matrix(c(1:12), nrow=2, byrow=T))
# outdat <- list()
# for(i in unique(a$library)){
	# df <- subset(a, a$library==i)
	# outdat[[i]] <- plotDist(df, main=i)
# }
# dev.off()

# # cells pass rd knee
# message(" - writing output ...")
# all <- do.call(rbind, outdat)
# message(" - number of cells = ",nrow(all))
# write.table(all, file=output, sep="\t", quote=F, row.names=T, col.names=T)
