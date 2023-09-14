# run Socrates on merged socrates object #

# libraries
library(Socrates)
library(harmony)
library(symphony)
library(igraph)
library(Matrix)

args <- commandArgs(trailingOnly=T)
if(length(args) != 6){stop("Rscript plot_marker_accessibility.R [rds] [meta] [feature_rate] [resolution] [PC_num] [prefix]")}

rds <- as.character(args[1])
meta <- as.character(args[2])
feature_rate <- as.numeric(args[3])
res <- as.numeric(args[4])
pcs <- as.numeric(args[5])
out <- as.character(args[6])



#module load R/4.1.0-foss-2019b

# functions --------------------------------------------------------------

# new TFIDF methods
tfidf <- function(obj,
                  frequencies=T,
                  log_scale_tf=T,
                  scale_factor=10000,
                  doL2=F,
                  slotName="residuals"){
  
  # set bmat
  bmat <- obj$counts
  
  # hidden functions
  .safe_tfidf       <- function(tf, idf,  block_size=2000e6){
    result = tryCatch({
      result = tf * idf
      result
    }, error = function(e) {
      options(DelayedArray.block.size=block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      
      tf = DelayedArray(tf)
      idf = as.matrix(idf)
      
      result = tf * idf
      result
    })
    return(result)
  }
  
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = .safe_tfidf(tf, idf)
  
  # do L2?
  if(doL2){
    l2norm <- function(x){x/sqrt(sum(x^2))}
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x=1/colNorm)
  }
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse=T)
  obj$norm_method <- "tfidf"
  
  # return
  return(obj)
}

# call clusters
callClusters  <- function(obj,
                          res=0.4,
                          k.near=30,
                          clustOB="svd",
                          cname="LouvainClusters",
                          min.reads=5e4,
                          m.clst=25,
                          threshold=5,
                          umap1="umap1",
                          umap2="umap2",
                          e.thresh=3,
                          cleanCluster=T,
                          cl.method=1,
                          svd_slotName="PCA",
                          umap_slotName="UMAP",
                          cluster_slotName="Clusters",
                          verbose=FALSE,
                          ...){
  
  # funcs
  filterSingle  <- function(pro,
                            k=50,
                            threshold=3,
                            type="umap",
                            m1="umap1",
                            m2="umap2"){
    
    # set column names
    if(type=="umap"){
      vars <- c(m1, m2)
    }else if(type=="pca"){
      vars <- colnames(pro)
    }
    
    # ensure that k is less than number of samples
    if(k > nrow(pro)){
      k <- nrow(pro)-1
    }
    
    # get nearest neighbors
    topk <- FNN::get.knn(pro[,vars], k=k)
    cell.dists <- as.matrix(topk$nn.dist)
    rownames(cell.dists) <- rownames(pro)
    colnames(cell.dists) <- paste0("k",seq(1:ncol(cell.dists)))
    aves <- apply(cell.dists, 1, mean)
    zscore <- as.numeric(scale(aves))
    names(zscore) <- rownames(pro)
    
    # thresholds
    p.zscore <- zscore[order(zscore, decreasing=T)]
    num.pass <- length(zscore[zscore < threshold])
    
    # filter
    prop.good <- zscore[zscore < threshold]
    ids <- names(prop.good)
    out <- pro[rownames(pro) %in% ids,]
    
    # return object
    return(out)
  }
  filtDistClst <- function(b,
                           umap1="umap1",
                           umap2="umap2",
                           threshold2=2){
    
    # iterate over each cluster
    clusts <- unique(b$seurat_clusters)
    out <- lapply(clusts, function(x){
      b.sub <- subset(b, b$seurat_clusters == x)
      b.umap <- b.sub[,c(umap1, umap2)]
      out.b.umap <- filterSingle(b.umap, k=25, threshold=threshold2, m1=umap1, m2=umap2)
      return(rownames(out.b.umap))
    })
    out <- do.call(c, out)
    b.out <- b[rownames(b) %in% as.character(out),]
    message("   * total number of cells surviving subcluster filtering = ", nrow(b.out))
    return(b.out)
    
  }
  
  # filter umap coordinates
  if(verbose){message(" - filtering outliers in UMAP manifold (z-score e.thresh = ", e.thresh, ") ...")}
  umap.original <- obj[[umap_slotName]]
  umap.filtered <- filterSingle(obj[[umap_slotName]], threshold=e.thresh, m1=umap1, m2=umap2)
  counts.filtered <- obj$counts[,rownames(umap.filtered)]
  meta.filtered <- obj$meta[rownames(umap.filtered),]
  pca.filtered <-  obj[[svd_slotName]][rownames(umap.filtered),]
  
  # run graph-based clustering
  if(verbose){message(" - creating seurat object for graph-based clustering ...")}
  
  # create Seurat object, find clusters
  sro <- SeuratObject::CreateSeuratObject(counts.filtered, min.cells=0, min.features=0)
  sro[["svd"]] <- CreateDimReducObject(embeddings = pca.filtered, key = "PC_", assay = DefaultAssay(sro))
  sro[["umap"]] <- CreateDimReducObject(embeddings=as.matrix(umap.filtered), key="UMAP_", assay=DefaultAssay(sro))
  sro <- AddMetaData(sro, meta.filtered)
  nn.eps.val <- 0
  n.starts <- 100
  
  sro <- FindNeighbors(sro, dims = 1:ncol(sro[[clustOB]]), reduction=clustOB, 
                       nn.eps=nn.eps.val, k.param=k.near, annoy.metric="euclidean")
  if(cl.method==4){
    graph.obj <- graph_from_adjacency_matrix(as(sro[[paste0(DefaultAssay(sro),"_snn")]], "dgCMatrix"),
                                             mode="undirected",
                                             weighted=T)
    
    graph.obj <- cluster_leiden(graph.obj, 
                                objective_function="modularity",
                                resolution_parameter=res)
    
    clustered <- membership(graph.obj)
    sro.meta <- data.frame(sro@meta.data)
    sro.meta$seurat_clusters <- factor(as.numeric(clustered[rownames(sro@meta.data)]))
    
  }else{
    sro <- FindClusters(sro, resolution=res, n.start=n.starts, algorithm=cl.method, ...)
    sro.meta <- data.frame(sro@meta.data)
    sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
  }  
  if(verbose){message(" - finished graph-based clustering ...")}
  
  # remove temp
  rm(umap.filtered)
  rm(counts.filtered)
  rm(meta.filtered)
  rm(pca.filtered)
  suppressMessages(gc())
  
  # remove outliers?
  if(cleanCluster){
    
    # verbose
    if(verbose){message("   * removing low quality clusters ...")}
    
    # prep
    sro.umap <- obj[[umap_slotName]][rownames(sro.meta),]
    colnames(sro.umap) <- c(umap1, umap2)
    sro.meta <- cbind(sro.meta, sro.umap)
    
    # filter by cluster size and number of cells
    agg.reads <- aggregate(sro.meta$nSites~sro.meta$seurat_clusters, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(sro.meta$seurat_clusters)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells >= m.clst & agg.reads$readDepth >= min.reads)
    sro.filt <- sro.meta[as.character(sro.meta$seurat_clusters) %in% as.character(agg.pass$clusters),]
    sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
    sro.filt$seurat_clusters <- as.numeric(factor(sro.filt$seurat_clusters))
    
    # remove outliers in the embedding
    if(verbose){message("   * filtering per-cluster outliers (z-score filtDistClst2 = ", threshold, ") ...")}
    sro.meta <- filtDistClst(sro.filt, umap1=umap1, umap2=umap2, threshold=threshold)
    sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
    
  }
  
  # filter by cluster size
  if(verbose){message(" - filtering clusters with low cell/read counts ...")}
  agg.reads <- aggregate(sro.meta$nSites~sro.meta$seurat_clusters, FUN=sum)
  colnames(agg.reads) <- c("clusters","readDepth")
  clust.cnts <- table(sro.meta$seurat_clusters)
  agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
  agg.pass <- subset(agg.reads, agg.reads$num_cells>=m.clst & agg.reads$readDepth>=min.reads)
  sro.filt <- sro.meta[sro.meta$seurat_clusters %in% agg.pass$clusters,]
  sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
  sro.filt$seurat_clusters <- as.numeric(as.factor(sro.filt$seurat_clusters))
  
  # rename output
  final.UMAP <- umap.original[rownames(sro.filt),]
  clusts <- sro.filt$seurat_clusters
  obj[[cluster_slotName]] <- obj$meta[rownames(sro.filt),]
  obj[[cluster_slotName]][,cname] <- factor(clusts)
  obj[[cluster_slotName]][,c(umap1)] <- final.UMAP[,c(umap1)]
  obj[[cluster_slotName]][,c(umap2)] <- final.UMAP[,c(umap2)]
  
  # return
  return(obj)
}


# load rds files and pre-processing --------------------------------------
soc.obj <- readRDS(rds)
meta.data <- read.table(meta)
rownames(meta.data) <- meta.data$cellID

# harmonize
soc.obj$meta <- meta.data
soc.obj$counts <- soc.obj$counts[,colnames(soc.obj$counts) %in% rownames(meta.data)]
soc.obj$meta <- soc.obj$meta[colnames(soc.obj$counts),]


# get per cell feature counts --------------------------------------------
cell.counts <- log10(Matrix::colSums(soc.obj$counts))  # count number of features with Tn5 insertions per cell
cell.counts.z <- as.numeric(scale(cell.counts)) # convert features counts into Z-scores
cell.counts.threshold <- max(c((10^cell.counts[cell.counts.z < -3]), 100)) # minimum feature counts (greater of 1 std or 1000)


# clean sparse counts matrix ---------------------------------------------
soc.obj <- cleanData(soc.obj, 
                     min.c=cell.counts.threshold,  # minimum number of accessible features per cell
                     min.t=0.0015, # minimum feature frequency across cells
                     max.t=0, # maximum feature frequency across cells
                     verbose=T)


# normalize with TFIDF ---------------------------------------------------
soc.obj <- tfidf(soc.obj, doL2=T)
#number.sites <- 50000 #ceiling(nrow(soc.obj$counts)*0.5)
number.sites <- ceiling(nrow(soc.obj$counts)*feature_rate)

# project with SVD -------------------------------------------------------
soc.obj <- reduceDims(soc.obj,
                      method="SVD", 
                      n.pcs=pcs,
                      cor.max=0.5,
                      num.var=number.sites,
                      verbose=T,
                      scaleVar=T,
                      doSTD=F,
                      doL1=F,
                      doL2=T,
                      refit_residuals=F)
					  
soc.obj <- projectUMAP(soc.obj, metric="cosine", k.near=15, svd_slotName="PCA", umap_slotName="UMAP")

#detect and remove potential doublets-----------------------------------
soc.obj <- detectDoublets(soc.obj, threads=8, nTrials=5, nSample=1000, rdMethod = "SVD", svd_slotName="PCA")
soc.obj <- filterDoublets(soc.obj, filterRatio=0.05, embedding="UMAP", umap_slotname="UMAP", svd_slotname="PCA", removeDoublets=T)

# extract feature loadings and var/mean of tfidf
ids <- rownames(soc.obj$PCA)
soc.obj$PCA.nonintegrated <- soc.obj$PCA

# remove batch effects with harmony --------------------------------------
soc.obj$PCA <- HarmonyMatrix(soc.obj$PCA, meta_data=soc.obj$meta, 
                         vars_use=c("genotype"), 
                         do_pca=F,
                         theta=c(2),
                         tau=c(5),
                         sigma=0.1,
                         lambda=c(0.1),
                         nclust=50,
                         max.iter.cluster=100,
                         max.iter.harmony=30)
colnames(soc.obj$PCA) <- paste0("PC_", 2:(ncol(soc.obj$PCA)+1))
rownames(soc.obj$PCA) <- ids


# reduce to 2-dimensions with UMAP ---------------------------------------
umap.modelout <- umap(soc.obj$PCA, 
                      metric="cosine",
                      n_neighbors=30,
                      a=1.95, b=0.75,
                      ret_model=T)

# update UMAP slot
soc.obj$UMAP <- umap.modelout$embedding
colnames(soc.obj$UMAP) <- c("umap1","umap2")
rownames(soc.obj$UMAP) <- rownames(soc.obj$PCA)


# identify clusters using neighborhood graph -----------------------------
soc.obj <- callClusters(soc.obj, 
                        res=res,
                        k.near=30,
                        verbose=T,
                        cleanCluster=F,
                        cl.method=4,
                        e.thresh=3,
                        threshold=3,
                        m.clst=100)


# plot cluster membership on UMAP embedding ------------------------------
pdf(paste0(out,".UMAP.clusters.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", cex=0.2)
dev.off()

pdf(paste0(out,".UMAP.nSites.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", column="log10nSites", cex=0.2)
dev.off()

pdf(paste0(out,".UMAP.tn5B_bc.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", column="tn5B_bc", cex=0.2)
dev.off()

pdf(paste0(out,".UMAP.genotype.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", column="genotype", cex=0.2)
dev.off()


# save data --------------------------------------------------------------
saveRDS(soc.obj, file=paste0(out,".socrates.processed.rds"))

# output text files
meta <- soc.obj$Clusters
rd <- soc.obj$PCA

# write data
write.table(meta, file=paste0(out, ".metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(rd, file=paste0(out, ".reduced_dimensions.txt"), quote=F, row.names=T, col.names=T, sep="\t")


