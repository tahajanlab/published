### source - https://github.com/sqjin/CellChat/issues/167
### Copied + pasted with no due diligence.
### Check these functions over if you want to be sure of the code.

computeNetSimilarity_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, thresh = NULL) {
     type <- match.arg(type)
     prob = methods::slot(object, slot.name)$prob
     if (is.null(k)) {
         if (dim(prob)[3] <= 25) {
             k <- ceiling(sqrt(dim(prob)[3]))
           } else {
               k <- ceiling(sqrt(dim(prob)[3])) + 1
             }
       }
     if (!is.null(thresh)) {
         prob[prob < quantile(c(prob[prob != 0]), thresh)] <- 0
       }
     if (type == "functional") {
         # compute the functional similarity
           D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
           S2 <- D_signalings; S3 <- D_signalings;
           for (i in 1:(dim(prob)[3]-1)) {
               for (j in (i1):dim(prob)[3]) {
                   Gi <- (prob[ , ,i] > 0)*1
                   Gj <- (prob[ , ,j] > 0)*1
                   S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
                 }
             }
           # define the similarity matrix
             S3[is.na(S3)] <- 0; S3 <- S3 + t(S3); diag(S3) <- 1
             # S_signalings <- S1 *S2
               S_signalings <- S3
             } else if (type == "structural") {
                 # compute the structure distance
                   D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
                   for (i in 1:(dim(prob)[3]-1)) {
                       for (j in (i1):dim(prob)[3]) {
                           Gi <- (prob[ , ,i] > 0)*1
                           Gj <- (prob[ , ,j] > 0)*1
                           D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
                         }
                     }
                   # define the structure similarity matrix
                     D_signalings[is.infinite(D_signalings)] <- 0
                     D_signalings[is.na(D_signalings)] <- 0
                     D_signalings <- D_signalings + t(D_signalings)
                     S_signalings <- 1-D_signalings
                   }
     # smooth the similarity matrix using SNN
       SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
       Similarity <- as.matrix(S_signalings*SNN)
       rownames(Similarity) <- dimnames(prob)[[3]]
       colnames(Similarity) <- dimnames(prob)[[3]]
       comparison <- "single"
       comparison.name <- paste(comparison, collapse = "-")
       if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
           methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
         }
       methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
       return(object)
}


computeNetSimilarityPairwise_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, thresh = NULL) {
         type <- match.arg(type)
         if (is.null(comparison)) {
             comparison <- 1:length(unique(object@meta$datasets))
           }
         cat("Compute signaling network similarity for datasets", as.character(comparison), '\n')
         comparison.name <- paste(comparison, collapse = "-")
         net <- list()
         signalingAll <- c()
         object.net.nameAll <- c()
         # 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))
           for (i in 1:length(comparison)) {
               object.net <- methods::slot(object, slot.name)[[comparison[i]]]
               object.net.name <- names(methods::slot(object, slot.name))[comparison[i]]
               object.net.nameAll <- c(object.net.nameAll, object.net.name)
               net[[i]] = object.net$prob
               signalingAll <- c(signalingAll, paste0(dimnames(net[[i]])[[3]], "--", object.net.name))
               # signalingAll <- c(signalingAll, dimnames(net[[i]])[[3]])
               }
         names(net) <- object.net.nameAll
         net.dim <- sapply(net, dim)[3,]
         nnet <- sum(net.dim)
         position <- cumsum(net.dim); position <- c(0,position)
         if (is.null(k)) {
             if (nnet <= 25) {
                 k <- ceiling(sqrt(nnet))
               } else {
                   k <- ceiling(sqrt(nnet)) + 1
                 }
           }
         if (!is.null(thresh)) {
             for (i in 1:length(net)) {
                 neti <- net[[i]]
                 neti[neti < quantile(c(neti[neti != 0]), thresh)] <- 0
                 net[[i]] <- neti
               }
           }
         if (type == "functional") {
             # compute the functional similarity
               S3 <- matrix(0, nrow = nnet, ncol = nnet)
               for (i in 1:nnet) {
                   for (j in 1:nnet) {
                       idx.i <- which(position - i >= 0)[1]
                       idx.j <- which(position - j >= 0)[1]
                       net.i <- net[[idx.i-1]]
                       net.j <- net[[idx.j-1]]
                       Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
                       Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
                       S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
                     }
                 }
               # define the similarity matrix
                 S3[is.na(S3)] <- 0;  diag(S3) <- 1
                 S_signalings <- S3
               } else if (type == "structural") {
                   # compute the structure distance
                     D_signalings <- matrix(0, nrow = nnet, ncol = nnet)
                     for (i in 1:nnet) {
                         for (j in 1:nnet) {
                             idx.i <- which(position - i >= 0)[1]
                             idx.j <- which(position - j >= 0)[1]
                             net.i <- net[[idx.i-1]]
                             net.j <- net[[idx.j-1]]
                             Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
                             Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
                             D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
                           }
                       }
                     # define the structure similarity matrix
                       D_signalings[is.infinite(D_signalings)] <- 0
                       D_signalings[is.na(D_signalings)] <- 0
                       S_signalings <- 1-D_signalings
                     }
         # smooth the similarity matrix using SNN
           SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
           Similarity <- as.matrix(S_signalings*SNN)
           rownames(Similarity) <- signalingAll
           colnames(Similarity) <- rownames(Similarity)
           if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
               methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
             }
           # methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
             methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
             return(object)
           }

### netEmbedding
netEmbedding_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, pathway.remove = NULL, k = NULL) {
           if (object@options$mode == "single") {
               comparison <- "single"
               cat("Manifold learning of the signaling networks for a single dataset", '\n')
             } else if (object@options$mode == "merged") {
                 if (is.null(comparison)) {
                     comparison <- 1:length(unique(object@meta$datasets))
                   }
                 cat("Manifold learning of the signaling networks for datasets", as.character(comparison), '\n')
               }
           comparison.name <- paste(comparison, collapse = "-")
           Similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
           if (is.null(pathway.remove)) {
               pathway.remove <- rownames(Similarity)[which(colSums(Similarity) == 1)]
             }
           if (length(pathway.remove) > 0) {
               pathway.remove.idx <- which(rownames(Similarity) %in% pathway.remove)
               Similarity <- Similarity[-pathway.remove.idx, -pathway.remove.idx]
             }
           if (is.null(k)) {
               k <- ceiling(sqrt(dim(Similarity)[1])) + 1
             }
           options(warn = -1)
           # dimension reduction
             Y <- runUMAP(Similarity)
             if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$dr)) {
                 methods::slot(object, slot.name)$similarity[[type]]$dr <- NULL
               }
             methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y
             return(object)
}



netClustering_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, methods = "kmeans", do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, k.eigen = NULL) {
     type <- match.arg(type)
     if (object@options$mode == "single") {
         comparison <- "single"
         cat("Classification learning of the signaling networks for a single dataset", '\n')
       } else if (object@options$mode == "merged") {
           if (is.null(comparison)) {
               comparison <- 1:length(unique(object@meta$datasets))
             }
           cat("Classification learning of the signaling networks for datasets", as.character(comparison), '\n')
         }
     comparison.name <- paste(comparison, collapse = "-")
     
       Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
       Y[is.na(Y)] <- 0
       data.use <- Y
       if (methods == "kmeans") {
           if (!is.null(k)) {
               clusters = kmeans(data.use,k,nstart=10)$cluster
             } else {
                 N <- nrow(data.use)
                 kRange <- seq(2,min(N-1, 10),by = 1)
                 if (do.parallel) {
                     future::plan("multiprocess", workers = nCores)
                     options(future.globals.maxSize = 1000 * 1024^2)
                   }
                 my.sapply <- ifelse(
                     test = future::nbrOfWorkers() == 1,
                     yes = pbapply::pbsapply,
                     no = future.apply::future_sapply
                   )
                 results = my.sapply(
                     X = 1:length(kRange),
                     FUN = function(x) {
                         idents <- kmeans(data.use,kRange[x],nstart=10)$cluster
                         clusIndex <- idents
                         #adjMat0 <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
                           adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
                           return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
                         },
                     simplify = FALSE
                   )
                 adjMat <- lapply(results, "[[", 1)
                 CM <- Reduce('+', adjMat)/length(kRange)
                 res <- computeEigengap(as.matrix(CM))
                 numCluster <- res$upper_bound
                 clusters = kmeans(data.use,numCluster,nstart=10)$cluster
                 if (do.plot) {
                     gg <- res$gg.obj
                     ggsave(filename= paste0("estimationNumCluster_",fig.id,"_",type,"_dataset_",comparison.name,".pdf"), plot=gg, width = 3.5, height = 3, units = 'in', dpi = 300)
                   }
               }
         } else if (methods == "spectral") {
             A <- as.matrix(data.use)
             D <- apply(A, 1, sum)
             L <- diag(D)-A                       # unnormalized version
             L <- diag(D^-0.5)%*%L%*% diag(D^-0.5) # normalized version
             evL <- eigen(L,symmetric=TRUE)  # evL$values is decreasing sorted when symmetric=TRUE
             # pick the first k first k eigenvectors (corresponding k smallest) as data points in spectral space
               plot(rev(evL$values)[1:30])
             Z <- evL$vectors[,(ncol(evL$vectors)-k.eigen1):ncol(evL$vectors)]
             clusters = kmeans(Z,k,nstart=20)$cluster
           }
       if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
           methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
         }
       methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
       return(object)
     }

