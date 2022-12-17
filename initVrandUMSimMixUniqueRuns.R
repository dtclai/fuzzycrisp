#####################
# call packages
########################

rm(list=ls())
library(MASS)
library(fpc)
library(plyr)
library(cluster)
library(vcd)
library(psych)
library(e1071)
library(permute)
library(Rtsne)
library(MixSim)


##############
#FUNCTIONS
##############

#FUNCTION:COMPUTE CENTERS OF CLUSTERS(PROTOTYPES)
#EQ(7)
computeClusterCenters<-function(alpha, beta, memMat,dataMat){
  centroids<-matrix(0,dim(memMat)[1],dim(dataMat)[2])
  for(i in seq(1,dim(memMat)[1],by=1)){
    #calculate weight for each centre
    term1 <-((memMat[i,]*memMat[i,])/(2*alpha[i]))+ (memMat[i,]/beta[i])
    #print (paste0("term1::", term1))
    #multiply weight with datapoint x and sum it for all k
    nume<-apply((term1*dataMat),2, sum)
    deno<-sum(term1)
    #print(paste0("nume=",nume))
    #print(paste0("enod=",deno))
    centroids[i,]<-(nume/deno)
  }
  #print(centroids)
  return(centroids)
}

#FUNCTION: CALCULATE DISTANCE TO CENTRES
calculateDist<-function(dataMat, centMat){
  diMatrix<-matrix(0,dim(centMat)[1],dim(dataMat)[1])
  for(i in seq(1,dim(centMat)[1])){
    for(k in seq(1,dim(dataMat)[1])){
      temp<-as.matrix(dataMat[k,]- centMat[i,])
      diMatrix[i,k]<-temp%*%t(temp)
    }
  }
  return (diMatrix)
}

#FUNCTION:COMPUTE ALPHA
#EQ(9)
computeAlpha<-function(memMat,dataMat, centMat, distMat){
  alpha<-rep(0,dim(memMat)[1])
  nume<-0
  for(i in seq(1,dim(memMat)[1],by=1)){
    nume[i] <- sqrt(sum((memMat[i,]*memMat[i,]) * distMat[i,]))
  }
  alpha<-nume/sum(nume)
  return(alpha)
}

#FUNCTION:COMPUTE BETA
#EQ(10)
computeBeta<-function(memMat,dataMat, centMat, distMat){
  beta<-rep(0,dim(memMat)[1])
  nume<-0
  for(i in seq(1,dim(memMat)[1],by=1)){
    nume[i] <- sqrt(sum(memMat[i,] * distMat[i,]))
    # nume[i] <- sqrt(sum((memMat[i,]*memMat[i,]) * distMat[i,]))
  }
  beta<-nume/sum(nume)
  return(beta)
}

#FUNCTION:COMPUTE OBJECTIVE FUNCTION VALUE
#EQ(1)
calculateObjVal <- function(alpha, beta, memMat, distMat){
  obv <- matrix(0, dim(memMat)[1],dim(memMat)[2])
  for(i in seq(dim(memMat)[1])){
    obv[i,]<-(((memMat[i,]*memMat[i,])/(2*alpha[i])) + (memMat[i,]/beta[i]))*distMat[i,]
  }
  return(sum(apply(obv,2,sum)))
}

createData<-function(clustnum, dimnum, avoverlap, samplenum,clustratio){
  #use this
  K <-clustnum
  p <-dimnum
  repeat{
    Q <- MixSim(BarOmega = avoverlap, K = K, p = p)
    if (Q$fail == 0) break
  }
  #Q$Pi shows cluster ratios
  Q$Pi <- clustratio
  print("clusteratio")
  print(Q$Pi)
  # simulate a dataset of size 300 and add 10 outliers simulated on (0,1)x(0,1)
  # A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 10, int = c(0, 1))
  A <- simdataset(n = samplenum, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.out = 0, int = c(0, 0))
  colors <- c("red", "green", "blue", "brown", "magenta", "cyan")
  plot(A$X, xlab = "x1", ylab = "x2", type = "n")
  for (k in 0:K){
    points(A$X[A$id == k, ], col = colors[k+1], pch = 19, cex = 0.5)
  }
  return (cbind(A$X, A$id))
}

#Permutate cluster assignments; to align them with other cluster solutions
uniqueperm2 <- function(d) {
  dat <- factor(d)
  N <- length(dat)
  n <- tabulate(dat)
  ng <- length(n)
  if(ng==1) return(d)
  a <- N-c(0,cumsum(n))[-(ng+1)]
  foo <- lapply(1:ng, function(i) matrix(combn(a[i],n[i]),nrow=n[i]))
  out <- matrix(NA, nrow=N, ncol=prod(sapply(foo, ncol)))
  xxx <- c(0,cumsum(sapply(foo, nrow)))
  xxx <- cbind(xxx[-length(xxx)]+1, xxx[-1])
  miss <- matrix(1:N,ncol=1)
  for(i in seq_len(length(foo)-1)) {
    l1 <- foo[[i]]
    nn <- ncol(miss)
    miss <- matrix(rep(miss, ncol(l1)), nrow=nrow(miss))
    k <- (rep(0:(ncol(miss)-1), each=nrow(l1)))*nrow(miss) + 
      l1[,rep(1:ncol(l1), each=nn)]
    out[xxx[i,1]:xxx[i,2],] <- matrix(miss[k], ncol=ncol(miss))
    miss <- matrix(miss[-k], ncol=ncol(miss))
  }
  k <- length(foo)
  out[xxx[k,1]:xxx[k,2],] <- miss
  out <- out[rank(as.numeric(dat), ties="first"),]
  foo <- cbind(as.vector(out), as.vector(col(out)))
  out[foo] <- d
  t(out)
}

#@param d - number of dimensions
#@param N - number of objects
#@param DM - data matrix
#@param CL - class labels

#datalist<-c("wobc","wdbc","appendicitis","sonar",
#            "iris","wine","cmc","seeds", "glass",
#             "segment","vehicle","dermatology","ecoli")

#datalist<-c("iris","appendicitis",
#            "yeast-10",
#            "wine-quality-white-7",
#            "wine-quality-white-5",
#            "colon-2",
#            "breast-tissue-6",
#            "breast-tissue-4",
#            "cardiotocography-10",
#            "cardiotocography-3",
#            "ionosphere")

# datalist<-c("wine","iris")
#datalist<-c("continuousdata\\simdata6","iris")


#################
# Main
#############

datafolder <- "D:\\Data\\"
path <- "D:\\FuzzyCrisp"
set.seed(5)
for(datai in c(19)){
  runNum<-10
  run<-1
  dataseq<-datai  #7-4; 8-12
  datasetname<-c(paste0("continuousdata\\simdata", dataseq))

cn=0
dn=0
avo=0
sn=0
cr=c(0,0)


if(dataseq==20){
  dn=2
  avo=0.005
  sn=50
  cr=c(0.7,0.2,0.1) 
  cn=length(cr)
}else
if(dataseq==19){
  dn=2
  avo=0.005
  sn=50
  cr=c(0.8,0.1,0.1) 
  cn=length(cr)
}else
if(dataseq==18){
  dn=2
  avo=0.005
  sn=50
  cr=c(0.7,0.3) 
  cn=length(cr)
}else
if(dataseq==17){
  dn=2
  avo=0.005
  sn=50
  cr=c(0.8,0.2) 
  cn=length(cr)
}else
if(dataseq==16){
  dn=2
  avo=0.000
  sn=50
  cr=c(0.6,0.2,0.1,0.1) 
  cn=length(cr)
}else
if(dataseq==15){
  cn=3
  dn=2
  avo=0.000
  sn=50
  cr=c(0.7,0.2,0.1)  
}else if(dataseq==14){
  cn=2
  dn=2
  avo=0.001
  sn=30
  cr=c(0.9,0.1) 
}else if(dataseq==13){
  cn=2
  dn=2
  avo=0.001
  sn=30
  cr=c(0.7,0.3)  
} else if(dataseq==12){
  cn=2
  dn=2
  avo=0.05
  sn=200
  cr=c(0.7,0.3)
} else if (dataseq==11){
  cn=2
  dn=2
  avo=0.01
  sn=200
  cr=c(0.7,0.3)  
}else if(dataseq==10){
  cn=2
  dn=2
  avo=0.0001
  sn=200
  cr=c(0.8,0.2)  
}else if(dataseq==9){
  cn=2
  dn=2
  avo=0.0005
  sn=200
  cr=c(0.9,0.1)  
}else if(dataseq==8){
  cn=2
  dn=2
  avo=0.01
  sn=200
  cr=c(0.9,0.1)  
}else if(dataseq==7){
  cn=2
  dn=2
  avo=0.05
  sn=200
  cr=c(0.9,0.1)  
}else if(dataseq==6){
  cn=2
  dn=2
  avo=0.05
  sn=200
  cr=c(0.5,0.5)  
}else if(dataseq==5){
  cn=2
  dn=2
  avo=0.05
  sn=200
  cr=c(0.2,0.8)  
}else{
  cn=2
  dn=2
  avo=0.01
  sn=200
  cr=c(0.2,0.8)  
}


# firstdata<-createData(clustnum=cn, dimnum=dn, avoverlap=avo, samplenum=sn,clustratio=cr)

randIndex<-matrix(0,runNum,3)
entropy<-matrix(0,runNum,3)
sindex<-matrix(0,runNum,3)
dunn2<-matrix(0,runNum,3)

runresult=NULL
classresult=matrix(0, nrow=sn, ncol=runNum+1)
classeval<-NULL
fuzzycrispDuration<-NULL

#same dataset for all runs
# myData <- createData(clustnum=cn, dimnum=dn, avoverlap=avo, samplenum=sn,clustratio=cr)


#for(run in start:runNum){
while(run<=runNum){
  
  # putting it here will create new synthetic data 30 times; each run with unique dataset
  myData <- createData(clustnum=cn, dimnum=dn, avoverlap=avo, samplenum=sn,clustratio=cr)

  colnames(myData)
  #str(myData)
  
  datasetname<-basename(datasetname)
  DM<- as.data.frame(scale(myData[,-dim(myData)[2]]))
  #DM<- myData[,-dim(myData)[2]]
  str(DM)
  if(datasetname=="segment") DM<-DM[,-3]
  
  CL<-myData[,dim(myData)[2]]
  clusnum<-max(CL)
  cl<-clusnum
  clcount<-table(CL)
  
  f<-dim(DM)[2]
  N<-dim(DM)[1]
  
  object<-dim(DM)[1]
  dimension<-dim(DM)[2]
  clusternum<-max(CL)
  cl<-clusternum
  a<-NULL
  
  png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"truth_tsne.png"))
  tsne_results <- Rtsne(DM, perplexity=5, check_duplicates = FALSE) # You can change the value of perplexity and see how the plot changes
  
  ## Generate the t_SNE plot
  #par(mfrow=c(1,2)) # To plot two images side-by-side
  plot(tsne_results$Y,  bg= CL, pch = 21, cex = 1.5)
  dev.off()
  
  #if(run<5)   sink(file(paste0(path,"\\20\\",datasetname,"\\fixedbeta\\0",run,"_consoleLog.txt")), append=TRUE)
  #run <-1
  alltime<-NULL
  fcStart <- Sys.time()
  print(paste0("==Run::0",run,":::Data::",basename(datasetname),"=========="))
  
  perf<-NULL
  perf$obj<-NULL
  perf$iter<-NULL
  result<- NULL
  
  iteration<-0
  threshold<-1
  max.iter <-100
  objval<-0

  #random initialisation of membership
  UM<-matrix(runif(N*cl, 0,1), clusnum, N)
  a<- rep(1/clusternum,clusternum)
  b <- rep(1/clusternum,clusternum)
  # b <- rep(runif(1,0,1),clusternum)
   # b <- rep(1,clusternum)
   
  
  initV<-NULL
  for(dgrp in 1:cl){
    initV<-rbind(initV,apply(DM[sample(which(CL==dgrp),2),],2,mean))
  }
  
  #variables for keeping the last best
  prevUM<-NULL
  prevCentre<-NULL
  prevAlpha<-NULL
  
  print("####Initial States###")
  print("####Initial Membership####")
  print(UM)
  print("####Initial alpha####")
  print(a)
  print("####Initial beta####")
  print(b)

#while(threshold>0.001){
#ignore threshold>0.001 at end of iteration 0 and 1
while(iteration <= max.iter && threshold>0.0001){
   
  #2.1 Determining cluster centre 
  if (iteration == 0){
    centre <- initV
  }else{
    centre <- computeClusterCenters(a, b, UM,DM)
  }
  distM <- calculateDist(DM,centre)

  if (iteration > 0){
  	a<-computeAlpha(UM,DM,centre,distM)
  	#b<-computeBeta(UM,DM,centre,distM)

  }
  print(paste0("====alpha===for iteration:===",iteration))
  print(a)
  print(paste0("====beta===for iteration:===",iteration))
  print(b)
  
  prevUM<-UM
  prevCentre<-centre
  prevAlpha<-a
  
  #2.3 Determining membership xk
  #for each datapoint xk
  for(k in 1:N) {
    xk <- DM[k,]  #get each datapoint
    v <- 0           #step 1 
    Hv <- seq(1, cl)  #step 1  set Hv to be 1 to c
    stopflag <- TRUE
    print(paste0("====k=",k,"===="))
    
    while(!is.na(Hv) && stopflag==TRUE){
    print(paste("Hv", Hv))
      
    #calculate distance from xk to centres  
    kvdist <- apply(centre, 1, FUN=function(cent) dist(rbind(xk, cent)))
    kvdist<-kvdist*kvdist
    print(paste("kvdist", kvdist))
    
    #calculate distance to Hv
    distH<- (1 + sum(a[Hv]/b[Hv]))/sum(a[Hv]/kvdist[Hv]) #step2, Eq 13
    print(paste("distH", distH))
    #str(distH)
    
    if(all((kvdist[Hv]/b[Hv])<distH)) {   #step 3
      #Hv <- Hv[kvdist < distH]
      notHv <- tryCatch((1:cl)[! 1:cl %in% Hv],
                        error= integer(0))

      UM[Hv,k] <- a[Hv] * ((distH/kvdist)[Hv] - (1/b)[Hv]) #step 3, Eq 15
      if(length(notHv)!=0) {
        UM[notHv,k]<- 0
        print(paste0("fuzzy crisp assigned"))
      }
      stopflag=FALSE #stop the loop if this condition true: all(kvdist<distH)
    } else{
      Hv <- na.omit(Hv[(kvdist/b) < distH]) #sometimes Hv has NA values
    }
    if(length(Hv)<1) {
      print("Problem::::Hv length is less than 1")
      break
    }
    v <- v+1#step 5
    print(paste0("==========v= ", v, "==k= ", k, "=============iteration=", iteration))
    }#  while(!is.na(Hv) && stopflag==TRUE){

    print(paste0("====k=",k,"==UM[,k]=", UM[,k]))
    print(paste("sum membership k=", sum(UM[,k])))
    print(paste("v=",v))
  #str(UM)
  }#end loop for each datapoint
  prevObjVal <- objval
  newobjval <- calculateObjVal (a, b, UM, distM)
  threshold<- abs(objval - newobjval)
  #threshold<- (objval - newobjval)
  perf$obj<-c(perf$obj,newobjval)
  perf$iter<-c(perf$iter,iteration)
  objval <- newobjval
  print(paste0("==========ObjVal=",newobjval,"=====================Threshold=",threshold))
  iteration<-iteration + 1
} #end while threshold  

print(paste0("====membership===for iteration:===",iteration))
print(UM)

#if(run<5) sink() #end sink console log

alltime[1]<-Sys.time() - fcStart
print(UM)
labels<-apply(UM,2,which.max)

if(length(unique(labels))!=max(CL)) next
d <- dist(DM, method = "euclidean") # distance matrix

result<-cluster.stats(d, labels, CL)
result$perf<-perf
result$UM <- as.matrix(UM)

# if(run<=3){
# write.table(centre, file=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fcCentres.txt"), sep="\t")
# write.table(a, file=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fcAlpha.txt"), sep="\t")
# write.table(UM, file=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fcMembership.txt"), sep="\t")
# 
# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fc_perf.png"))
# plot(perf$iter, perf$obj, xlab = "Iterations", ylab="Objective function") #color is truth value
# dev.off()
# }
result$centre <- centre
result$a <- a
result$b<- b

kStart <- Sys.time()
# kclus <- tryCatch(expr = {kmeans(DM,initV, iter.max=max.iter)},error= function(e) NA)
kclus <- kmeans(DM,initV, iter.max=max.iter)

# if(is.na(kclus)) next
alltime[2] <- Sys.time() - kStart


fStart <- Sys.time()  
fclus <- cmeans(DM,initV, iter.max=max.iter)
alltime[3] <- Sys.time() - fStart
fuzzycrispDuration<- rbind(fuzzycrispDuration, alltime)
if(run<=3){
  write.table(t(fclus$membership), file=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fMembership.txt"), sep="\t")
  write.table(t(c(iteration, kclus$iter,fclus$iter)), 
              file=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\","Iterations.txt"), 
              sep="\t", append=TRUE, col.names = FALSE)
}
kresult<-cluster.stats(d,kclus$cluster, CL)
fresult<-cluster.stats(d,fclus$cluster, CL)
randIndex[run,]<-c(result$corrected.rand,kresult$corrected.rand, fresult$corrected.rand)
entropy[run,]<-c(result$entropy ,kresult$entropy, fresult$entropy)
sindex[run,]<-c(result$sindex ,kresult$sindex, fresult$sindex)
dunn2[run,]<-c(result$dunn2 ,kresult$dunn2, fresult$dunn2)

sink(paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_result.txt"))
print("===========Fuzzy-crisp results===========")
print(result)
print("===========Kmeans results===========")
print(kresult)
print("===========FCM results===========")
print(fresult)
sink()

#Draw fc plots
# Centroid Plot against 1st 2 discriminant functions
vcol <- c("blue","green","purple","red","yellow", "black", 
          "pink", "cyan", "brown", "gray")
if(run<=10){
# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_truth_fpc.png"))
# plotcluster(DM, CL, col=vcol[CL]) #color is truth value
# dev.off()

png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fc_fpc.png"))
plotcluster(DM, labels, col=vcol[CL]) #color is truth value
dev.off()

# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fc_cluster.png"))
# clusplot(DM, labels, color=TRUE, shade=TRUE,
#          labels=2, lines=0, col.p=CL)
# dev.off()

# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_fc_tsne.png"))
# plot(tsne_results$Y,  bg= labels, pch = 21, cex = 1.5)
# dev.off()

# Draw kclus plots
png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_k_fpc.png"))
plotcluster(DM, as.numeric(kclus$cluster), col=vcol[CL])
dev.off()

# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_k_cluster.png"))
# clusplot(DM, kclus$cluster, color=TRUE, shade=TRUE,
#          labels=2, lines=0, col.p=CL)
# dev.off()

# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_k_tsne.png"))
# plot(tsne_results$Y,  bg= kclus$cluster, pch = 21, cex = 1.5)
# dev.off()

# Draw fclus plots
png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_f_fpc.png"))
vcol <- c("blue","green","purple","red","yellow", "black", 
          "pink", "cyan", "brown", "gray")
plotcluster(DM, as.numeric(fclus$cluster), col=vcol[CL])
dev.off()

# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_f_cluster.png"))
# clusplot(DM, fclus$cluster, color=TRUE, shade=TRUE,
#          labels=2, lines=0, col.p=CL)
# dev.off()

# png(filename=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\0",run,"_f_tsne.png"))
# plot(tsne_results$Y,  bg= fclus$cluster, pch = 21, cex = 1.5)
# dev.off()
}
sink(paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\index.txt"))
print("===========RandIndex Fuzzy Crisp VS KMeans VS FCM===========")
print(randIndex)
print(c("Mean",apply(randIndex,2,
                     function(x){round(mean(x),digits=3)})))
print(c("StDev",apply(randIndex,2,
                      function(x){round(sd(x),digits=3)})))
print(c("Best",apply(randIndex,2,
                     function(x){round(max(x),digits=3)})))
print(c("Min",apply(randIndex,2,
                    function(x){round(min(x),digits=3)})))


print("===========Entropy Fuzzy Crisp VS KMeans VS FCM===========")
print(entropy)
print(c("Mean",apply(entropy,2,
                     function(x){round(mean(x),digits=3)})))
print(c("StDev",apply(entropy,2,
                      function(x){round(sd(x),digits=3)})))
print(c("Best",apply(entropy,2,
                     function(x){round(max(x),digits=3)})))
print(c("Min",apply(entropy,2,
                    function(x){round(min(x),digits=3)})))

print("===========Sindex Fuzzy Crisp VS KMeans VS FCM===========")
print(sindex)
print(c("Mean",apply(sindex,2,
                     function(x){round(mean(x),digits=3)})))
print(c("StDev",apply(sindex,2,
                      function(x){round(sd(x),digits=3)})))
print(c("Best",apply(sindex,2,
                     function(x){round(max(x),digits=3)})))
print(c("Min",apply(sindex,2,
                    function(x){round(min(x),digits=3)})))

print("===========Dunn2 Fuzzy Crisp VS KMeans VS FCM===========")
print(dunn2)
print(c("Mean",apply(dunn2,2,
                     function(x){round(mean(x),digits=3)})))
print(c("StDev",apply(dunn2,2,
                      function(x){round(sd(x),digits=3)})))
print(c("Best",apply(dunn2,2,
                     function(x){round(max(x),digits=3)})))
print(c("Min",apply(dunn2,2,
                    function(x){round(min(x),digits=3)})))
sink()

#evaluate cluster solutions
runresult[[run]]<-cluster.stats(d = dist(DM), clustering=as.numeric(labels), alt.clustering = CL)
classresult[,run]<-labels
actual = CL # actual labels
predicted = labels # predicted labels
cm = as.matrix(table(Actual = actual, Predicted = predicted)) # create the confusion matrix
cm

conMatrix<-NULL
n<-0
nc<-0
diag<-NULL
accuracy<-0
bestPredicted<-NULL
bestAccuracy<-0

#realignment of labels    
permuteMat <- uniqueperm2(1:clusnum)

for(perj in seq(1,dim(permuteMat)[1])){
  #print(predicted)
  newpredicted<-predicted
  #print(permuteMat[perj,])
  #levels(newpredicted)= permuteMat [perj,]
  #print(levels(newpredicted))
  #print(newpredicted)
  newnewpredicted <- factor(newpredicted, levels = permuteMat [perj,])
  
  conMatrix = table(actual, newnewpredicted) # create the confusion matrix
  #print(conMatrix)
  n = sum(conMatrix) # number of instances
  nc = nrow(conMatrix) # number of classes
  diag = diag(conMatrix) # number of correctly classified instances per class 
  #print(diag)
  
  accuracy = sum(diag) / n
  print(accuracy)
  
  if(perj==1||accuracy>bestAccuracy){
    bestPredicted=newnewpredicted
    bestAccuracy=accuracy
  }
}#end for permutation

print(bestAccuracy)
predicted=bestPredicted
conMatrix = table(actual, predicted) # create the confusion matrix
conMatrix

n = sum(conMatrix) # number of instances
nc = nrow(conMatrix) # number of classes
diag = diag(conMatrix) # number of correctly classified instances per class 

accuracy = sum(diag) / n
#print(levels(predicted)) 
#print(accuracy)

rowsums = apply(conMatrix, 1, sum) # number of instances per class
colsums = apply(conMatrix, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes

precision = na.omit(diag / colsums)
recall = na.omit(diag / rowsums) 

f1 = sum(na.omit(2 * precision * recall / (precision + recall))) /clusnum

if(length(precision) == length(recall)){  
  expAccuracy = sum(p*q)
  kappa = (accuracy - expAccuracy) / (1 - expAccuracy)
  kappa
  
  #Kappa is always less than or equal to 1. 
  #A value of 1 implies perfect agreement and values less than 1 imply less than perfect agreement.
  #In rare situations, Kappa can be negative. 
  #This is a sign that the two observers agreed less than would be expected just by chance.
  
  Kappa(conMatrix)
  cohen.kappa(conMatrix)
  cohenKappa<-Kappa(conMatrix)$Unweighted[1]
}
else{
  f1=NA
  kappa=NA
  cohenKappa=NA
}

if(is.nan(f1)) f1=NA

print(labels)

classeval[[run]]<-list(cm, accuracy, 
                       mean(precision),
                       mean(recall),
                       mean(f1),
                       Kappa(conMatrix),
                       cohen.kappa(conMatrix),
                       cohenKappa,
                       data.frame(precision, recall, f1))
print("==Accuracy==")
print(accuracy)
print("==f1==")
print(f1)

run=run+1

print(paste0("dataset::: ",basename(datasetname)))
print(paste0("run::: 0",run))
}#end for run
runresult[[run]]<-cluster.stats(d = dist(DM), CL)
classresult[,run]<-CL
runresult.df = as.data.frame(do.call(rbind, runresult))
classeval.df = as.data.frame(do.call(rbind, list(classeval)))

foldername=paste0(path,"\\20\\",basename(datasetname),"\\fixedbeta\\")
#write.table(as.matrix(spiterFrame),paste0(foldername,"spiter_",filename), sep="\t")
#write.table(as.matrix(iterFrame),paste0(foldername,"iter_",filename), sep="\t")
write.table(t(as.matrix(runresult.df)),paste0(foldername,"runresults.txt"), sep="\t")
write.table(t(as.matrix(classresult)),paste0(foldername,"classresults.txt"), sep="\t")
write.table(t(as.matrix(classeval.df)),paste0(foldername,"classevalresults.txt"), sep="\t")
write.table(t(as.matrix(fuzzycrispDuration)),paste0(foldername,"fctime.txt"), sep="\t")
# } #end datalist

print(datasetname)
}#end datai loop



