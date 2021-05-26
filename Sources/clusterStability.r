clusterStability <- function(data=NULL, clustermethod=NULL, randomTests = 20, trainFraction = 0.5,pac.thr=0.1, ...)
{
  clusterLabels <- list();
  randomSamples <- list();
  numberofClusters <- 0;
  testCounts <- numeric(nrow(data))
  randomSeeds <- sample(randomTests);

  for (i in 1:randomTests)
  {
    randomSamples[[i]] <- sample(nrow(data),trainFraction*nrow(data));
    mod1 <- clustermethod(data[randomSamples[[i]],],...);
    clusterLabels[[i]] <- predict(mod1,data); #data[-randomSamples[[i]],]
    names(clusterLabels[[i]]$classification) <- rownames(data)
    plot(data[,1:2],col = clusterLabels[[i]]$classification,main=sprintf("%d",i));
    numberofClusters <- numberofClusters + length(table(clusterLabels[[i]]$classification))
    testCounts[-randomSamples[[i]]] <- testCounts[-randomSamples[[i]]] + 1;
    set.seed(randomSeeds[i]);
  }

  numberofClusters <- numberofClusters/randomTests;
  cat("Done Testing:")
  randIndex <- numeric();
  jaccIndex <- numeric();
  meanJaccard <- numeric();
  jaccardpoint <- numeric(nrow(data));
  jaccardpointcount <- numeric(nrow(data));
  trainrandIndex <- numeric();
  trainjaccIndex <- numeric();
  trainmeanJaccard <- numeric();
  trainjaccardpoint <- numeric(nrow(data));
  trainjaccardpointcount <- numeric(nrow(data));
  for (i in 1:(randomTests - 1))
  {
    for (j in (i + 1):randomTests)
    {
      outsamples <- unique(c(randomSamples[[i]],randomSamples[[j]]))
      if ((nrow(data) - length(outsamples)) > 10)
      {
        randIndex <- c(randIndex,adjustedRandIndex(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]));
        jaccard <- jaccardMatrix(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]);
        #cat(jaccard$balancedMeanJaccard)
        jaccIndex <- c(jaccIndex,jaccard$balancedMeanJaccard);
        meanJaccard <- c(meanJaccard,mean(jaccard$elementJaccard));
        jaccardpoint[-outsamples] <- jaccardpoint[-outsamples] + jaccard$elementJaccard;
        jaccardpointcount[-outsamples] <- jaccardpointcount[-outsamples] + 1;
      }
      insamples <- randomSamples[[i]][randomSamples[[i]] %in% randomSamples[[j]]]
      if ((nrow(data) - length(insamples)) > 10)
      {
        trainrandIndex <- c(trainrandIndex,adjustedRandIndex(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]));
        trainjaccard <- jaccardMatrix(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]);
        trainjaccIndex <- c(trainjaccIndex,trainjaccard$balancedMeanJaccard);
        trainmeanJaccard <- c(trainmeanJaccard,mean(trainjaccard$elementJaccard));
        trainjaccardpoint[insamples] <- trainjaccardpoint[insamples] + trainjaccard$elementJaccard;
        trainjaccardpointcount[insamples] <- trainjaccardpointcount[insamples] + 1;
      }
    }
  }
  cat("After Jacckard:")
  jaccardpoint[jaccardpointcount > 0] <- jaccardpoint[jaccardpointcount > 0]/jaccardpointcount[jaccardpointcount > 0];
   names(jaccardpoint) <- rownames(data);
  trainjaccardpoint[trainjaccardpointcount > 0] <- trainjaccardpoint[trainjaccardpointcount > 0]/trainjaccardpointcount[trainjaccardpointcount > 0];
  names(trainjaccardpoint) <- rownames(data);
  cat(nrow(data))
  testConsesus <- matrix(0,nrow = nrow(data), ncol = nrow(data))
  colnames(testConsesus) <- rownames(data)
  rownames(testConsesus) <- rownames(data)
  countMat <- testConsesus;
  dataConcensus <- testConsesus
  totwts <- 0;
  for (i in 1:randomTests)
  {
    cat("randomTests: ")
    cat(i)
    cat('\n')
    testset <- rownames(data[-randomSamples[[i]],])
    aclassLabels <- clusterLabels[[i]]$classification;
    nclus <- length(table(aclassLabels))
    wts <- (1.0-0.99*(nclus < 2))/(1.0+abs(nclus-numberofClusters));
    classLabels <- aclassLabels[testset];
    btestset <- rownames(data) %in% testset;
    for (id in testset)
    {
      testConsesus[id,btestset] <- testConsesus[id,btestset] + wts*(classLabels == aclassLabels[id]);
      countMat[id,btestset] <- countMat[id,btestset] + wts;
    }
    cat("End for one: ")
    cat(i)
    cat('\n')
    classLabels <- clusterLabels[[i]]$classification;
    for (id in 1:nrow(data))
    {
      if(id%%1000==0)
      {
        cat("for two: ")
        cat(id)
        cat('\n')
      }
      dataConcensus[id,] <- dataConcensus[id,] + wts*(classLabels == classLabels[id]);
    }
    totwts <- totwts + wts;
  }
  cat("After Counting.")
  testConsesus[countMat > 0] <- testConsesus[countMat > 0]/countMat[countMat > 0];
  dataConcensus <- dataConcensus/totwts;
  pac <- sum(testConsesus[(testConsesus > pac.thr) & (testConsesus < (1.0 - pac.thr))])/nrow(data)/nrow(data);


  result <- list(randIndex = randIndex,jaccIndex = jaccIndex,meanJaccard = meanJaccard,randomSamples = randomSamples,
                 clusterLabels=clusterLabels, jaccardpoint=jaccardpoint,averageNumberofClusters=numberofClusters,
                 testConsesus=testConsesus,trainRandIndex = trainrandIndex,trainJaccIndex = trainjaccIndex,trainMeanJaccard = trainmeanJaccard,
                 trainJaccardpoint=trainjaccardpoint,PAC=pac,dataConcensus=dataConcensus);
  class(result) <- "ClusterStability"
  return(result);
}


plot.ClusterStability <- function(object,...)
{
  plot(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
  boxplot(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
}

summary.ClusterStability <- function(object,...)
{
  summary(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
}
