# Label the subjects that shere the same connectivity

# Label the subjects that shere the same connectivity
getConsensusCluster <- function(object,who="training",thr=seq(0.70,0.30,-0.05))
{
  
  orgnames <-  rownames(object$dataConcensus);
  if (who != "training")
  {
    orgnames <-  rownames(object$testConsesus);
    pointJaccard <- object$jaccardpoint;
    names(pointJaccard) <- orgnames;
    concensusMat <- object$testConsesus[order(-pointJaccard),]
  }
  else
  {
    pointJaccard <- 0.99*object$trainJaccardpoint + 0.01*object$jaccardpoint;
    names(pointJaccard) <- orgnames;
    concensusMat <- object$dataConcensus[order(-pointJaccard),]
  }
  concensusMat <- concensusMat[,order(-pointJaccard)]
  classID <- numeric(nrow(concensusMat));
  pointJaccard <- pointJaccard[order(-pointJaccard)];
  names(classID) <-  names(pointJaccard);
  npoints <- length(pointJaccard)
  label <- 1;
  for (lthr in thr)
  {
    totlabeled <- sum(classID > 0);
    if (totlabeled < npoints)
    {
      added <- 1;
      while (added > 0)
      {
        added <- 0;
        for (i in 1:npoints)
        {
          if (classID[i] == 0)
          {
            maxConnectedLabel <- label;
            wcon <- concensusMat[i,];
            consensA <- (wcon > lthr) & (classID > 0)
            consensB <- (wcon >= lthr) & (classID == 0)
            SconA <- sum(pointJaccard[consensA]*wcon[consensA]);
            SconB <- sum(pointJaccard[consensB]*wcon[consensB]) - pointJaccard[i]*wcon[i];
            #              cat("A:",SconA,"B:",SconB,"P:",pointJaccard[i],"\n")
            
            if ( (SconB > 0.0075*npoints) || (SconA > 0.01*npoints) )
            {
              if (SconB > 0.75*SconA)
              {
                classID[consensB] <- label;
                added <- 1;
                label <- label + 1;
              }
              else
              {
                if (SconB >= 0)
                {
                  if (SconA > 0)
                  {
                    tb <- table(classID[consensA])
                    smp <- tb;
                    for (nt in names(tb))
                    {
                      ptss <- consensA & (classID == as.numeric(nt));
                      smp[nt] <- sum(pointJaccard[ptss]*wcon[ptss]);
                    }
                    maxConnectedLabel <- names(which.max(smp))[1]
                    if ((0.75*smp[maxConnectedLabel]) <= SconB )
                    {
                      maxConnectedLabel <- label;
                    }
                    maxConnectedLabel <- as.numeric(maxConnectedLabel);
                    #                      print(c(label,smp,SconB,maxConnectedLabel))
                  }
                  classID[consensB] <- maxConnectedLabel;
                  added <- 1;
                  if (maxConnectedLabel == label)
                  {
                    label <- label + 1;
                  }
                }
              }
            }
          }
        }
      }
      totlabeled <- sum(classID > 0);
      cat(maxConnectedLabel,":",sprintf("(%5.3f)",lthr),": ",totlabeled,": ")
    }
  }
  classID <- classID[orgnames];
  tb <- table(classID)
  tb <- tb[order(-tb)];
  classID <- match(as.character(classID),names(tb))
  names(classID) <- orgnames;
  
  attr(classID,"who") <- who;
  
  return (classID);
}

relaxConsensusCluster <- function(object,classID,ww=19,loops=21)
{
  if (attr(classID,"who") != "training")
  {
    concensusMat <- object$testConsesus;
  }
  else
  {
    concensusMat <- object$dataConcensus;
  }
  cat(sum(names(classID) != rownames(concensusMat) ),"\n")
  
  orgnames <-  names(classID);
  cm <- as.data.frame(cor(concensusMat,method="pearson"));
  orderlist <- as.data.frame(matrix(0,nrow=ww,ncol=nrow(cm)));
  colnames(cm) <- orgnames;
  rownames(cm) <- orgnames;
  relaxclass <- classID;
  
  lp <- 0;
  changes <- 1;
  wws <- 3;
  while ((lp < loops) && (changes>0))
  {
    for (i in orgnames)
    {
      if (lp == 0)
      {
        orderlist[,i] <- order(-cm[,i])[1:ww];
      }
      tb <- table(classID[orderlist[,i]][1:wws]);
      tb <- tb[order(-tb)];
      if (tb[1] > 1)
      {
        relaxclass[i] <- as.integer(names(tb)[1]);
      }
    }
    lp <- lp + 1;
    changes <- sum(relaxclass != classID);
    cat(changes,":")
    classID <- relaxclass;
    if (wws < ww)
    {
      wws <- wws + 2;
      changes <- 1;
    }
  }
  cat("\n")
  return (relaxclass);
}