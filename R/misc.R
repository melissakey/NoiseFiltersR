weightedC45 <- function(data,weights,formu){
  N <- 1/min(weights[which(weights>1e-5)])
  dataTrain <- data[rep(seq_len(nrow(data)),times=round(weights*N)),]
  return(RWeka::J48(formu,dataTrain))
}

# network functions
distt <- function(x,y){
  class <- sapply(x,class)
  if("factor"%in%class){
    out <- sum((x[1,class!="factor"]-y[1,class!="factor"])^2)+sum(x[1,class=="factor"]!=y[1,class=="factor"])
  }
  else{
    out <- sum((x-y)^2)
  }
  out <- sqrt(out)
  return(out)
}

isNeighbor <- function(distMat, fun) {
  mat <- matrix(NA, nrow = nrow(distMat), ncol = ncol(distMat))
  sapply(1:(nrow(distMat)-1), function(i){
    sapply((i + 1):nrow(distMat), function(j) {
      mat[i,j] <<- mat[j,i] <<- all(distMat[i,j] <= fun(distMat[-c(i,j),i],distMat[-c(i,j),j]) + 1e-10)
    })
  })
  mat
  # mat <- mat + t(mat)
}


# DROP functions
distEnemy <- function(i,data,classColumn){
  isEnemy <- data[,classColumn]!=data[i,classColumn]
  dist <- kknn::kknn(as.formula(paste(names(data)[classColumn],"~.",sep = "")), train = data[isEnemy,], test = data[i,], k = 1)$D
  return(dist)
}




### PRISM functions
cld <- function(i,data,classColumn){
  classes <- unique(data[,classColumn])
  probs <- sapply(classes,function(cl){
    thisClass <- data[,classColumn]==cl
    sapply(setdiff(1:ncol(data),classColumn),function(att){
      if(is.factor(data[,att]) | is.logical(data[,att])){
        sum(data[thisClass,att]==data[i,att])/sum(thisClass)
      }
      else{
        m <- mean(data[thisClass,att])
        d <- stats::sd(data[thisClass,att])
        stats::dnorm(data[i,att], mean = m, sd = d)
      }
    })
  })
  probsPerClass <- apply(probs,2,prod)
  classIdx <- which(classes==data[i,classColumn])
  probsPerClass[classIdx]-max(probsPerClass[-classIdx])
}

dn <- function(i,data,classColumn){
  class <- as.character(data[i,classColumn])
  form <- as.formula(paste(names(data)[classColumn],"~.",sep=""))
  nn <- kknn::kknn(formula = form,
    train = data[-i,],
    test = data[i,],
    k = 17,
    kernel = "rectangular")$CL
  dns <- sapply(1:17,function(i){
    sum(nn[1:i]!=class)/i
  })
  mean(dns)
}

dcp <- function(i,data,classColumn){
  form <- as.formula(paste(names(data)[classColumn],"~.",sep=""))
  tree <- RWeka::J48(form,data)
  probs <- predict(tree,data[i,],type="probability")
  probs[1,colnames(probs)==data[i,classColumn]]
}

ds <- function(i,data,classColumn){
  equalities <- apply(data,1,function(x){
    all(x==data[i,])
  })
  sum(equalities)-1
}

# INFFC functions
FusionClassifiers <- function(data, trainingIndexes, majThreshold, returnNoisy=FALSE){
  predC45 <- predict(RWeka::J48(formula = class~., data = data[trainingIndexes,]),data)
  pred3NN <- sapply(1:nrow(data),function(i){kknn::kknn(class~.,train=data[setdiff(trainingIndexes,i),],test=data[i,], k=3)$fitted.values})
  invisible(utils::capture.output(predLOG <- predict(nnet::multinom(class~., data[trainingIndexes,]), data))) #To avoid some messages getting displayed in the console
  votes <- (predC45 != data$class) + (pred3NN != data$class) + (predLOG != data$class)
  if(returnNoisy){
    return(which(votes >= majThreshold))
  }else{
    return(which(votes < majThreshold))
  }
}

NoiseScore <- function(data,NoisyIndexes,k,indexToScore){
  neighborsIndexes <- kknn::kknn(class~., train = data[-indexToScore,], test = data[indexToScore, ], k = k)$C
  sum(sapply(neighborsIndexes,function(i){Confidence(data,NoisyIndexes,k,i)*Clean(data,NoisyIndexes,k,i)*ifelse(data[i,]$class==data[indexToScore,]$class,-1,1)}))/k
}

Confidence <- function(data,NoisyIndexes,k,index){
  t <- sum(sapply(NoisyIndexes,function(i){index %in% kknn::kknn(class~., train = data[-i,], test = data[i,], k = k)$C}))
  1/sqrt(1+t^2)
}

Clean <- function(data,NoisyIndexes,k,index){
  neighborsIndexes <- kknn::kknn(class~., train = data[-index,], test = data[index,], k = k)$C
  n <- sum(sapply(neighborsIndexes,function(i){i %in% NoisyIndexes}))
  (k+ifelse(index %in% NoisyIndexes,1,-1)*(n-k))/(2*k)
}
