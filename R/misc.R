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


# saturation filter functions

twoClassSaturationFilter <- function(data, noiseThreshold = NULL){
  initialInstances <- nrow(data)
  yes <- levels(data$class)[1]
  no <- levels(data$class)[2]
  
  if(is.null(noiseThreshold)){
    noiseThreshold <- getThreshold(min( sum(data$class==yes),sum(data$class==no) ))
  }
  
  cove <- coverage3array(data)
  
  indexesToRemove <- integer(0)
  KeepOn <- TRUE
  count <- 1
  while(KeepOn){
    #cat("Round:",count,"\n")
    if(dim(cove)[1]==0 || dim(cove)[2]==0){
      stop("a whole class was filtered!")
    }
    weightsPos <- numeric(dim(cove)[1])
    weightsNeg <- numeric(dim(cove)[2])
    literals <- optimLiterals(cove)
    for(k in literals){
      positive <- integer(0)
      negative <- integer(0)
      indexes <- arrayInd(which(cove[,,k] & apply(cove[,,literals],c(1,2),sum)==1),dim(cove)[c(1,2)])
      weightsPos[unique(indexes[,1])] <- weightsPos[unique(indexes[,1])]+1/length(unique(indexes[,1]))
      weightsNeg[unique(indexes[,2])] <- weightsNeg[unique(indexes[,2])]+1/length(unique(indexes[,2]))
    }
    if(max(c(weightsNeg,weightsPos))>noiseThreshold){
      count <- count+1
      if(max(weightsNeg)>max(weightsPos)||(max(weightsNeg)==max(weightsPos) && dim(cove)[1]<dim(cove)[2])){
        indexesToRemove <- c(indexesToRemove,which(data$class[setdiff(1:initialInstances,indexesToRemove)]==no)[which.max(weightsNeg)])
        cove <- cove[,- which.max(weightsNeg),]
      }
      else{
        indexesToRemove <- c(indexesToRemove,which(data$class[setdiff(1:initialInstances,indexesToRemove)]==yes)[which.max(weightsPos)])
        cove <- cove[- which.max(weightsPos),,]
      }
    }else{
      KeepOn <- FALSE
    }
  }
  indexesToRemove
}

coverage3array <- function(data){
  yes <- levels(data$class)[1]
  no <- levels(data$class)[2]
  posIndexes <- which(data$class==yes)
  negIndexes <- which(data$class==no)
  numLiterals <- 0
  out <- logical(0)
  for(col in 1:ncol(data)){
    if(names(data)[col]=="class"){
      next
    }
    if(class(data[,col])=="factor" || class(data[,col])=="integer"){
      #literals created from 'positive' instances
      literals <- unique(data[posIndexes,col])
      numLiterals <- numLiterals+length(literals)
      out <- c(out,sapply(literals,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]==l && data[iNeg,col]!=l})})}))
      #literals created from 'negative' instances
      literals <- unique(data[negIndexes,col])
      numLiterals <- numLiterals+length(literals)
      out <- c(out,sapply(literals,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]!=l && data[iNeg,col]==l})})}))
    }
    if(class(data[,col])=="numeric" || class(data[,col])=="integer"){
      literalsLQ <- numeric(0)
      literalsG <- numeric(0)
      v <- unique(data[posIndexes,col])
      names(v) <- rep(yes,length(v))
      w <- unique(data[negIndexes,col])
      names(w) <- rep(no,length(w))
      values <- ordering(c(v,w))
      currentClass <- names(values)[1]
      for(i in 1:(length(values)-1)){
        if(names(values)[i+1]!=currentClass){
          if(currentClass==yes){
            literalsLQ <- c(literalsLQ,mean(values[c(i,i+1)]))
          }
          else{
            literalsG <- c(literalsG,mean(values[c(i,i+1)]))
          }
          currentClass <- names(values)[i+1]
        }
      }
      #                   if(data[1,"class"]==yes){
      #                         cat("Column",col,":\n")
      #                         cat("Values:\n")
      #                         print(values)
      #                         cat("LiteralsLQ:\n",literalsLQ,"\n")
      #                         cat("LiteralsG:\n",literalsG,"\n")
      #                   }
      numLiterals <- numLiterals+length(literalsLQ)+length(literalsG)
      out <- c(out,sapply(literalsLQ,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]<=l && data[iNeg,col]>=l && data[iPos,col]!=data[iNeg,col]})})}))
      out <- c(out,sapply(literalsG,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]>=l && data[iNeg,col]<=l && data[iPos,col]!=data[iNeg,col]})})}))
    }
  }
  if(class(out)=="list"){
    out <- unlist(out)
  }
  dim(out) <- c(length(posIndexes),length(negIndexes),numLiterals)
  out
}

optimLiterals <- function(cove){
  coveCopy <- cove
  out <- integer(0)
  subround <- 1
  while(!all(is.na(cove))){
    #cat("Subround:",subround,"\n")
    #print(str(apply(cove,c(1,2),sum)))
    #print(all(is.na(apply(cove,c(1,2),sum))))
    #summing <- apply(cove,c(1,2),sum)
    #print(str(summing[which(!is.na(summing))]))
    weights <- 1/apply(cove,c(1,2),sum)
    #PRINT CLASS WEIGHTS TO MAKE SURE IT IS A MATRIX
    #print(str(weights))
    #print(all(is.na(weights)))
    #print(summary(as.vector(weights)))
    indMax <- arrayInd(which.max(weights),dim(cove)[c(1,2)])
    #print(weights[indMax])
    litToTest <- which(cove[indMax[1],indMax[2],])
    #print(length(litToTest))
    weightsLit <- numeric(length(litToTest))
    for(k in seq_along(litToTest)){
      arrIndexes <- arrayInd(which(!is.na(cove[,,litToTest[k]]) & cove[,,litToTest[k]]),dim(cove)[c(1,2)])
      #print(summary(apply(arrIndexes,2,function(ind){weights[ind]})))
      #print(class(arrIndexes))
      #apply(arrIndexes,2,function(ind){weights[ind]})
      weightsLit[k] <- sum(apply(arrIndexes,1,function(ind){weights[ind[1],ind[2]]}))
    }
    #print(weightsLit)
    litToKeep <- litToTest[which.max(weightsLit)]
    #print(weightsLit[which.max(weightsLit)])
    #print(str(weightsLit))
    out <- c(out,litToKeep)
    covered <- arrayInd(which(!is.na(cove[,,litToKeep]) & cove[,,litToKeep]),dim(cove)[c(1,2)])
    #print(length(which(!is.na(cove[,,litToKeep]) & cove[,,litToKeep])))
    for(i in 1:nrow(covered)){
      cove[covered[i,1],covered[i,2],] <- NA
    }
    cove <- cove[,,-litToKeep]
    #print(any(apply(cove,c(1,2),function(v){any(is.na(v)) && !all(is.na(v))})))
  }
  for(k in out){
    if(!any(coveCopy[,,k] & apply(coveCopy[,,out],c(1,2),sum)==1)){
      out <- setdiff(out,k)
    }
  }
  out
}

getThreshold <- function(x){
  if(x<=50){
    out <- 1.5
  }else if(50<x && x<=100){
    out <- 1
  }else if(100<x && x<=200){
    out <- 0.75
  }else{
    out <- 0.5
  }
  out
}

ordering <- function(values){
  out <- values[order(values)]
  for(i in 2:(length(out)-1)){
    if(out[i]==out[i+1] && names(out)[i+1]==names(out)[i-1]){
      names(out)[i+1] <- names(out)[i]
      names(out)[i] <- names(out)[i-1]
    }
  }
  #       if(out[1]==out[2] && names(out)[1]==names(out)[3]){
  #             names(out)[1] <- names(out)[2]
  #             names(out)[2] <- names(out)[3]
  #       }
  out
}

