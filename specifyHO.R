# packages required to use this function
# install.packages(c('combinat','calculus'))


# This function converts a lavaan model containing composites specified via the '<~' command
# into a lavaan model in which composites are specified by means of the Henseler-Ogasawara (H-O) specification.
# Note, if you want to apply the sem or cfa function of the lavaan package to the model, 
# make sure that orthogonal is set to TRUE; otherwise, not all necessary constraints are imposed. 
# 
# The function provides the following argument:
# .model: lavaan model; this model should not include any additional specifications such as *value, *NA or comments (#). 
# Moreover, composites need to be specified in one line, i.e., it is currently not possible to specify a composite c1 that is formed 
# by x1 and x2 as c1<~ x1 and c1<~ x2.
# # .typeHO: the type of H-O specification:
#   -'refined': Refined H-O specification
#   -'normal': Normal H-O specification as presented in Schuberth (in press)
# 
# .order_indicators: determines the order of the indicators:
#   -'exact': use the original order of the lavaan model provided by the user
#   -'random': randomly orders the indicators; this allows to try different specifications, e.g., in case of convergence issues. 
# 
# .determine_weights: determines whether the weights should be calculated.
# 
# 
#  .print_to_console: determines whether the model is printed to the console. 
# 
#  .seed: can be used to set a seed to fix the order of the indicator names in case of .order_indicators = 'random'.

specifyHO <- function(.model = NULL,
                      .typeHO = c('normal','refined'), 
                      .order_indicators=c('exact','random'),
                      .determine_weights=TRUE,
                      .print_to_console=FALSE,
                      .seed=NULL){
  
  
  # Load required packages ----
  require(combinat);require(calculus)
  
  # match arguments ----
  .typeHO <- match.arg(.typeHO)
  .order_indicators <- match.arg(.order_indicators)
  
  
  eachline <- strsplit(x=.model,split='\\n')[[1]]
  positionEmergent <- grep(x=eachline,pattern = '<~')
  
  
  # Check whether a composite has been specified over various lines  
  temp=lapply(eachline,function(x){
    temp <- strsplit(x,split="<~")[[1]]
  })
  
  temp1=sapply(temp,function(x){
    if(length(x)==2){
      x[1]
    }else{
      NA
    }
  })
  
  temp1=temp1[!is.na(temp1)]
  
  if(sum(duplicated(temp1))!=0){
    stop("Please specify your composites in a single line.")
  }
  
  # loop over all lines that contain a composite
  for(Line in positionEmergent){
    
    # split into composite and indicators
    nameEmergentAndIndicators <- strsplit(x=eachline[Line],split='<~')[[1]]
    
    # remove potential spaces
    nameEmergentAndIndicators <- gsub(" ", "", nameEmergentAndIndicators, fixed = TRUE)
    
    nameEmergent <- nameEmergentAndIndicators[1]
    
    # label variance of the emergent variable
    varEmer <- paste0(nameEmergent,'~~',paste0('v',nameEmergent),'*',nameEmergent)
    
    
    nameIndicatorstemp <- strsplit(x=nameEmergentAndIndicators[2],split='+',fixed=TRUE)[[1]]
    nameIndicatorstemp1 <- strsplit(x=nameIndicatorstemp,split='*',fixed=TRUE)
    
    nameIndicators <- sapply(nameIndicatorstemp1,function(x){
      if(length(x)==1){
        x
      }else if(length(x)==2){
        x[2]
      } else {
        stop("Something went wrong.")
      }
    })
    
    # Store fixed weight values
    WeightValues <- sapply(nameIndicatorstemp1,function(x){
      if(length(x)==1){
        NA
      }else if(length(x)==2){
        x[1]
      } else {
        stop("Something went wrong.")
      }
    })
    
    names(WeightValues) <- nameIndicators
    
    ThereArePresetWeights<- length(WeightValues[!is.na(WeightValues)])!=0
    
    # If there are fixed weights, use refined HO
    if(ThereArePresetWeights & .typeHO!="refined"){
      .typeHO = "refined"
      warning("Since weights are preset, `.typeHO` is set to `refined`.")
    }
    
    if(.order_indicators == 'random'){#else the indicator names are used in order as provided by the user
      # Leads to problems if there are too many indicators
      # allindicatorCombinations <- combinat::permn(nameIndicators)
      # nameIndicators <- allindicatorCombinations[[sample(x = 1:length(allindicatorCombinations), size = 1)]]
      if(!is.null(.seed)){
        set.seed(.seed)
      }
      nameIndicators <- nameIndicators[sample.int(length(nameIndicators))]
    }
    
    # reorder weight values
    WeightValues <- WeightValues[nameIndicators]
    
    # If there are no fixed weights
    if(!ThereArePresetWeights){
      # Specify equation for the emergent variable
      tempEmer <- paste0(nameEmergent,'=~',paste0(nameIndicators,collapse = '+'),
                         if(length(nameIndicators)>1){paste0('+',paste0('start(0)','*',nameIndicators[-1],collapse = '+'))})
      
      labLoadEmer <- c(1,paste0('l',nameEmergent,2:length(nameIndicators)))
      tempEmer <- paste0(tempEmer,'+', paste0(labLoadEmer,'*',nameIndicators,collapse = '+'))
      
    }else if(ThereArePresetWeights){ # if there are fixed weights use phantom variable approach
      
      namePhantom <- paste0('p',nameIndicators)
      names(namePhantom) <- nameIndicators
      
      tempEmer <- paste0(nameEmergent,'=~',paste0(namePhantom,collapse = '+'))
      tempEmer <- paste0(tempEmer,"+NA*",namePhantom[1])
      labLoadEmer <- paste0('l',nameEmergent,1:length(namePhantom))                   
      
      tempEmer <- paste0(tempEmer,'+', paste0(labLoadEmer,'*',namePhantom,collapse = '+'))                   
      
      # Add constraint on the sum of the emergent variable loadings
      ConstraintEmer <- paste0(paste0('l',nameEmergent,1:length(nameIndicators),collapse="+"),"==1")
      
      ErrcovPhan <-  paste0(namePhantom,'~~0*',namePhantom,collapse='\n')
      
      # Specify relations between indicators and phantom variables
      relPhantomInd = list()
      for(j in nameIndicators){
        relPhantomInd[j] <- paste0(namePhantom[j],"=~",
                                   if(!is.na(WeightValues[j])){paste0("(1/",WeightValues[j],")")
                                   }else{
                                     WeightValues[j] 
                                   },"*",j)
        
        # label the relationship between the phantom variables and the indicators
        relPhantomInd[j] = paste0(relPhantomInd[j],if(is.na(WeightValues[j])){paste0("+",paste0("l",j),"*",j)})
      }
      
      RelPhantom <- paste0(unlist(relPhantomInd),collapse = '\n')    
    }
    
    # Constrain measurement error variances to zero 
    tempErrcov <- paste0(nameIndicators,'~~0*',nameIndicators,collapse='\n')
    
    if(length(nameIndicators)>1){
      
      nameExcrescents <- paste0('e',Line,1:(length(nameIndicators)-1))
      
      # Label variances of the excrescent variables
      varExcr <- paste0(nameExcrescents, '~~', paste0('v',nameExcrescents),'*',nameExcrescents,collapse='\n')
      
      # Specify excrescent variables and their relations with the indicators/phantom variables
      # Loading pattern of the excrescent variables
      if(!ThereArePresetWeights){
        Loadingmatrix <- as.data.frame(matrix(0,ncol=length(nameIndicators),nrow=length(nameExcrescents),dimnames=list(nameExcrescents,nameIndicators)))  
      } else { #if fixed weights are used, use always refined HO
        Loadingmatrix <- as.data.frame(matrix(0,ncol=length(namePhantom),nrow=length(nameExcrescents),dimnames=list(nameExcrescents,namePhantom)))  
      }
      
      if(.typeHO == 'refined'){
        # fill loading matrix
        for(j in 1:nrow(Loadingmatrix)){
          if(!ThereArePresetWeights){
            Loadingmatrix[j,] <- c(rep(0,times=(j-1)),'NA',1,rep(0,ncol(Loadingmatrix)-(j+1))) 
          } else if(ThereArePresetWeights) { 
            Loadingmatrix[j,] <- c(rep(0,times=(j-1)),-1,1,rep(0,ncol(Loadingmatrix)-(j+1))) 
          }
        }
      }
      if(.typeHO == 'normal'){
        for(j in 1:nrow(Loadingmatrix)){
          Loadingmatrix[j,] <- c(rep(0,times=(j-1)),'NA',1,rep('NA',ncol(Loadingmatrix)-(j+1))) 
        } 
      }
      
      tempExcr <- list()
      for(j in nameExcrescents){
        temp = paste0(j,'=~', paste0(paste0(Loadingmatrix[j,],'*'),colnames(Loadingmatrix),collapse = '+'))
        if(!ThereArePresetWeights){          
          tempExcr[[j]] <- paste0(temp ,'+',
                                  paste0(if(.typeHO == 'refined'){
                                    'start(-1)'
                                  } else if(.typeHO == 'normal'){
                                    'start(0)' 
                                  },'*',colnames(Loadingmatrix[j,])[which(Loadingmatrix[j,]=='NA')],collapse = '+'))  
        } else {
          tempExcr[[j]] <- temp
        }
      }
      
      LoadingNamesmatrix<-Loadingmatrix 
      temp<-which(Loadingmatrix=='NA',arr.ind = T)
      LoadingNamesmatrix[temp] <- paste0('l',nameEmergent,apply(temp,1,paste0,collapse=''))
      
      for(j in nameExcrescents){
        # add parameter label
        if(sum(Loadingmatrix[j,]=="NA")>0){
          tempExcr[[j]] <- paste0(tempExcr[[j]],'+',paste0(LoadingNamesmatrix[j,][which(Loadingmatrix[j,]=="NA")],
                                                           '*',colnames(LoadingNamesmatrix[j,])[which(Loadingmatrix[j,]=="NA")],collapse='+'))
        }
      }
      
      
      # Specify covariances between the excrescent variables
      tempExcrCov <- list()
      if(.typeHO == 'refined'){
        ExcrCov <- as.data.frame(matrix(0,ncol=length(nameExcrescents),nrow=length(nameExcrescents),dimnames=list(nameExcrescents,nameExcrescents))) 
        
        # Fill covariance matrix of the excrescent variables to allow for covariances between the excrescent variables
        ExcrCov[upper.tri(ExcrCov)]<-1
        
        temp<-which(ExcrCov=='1',arr.ind = T)
        
        # This might create issue if the first row is empty
        for(j in unique(temp[,'row'])){
          tempExcrCov[[j]] <- paste0(rownames(ExcrCov)[j],'~~',paste0(rownames(ExcrCov)[temp[temp[,'row']==j,'col']],collapse = '+'))
          
          ExcrCovLabtemp <- list()
          for(i in temp[temp[,'row']==j,'col']){
            ExcrCovLabtemp[[i]] <- paste0('r',rownames(ExcrCov)[j],colnames(ExcrCov)[i],'*',colnames(ExcrCov)[i])
          }
          # Remove NULL from list 
          idx <- sapply(ExcrCovLabtemp,is.null)
          ExcrCovLabtemp <- ExcrCovLabtemp[idx==FALSE]
          
          # Extend the original syntax by the labels
          tempExcrCov[[j]] <- paste0(tempExcrCov[[j]],'+',paste0(ExcrCovLabtemp,collapse='+'))
          
        }
        
      }
      

      # Calculate weights
      if(.determine_weights == TRUE){
        
        # Create loading matrix of the emergent and excrescent variables
        mL <- rbind(labLoadEmer,
                    as.matrix(LoadingNamesmatrix))
        
        outW <- mxinv(t(mL))
        
        if(ThereArePresetWeights){
          weightValuestemp = WeightValues
          for(i in 1:length(WeightValues)){
            if(is.na(WeightValues[i])){
              weightValuestemp[i] <- paste0("l",nameIndicators[i])
            }else if(!is.na(WeightValues[i])){
              weightValuestemp[i] <- paste0("1/",WeightValues[i])
            }
          }
          
          # construct the loading matrix between the indicators and the phantom variables
          mLL <- matrix(0,nrow=length(nameIndicators),ncol=length(namePhantom),dimnames=list(nameIndicators,namePhantom))
          diag(mLL) <-weightValuestemp
          
          outW <- mx(outW,mxinv(mLL))
        }
        
        Wspec <- paste0('w',nameIndicators,':=',outW[1,],collapse='\n' )
        

        # determine variances of the indicators
        vcvemerexcr <- matrix(0,nrow=ncol(mL),ncol=ncol(mL))
        if(.typeHO=='normal'){
          diag(vcvemerexcr) <- c(paste0('v',nameEmergent),paste0('v',nameExcrescents))
        }
        
        if(.typeHO=='refined'){
          # determine the vcv of the emergent and excrescent variables
          vcvemerexcr <- as.matrix(Matrix::bdiag(1,as.matrix(ExcrCov)))
          dimnames(vcvemerexcr) <- list(c(nameEmergent,nameExcrescents),c(nameEmergent,nameExcrescents))
          
          ExcrCovMat <- ExcrCov
          temp <- which(ExcrCovMat==1,arr.ind=T)
          ExcrCovMat[temp] <- paste0('r',rownames(ExcrCov)[temp[,'row']],colnames(ExcrCov)[temp[,'col']])
          # make symmetric
          ExcrCovMat[lower.tri(ExcrCovMat)]<-t(ExcrCovMat)[lower.tri(ExcrCovMat)]
          diag(ExcrCovMat) <- paste0('v',nameExcrescents)
          vcvemerexcr[nameExcrescents,nameExcrescents] <- as.matrix(ExcrCovMat)
          vcvemerexcr[nameEmergent,nameEmergent]  <- paste0('v',nameEmergent)
        }
        
        
        if(!ThereArePresetWeights){
        vcvInd <- mx(mx(t(mL),vcvemerexcr),mL)
        } else if(ThereArePresetWeights){
        vcvInd <- mx(mx(mx(mx(mLL,t(mL)),vcvemerexcr),mL),mLL)
        }
        varInd <- paste0('v',nameIndicators,':=',diag(vcvInd),collapse='\n')
        
        
        # calculate standardized weights
        SDInd<-paste0('sqrt(',paste0('v',nameIndicators),')')
        SDIndMatTemp <- matrix(0,nrow=length(SDInd),ncol=length(SDInd))
        diag(SDIndMatTemp) <- SDInd
        SDEmerMatTemp <- matrix(0,nrow=length(SDInd),ncol=length(SDInd))
        diag(SDEmerMatTemp) <- paste0('1/sqrt(',paste0('v',nameEmergent),')')
        
        wstd <- mx(mx(paste0('w',nameIndicators),SDIndMatTemp),SDEmerMatTemp)
        
        wspecstd <- paste0('wstd',nameIndicators,':=', wstd,collapse='\n')
        
      }
      
      # Put everything together
      eachline[Line] <- paste0(tempEmer,'\n',
                               varEmer,'\n',
                               paste0(unlist(tempExcr),collapse='\n'),'\n',
                               varExcr,'\n',
                               paste0(unlist(tempExcrCov),collapse='\n'),'\n',
                               tempErrcov,
                               if(ThereArePresetWeights){
                                 paste0("\n",RelPhantom,"\n",ErrcovPhan,"\n",ConstraintEmer)
                               },
                               if(.determine_weights == TRUE){
                                 paste0('\n',Wspec,'\n',varInd, '\n',wspecstd)})
      
    }else{ #single-indicator composite
      eachline[Line] <- paste0(tempEmer,'\n',tempErrcov)
    }
    
  }
  
  
  # Put all lines together
  out <- paste0(eachline,collapse = '\n')
  
  if(.print_to_console == TRUE){
    cat(paste0(eachline,collapse = '\n'))
  }
  
  out
}
