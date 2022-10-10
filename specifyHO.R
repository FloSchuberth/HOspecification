# packages required to use this function
install.packages('combinat')
# This function converts a lavaan model containing composites specified via the '<~' command
# into a lavaan model in which composites are specified by means of the Henseler-Ogasawara (H-O) specification.
# Note, if you want to apply the sem or cfa function of the lavaan package to the model, 
# make sure that orthogonal is set to TRUE; otherwise, not all necessary constraints are imposed. 
# 
# The function provides the following argument:
# .model: lavaan model; this model should not include any additional specifications such as *value, *NA or comments (#). 
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

specifyHO <- function(.model = NULL,
                      .typeHO = c('normal','refined'), 
                      .order_indicators=c('exact','random'),
                      .determine_weights=TRUE,
                      .print_to_console=FALSE){
  
  
  # Load required packages ----
  require(combinat);require(calculus)
  
  # match arguments ----
  .typeHO <- match.arg(.typeHO)
  .order_indicators <- match.arg(.order_indicators)

  eachline <- strsplit(x=.model,split='\\n')[[1]]
  positionEmergent <- grep(x=eachline,pattern = '<~')
  
  # loop over all lines that contain a composite
  for(Line in positionEmergent){
    
    # split into composite and indicators
    nameEmergentAndIndicators <- strsplit(x=eachline[Line],split='<~')[[1]]
    
    # remove potential spaces
    nameEmergentAndIndicators <- gsub(" ", "", nameEmergentAndIndicators, fixed = TRUE)
    
    nameEmergent <- nameEmergentAndIndicators[1]
    nameIndicators <- strsplit(x=nameEmergentAndIndicators[2],split='+',fixed=TRUE)[[1]]

    if(.order_indicators == 'random'){#else the indicator names are used in order as provided by the user
      allindicatorCombinations <- combinat::permn(nameIndicators)
      nameIndicators <- allindicatorCombinations[[sample(x = 1:length(allindicatorCombinations), size = 1)]]
    }
    
    # Specify equation for the emergent variable
    tempEmer <- paste0(nameEmergent,'=~',paste0(nameIndicators,collapse = '+'),if(length(nameIndicators)>1){paste0('+',paste0('start(0)','*',nameIndicators[-1],collapse = '+'))})

    # Fix all random measurement errors to zero
    tempErrcov <- paste0(nameIndicators,'~~0*',nameIndicators,collapse='\n')
    
    if(length(nameIndicators)>1){
      
      # excrescentNames <- paste0(nameEmergent,'e',1:(length(nameIndicators)-1))
      excrescentNames <- paste0('e',Line,1:(length(nameIndicators)-1))
      
      # Loading pattern of the excrescent variables
      Loadingmatrix <- as.data.frame(matrix(0,ncol=length(nameIndicators),nrow=length(excrescentNames),dimnames=list(excrescentNames,nameIndicators)))  
      if(.typeHO == 'refined'){
        # fill loading matrix
        for(j in 1:nrow(Loadingmatrix)){
          Loadingmatrix[j,] <- c(rep(0,times=(j-1)),'NA',1,rep(0,ncol(Loadingmatrix)-(j+1))) 
        }
      }
      if(.typeHO == 'normal'){
        for(j in 1:nrow(Loadingmatrix)){
          Loadingmatrix[j,] <- c(rep(0,times=(j-1)),'NA',1,rep('NA',ncol(Loadingmatrix)-(j+1))) 
        } 
      }
      

      # specify excrescent variables
      tempExcr <- list()
      for(j in excrescentNames){
        tempExcr[[j]] <- paste0(paste0(j,'=~', paste0(paste0(Loadingmatrix[j,],'*'),nameIndicators,collapse = '+')),'+',
                                paste0(if(.typeHO == 'refined'){
                                  'start(-1)'
                                } else if(.typeHO == 'normal'){
                                  'start(0)' 
                                },'*',colnames(Loadingmatrix[j,])[which(Loadingmatrix[j,]=='NA')],collapse = '+'))
        
      }
      
      # Specify covariances among the excrescent variables
      tempExcrCov <- list()
      if(.typeHO == 'refined'){
        
        ExcrCov <- as.data.frame(matrix(0,ncol=length(excrescentNames),nrow=length(excrescentNames),dimnames=list(excrescentNames,excrescentNames))) 
        
        # Fill covariance matrix of the excrescent variables
        ExcrCov[upper.tri(ExcrCov)]<-1

        temp<-which(ExcrCov=='1',arr.ind = T)
        
        # This might create issue if the first row is empty
        for(j in unique(temp[,'row'])){
          tempExcrCov[[j]] <- paste0(rownames(ExcrCov)[j],'~~',paste0(rownames(ExcrCov)[temp[temp[,'row']==j,'col']],collapse = '+'))
        }
      }

      

# Add labels for the existing equations
      if(.determine_weights == TRUE){
        
        # Label loadings of the emergent variable
        labLoadEmer <- paste0('l',nameEmergent,2:length(nameIndicators))
        tempEmer <- paste0(tempEmer,'+', paste0(labLoadEmer,'*',nameIndicators[2:length(nameIndicators)],collapse = '+'))
        
        # Label loadings of the excrescent variables
        LoadingNamesmatrix<-Loadingmatrix 
        # LoadingNamesmatrix[Loadingmatrix==1] <- 0 
        temp<-which(Loadingmatrix=='NA',arr.ind = T)
        LoadingNamesmatrix[temp] <- paste0('l',nameEmergent,apply(temp,1,paste0,collapse=''))
        
        # ExcrLabtemp=list()
        for(j in excrescentNames){
          # add parameter label
          tempExcr[[j]] <- paste0(tempExcr[[j]],'+',paste0(LoadingNamesmatrix[j,][which(LoadingNamesmatrix[j,]!=0&LoadingNamesmatrix[j,]!=1)],
                                                           '*',colnames(LoadingNamesmatrix[j,])[which(LoadingNamesmatrix[j,]!=0&LoadingNamesmatrix[j,]!=1)],collapse='+'))
        }

        
        # Label the covariances among the excrescent variables
        if(.typeHO=='refined'){        
        # Specify covariances among the excrescent variables if refined H-O is used. 
        temp <- which(ExcrCov=='1',arr.ind = T)

        for(i in unique(temp[,'row'])){
          
          ExcrCovLabtemp <- list()
          for(ii in temp[temp[,'row']==i,'col']){
            ExcrCovLabtemp[[ii]] <- paste0('r',rownames(ExcrCov)[i],colnames(ExcrCov)[ii],'*',colnames(ExcrCov)[ii])
          }
      # Remove NULL from list 
          idx <- sapply(ExcrCovLabtemp,is.null)
          ExcrCovLabtemp <- ExcrCovLabtemp[idx==FALSE]

          # Extend the original syntax by the labels
          tempExcrCov[[i]] <- paste0(tempExcrCov[[i]],'+',paste0(ExcrCovLabtemp,collapse='+'))
        }
        }
        
        # calculate weights
        # add loadings of emergent variable
        mL <- rbind(c(1,labLoadEmer),
                      as.matrix(LoadingNamesmatrix))
        
        outW <- mxinv(t(mL))

        Wspec <- paste0('w',nameIndicators,':=',outW[1,],collapse='\n' )
        
        
        # label variance of the emergent variable
        varEmer <- paste0(nameEmergent,'~~',paste0('v',nameEmergent),'*',nameEmergent)
        
        # Label variances of the excrescent variables
        varExcr <- paste0(excrescentNames, '~~', paste0('v',excrescentNames),'*',excrescentNames,collapse='\n')
        
        # determine variances of the indicators
        vcvemerexcr <- matrix(0,nrow=ncol(mL),ncol=ncol(mL))
        if(.typeHO=='normal'){
          diag(vcvemerexcr) <- c(paste0('v',nameEmergent),paste0('v',excrescentNames))
        }
        
        if(.typeHO=='refined'){
          # determine the vcv of the emergent and excrescent variables
          vcvemerexcr <- as.matrix(Matrix::bdiag(1,as.matrix(ExcrCov)))
          dimnames(vcvemerexcr) <- list(c(nameEmergent,excrescentNames),c(nameEmergent,excrescentNames))
          
          ExcrCovMat <- ExcrCov
          temp <- which(ExcrCovMat==1,arr.ind=T)
          ExcrCovMat[temp] <- paste0('r',rownames(ExcrCov)[temp[,'row']],colnames(ExcrCov)[temp[,'col']])
          # make symmetric
          ExcrCovMat[lower.tri(ExcrCovMat)]<-t(ExcrCovMat)[lower.tri(ExcrCovMat)]
          diag(ExcrCovMat) <- paste0('v',excrescentNames)
          vcvemerexcr[excrescentNames,excrescentNames] <- as.matrix(ExcrCovMat)
          vcvemerexcr[nameEmergent,nameEmergent]  <- paste0('v',nameEmergent)
          }
        vcvInd <- mx(mx(t(mL),vcvemerexcr),mL)
        varInd <- paste0('v',colnames(mL),':=',diag(vcvInd),collapse='\n')
        
        
        # calculate standardized weights
        SDInd<-paste0('sqrt(',paste0('v',colnames(mL)),')')
        SDIndMatTemp <- matrix(0,nrow=length(SDInd),ncol=length(SDInd))
        diag(SDIndMatTemp) <- SDInd
        SDEmerMatTemp <- matrix(0,nrow=length(SDInd),ncol=length(SDInd))
        diag(SDEmerMatTemp) <- paste0('1/sqrt(',paste0('v',nameEmergent),')')
        
        wstd <- mx(mx(paste0('w',nameIndicators),SDIndMatTemp),SDEmerMatTemp)
        
        wspecstd <- paste0('wstd',nameIndicators,':=', wstd,collapse='\n')
        
      }
      
      # Put everything together
      eachline[Line] <- paste0(tempEmer,'\n',paste0(unlist(tempExcr),collapse='\n'),'\n',paste0(unlist(tempExcrCov),collapse='\n'),'\n',tempErrcov,if(.determine_weights == TRUE){
        paste0('\n',Wspec,'\n',varEmer, '\n',varExcr,'\n',varInd, '\n',wspecstd)})
      
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
 