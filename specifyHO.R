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
#  .print_to_console: determines whether the model is printed to the console or not. 

specifyHO <- function(.model = NULL,
                      .typeHO = c('normal','refined'), 
                      .order_indicators=c('exact','random'),
                      .print_to_console=FALSE){
  
  
  # Load required packages
  require(combinat)
  
  # match arguments
  .typeHO <- match.arg(.typeHO)
  .order_indicators <- match.arg(.order_indicators)

  eachline <- strsplit(x=.model,split='\\n')[[1]]
  positionEmergent <- grep(x=eachline,pattern = '<~')
  
  # loop over all lines that contain a composite
  for(i in positionEmergent){
    
    # split into composite and indicators
    nameEmergentAndIndicators <- strsplit(x=eachline[i],split='<~')[[1]]
    # remove potential spaces
    nameEmergentAndIndicators <- gsub(" ", "", nameEmergentAndIndicators, fixed = TRUE)
    
    nameEmergent <- nameEmergentAndIndicators[1]
    nameIndicators <- strsplit(x=nameEmergentAndIndicators[2],split='+',fixed=TRUE)[[1]]
    
    # if(length(nameIndicators)==1){
    #   stop('This function does not work for single-indicator composites')
    # }
    # 
    
    if(.order_indicators == 'random'){#else the indicator names are used in order as provided by the user
      allindicatorCombinations <- combinat::permn(nameIndicators)
      nameIndicators <- allindicatorCombinations[[sample(x = 1:length(allindicatorCombinations), size = 1)]]
    }
    
    # Specify equation for the emergent variable
    tempEmer <- paste0(nameEmergent,'=~',paste0(nameIndicators,collapse = '+'),'+',paste0('start(0)','*',nameIndicators[-1],collapse = '+'))
    
    # Fix all random measurement errors to zero
    tempErrcov <- paste0(nameIndicators,'~~0*',nameIndicators,collapse='\n')
    
    if(length(nameIndicators)>1){
      
      excrescentNames <- paste0(nameEmergent,'nu',1:(length(nameIndicators)-1))
      
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
      
      # Specify covariances among the excrescent variables if refined H-O is used. 
      tempExcrCov <- list()
      if(.typeHO == 'refined'){
        if(length(excrescentNames)>1){
          for(j in 1:(length(excrescentNames)-1)){
            tempExcrCov[[j]] <- paste0(excrescentNames[j],'~~', paste0(excrescentNames[-(1:j)],collapse = '+'))
          }
        }
      }
      
      # Put everything together
      eachline[i] <- paste0(tempEmer,'\n',paste0(unlist(tempExcr),collapse='\n'),'\n',paste0(unlist(tempExcrCov),collapse='\n'),'\n',tempErrcov)
      
    }else{ #single-indicator composite
      eachline[i] <- paste0(tempEmer,'\n',tempErrcov)
    }
    
  }
  
  # Put all lines together
  out <- paste0(eachline,collapse = '\n')

  if(.print_to_console == TRUE){
    cat(paste0(eachline,collapse = '\n'))
  }
  
  out
}
 