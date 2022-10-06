.model='
 PSS <~PSS1+PSS2+ PSS3+PSS4+PSS5
GRS=~GRS1+GRS2+GRS3+GRS4+GRS5
SRS<~SRS1+SRS2+SRS3+SRS4+SRS5

PSS~~GRS+SRS
GRS~~SRS
'

# lav_model=modelcfa


# Function should be used with orthogonal = TRUE
HOspecification <- function(.model = NULL,
                         .type = c('normal','refined'), 
                         .model_specification=c('exact','random')){

  
  # Load packages
  require(combinat)
  
  # match arguments
  .type <- match.arg(.type)
  .model_specification <- match.arg(.model_specification)
  
  
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

  if(length(nameIndicators)==1){
    stop('This function does not work for single-indicator composites')
  }
  
  
  allindicatorCombinations <- combinat::permn(nameIndicators)
  
  if(.model_specification == 'random'){#else the indicator names are used in order as provided by the user
    nameIndicators <- allindicatorCombinations[[sample(x = 1:length(allindicatorCombinations), size = 1)]]
  }

  # Specify equation for the emergent variable
  tempEmer <- paste0(nameEmergent,'=~',paste0(nameIndicators,collapse = '+'),'+',paste0('start(0)','*',nameIndicators[-1],collapse = '+'))

  # Fix all random measurement errors to zero
  tempErrcov <- paste0(nameIndicators,'~~0*',nameIndicators,collapse='\n')
  
  excrescentNames <- paste0(nameEmergent,'nu',1:(length(nameIndicators)-1))
  
  # Loading pattern of the excrescent variables
  Loadingmatrix <- as.data.frame(matrix(0,ncol=length(nameIndicators),nrow=length(excrescentNames),dimnames=list(excrescentNames,nameIndicators)))  
  if(.type == 'refined'){
    # fill loading matrix
    for(j in 1:nrow(Loadingmatrix)){
      Loadingmatrix[j,] <- c(rep(0,times=(j-1)),'NA',1,rep(0,ncol(Loadingmatrix)-(j+1))) 
    }
  }
  if(.type == 'normal'){
    for(j in 1:nrow(Loadingmatrix)){
      Loadingmatrix[j,] <- c(rep(0,times=(j-1)),'NA',1,rep('NA',ncol(Loadingmatrix)-(j+1))) 
    } 
  }
  
  # specify excrescent variables
  tempExcr <- list()
  for(j in excrescentNames){
    tempExcr[[j]] <- paste0(paste0(j,'=~', paste0(paste0(Loadingmatrix[j,],'*'),nameIndicators,collapse = '+')),'+',
                            paste0(if(.type == 'refined'){
                              'start(-1)'
                            } else if(.type == 'normal'){
                              'start(0)' 
                            },'*',colnames(Loadingmatrix[j,])[which(Loadingmatrix[j,]=='NA')],collapse = '+'))
    
  }
  
  # Specify covariances among the excrescent variables if refined H-O is used. 
  tempExcrCov=list()
  if(.type == 'refined'){
    if(length(excrescentNames)>1){
      for(j in 1:(length(excrescentNames)-1)){
        tempExcrCov[[j]]=paste0(excrescentNames[j],'~~', paste0(excrescentNames[-(1:j)],collapse = '+'))
      }
    }
  }
  
  # Put everything together
  eachline[i] = paste0(tempEmer,'\n',paste0(unlist(tempExcr),collapse='\n'),'\n',paste0(unlist(tempExcrCov),collapse='\n'),'\n',tempErrcov)

}

mod=paste0(eachline,collapse = '\n')

return(mod)
}
 
model=HOspecification(.model,.type='normal',.model_specification = 'random')

out= sem(model = model,data=dataset,orthogonal=T)
summary(out)
