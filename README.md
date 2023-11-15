# Function that facilitates the H-O specification

To load the specifyHO, you can use the following syntax:
``` r
devtools::source_url("https://github.com/FloSchuberth/HOspecification/blob/main/specifyHO.R?raw=TRUE")
```

In a next step, the model is specified in lavaan syntax using the `<~` operator to specify emergent variables. 
It is noteworthy that the model syntax for the emergent variables must not contain any parameter labels or constraints.  
An example can be found in the following:
``` r
model='
# Specify the emergent variables
PSS <~ PSS1 + PSS2 + PSS3 + PSS4 + PSS5
GRS <~ GRS1 + GRS2 + GRS3 + GRS4 + GRS5
SRS <~ SRS1 + SRS2 + SRS3 + SRS4 + SRS5

# Specify the relationships among the emergent variables
GRS ~~ PSS + SRS + Reader
PSS ~~ SRS + Reader
SRS ~~ Reader

'
``` 
To obtain the H-O specification of this model, the `specifyHO` function can be used. 
``` r
modelHO = specifyHO(.model=model,
                    .typeHO='refined',
                    .order_indicators = 'exact',
                    .print_to_console = T,
                    .determine_weights = T,
                    .seed = NULL)
modelHO
```
Here the resulting model syntax:
```
# Specify the emergent variables
PSS=~PSS1+PSS2+PSS3+PSS4+PSS5+start(0)*PSS2+start(0)*PSS3+start(0)*PSS4+start(0)*PSS5+lPSS2*PSS2+lPSS3*PSS3+lPSS4*PSS4+lPSS5*PSS5
e31=~NA*PSS1+1*PSS2+0*PSS3+0*PSS4+0*PSS5+start(-1)*PSS1+lPSS11*PSS1
e32=~0*PSS1+NA*PSS2+1*PSS3+0*PSS4+0*PSS5+start(-1)*PSS2+lPSS22*PSS2
e33=~0*PSS1+0*PSS2+NA*PSS3+1*PSS4+0*PSS5+start(-1)*PSS3+lPSS33*PSS3
e34=~0*PSS1+0*PSS2+0*PSS3+NA*PSS4+1*PSS5+start(-1)*PSS4+lPSS44*PSS4
e31~~e32+e33+e34+re31e32*e32+re31e33*e33+re31e34*e34
e32~~e33+e34+re32e33*e33+re32e34*e34
e33~~e34+re33e34*e34
PSS1~~0*PSS1
PSS2~~0*PSS2
PSS3~~0*PSS3
PSS4~~0*PSS4
PSS5~~0*PSS5
wPSS1:=((1)*((1)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lPSS2)*((lPSS11)*((1)*((1)*((1))))) + (lPSS3)*((lPSS11)*((lPSS22)*((1)*((1))))) + -(lPSS4)*((lPSS11)*((lPSS22)*((lPSS33)*((1))))) + (lPSS5)*((lPSS11)*((lPSS22)*((lPSS33)*((lPSS44))))))
wPSS2:=-((lPSS11)*((1)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lPSS2)*((lPSS11)*((1)*((1)*((1))))) + (lPSS3)*((lPSS11)*((lPSS22)*((1)*((1))))) + -(lPSS4)*((lPSS11)*((lPSS22)*((lPSS33)*((1))))) + (lPSS5)*((lPSS11)*((lPSS22)*((lPSS33)*((lPSS44))))))
wPSS3:=((lPSS11)*((lPSS22)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lPSS2)*((lPSS11)*((1)*((1)*((1))))) + (lPSS3)*((lPSS11)*((lPSS22)*((1)*((1))))) + -(lPSS4)*((lPSS11)*((lPSS22)*((lPSS33)*((1))))) + (lPSS5)*((lPSS11)*((lPSS22)*((lPSS33)*((lPSS44))))))
wPSS4:=-((lPSS11)*((lPSS22)*((lPSS33)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lPSS2)*((lPSS11)*((1)*((1)*((1))))) + (lPSS3)*((lPSS11)*((lPSS22)*((1)*((1))))) + -(lPSS4)*((lPSS11)*((lPSS22)*((lPSS33)*((1))))) + (lPSS5)*((lPSS11)*((lPSS22)*((lPSS33)*((lPSS44))))))
wPSS5:=((lPSS11)*((lPSS22)*((lPSS33)*((lPSS44))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lPSS2)*((lPSS11)*((1)*((1)*((1))))) + (lPSS3)*((lPSS11)*((lPSS22)*((1)*((1))))) + -(lPSS4)*((lPSS11)*((lPSS22)*((lPSS33)*((1))))) + (lPSS5)*((lPSS11)*((lPSS22)*((lPSS33)*((lPSS44))))))
PSS~~vPSS*PSS
e31~~ve31*e31
e32~~ve32*e32
e33~~ve33*e33
e34~~ve34*e34
vPSS1:=((1) * (vPSS)) * (1) + ((lPSS11) * (ve31)) * (lPSS11)
vPSS2:=((lPSS2) * (vPSS)) * (lPSS2) + ((1) * (ve31) + (lPSS22) * (re31e32)) * (1) + ((1) * (re31e32) + (lPSS22) * (ve32)) * (lPSS22)
vPSS3:=((lPSS3) * (vPSS)) * (lPSS3) + ((1) * (ve32) + (lPSS33) * (re32e33)) * (1) + ((1) * (re32e33) + (lPSS33) * (ve33)) * (lPSS33)
vPSS4:=((lPSS4) * (vPSS)) * (lPSS4) + ((1) * (ve33) + (lPSS44) * (re33e34)) * (1) + ((1) * (re33e34) + (lPSS44) * (ve34)) * (lPSS44)
vPSS5:=((lPSS5) * (vPSS)) * (lPSS5) + ((1) * (ve34)) * (1)
wstdPSS1:=((wPSS1) * (sqrt(vPSS1))) * (1/sqrt(vPSS))
wstdPSS2:=((wPSS2) * (sqrt(vPSS2))) * (1/sqrt(vPSS))
wstdPSS3:=((wPSS3) * (sqrt(vPSS3))) * (1/sqrt(vPSS))
wstdPSS4:=((wPSS4) * (sqrt(vPSS4))) * (1/sqrt(vPSS))
wstdPSS5:=((wPSS5) * (sqrt(vPSS5))) * (1/sqrt(vPSS))
GRS=~GRS1+GRS2+GRS3+GRS4+GRS5+start(0)*GRS2+start(0)*GRS3+start(0)*GRS4+start(0)*GRS5+lGRS2*GRS2+lGRS3*GRS3+lGRS4*GRS4+lGRS5*GRS5
e41=~NA*GRS1+1*GRS2+0*GRS3+0*GRS4+0*GRS5+start(-1)*GRS1+lGRS11*GRS1
e42=~0*GRS1+NA*GRS2+1*GRS3+0*GRS4+0*GRS5+start(-1)*GRS2+lGRS22*GRS2
e43=~0*GRS1+0*GRS2+NA*GRS3+1*GRS4+0*GRS5+start(-1)*GRS3+lGRS33*GRS3
e44=~0*GRS1+0*GRS2+0*GRS3+NA*GRS4+1*GRS5+start(-1)*GRS4+lGRS44*GRS4
e41~~e42+e43+e44+re41e42*e42+re41e43*e43+re41e44*e44
e42~~e43+e44+re42e43*e43+re42e44*e44
e43~~e44+re43e44*e44
GRS1~~0*GRS1
GRS2~~0*GRS2
GRS3~~0*GRS3
GRS4~~0*GRS4
GRS5~~0*GRS5
wGRS1:=((1)*((1)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lGRS2)*((lGRS11)*((1)*((1)*((1))))) + (lGRS3)*((lGRS11)*((lGRS22)*((1)*((1))))) + -(lGRS4)*((lGRS11)*((lGRS22)*((lGRS33)*((1))))) + (lGRS5)*((lGRS11)*((lGRS22)*((lGRS33)*((lGRS44))))))
wGRS2:=-((lGRS11)*((1)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lGRS2)*((lGRS11)*((1)*((1)*((1))))) + (lGRS3)*((lGRS11)*((lGRS22)*((1)*((1))))) + -(lGRS4)*((lGRS11)*((lGRS22)*((lGRS33)*((1))))) + (lGRS5)*((lGRS11)*((lGRS22)*((lGRS33)*((lGRS44))))))
wGRS3:=((lGRS11)*((lGRS22)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lGRS2)*((lGRS11)*((1)*((1)*((1))))) + (lGRS3)*((lGRS11)*((lGRS22)*((1)*((1))))) + -(lGRS4)*((lGRS11)*((lGRS22)*((lGRS33)*((1))))) + (lGRS5)*((lGRS11)*((lGRS22)*((lGRS33)*((lGRS44))))))
wGRS4:=-((lGRS11)*((lGRS22)*((lGRS33)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lGRS2)*((lGRS11)*((1)*((1)*((1))))) + (lGRS3)*((lGRS11)*((lGRS22)*((1)*((1))))) + -(lGRS4)*((lGRS11)*((lGRS22)*((lGRS33)*((1))))) + (lGRS5)*((lGRS11)*((lGRS22)*((lGRS33)*((lGRS44))))))
wGRS5:=((lGRS11)*((lGRS22)*((lGRS33)*((lGRS44))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lGRS2)*((lGRS11)*((1)*((1)*((1))))) + (lGRS3)*((lGRS11)*((lGRS22)*((1)*((1))))) + -(lGRS4)*((lGRS11)*((lGRS22)*((lGRS33)*((1))))) + (lGRS5)*((lGRS11)*((lGRS22)*((lGRS33)*((lGRS44))))))
GRS~~vGRS*GRS
e41~~ve41*e41
e42~~ve42*e42
e43~~ve43*e43
e44~~ve44*e44
vGRS1:=((1) * (vGRS)) * (1) + ((lGRS11) * (ve41)) * (lGRS11)
vGRS2:=((lGRS2) * (vGRS)) * (lGRS2) + ((1) * (ve41) + (lGRS22) * (re41e42)) * (1) + ((1) * (re41e42) + (lGRS22) * (ve42)) * (lGRS22)
vGRS3:=((lGRS3) * (vGRS)) * (lGRS3) + ((1) * (ve42) + (lGRS33) * (re42e43)) * (1) + ((1) * (re42e43) + (lGRS33) * (ve43)) * (lGRS33)
vGRS4:=((lGRS4) * (vGRS)) * (lGRS4) + ((1) * (ve43) + (lGRS44) * (re43e44)) * (1) + ((1) * (re43e44) + (lGRS44) * (ve44)) * (lGRS44)
vGRS5:=((lGRS5) * (vGRS)) * (lGRS5) + ((1) * (ve44)) * (1)
wstdGRS1:=((wGRS1) * (sqrt(vGRS1))) * (1/sqrt(vGRS))
wstdGRS2:=((wGRS2) * (sqrt(vGRS2))) * (1/sqrt(vGRS))
wstdGRS3:=((wGRS3) * (sqrt(vGRS3))) * (1/sqrt(vGRS))
wstdGRS4:=((wGRS4) * (sqrt(vGRS4))) * (1/sqrt(vGRS))
wstdGRS5:=((wGRS5) * (sqrt(vGRS5))) * (1/sqrt(vGRS))
SRS=~SRS1+SRS2+SRS3+SRS4+SRS5+start(0)*SRS2+start(0)*SRS3+start(0)*SRS4+start(0)*SRS5+lSRS2*SRS2+lSRS3*SRS3+lSRS4*SRS4+lSRS5*SRS5
e51=~NA*SRS1+1*SRS2+0*SRS3+0*SRS4+0*SRS5+start(-1)*SRS1+lSRS11*SRS1
e52=~0*SRS1+NA*SRS2+1*SRS3+0*SRS4+0*SRS5+start(-1)*SRS2+lSRS22*SRS2
e53=~0*SRS1+0*SRS2+NA*SRS3+1*SRS4+0*SRS5+start(-1)*SRS3+lSRS33*SRS3
e54=~0*SRS1+0*SRS2+0*SRS3+NA*SRS4+1*SRS5+start(-1)*SRS4+lSRS44*SRS4
e51~~e52+e53+e54+re51e52*e52+re51e53*e53+re51e54*e54
e52~~e53+e54+re52e53*e53+re52e54*e54
e53~~e54+re53e54*e54
SRS1~~0*SRS1
SRS2~~0*SRS2
SRS3~~0*SRS3
SRS4~~0*SRS4
SRS5~~0*SRS5
wSRS1:=((1)*((1)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lSRS2)*((lSRS11)*((1)*((1)*((1))))) + (lSRS3)*((lSRS11)*((lSRS22)*((1)*((1))))) + -(lSRS4)*((lSRS11)*((lSRS22)*((lSRS33)*((1))))) + (lSRS5)*((lSRS11)*((lSRS22)*((lSRS33)*((lSRS44))))))
wSRS2:=-((lSRS11)*((1)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lSRS2)*((lSRS11)*((1)*((1)*((1))))) + (lSRS3)*((lSRS11)*((lSRS22)*((1)*((1))))) + -(lSRS4)*((lSRS11)*((lSRS22)*((lSRS33)*((1))))) + (lSRS5)*((lSRS11)*((lSRS22)*((lSRS33)*((lSRS44))))))
wSRS3:=((lSRS11)*((lSRS22)*((1)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lSRS2)*((lSRS11)*((1)*((1)*((1))))) + (lSRS3)*((lSRS11)*((lSRS22)*((1)*((1))))) + -(lSRS4)*((lSRS11)*((lSRS22)*((lSRS33)*((1))))) + (lSRS5)*((lSRS11)*((lSRS22)*((lSRS33)*((lSRS44))))))
wSRS4:=-((lSRS11)*((lSRS22)*((lSRS33)*((1))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lSRS2)*((lSRS11)*((1)*((1)*((1))))) + (lSRS3)*((lSRS11)*((lSRS22)*((1)*((1))))) + -(lSRS4)*((lSRS11)*((lSRS22)*((lSRS33)*((1))))) + (lSRS5)*((lSRS11)*((lSRS22)*((lSRS33)*((lSRS44))))))
wSRS5:=((lSRS11)*((lSRS22)*((lSRS33)*((lSRS44))))) / ((1)*((1)*((1)*((1)*((1))))) + -(lSRS2)*((lSRS11)*((1)*((1)*((1))))) + (lSRS3)*((lSRS11)*((lSRS22)*((1)*((1))))) + -(lSRS4)*((lSRS11)*((lSRS22)*((lSRS33)*((1))))) + (lSRS5)*((lSRS11)*((lSRS22)*((lSRS33)*((lSRS44))))))
SRS~~vSRS*SRS
e51~~ve51*e51
e52~~ve52*e52
e53~~ve53*e53
e54~~ve54*e54
vSRS1:=((1) * (vSRS)) * (1) + ((lSRS11) * (ve51)) * (lSRS11)
vSRS2:=((lSRS2) * (vSRS)) * (lSRS2) + ((1) * (ve51) + (lSRS22) * (re51e52)) * (1) + ((1) * (re51e52) + (lSRS22) * (ve52)) * (lSRS22)
vSRS3:=((lSRS3) * (vSRS)) * (lSRS3) + ((1) * (ve52) + (lSRS33) * (re52e53)) * (1) + ((1) * (re52e53) + (lSRS33) * (ve53)) * (lSRS33)
vSRS4:=((lSRS4) * (vSRS)) * (lSRS4) + ((1) * (ve53) + (lSRS44) * (re53e54)) * (1) + ((1) * (re53e54) + (lSRS44) * (ve54)) * (lSRS44)
vSRS5:=((lSRS5) * (vSRS)) * (lSRS5) + ((1) * (ve54)) * (1)
wstdSRS1:=((wSRS1) * (sqrt(vSRS1))) * (1/sqrt(vSRS))
wstdSRS2:=((wSRS2) * (sqrt(vSRS2))) * (1/sqrt(vSRS))
wstdSRS3:=((wSRS3) * (sqrt(vSRS3))) * (1/sqrt(vSRS))
wstdSRS4:=((wSRS4) * (sqrt(vSRS4))) * (1/sqrt(vSRS))
wstdSRS5:=((wSRS5) * (sqrt(vSRS5))) * (1/sqrt(vSRS))

# Specify the relationships among the emergent variables
GRS ~~ PSS + SRS + Reader
PSS ~~ SRS + Reader
SRS ~~ Reader
```

As one can see, the output of the `specifyHO` function is a character string that contains the model which can be directly used as input for the `sem()` or `cfa()` function of lavaan. Hereby it is important that the `orthogonal` argument is set to `TRUE` to ensure the correct parameterization:
``` r
library(lavaan)
out = sem(model = modelHO, data = dataset, orthogonal = TRUE)
summary(out, standardized = TRUE)
```

