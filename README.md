# News

Update 01-10-2025:

- Implement the blended and phantom-variable H-O specification
- Remove the standardized weights (these cannot be calculated when a composite is a dependent variable)
- Remove the option to fix weights (will be implemented in the future)
- The original version of the specifyHO() function can still be access as specifyHO_deprecated (originalHOspecify.R)

# Function that facilitates using the H-O specification

To load the specifyHO, you can use the following syntax:
``` r
devtools::source_url("https://github.com/FloSchuberth/HOspecification/blob/main/specifyHO.R?raw=TRUE")
```

In a next step, the model is specified in lavaan syntax using the `<~` operator to specify emergent variables. 
It is noteworthy that the model syntax for the emergent variables must not contain any parameter labels and 
emergent variables must not be specified in various lines.  
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
cat(modelHO)
```
Here the resulting model syntax:
```
# Specify the emergent variables
PSS=~PSS1+PSS2+PSS3+PSS4+PSS5+1*PSS1+lPSS2*PSS2+lPSS3*PSS3+lPSS4*PSS4+lPSS5*PSS5+start(0)*PSS2+start(0)*PSS3+start(0)*PSS4+start(0)*PSS5
PSS~~vPSS*PSS
e31=~NA*PSS1+1*PSS2+0*PSS3+0*PSS4+0*PSS5+start(-1)*PSS1+lPSS11*PSS1
e32=~0*PSS1+NA*PSS2+1*PSS3+0*PSS4+0*PSS5+start(-1)*PSS2+lPSS22*PSS2
e33=~0*PSS1+0*PSS2+NA*PSS3+1*PSS4+0*PSS5+start(-1)*PSS3+lPSS33*PSS3
e34=~0*PSS1+0*PSS2+0*PSS3+NA*PSS4+1*PSS5+start(-1)*PSS4+lPSS44*PSS4
e31~~ve31*e31
e32~~ve32*e32
e33~~ve33*e33
e34~~ve34*e34
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
GRS=~GRS1+GRS2+GRS3+GRS4+GRS5+1*GRS1+lGRS2*GRS2+lGRS3*GRS3+lGRS4*GRS4+lGRS5*GRS5+start(0)*GRS2+start(0)*GRS3+start(0)*GRS4+start(0)*GRS5
GRS~~vGRS*GRS
e41=~NA*GRS1+1*GRS2+0*GRS3+0*GRS4+0*GRS5+start(-1)*GRS1+lGRS11*GRS1
e42=~0*GRS1+NA*GRS2+1*GRS3+0*GRS4+0*GRS5+start(-1)*GRS2+lGRS22*GRS2
e43=~0*GRS1+0*GRS2+NA*GRS3+1*GRS4+0*GRS5+start(-1)*GRS3+lGRS33*GRS3
e44=~0*GRS1+0*GRS2+0*GRS3+NA*GRS4+1*GRS5+start(-1)*GRS4+lGRS44*GRS4
e41~~ve41*e41
e42~~ve42*e42
e43~~ve43*e43
e44~~ve44*e44
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
SRS=~SRS1+SRS2+SRS3+SRS4+SRS5+1*SRS1+lSRS2*SRS2+lSRS3*SRS3+lSRS4*SRS4+lSRS5*SRS5+start(0)*SRS2+start(0)*SRS3+start(0)*SRS4+start(0)*SRS5
SRS~~vSRS*SRS
e51=~NA*SRS1+1*SRS2+0*SRS3+0*SRS4+0*SRS5+start(-1)*SRS1+lSRS11*SRS1
e52=~0*SRS1+NA*SRS2+1*SRS3+0*SRS4+0*SRS5+start(-1)*SRS2+lSRS22*SRS2
e53=~0*SRS1+0*SRS2+NA*SRS3+1*SRS4+0*SRS5+start(-1)*SRS3+lSRS33*SRS3
e54=~0*SRS1+0*SRS2+0*SRS3+NA*SRS4+1*SRS5+start(-1)*SRS4+lSRS44*SRS4
e51~~ve51*e51
e52~~ve52*e52
e53~~ve53*e53
e54~~ve54*e54
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

