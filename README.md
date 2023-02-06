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
PSS <~ PSS1+PSS2 +PSS3+ PSS4+ PSS5
GRS<~GRS1+ GRS2+GRS3+GRS4+GRS5
SRS<~SRS1+SRS2+SRS3+SRS4+SRS5
PSS~~GRS + SRS
GRS~~SRS
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
The output of the `specifyHO` function is a character string that contains the model and that can be directly used as input for the `sem()` or `cfa()` function of lavaan. Hereby it is important that the `orthogonal` argument is set to `TRUE` to ensure the correct parameterization:
``` r
library(lavaan)
out = sem(model = modelHO, data = dataset, orthogonal = TRUE)
summary(out, standardized = TRUE)
```
