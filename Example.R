# Load the function from GitHub
devtools::source_url("https://github.com/FloSchuberth/HOspecification/blob/main/specifyHO.R?raw=TRUE")

model='
 PSS <~ PSS1+PSS2 +PSS3+ PSS4+ PSS5
GRS<~GRS1+ GRS2+GRS3+GRS4+GRS5
SRS<~SRS1+SRS2+SRS3+SRS4+SRS5

PSS~~GRS + SRS
GRS~~SRS
'

modelHO=specifyHO(.model=model,
                .typeHO='refined',
                .order_indicators = 'exact',
                .print_to_console = T,
                .determine_weights = T)

modelHO
library(lavaan)
out= sem(model = modelHO,data=dataset,orthogonal=T)
summary(out,standardized=T)

