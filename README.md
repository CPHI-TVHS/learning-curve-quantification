## Getting started:
In this reporsitory, we have provided an example to execute a learning curve detection and quantification pipeline. To get started please explore the file `execute.R`. A sample dataset named `sample.csv` has been provided in `/data` and the necessary functions are located in file `functions.R`.

## Running the pipeline:  
Specify the data and arguments then run the script `source(execute.R)`. The main function `PLCAnalysis()` is executed and and results R object is generated.  

## Results:  
The results R object contains a number of outputs. For example;  

### Learning effects level detected  
`output$results$learningDetection$level`  

### Learning curve estimated  
`output$results$lcQuantification$formEstimated`  

### Name of the functional form identified  
`output$results$lcQuantification$formEstimatedName`  

### Parameters of the estimated form  
`output$results$lcQuantification$estimatedFormDetails$selected`  

### Learning curve plot  
`x<-output$results$lcQuantification$estimatedLearningCurve$CaseOrder`  
`y<-output$results$lcQuantification$estimatedLearningCurve$estLearningRisk`  
`base::plot(x, y, xlab="Case order (Experience)", ylab="Adverse event risk (%)", main="Estimated learning curve")`  
