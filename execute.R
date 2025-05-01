datasetNM<-"sample"
wdir<-"~/learning-curve-framework"
setwd(wdir)
source("functions.R")

df<-read.csv(file.path(wdir, "data", str_glue("{datasetNM}.csv")))
covariates <- colnames(df)[startsWith(colnames(df), 'Pt_')]

output<-PLCAnalysis(
        #==========================================
        # Part 1: Data management
        #==========================================
        data=df,
        datasetIdentifier=datasetNM,
        caseIDFieldNM="Patient",
        caseDateFieldNM="ProcDate",
        outcomeFieldNM="Outcome_Final",
        orderFieldNM="CaseOrder_All",
        operatorFieldNM="Operator",
        covariateFieldNMs = covariates,
        #==========================================
        # Part 2: Exposure 
        #==========================================    
        exposureFieldNM="Device",
        exposureOfInterestNM="B",
        exposureOfInterestOperatorCaseSeriesFieldNM="CaseOrder_Op_DevB",
        #==========================================
        # Part 3: IPTW
        #==========================================    
        allowIPTWWeighting=FALSE,
        normalizeIPTWWeights=FALSE,
        useGeneralCovariateFieldNMs=TRUE,
        iptwCovariateFieldNMs=c(), 
        #==========================================
        # Part 4: GAM Modelling
        #==========================================    
        allowOperatorClustering=FALSE,
        allowGAMBasedAVS=FALSE,
        #==========================================
        # Part 5: Learning effect detection
        #==========================================        
        learningUnadjustedSignalRequired=TRUE,
        leDetectionAlpha=0.05,
        #==========================================
        # Part 6: Learning curve estimation
        #==========================================    
        lcEstimationAlpha=0.05,
        bootstraps=1,
        forceEstimationAsymptoteToZero=TRUE,
        #==========================================
        # Part 7: Prospective analysis
        #==========================================    
        allowProspectiveAnalysis=FALSE,
        allowAlphaSpending=TRUE,
        alphaSpendingMethod="asKD", # asKD, asOF
        alphaToSpend=0.05,
        sequencerFieldNM="ProcDate",
        exposurePeriodFieldNM="OutcomePeriod",
        outcomePeriodFieldNM="OutcomePeriod"
)

# Learning effects level detected
output$results$learningDetection$level

# Learning curve estimated
output$results$lcQuantification$formEstimated

# Name of the functional form identified
output$results$lcQuantification$formEstimatedName

# Parameters of the estimated form
output$results$lcQuantification$estimatedFormDetails$selected

# Learning curve plot
x<-output$results$lcQuantification$estimatedLearningCurve$CaseOrder
y<-output$results$lcQuantification$estimatedLearningCurve$estLearningRisk
base::plot(x, y, xlab="Case order (Experience)", ylab="Adverse event risk (%)", main="Estimated learning curve")








