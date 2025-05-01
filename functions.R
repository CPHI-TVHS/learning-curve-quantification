suppressPackageStartupMessages(library(tidyverse, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(gbm, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(mgcv, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(parameters, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(lmtest, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(minpack.lm, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(plotly, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(rpact, quietly=T, warn.conflicts=F))
suppressPackageStartupMessages(library(optparse, quietly=T, warn.conflicts=F))

parametric.data.management<-function(data, caseIDFieldNM, caseDateFieldNM, outcomeFieldNM, exposureFieldNM, exposureOfInterestNM, orderFieldNM, operatorFieldNM, exposureOfInterestOperatorCaseSeriesFieldNM, covariateFieldNMs, datasetIdentifier, periodIdentifier){
    # Data management for the parametric pipeline-returns simple metrics and a clean/formated dataset for further procesess.
    data <-as.data.frame(data)
    totalCases<-dim(data)[1]
    if (totalCases==0){
        # Stop if dataframe is empty
        output<-list(
            datasetIdentifier=datasetIdentifier, 
            periodIdentifier=periodIdentifier,
            exposureOfInterestNM=exposureOfInterestNM,
            message<-"Data provided has no records...",
            totalCases=totalCases,
            totalOutcomeCount=0,
            exposureOfInterestCount=0,
            exposureOfInterestOutcomeCount=0,
            covariatesCount=length(covariateFieldNMs),
            data=data,
            passed=F
        )
    } else {
        # We have data, lets go
        data$caseID<-data[[caseIDFieldNM]]
        data$caseDate<-data[[caseDateFieldNM]]
        data$outcome <- data[[outcomeFieldNM]]
        data$exposureNM <- data[[exposureFieldNM]]
        data$orderID <- data[[orderFieldNM]]
        data$operatorNM <- data[[operatorFieldNM]]
        data$operatorSeries <- data[[exposureOfInterestOperatorCaseSeriesFieldNM]]
        keyFields<-c('caseID', 'caseDate', 'exposureNM', 'exposure', 'orderID', 'operatorNM', 'operatorSeries', 'outcome', covariateFieldNMs)
        data<-data %>% 
            dplyr::arrange(orderID) %>% 
            dplyr::mutate(
                operatorNM=as.factor(operatorNM),
                exposure=ifelse(exposureNM==exposureOfInterestNM, 1, 0),
                caseDate=as.Date(caseDate, format="%Y-%m-%d"),
                operatorSeries=as.integer(operatorSeries)
            ) %>% 
            dplyr::select(all_of(keyFields)) %>% 
            dplyr::arrange(orderID)
        # Generate counts
        data.a<-data %>% 
            dplyr::filter(exposure==1)
        data.b<-data.a %>% 
            dplyr::filter(outcome==1)
        data.c<-data %>% 
            dplyr::filter(outcome==1)
        # Detect any covariates that are characters and convert them to factors
        for (i in covariateFieldNMs){
            if (is.character(data[[i]])==T){
                data[[i]]<-as.factor(data[[i]])
            }
        }
        output<-list(
            datasetIdentifier=datasetIdentifier, 
            periodIdentifier=periodIdentifier,
            message<-stringr::str_glue("Data provided has {totalCases} records..."),
            exposureOfInterestNM=exposureOfInterestNM,
            totalCases=totalCases,
            totalOutcomeCount=dim(data.c)[1],
            exposureOfInterestCount=dim(data.a)[1],
            exposureOfInterestOutcomeCount=dim(data.b)[1],
            covariatesCount=length(covariateFieldNMs),
            data=data,
            passed=T
        )
    } 
    return(output)
}

models.generate.formular.string<-function(target, features){
    f <- as.formula(
        paste(target, 
              paste(features, collapse = " + "), 
              sep = " ~ "))
    #print(f)
    return(f)    
}

generate.iptw.weights<-function(data, allowIPTWWeighting, normalizeIPTWWeights, useGeneralCovariateFieldNMs, covariateFieldNMs, iptwCovariateFieldNMs){
    if (allowIPTWWeighting==TRUE){
        if (useGeneralCovariateFieldNMs==T){
            covariates<-covariateFieldNMs
        } else {
            covariates<-iptwCovariateFieldNMs
        }
        ft<-models.generate.formular.string('exposure', covariates)
        model.t<-gbm(ft, distribution = "bernoulli", data=data, n.trees = 1000)
        # summary(model.t)
        data$pscores<-predict(model.t, type='response')
        # summary(data$pscores)
        data$weights <- ifelse(data$exposure==1, 1/data$pscores, 1/(1-data$pscores))
        # summary(data$weights)
        #--------------------------------------------------------------------------------------------------
        # # STABILIZATION OF WEIGHTS
        # if (config$stabilizeTreatmentCofoundingWeights==TRUE){
        #     a<-dim(subset(data, exposure==0))[1]/dim(data)[1]
        #     b<-dim(subset(data, exposure==1))[1]/dim(data)[1]
        #     data$weights.stablized<-ifelse(data$exposure==0, data$weights*a, NA)
        #     data$weights.stablized<-ifelse(data$exposure==1, data$weights*b, data$weights.stablized)
        #     data$weights<-data$weights.stablized
        # } 
        # # TRIMMING WEIGHTS
        # if (config$trimTreatmentCofoundingWeights==TRUE){
        #     Q05<- stats::quantile(data$weights, c(config$weightTrimLowerCutoff))
        #     Q95 <- stats::quantile(data$weights, c(config$weightTrimUpperCutoff))
        #     data$lower_cutoff<-Q05[[1]]
        #     data$upper_cutoff<-Q95[[1]]
        #     #set trimmed equal to raw scores
        #     data$trimmed_weights<-0
        #     #apply limits
        #     data$trimmed_weights <- ifelse(data$weights < data$lower_cutoff, data$lower_cutoff, data$weights)
        #     data$trimmed_weights <- ifelse(data$weights > data$upper_cutoff, data$upper_cutoff, data$weights)
        #     data$weights<-data$trimmed_weights
        # }
        # NORMALIZE WEIGHTS FOR GAM
        if (normalizeIPTWWeights==TRUE){
            data$weights<-data$weights/mean(data$weights)
        }
        #summary(data$weights)
    } else {
        data$weights<-NA
    }
    return(data)
}

signal_detection_gam<-function(target, features, data, method="REML", allowAVS=F, allowWeights=F, alphaLevel=0.05){
    fm<-models.generate.formular.string(target, features)
    smodel<-mgcv::bam(fm, data=data, method=method, family = binomial('logit'), select=allowAVS)
    if (allowWeights==T){
        smodel<-mgcv::bam(fm, data=data, method=method, family = binomial('logit'), select=allowAVS, weights=data$weights)
    }
    data$estimated_baseline_probabilities<-predict(smodel, type='response')
    data<-data %>%
        dplyr::select(caseID, outcome, estimated_baseline_probabilities)
    mps<-parameters::model_parameters(smodel, exponentiate = TRUE, pretty_names = TRUE, summary=TRUE)
    sart<-as.data.frame(mps)
    sres<-subset(sart, Parameter=='exposure')
    od<-round(as.numeric(sres['Coefficient']), 4)
    odl<-round(as.numeric(sres['CI_low']), 4)
    odu<-round(as.numeric(sres['CI_high']), 4)
    pv<-round(as.numeric(sres['p']), 8) 
    yes<-ifelse(odl>1 & pv<=alphaLevel, TRUE, FALSE)
    auc<-pROC::auc(data$outcome, data$estimated_baseline_probabilities)
    return(list(
        form=fm,
        model=smodel,
        model.type='gam',
        estimates=sart,
        hasSignal=yes,
        signal=od,
        signalLower=odl,
        signalUpper=odu,
        pvalue=pv,
        data=data,
        AIC=AIC(smodel),
        AUC=as.numeric(auc),
        alphaLevel=alphaLevel
    ))
}

learning_detection<-function(liuModel, liaModel, leDetectionAlpha){
    lrt<-lmtest::lrtest(liuModel, liaModel)
    lrt.df<-as.data.frame(lrt)
    lrt.pvalue<-lrt.df$`Pr(>Chisq)`[[2]]
    lrt.Chisq<-lrt.df$Chisq[[2]]
    det<-list(level='none')
    if (lrt.pvalue<leDetectionAlpha){
        det$level<-'operator'
    }
    return(
        list(
            lrTestPvalue=round(lrt.pvalue, 5),
            lrTestAlphaUsed=leDetectionAlpha,
            lrTestChisq=round(lrt.Chisq, 2),
            level=det$level
        )
    )
}

lc_estimate_form<-function(data, sequencer, learningRiskFieldNM, alpha=0.01, prefer.rse=TRUE){
    raw<-NA
    data$X<-data[[sequencer]]
    data$Y<-data[[learningRiskFieldNM]]
    data<-data %>% dplyr::select(X, Y)
    #print(head(data))
    wc<-FALSE
    pc<-FALSE
    ec<-FALSE
    lc<-FALSE
    rc<-FALSE
    model<-NA
    #------------------------------------------------------------------------------- 
    a<-tryCatch({
        weibull<-function(x, a, b, c, d){a+b*exp(-1*c*x^d)}
        wf<-as.formula(Y ~ weibull(X, a, b, c, d))
        wm <- nlsLM(wf, data=data, start = list(a=1, b=1, c=1, d=1), control = nls.lm.control(maxiter = 1000))
        model<-wm
        #print(summary(wm))
        wc<-wm$convInfo$isConv
        #print(wc)
    }, warning = function(w){
        #code that handles the warnings
        #print(w)
    }, error = function(e){
        #code that handles the errors
        #print(e)
    }, finally = function(f){
        #clean-up code
    }) 
    #print(a)
    #-------------------------------------------------------------------------------
    a<-tryCatch({
        power<-function(x, a, b, c){a+b*x^(-1 * c)}
        pf<-as.formula(Y ~ power(X, a, b, c))
        pm <- nlsLM(pf, data=data, start = list(a=1, b=1, c=1), control = nls.lm.control(maxiter = 1000))
        pc<-pm$convInfo$isConv
        model<-pm
        #print(summary(pm))        
    }, warning = function(w){
        #code that handles the warnings
    }, error = function(e){
        #code that handles the errors
    }, finally = function(f){
        #clean-up code
    })    
    #-------------------------------------------------------------------------------
    a<-tryCatch({
        expo<-function(x, a, b, c){a+b*base::exp(-1 * c*x)}
        ef<-as.formula(Y ~ expo(X, a, b, c))
        em <- nlsLM(ef, data=data, start = list(a=1, b=1, c=1), control = nls.lm.control(maxiter = 1000))
        ec<-em$convInfo$isConv
        model<-em
        #print(summary(em))        
    }, warning = function(w){
        #code that handles the warnings
    }, error = function(e){
        #code that handles the errors
    }, finally = function(f){
        #clean-up code
    })    
    #-------------------------------------------------------------------------------   
    # a<-tryCatch({
    #     rec<-function(x, a, b){a+(b/x)}
    #     rf<-as.formula(Y ~ rec(sequencer, a, b))
    #     rm <- nlsLM(rf, data=data, start = list(a=1, b=1), control = nls.lm.control(maxiter = 1000))
    #     rc<-rm$convInfo$isConv
    # model<-rm
    #     #print(summary(rm))        
    # }, warning = function(w){
    #     #code that handles the warnings
    # }, error = function(e){
    #     #code that handles the errors
    # }, finally = function(f){
    #     #clean-up code
    # })    
    #-------------------------------------------------------------------------------
    # a<-tryCatch({
    #     logistic<-function(x, a, b, c, d){a+b/(1+c*base::exp(d*x))}
    #     lf<-as.formula(Y ~ logistic(sequencer, a, b, c, d))
    #     lm <- nlsLM(lf, data=data, start = list(a=1, b=1, c=1, d=1), control = nls.lm.control(maxiter = 1000))
    #     lc<-lm$convInfo$isConv
    #     model<-lm
    #     #print(summary(lm))        
    # }, warning = function(w){
    #     #code that handles the warnings
    # }, error = function(e){
    #     #code that handles the errors
    # }, finally = function(f){
    #     #clean-up code
    # })    
    #-------------------------------------------------------------------------------    
    sdf0<-data.frame(var=as.character(), Estimate=as.numeric(), `Std. Error`=as.numeric(), `Pr(>|t|)`=as.numeric(), lci=as.numeric(), uci=as.numeric(), ci=as.numeric(), form=as.character(), converged=as.logical(), iterations=as.integer())
    forms<-list()
    cv<-list()
    ses<-list()
    if (pc==TRUE){
        cv$power<-TRUE
        x<-summary(pm)
        ses$power<-x$sigma
        sdf<-as.data.frame(x$coefficients)
        sdf$var<-row.names(sdf)
        sdf<-sdf[c('var', 'Estimate', 'Std. Error', 'Pr(>|t|)')]
        sdf$lci<-sdf$Estimate-(sdf$`Std. Error`*1.96)
        sdf$uci<-sdf$Estimate+(sdf$`Std. Error`*1.96)
        sdf$ci<-sdf$uci-sdf$lci
        sdf$form<-'power'
        sdf$converged<-pm$convInfo$isConv
        sdf$iterations<-pm$convInfo$finIter
        sdf$rse<-sum(resid(pm)^2)# x$sigma
        sdf0<-base::rbind(sdf0, sdf)
        raw<-sdf0
        forms$power<-list(
            form='a+b*x^(-1 * c)', 
            name='Power', 
            b=subset(sdf, var=='b')$Estimate,
            bll=subset(sdf, var=='b')$Estimate-subset(sdf, var=='b')$`Std. Error`,
            bul=subset(sdf, var=='b')$Estimate+subset(sdf, var=='b')$`Std. Error`,
            c=subset(sdf, var=='c')$Estimate,
            cll=subset(sdf, var=='c')$Estimate-subset(sdf, var=='c')$`Std. Error`,
            cul=subset(sdf, var=='c')$Estimate+subset(sdf, var=='c')$`Std. Error`,
            d=as.character(NA),
            dll=as.character(NA),
            dul=as.character(NA)
        )
    }
    if (wc==TRUE){
        cv$weibull<-TRUE
        x<-summary(wm)
        ses$weibull<-x$sigma
        sdf<-as.data.frame(x$coefficients)
        sdf$var<-row.names(sdf)
        sdf<-sdf[c('var', 'Estimate', 'Std. Error', 'Pr(>|t|)')]
        sdf$lci<-sdf$Estimate-(sdf$`Std. Error`*1.96)
        sdf$uci<-sdf$Estimate+(sdf$`Std. Error`*1.96)
        sdf$ci<-sdf$uci-sdf$lci
        sdf$form<-'weibull'
        sdf$converged<-wm$convInfo$isConv
        sdf$iterations<-wm$convInfo$finIter  
        sdf$rse<-sum(resid(wm)^2)# x$sigma
        sdf0<-base::rbind(sdf0, sdf)
        raw<-sdf0
        forms$weibull<-list(
            form='a+b exp(-cx^d)', 
            name='Weibull',
            b=subset(sdf, var=='b')$Estimate, 
            bll=subset(sdf, var=='b')$Estimate-subset(sdf, var=='b')$`Std. Error`,
            bul=subset(sdf, var=='b')$Estimate+subset(sdf, var=='b')$`Std. Error`,        
            c=subset(sdf, var=='c')$Estimate, 
            cll=subset(sdf, var=='c')$Estimate-subset(sdf, var=='c')$`Std. Error`,
            cul=subset(sdf, var=='c')$Estimate+subset(sdf, var=='c')$`Std. Error`,        
            d=subset(sdf, var=='d')$Estimate,
            dll=subset(sdf, var=='d')$Estimate-subset(sdf, var=='d')$`Std. Error`,
            dul=subset(sdf, var=='d')$Estimate+subset(sdf, var=='d')$`Std. Error`        
        )    
    }
    if (rc==TRUE){
        cv$reciprocal<-TRUE
        x<-summary(rm)
        ses$reciprocal<-x$sigma
        sdf<-as.data.frame(x$coefficients)
        sdf$var<-row.names(sdf)
        sdf<-sdf[c('var', 'Estimate', 'Std. Error', 'Pr(>|t|)')]
        sdf$lci<-sdf$Estimate-(sdf$`Std. Error`*1.96)
        sdf$uci<-sdf$Estimate+(sdf$`Std. Error`*1.96)
        sdf$ci<-sdf$uci-sdf$lci
        sdf$form<-'reciprocal'
        sdf$converged<-rm$convInfo$isConv
        sdf$iterations<-rm$convInfo$finIter 
        sdf$rse<-sum(resid(rm)^2)# x$sigma
        sdf0<-base::rbind(sdf0, sdf)
        raw<-sdf0
        forms$reciprocal<-list(
            form='a+b/x',
            name='Reciprocal',
            b=subset(sdf, var=='b')$Estimate,
            bll=subset(sdf, var=='b')$Estimate-subset(sdf, var=='b')$`Std. Error`,
            bul=subset(sdf, var=='b')$Estimate+subset(sdf, var=='b')$`Std. Error`,        
            c=as.character(NA), cll=as.character(NA), cul=as.character(NA),
            d=as.character(NA), dll=as.character(NA), dul=as.character(NA)
        )    
    }
    if (ec==TRUE){
        cv$exponential<-TRUE
        x<-summary(em)
        ses$exponential<-x$sigma
        sdf<-as.data.frame(x$coefficients)
        sdf$var<-row.names(sdf)
        sdf<-sdf[c('var', 'Estimate', 'Std. Error', 'Pr(>|t|)')]
        sdf$lci<-sdf$Estimate-(sdf$`Std. Error`*1.96)
        sdf$uci<-sdf$Estimate+(sdf$`Std. Error`*1.96)
        sdf$ci<-sdf$uci-sdf$lci
        sdf$form<-'exponential'
        sdf$converged<-em$convInfo$isConv
        sdf$iterations<-em$convInfo$finIter 
        sdf$rse<-sum(resid(em)^2)# x$sigma
        sdf0<-base::rbind(sdf0, sdf)
        raw<-sdf0
        forms$exponential<-list(
            form='a+b exp(-cx)', 
            name='Exponential', 
            b=subset(sdf, var=='b')$Estimate,
            bll=subset(sdf, var=='b')$Estimate-subset(sdf, var=='b')$`Std. Error`,
            bul=subset(sdf, var=='b')$Estimate+subset(sdf, var=='b')$`Std. Error`,
            c=subset(sdf, var=='c')$Estimate,
            cll=subset(sdf, var=='c')$Estimate-subset(sdf, var=='c')$`Std. Error`,
            cul=subset(sdf, var=='c')$Estimate+subset(sdf, var=='c')$`Std. Error`,
            d=as.character(NA),
            dll=as.character(NA),
            dul=as.character(NA)
        )
    }
    if (lc==TRUE){
        cv$logistic<-TRUE
        x<-summary(lm)
        ses$logistic<-x$sigma
        sdf<-as.data.frame(x$coefficients)
        sdf$var<-row.names(sdf)
        sdf<-sdf[c('var', 'Estimate', 'Std. Error', 'Pr(>|t|)')]
        sdf$lci<-sdf$Estimate-(sdf$`Std. Error`*1.96)
        sdf$uci<-sdf$Estimate+(sdf$`Std. Error`*1.96)
        sdf$ci<-sdf$uci-sdf$lci
        sdf$form<-'logistic'
        sdf$converged<-lm$convInfo$isConv
        sdf$iterations<-lm$convInfo$finIter
        sdf$rse<-sum(resid(lm)^2)# x$sigma
        #print(sdf)
        sdf0<-base::rbind(sdf0, sdf)
        raw<-sdf0
        forms$logistic<-list(
            form='a+b/(1 + c exp(dx))', 
            name='Logistic',
            b=subset(sdf, var=='b')$Estimate, 
            bll=subset(sdf, var=='b')$Estimate-subset(sdf, var=='b')$`Std. Error`,
            bul=subset(sdf, var=='b')$Estimate+subset(sdf, var=='b')$`Std. Error`,        
            c=subset(sdf, var=='c')$Estimate, 
            cll=subset(sdf, var=='c')$Estimate-subset(sdf, var=='c')$`Std. Error`,
            cul=subset(sdf, var=='c')$Estimate+subset(sdf, var=='c')$`Std. Error`,        
            d=subset(sdf, var=='d')$Estimate,
            dll=subset(sdf, var=='d')$Estimate-subset(sdf, var=='d')$`Std. Error`,
            dul=subset(sdf, var=='d')$Estimate+subset(sdf, var=='d')$`Std. Error`        
        )    
    } 
    if (length(cv)>0){
        sdf0<-data.columns.rename(sdf0, c('Pr(>|t|)'), c('pvalue'))
        df<-sdf0 %>% dplyr::mutate(pvalue=base::round(pvalue, 6))
        df$pvalue_passed<-base::ifelse(df$pvalue<alpha, 1, 0)
        df$isParameter<-base::ifelse(df$var %in% c('b', 'c', 'd'), 1, 0)
        df<-base::subset(df, converged==TRUE & isParameter==1)
        df2<-df %>%
            dplyr::group_by(form) %>%
            dplyr::summarise(form_passed=sum(pvalue_passed)/sum(isParameter)) %>%
            dplyr::mutate(form_passed_highest=ifelse(form_passed==max(form_passed), 1, 0)) %>%
            dplyr::select(form, form_passed_highest, form_passed)
        df<-base::merge(x=df, y=df2, by='form', all.x=TRUE)
        df$rse<-base::round(df$rse, 2)
        #flag negative parameters
        dfn<-df
        dfn$npara<-ifelse(dfn$Estimate<0, 1, 0)
        dfn <- dfn %>% 
            dplyr::group_by(form) %>%
            dplyr::summarise(flag=sum(npara)) %>%
            dplyr::mutate(npflag=ifelse(flag>0, 1, 0)) %>%
            dplyr::select(form, npflag)
        df<-base::merge(df, dfn, by='form')
        sdf01<-df         
    }
    
    if (length(cv)==0){
        return(list(model=model, raw=raw, all=sdf0, forms=forms, convs=cv, ses=ses, estimated=0, selected=list(name=NA, b=NA, bll=NA, bul=NA, c=NA, cll=NA, cul=NA, d=NA, dll=NA, dul=NA)))
    } else if (length(cv)==1){
        df3<-base::subset(df, form_passed==1.0 & npflag==0)
        if (length(base::unique(df3$form))==1){
            # form with most optimized parameters is selected... 
            return(list(model=model, raw=raw, all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df3$form)[1]]]))
        } else {
            return(list(model=model, raw=raw, all=sdf0, forms=forms, convs=cv, ses=ses, estimated=0, selected=list(name=NA, b=NA, bll=NA, bul=NA, c=NA, cll=NA, cul=NA, d=NA, dll=NA, dul=NA)))
        }
        return(list(model=model, raw=raw, all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[names(cv)[1]]]))
    } else if (length(cv)>1){
        df3<-base::subset(df, form_passed==1.0 & npflag==0)
        if (length(base::unique(df3$form))==1){
            # form with most optimized parameters is selected... 
            return(list(model=model, raw=raw,all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df3$form)[1]]]))
        } else {
            # if tied, selected the one with the lowest iterations if prefer rse=FALSE
            if (prefer.rse==F){
                #print('More than 1 with 0.5...')
                minIt<-base::min(df3$iterations)
                df3$hasMinIt<-base::ifelse(df3$iterations==minIt, 1, 0)
                df5<-base::subset(df3, hasMinIt==1  & npflag==0)
                if (length(base::unique(df5$form))==1){
                    # form with the 
                    return(list(model=model, raw=raw, all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df5$form)[1]]]))
                } else {
                    # select form with most parameters
                    #print('More than 1 has min iterations...')
                    df6<-df5 %>%
                        dplyr::group_by(form) %>%
                        dplyr::summarise(mostParams=sum(isParameter)) 
                    df5<-base::merge(x=df5, y=df6, by='form', all.x=TRUE)
                    df7<-base::subset(df5, mostParams==1  & npflag==0)
                    if (length(base::unique(df7$form))==1){
                        return(list(model=model, raw=raw,all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df7$form)[1]]]))
                    } else {
                        #print('FAILED to identify the form...')
                        return(list(model=model, raw=raw, all=sdf01, forms=forms, convs=cv, ses=ses, estimated=0, selected=list(name=NA, b=NA, bll=NA, bul=NA, c=NA, cll=NA, cul=NA, d=NA, dll=NA, dul=NA)))
                    }                
                }                
            } else {
                # give preference to form with the lowest rse and break tie with less more parameters?
                minRse<-base::min(df3$rse)
                df3$hasMinRse<-base::ifelse(df3$rse==minRse, 1, 0)
                df5<-base::subset(df3, hasMinRse==1 & npflag==0)
                if (length(base::unique(df5$form))==1){
                    # form with the 
                    return(list(model=model, raw=raw, all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df5$form)[1]]]))
                } else {
                    # select form with most parameters
                    df6<-df5 %>%
                        dplyr::group_by(form) %>%
                        dplyr::summarise(mostParams=sum(isParameter)) 
                    df5<-base::merge(x=df5, y=df6, by='form', all.x=TRUE)
                    df7<-base::subset(df5, mostParams==1 & npflag==0)
                    if (length(base::unique(df7$form))==1){
                        return( list(model=model, raw=raw,all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df7$form)[1]]]))
                    } else {
                        # select least iterations
                        minIt<-base::min(df3$iterations)
                        df3$hasMinIt<-base::ifelse(df3$iterations==minIt, 1, 0)
                        df5<-base::subset(df3, hasMinIt==1 & npflag==0)
                        if (length(base::unique(df5$form))==1){
                            # form with the 
                            return( list(model=model, raw=raw,all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df5$form)[1]]]))
                        } else {
                            # select form with most parameters among those with the least iterations
                            #print('More than 1 has min iterations...')
                            df6<-df5 %>%
                                dplyr::group_by(form) %>%
                                dplyr::summarise(mostParams=sum(isParameter)) 
                            df5<-base::merge(x=df5, y=df6, by='form', all.x=TRUE)
                            df7<-base::subset(df5, mostParams==1 & npflag==0)
                            if (length(base::unique(df7$form))==1){
                                return( list(model=model, raw=raw,all=sdf01, forms=forms, convs=cv, ses=ses, estimated=1, selected=forms[[base::unique(df7$form)[1]]]))
                            } else {
                                #print('FAILED to identify the form...')
                                return( list(model=model, raw=raw,all=sdf01, forms=forms, convs=cv, ses=ses, estimated=0, selected=list(name=NA, b=NA, bll=NA, bul=NA, c=NA, cll=NA, cul=NA, d=NA, dll=NA, dul=NA)))
                            }                
                        }
                    }                
                }                
            }
            
        }
    }
}

data.columns.rename<-function(data, names.old, names.new){
    #rename columns in old with new
    pos<-1
    for (i in names.old){
        names(data)[names(data) == i] <- names.new[pos]
        pos<-pos+1
    }
    return(data)
}

lc_apply_estimates<-function(model, data, selected, x.field, mse, y.field, n.boots=100, forceEstimationAsymptoteToZero){
    data$X<-data[[x.field]]
    data$Y<-data[[y.field]]
    ldf<-data %>% 
        dplyr::select(X, Y) %>%
        dplyr::arrange(X)
    # use bootstrapping method
    res<-nls_bootstrapping(model, n.boots, ldf, selected$name, forceEstimationAsymptoteToZero)
    pdf<-res$preds %>%
        dplyr::mutate(
            estLearningRisk=.fitted
        ) %>%
        dplyr::rename(CaseOrder=X, Est=.fitted, LCL=lwr_CI, UCL=upr_CI) %>%
        dplyr::mutate(Est=round(Est*100, 4), LCL=round(LCL*100, 4), UCL=round(UCL*100, 4)) %>%
        dplyr::select(CaseOrder, Est, LCL, UCL, estLearningRisk) %>%
        dplyr::arrange(CaseOrder)
    return(list(preds=pdf, fig=res$fig, params=res$params))
}

lc_quantification<-function(data, lcEstimationAlpha, learningRiskFieldNM, bootstraps, forceEstimationAsymptoteToZero){
    report<-list()
    estimatedForm<-lc_estimate_form(data=data, sequencer='operatorSeries', learningRiskFieldNM, alpha=lcEstimationAlpha, prefer.rse=TRUE) 
    report$estimatedFormDetails<-estimatedForm
    if (is.na(estimatedForm$selected$name)==FALSE){
        # E) IF THERE'S A FORM SELECTED, ESTIMATE THE LEARNING RISK BASED ON THE FORM SELECTED
        edt<-estimatedForm$raw
        report$parametersEstimates<-edt
        report$formEstimated<-TRUE
        report$formEstimatedName<-estimatedForm$selected$name
        report$formParameterB<-estimatedForm$selected$b
        report$formParameterBLCI<-estimatedForm$selected$bll
        report$formParameterBUCI<-estimatedForm$selected$bul
        report$formParameterC<-estimatedForm$selected$c
        report$formParameterCLCI<-estimatedForm$selected$cll
        report$formParameterCUCI<-estimatedForm$selected$cul
        report$formParameterD<-estimatedForm$selected$d
        report$formParameterDLCI<-estimatedForm$selected$dll
        report$formParameterDUCI<-estimatedForm$selected$dul
        x<-estimatedForm$raw %>%
            dplyr::filter(form==tolower(estimatedForm$selected$name) & var=='a') %>%
            dplyr::select('Estimate')
        asymp<-x$Estimate
        data$asymp<-asymp
        report$asymp<-asymp
        estimated<-lc_apply_estimates(model=estimatedForm$model, data=data, selected=estimatedForm$selected, x.field='operatorSeries', mse=estimatedForm$ses, y.field='learningRiskOP', n.boots=bootstraps, forceEstimationAsymptoteToZero)
        report$estimated<-estimated$preds
        estimated$estLearningRiskWithAsymp<-estimated$preds$estLearningRisk+asymp
        estimated3<-estimated$preds %>%
            dplyr::arrange(CaseOrder)
        estimated2<-estimated3 %>%
            dplyr::distinct(CaseOrder, estLearningRisk)
        report$estimatedLearningCurve<-estimated2
        report$estimatedLearningCurvePlot<-estimated$fig
    } else {
        # No functional form is selected
        report$formEstimated<-FALSE
        print('No functional form selected...')
    }
    return(report)
}  

nls_bootstrapping<-function(model, boots, data, form, zero_asymptote=TRUE, use.blocks=FALSE, get.cis=FALSE){
    # data must have fields-X and Y
    weibull<-function(x, a, b, c, d){a+b*base::exp(-c*x^d)}
    power<-function(x, a, b, c){a+b*x^(-c)}
    expo<-function(x, a, b, c){a+b*base::exp(-c*x)}
    logistic<-function(x, a, b, c, d){a+b/(1+c*base::exp(d*x))}
    ref<-list(
        Weibull=list(form=as.formula(Y ~ weibull(X, a, b, c, d)), np=4, starts=list(a=0, b=9.79e-08, c=4.43, d=0.239)),
        Power=list(form=as.formula(Y~power(X, a, b, c)), np=3, starts=list(a=1, b=1, c=1) ),
        Exponential=list(form=as.formula(Y~expo(X, a, b, c)), np=3, starts=list(a=1, b=1, c=1) ),
        Logistic=list(form=as.formula(Y~logistic(X, a, b, c, d)), np=4, starts=list(a=1, b=1, c=1, d=1) )
    )
    no_paras<-ref[[form]][['np']]
    fform<-ref[[form]][['form']]
    ss<-ref[[form]][['starts']] 
    # new data frame of predictions
    ds<-base::subset(data, select=c('X', 'Y'))
    new_preds <- ds %>%
        do(., data.frame( X = seq(min(.$X), max(.$X), by=1), stringsAsFactors = FALSE))
    new_preds0<-new_preds
    fit0<- model # nlsLM(fform, data=ds, start = ss, control = nls.lm.control(maxiter = 1000))
    if (zero_asymptote==TRUE){
        fa0<-fit0$m$getPars()
        fa0[['a']]<-0
        fit0$m$setPars(fa0)        
    } 
    new_preds0$.fitted<-predict(fit0, newdata=new_preds0)
    pdds<-data.frame(X=as.integer(), .fitted=as.numeric())
    fits<-list()
    params<-data.frame(boot_id=as.integer(), form=as.character(), term=as.character(), estimate=as.numeric(), std.error=as.numeric(), statistic=as.numeric(), p.value=as.numeric())
    blocks<-base::seq(1, length(unique(ds$X)), 1)
    boot.seq<-base::seq(1, boots, 1)
    print(str_glue('Bootstrapping on: {boots}'))
    if (get.cis==TRUE){
        for (k in boot.seq){
            #print(str_glue('Bootstrapping: {k} of {boots}'))
            if (use.blocks==TRUE){
                #create a boot dataset
                res<-data.frame(X=as.numeric(), Y=as.numeric())           
                for (i in blocks){
                    dsi<-base::subset(ds, X==i)
                    rns<-base::rownames(dsi)
                    idx<-base::sample(rns, size=length(rns), replace =T)
                    dss<-dsi[idx, ]
                    res<-base::rbind(res, dss)           
                }            
            } else {
                rns<-base::rownames(ds)
                idx<-base::sample(rns, size=length(rns), replace =T)
                res<-ds[idx, ]
            }
            # fit on the bootstap
            a<-tryCatch({
                fit <- nlsLM(fform, data=res, start = ss, control = nls.lm.control(maxiter = 1000))
                wc<-fit$convInfo$isConv
                #print(str_glue('Converged: {wc}'))
                if (wc==TRUE){
                    #print(summary(fit))
                    if (zero_asymptote==TRUE){
                        fa<-fit$m$getPars()
                        fa[['a']]<-0
                        fit$m$setPars(fa)        
                    }                
                    #fits<-append(fits, fit)
                    pms <- broom::tidy(fit)
                    pms$form<-form
                    pms$boot_id<-k
                    #print(pms)
                    params<-base::rbind(params, pms)
                    new_preds$.fitted<-predict(fit, newdata=new_preds)
                    pdds<-rbind(pdds, new_preds)                
                }
            }, warning = function(w){
                #code that handles the warnings
                #print(w)
            }, error = function(e){
                #code that handles the errors
                #print(e)
            }, finally = function(f){
                #clean-up code
            })         
        }
        pdds.ci<- pdds %>%
            dplyr::group_by(X) %>%
            dplyr::summarise(
                #est=mean(.fitted),
                #est_median=quantile(.fitted, 0.5),
                lwr_CI = quantile(.fitted, 0.025),
                upr_CI = quantile(.fitted, 0.975))        
        new_preds0<-base::merge(x=new_preds0, y=pdds.ci, by='X', all.x=T)
    } else {
        new_preds0$lwr_CI<-new_preds0$.fitted
        new_preds0$upr_CI<-new_preds0$.fitted
    }
    new_preds0<-as.data.frame(new_preds0)
    fig<-plot_ly(data = as.data.frame(new_preds0), x = ~X )
    fig <- fig %>% add_trace(y = ~.fitted, name = 'Estimated Learning Curve', mode = 'lines', type='scatter', line = list(color = 'rgb(0, 0, 255)', size = 10))
    fig <- fig %>% add_trace(x=ds$X, y=ds$Y, name = 'Predicted Learning Risk', type='scatter', mode = 'markers', marker = list(color = 'rgb(17, 157, 255)', size = 1, opacity = 0.5))
    fig <- fig %>% add_trace(y = ~lwr_CI, mode = 'lines', type='scatter', line=list(color=toRGB('cornflowerblue'), width=0.1), showlegend = F)
    fig <- fig %>% add_trace(y = ~upr_CI, mode = 'lines', type='scatter', line=list(color=toRGB('cornflowerblue'), width=0.1), fill = 'tonexty', fillcolor = toRGB("cornflowerblue", 0.3), showlegend = F)
    fig <- fig %>% layout(title = 'Outcomes by experience', xaxis = list(title = 'Case Number'), yaxis = list(title = 'Learning Adjustment (increased Risk)')) 
    #print(fig)
    return(list(
        fig=fig,
        preds=as.data.frame(new_preds0),
        params=as.data.frame(params)
    ))
}

obf.alpha.spending<-function(periods, events, totals, alpha=0.05){
    a<-1/length(periods)
    b<-(alpha)/2
    t <- seq(a, 1, length=length(periods))
    ASF <- c(1,1)
    xAlpha <- c(b, b)
    quarters <- periods
    n1 <- totals
    k1 <- events
    scoreciMEM <- function (names, x, n, zalpha, zalphaspend){
        phat <- x/n
        bound <- (zalpha * ((phat * (1 - phat) + (zalpha^2)/(4 * n))/n)^(1/2))/(1 + (zalpha^2)/n)
        boundspend <- (zalphaspend * ((phat * (1 - phat) + (zalphaspend^2)/(4 * n))/n)^(1/2))/(1 + (zalphaspend^2)/n)
        midpnt <- (phat + (zalpha^2)/(2 * n))/(1 + (zalpha^2)/n)
        midpntspend <- (phat + (zalphaspend^2)/(2 * n))/(1 + (zalphaspend^2)/n)
        uplim <- midpnt + bound
        uplimspend <- midpntspend + boundspend
        lowlim <- midpnt - bound
        lowlimspend <- midpntspend - boundspend
        cint <- c(names,x,n,lowlimspend, lowlim, phat, uplim, uplimspend,zalpha,zalphaspend)
        return(cint)
    }
    alpha1 <- bounds(t,n1,iuse=ASF,alpha=xAlpha)
    result1 <- t(mapply(scoreciMEM,quarters,k1,n1,1.96,alpha1$upper.bounds))
    colnames(result1) <- c("Name","k","n","lower_ci_spend0","lowerCI","p","upperCI","upper_ci_spend0","ciz_score0","alpha_ciz_score0")
    output<-as.data.frame(result1)
    output$lower_ci_spend<- format(round(as.numeric(as.character(output$lower_ci_spend0)), digits=4), nsmall=4)
    output$upper_ci_spend<- format(round(as.numeric(as.character(output$upper_ci_spend0)), digits=4), nsmall=4)
    output$ciz_score<- format(round(as.numeric(as.character(output$ciz_score0)), digits=4), nsmall=4)
    output$alpha_ciz_score<- format(round(as.numeric(as.character(output$alpha_ciz_score0)), digits=4), nsmall=4)
    output$period<- rownames(output)    
    output<-output[c("period", "lower_ci_spend", "upper_ci_spend","ciz_score","alpha_ciz_score")]
    return(output)
}

parametric_lc_analysis<-function(
    # Hyperparameters: Data management
    data,
    datasetIdentifier,
    periodIdentifier,
    caseIDFieldNM,
    caseDateFieldNM,
    outcomeFieldNM,
    orderFieldNM,
    operatorFieldNM,
    covariateFieldNMs,
    # Hyperparameters: exposures
    exposureFieldNM,
    exposureOfInterestNM,
    exposureOfInterestOperatorCaseSeriesFieldNM,
    # Hyperparameters: iptw
    allowIPTWWeighting=FALSE,
    normalizeIPTWWeights=FALSE,
    useGeneralCovariateFieldNMs=TRUE,
    iptwCovariateFieldNMs=c(),
    # Modelling
    allowOperatorClustering=FALSE,
    allowGAMBasedAVS=FALSE,
    # Learning effect detection
    learningUnadjustedSignalRequired=TRUE,
    leDetectionAlpha=0.05,
    # Leaning curve estimation
    lcEstimationAlpha=0.05,
    bootstraps=1,
    forceEstimationAsymptoteToZero=TRUE,
    # Intrinsic Device Signal Estimation
    deviceSignalEstimationAlpha=0.05
    
){
    # Perform a parametric learning curve analysis
    # Manage the data
    dm<-parametric.data.management(data, caseIDFieldNM, caseDateFieldNM, outcomeFieldNM, exposureFieldNM, exposureOfInterestNM, orderFieldNM, operatorFieldNM, exposureOfInterestOperatorCaseSeriesFieldNM, covariateFieldNMs, datasetIdentifier, periodIdentifier)
    output<-list(inputData=dm)
    ##Only proceed if data management passes
    if (dm$passed==TRUE){
        # Run IPTW if the user allows treatment weighting, weights will be added to the dataset or NA if not
        data.iptw<-generate.iptw.weights(dm$data, allowIPTWWeighting, normalizeIPTWWeights, useGeneralCovariateFieldNMs, covariateFieldNMs, iptwCovariateFieldNMs)

        # If operator clustering is allowed, create a string to be added to the GAM model formula
        clus<-c()
        if (allowOperatorClustering==T){
            clus<-'s(operatorNM, bs="re")'
        }
        # Let us check if there is a raw signal (no learning adjustment)
        liuFeatures<-c("exposure", covariateFieldNMs, clus)
        data.iptw<- data.iptw %>%
            dplyr::arrange(exposure, operatorNM, operatorSeries) %>%
            as.data.frame()
    
        signal.liu<-signal_detection_gam(target="outcome", features=liuFeatures, data=data.iptw, method="REML", allowAVS=allowGAMBasedAVS, allowWeights=allowIPTWWeighting, alphaLevel=deviceSignalEstimationAlpha)
        output$learningUnadjustedSignal=signal.liu
        
        if (learningUnadjustedSignalRequired==TRUE & signal.liu$hasSignal==FALSE){
            # Let us check if a raw signal exists before we attempt to detect learning effect presence
            output$learningUnadjustedSignalRequired<-TRUE
            output$signal.liu$hasSignal<-FALSE
            output$comment<-"Analysis stopped, learning unadjusted signal not detected..."
            print(output$comment)
            return(output)
        } else {
            # Proceed
            # Learning unadjusted model
            liu<-signal.liu$model

            # Learning adjusted model
            # For cases that don't have the exposure of interest, let us set their operatorSeries value to equal to the max operatorSeries of the exposure of interest
            eoi<-data.iptw %>%
                dplyr::filter(exposure==1) %>%
                dplyr::summarize(
                    maxOS=max(operatorSeries)
                )
            data.iptw$operatorSeries<-ifelse(data.iptw$exposure==0, eoi$maxOS, data.iptw$operatorSeries)
            liaFeatures<-c("exposure", covariateFieldNMs, clus, "s(operatorSeries, bs='ad')") #"s(operatorSeries, k=50, bs='ad', m=5)"
            signal.lia<-signal_detection_gam(target="outcome", features=liaFeatures, data=data.iptw, method="REML", allowAVS=allowGAMBasedAVS, allowWeights=allowIPTWWeighting, alphaLevel=deviceSignalEstimationAlpha)
            output$learningAdjustedSignal=signal.lia
            lia<-signal.lia$model
            output$GAMplot<-""#plot.gam(lia)
            #summary(lia) plot.gam(lia)

            # DETECT LEARNING USING THE LR TEST
            det<-learning_detection(liu, lia, leDetectionAlpha)
            output$learningDetection<-det
            if (det$level!="operator"){
                # No operator learning detected, stop
                output$comment<-"Analysis stopped, operator learning not detected..."
                print(output$comment)
                return(output)
            } else {
                # Proceed now that operator learning effect has been detected
                # ESTIMATE LEARNING RISK
                doi0<-data.iptw %>%
                    dplyr::filter(exposure==1) %>%
                    dplyr::mutate(operatorSeries2=operatorSeries)
                doi0$overallRisk<-predict(lia, type = 'response', newdata=doi0) #summary(doi0$overallRisk)
                doi0$operatorSeries<-eoi$maxOS
                doi0$patientDeviceRisk<-predict(lia, type = 'response', newdata=doi0)
                doi0<-doi0 %>%
                    dplyr::mutate(
                        operatorSeries=operatorSeries2,
                        learningRiskOP=overallRisk-patientDeviceRisk
                    ) %>%
                    dplyr::select(caseID, operatorNM, operatorSeries, overallRisk, patientDeviceRisk, learningRiskOP)
                # ESTIMATE LEARNING FORM
                lc_est<-lc_quantification(data=doi0, lcEstimationAlpha, learningRiskFieldNM='learningRiskOP', bootstraps, forceEstimationAsymptoteToZero)
                output$lcQuantification<-lc_est
                output$lcQuantification$dataUsed<-doi0
                return(output)
            }
        }
    }
    else {
        print(dm$message)
        return(output)
    }
}

derive_alpha_spending_alphas<-function(df, alpha=0.05){
    nAnalyses<-length(unique(df$exposurePeriodFieldNM))
    time <- seq(0.1, 1,length=nAnalyses)
    obf.bd <- ldBounds(t=time, alpha=alpha)
    summary(obf.bd)
    x<-summary(obf.bd)
    x$oalpha
    x$spending
    y<-as.data.frame(x$bounds) %>% 
        dplyr::mutate(
            useAlpha=round(`Nominal Alpha`, 3),
            useAlphaExit=round(`Exit pr.`, 3),
            period=row_number()
        ) %>% 
        dplyr::relocate(period)
    return(y)
}

alphaSpending<-function(alpha=0.05, total, periodicalCount, type="asKD"){
    design <- getDesignGroupSequential(
        sided = 2, 
        alpha = alpha, 
        beta = 0.2,
        informationRates = periodicalCount/total, 
        typeOfDesign = type,
        gammaA = 2.5
    )
    return(list(
        summary=summary(design),
        stageLevels=design$stageLevels,
        alphaSpent=design$alphaSpent
    ))
}

prospective_lc_analysis<-function(){
    # Performs sequential prospective lc analyses
    
}

PLCAnalysis<-function(
        #=====================
        # Part 1: Data management
        #=====================
        data,
        datasetIdentifier="Dataset1",
        caseIDFieldNM,
        caseDateFieldNM,
        outcomeFieldNM,
        orderFieldNM,
        operatorFieldNM,
        covariateFieldNMs,
        #=====================
        # Part 2: Exposure
        #=====================    
        exposureFieldNM,
        exposureOfInterestNM,
        exposureOfInterestOperatorCaseSeriesFieldNM,
        #=====================
        # Part 3: IPTW
        #=====================    
        allowIPTWWeighting=FALSE,
        normalizeIPTWWeights=FALSE,
        useGeneralCovariateFieldNMs=TRUE,
        iptwCovariateFieldNMs=c(),
        #=====================
        # Part 4: GAM Modelling
        #=====================    
        allowOperatorClustering=FALSE,
        allowGAMBasedAVS=FALSE,
        #=====================
        # Part 5: Learning effect detection
        #=====================
        learningUnadjustedSignalRequired=TRUE,
        leDetectionAlpha=0.05,
        #=====================
        # Part 6: Learning curve estimation
        #=====================    
        lcEstimationAlpha=0.05,
        bootstraps=0,
        forceEstimationAsymptoteToZero=TRUE,
        #==========================================
        # Part 7: Intrinsic Device Signal Estimation
        #==========================================    
        deviceSignalEstimationAlpha=0.05,
        #==========================================
        # Part 8: Prospective analysis
        #==========================================   
        allowProspectiveAnalysis=FALSE,
        allowAlphaSpending=TRUE,
        alphaSpendingMethod="asKD",
        alphaToSpend=0.05,
        sequencerFieldNM,
        exposurePeriodFieldNM,
        outcomePeriodFieldNM
){
    if (allowProspectiveAnalysis==FALSE){
        # Calls parametric_lc_analysis once or sequentially if prospective analysis requested.
        results<-parametric_lc_analysis(
            # Hyperparameters: Data management
            data,
            datasetIdentifier,
            "period-001",
            caseIDFieldNM,
            caseDateFieldNM,
            outcomeFieldNM,
            orderFieldNM,
            operatorFieldNM,
            covariateFieldNMs,
            # Hyperparameters: exposures
            exposureFieldNM,
            exposureOfInterestNM,
            exposureOfInterestOperatorCaseSeriesFieldNM,
            # Hyperparameters: iptw
            allowIPTWWeighting,
            normalizeIPTWWeights,
            useGeneralCovariateFieldNMs,
            iptwCovariateFieldNMs,
            # Modelling
            allowOperatorClustering,
            allowGAMBasedAVS,
            # Learning effect detection
            learningUnadjustedSignalRequired,
            leDetectionAlpha,
            # Leaning curve estimation
            lcEstimationAlpha,
            bootstraps,
            forceEstimationAsymptoteToZero,
            # Intrinsic Device Signal Estimation
            deviceSignalEstimationAlpha
        ) 
        return(list(
            allowProspectiveAnalysis=FALSE,
            results=results
        ))
    } else {
        # Prospective processing
        data$exposurePeriodFieldNM<-data[[exposurePeriodFieldNM]] # A period under which the exposure occured.
        data$outcomePeriodFieldNM<-data[[outcomePeriodFieldNM]] # A period under which the outcome occured.
        data$outcomeFieldNM<-data[[outcomeFieldNM]] # Name of outcome field
        
        # Perioidical counts
        pcs<-data %>% 
            dplyr::group_by(exposurePeriodFieldNM) %>% 
            dplyr::summarise(
                N=n()
            ) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(
                exposurePeriodFieldNM=as.integer(exposurePeriodFieldNM)
            ) %>% 
            dplyr::arrange(exposurePeriodFieldNM) %>% 
            dplyr::mutate(
                csN=cumsum(N)
            )
        print(pcs)
        # Determine how many analyses are needed
        asf<-alphaSpending(alpha=alphaToSpend, total=dim(data)[1], periodicalCount=pcs$csN, type=alphaSpendingMethod)
        alphas<-list(alpha=asf$stageLevels)
        alphaDf<-data.frame(alphas)
        alphaDf$id<-row.names(alphaDf) 
        alphaDf<-alphaDf %>% 
            dplyr::mutate(alpha=round(alpha*2, 4)) %>% 
            dplyr::mutate(id=as.integer(id)) %>% 
            dplyr::arrange(id)
        print(alphaDf)
        # If allowAlphaSpending=FALSE, Override asf and set all alphas to alphaToSpend
        if (allowAlphaSpending==FALSE){
            alphaDf<-alphaDf %>% 
                dplyr::mutate(alpha=alphaToSpend)
            print(alphaDf)
        }

        # Process each period
        periods<-alphaDf$id
        combo<-list()
        for (i in periods){
            print("-------------------------------------------------------------------------------")
            print(paste("Analysis period=", i))
            tag<-paste("period", i, sep="_")
            # get alpha
            iasf<-alphaDf %>% 
                dplyr::filter(id==i)
            useAlpha<-iasf$alpha
            
            # Manage data subset
            # Manage periodic data
            dfi<-data %>% 
                dplyr::filter(exposurePeriodFieldNM<=i) %>% 
                dplyr::mutate(
                    outcomeFieldNM=ifelse(outcomePeriodFieldNM>i, 0, outcomeFieldNM) # If outcome is in a future period set to 0
                )
            # Analyse subset
            results<-parametric_lc_analysis(
                # Hyperparameters: Data management
                datasetIdentifier=datasetIdentifier,
                periodIdentifier=tag,
                data=dfi,
                caseIDFieldNM,
                caseDateFieldNM,
                outcomeFieldNM='outcomeFieldNM',
                orderFieldNM,
                operatorFieldNM,
                covariateFieldNMs,
                # Hyperparameters: exposures
                exposureFieldNM,
                exposureOfInterestNM,
                exposureOfInterestOperatorCaseSeriesFieldNM,
                # Hyperparameters: iptw
                allowIPTWWeighting,
                normalizeIPTWWeights,
                useGeneralCovariateFieldNMs,
                iptwCovariateFieldNMs,
                # Modelling
                allowOperatorClustering,
                allowGAMBasedAVS,
                # Learning effect detection
                learningUnadjustedSignalRequired,
                leDetectionAlpha=useAlpha, # useAlpha as an experiement but too stringent for detection.leDetectionAlpha
                # Leaning curve estimation
                lcEstimationAlpha,
                bootstraps,
                forceEstimationAsymptoteToZero,        
                # Intrinsic Device Signal Estimation
                deviceSignalEstimationAlpha=useAlpha
            )
            combo[[tag]]<-results
        }
        return(list(
            allowProspectiveAnalysis=TRUE,
            results=combo
        ))
    }
}

