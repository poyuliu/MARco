
#required packages
require("randomForest")
require("plyr") # for the "arrange" function
require("rfUtilities") # to test model significance
require("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
require("ROCR")
require("pROC")
require("doMC")

#Removing Rare Features
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}



#Transforming Your Data
otu_trans_scale <- function(otu_table_rare_removed){
  otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
  otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)
  otu_table_scaled <- as.data.frame(otu_table_scaled)
  otu_table_scaled
}


#Running Model
#classification or regression
runRF <- function(otu_table_scaled,Y,metadata,model=c("classification","regression"),ntree=501,cores=4,plotALL=FALSE,plotTOP10=FALSE){
  require("randomForest")
  require("plyr") # for the "arrange" function
  require("rfUtilities") # to test model significance
  require("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
  try(if(model!="classification" & model!="regression") stop("Must select an RF model type: 'classification' or 'regression'!"))
  otu_table_scaled_data <- data.frame(t(otu_table_scaled))  
  attach(metadata)
  otu_table_scaled_data$meta <- Y
  
  registerDoMC(cores) #number of cores on the machine
  RF_model <- foreach(y=seq(10), .combine=combine ) %dopar% {
    set.seed(151)
    RF_model_parallel <- randomForest(x=otu_table_scaled_data[,1:(ncol(otu_table_scaled_data)-1)],
                                      y=otu_table_scaled_data[ , ncol(otu_table_scaled_data)],
                                      ntree=ntree, importance=TRUE, proximities=TRUE)
  }
  
  #Permutation Test
  RF_model_sig <- rf.significance(x=RF_model, xdata=otu_table_scaled_data[,1:(ncol(otu_table_scaled_data)-1)],nperm=1000 , ntree=ntree )  
  
  #Accuracy Estimated by Cross-validation
  fit_control <- trainControl( method = "LOOCV" )
  RF_model_loocv <- train(otu_table_scaled_data[,1:(ncol(otu_table_scaled_data)-1)],
                          y=otu_table_scaled_data[, ncol(otu_table_scaled_data)],
                          method="rf", ntree=ntree , tuneGrid=data.frame( mtry=RF_model$mtry ), trControl=fit_control )
  
  #Identifying Important Features
  model <- match.arg(model,)
  if(model=="classification"){
    RF_model_classify_imp <- as.data.frame( RF_model$importance )
    RF_model_classify_imp$features <- rownames( RF_model_classify_imp )
    RF_model_classify_imp_sorted <- arrange( RF_model_classify_imp  , desc(MeanDecreaseAccuracy)  )
    if(plotALL==TRUE & plotTOP10==TRUE){
      par(mfrow=c(1,2))
      barplot(RF_model_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")
      barplot(RF_model_classify_imp_sorted[1:10,"MeanDecreaseAccuracy"], names.arg=RF_model_classify_imp_sorted[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=2, main="Classification RF")  
    } else if(plotALL==TRUE & plotTOP10==FALSE){
      barplot(RF_model_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")
    } else if(plotALL==FALSE & plotTOP10==TRUE){
      barplot(RF_model_classify_imp_sorted[1:10,"MeanDecreaseAccuracy"], names.arg=RF_model_classify_imp_sorted[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=2, main="Classification RF")
    }
    return(list(RFmodel=RF_model, PermuTest=RF_model_sig, LOOCV=RF_model_loocv, SortFeature=RF_model_classify_imp_sorted))
    
  } else if(model=="regression"){
    RF_model_regress_imp <- as.data.frame( RF_model$importance )
    RF_model_regress_imp$features <- rownames( RF_model_regress_imp )
    RF_model_regress_imp_sorted <- arrange( RF_model_regress_imp  , desc(`%IncMSE`)  )
    if(plotALL==TRUE & plotTOP10==TRUE){
      par(mfrow=c(1,2))
      barplot(RF_model_regress_imp_sorted$`%IncMSE`, ylab="% Increase in Mean Squared Error (Variable Importance)", main="RF Regression Variable Importance Distribution")
      barplot(RF_model_regress_imp_sorted[1:10,"%IncMSE"], names.arg=RF_model_regress_imp_sorted[1:10,"features"] , ylab="% Increase in Mean Squared Error (Variable Importance)", las=2, main="Regression RF")  
    } else if(plotALL==TRUE & plotTOP10==FALSE){
      barplot(RF_model_regress_imp_sorted$`%IncMSE`, ylab="% Increase in Mean Squared Error (Variable Importance)", main="RF Regression Variable Importance Distribution")
    } else if(plotALL==FALSE & plotTOP10==TRUE){
      barplot(RF_model_regress_imp_sorted[1:10,"%IncMSE"], names.arg=RF_model_regress_imp_sorted[1:10,"features"] , ylab="% Increase in Mean Squared Error (Variable Importance)", las=2, main="Regression RF")  
    }
    return(list(RFmodel=RF_model, PermuTest=RF_model_sig, LOOCV=RF_model_loocv, SortFeature=RF_model_regress_imp_sorted))
  }
  
}



RF2ROC <- function(RFmodel,plot=FALSE,col="#003366"){
  require("ROCR")
  require("pROC")
  RFmodel <- RFmodel$RFmodel
  Importance <- importance(RFmodel,type = 2)
  predictions=as.vector(RFmodel$votes[,2])
  pred=prediction(predictions,as.numeric(RFmodel$y)-1)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve
  ci95 <- ci.auc(RFmodel$y,RFmodel$votes[,2])
  diagnosis <- table(RFmodel$predicted,RFmodel$y)
  accuracy <- sum(diag(diagnosis))/sum(diagnosis)
  cutff <- which.min(abs(perf_ROC@alpha.values[[1]] - accuracy)) #best cutoff point
  sensitivity <- perf_ROC@y.values[[1]][cutff]
  specificity <- 1-perf_ROC@x.values[[1]][cutff]
  if(plot==TRUE){
    par(las=1,lwd=1.5,mar=c(5.1,4.1,4.1,10.1))
    plot(perf_ROC,type="n",ylab="Sensitivity",xlab="1-Specificity")
    abline(a = 0,b = 1,lwd=1.5,lty=3,col="gray50")
    lines(perf_ROC@x.values[[1]],perf_ROC@y.values[[1]],col=col)
    points(1-specificity,sensitivity,pch=17,col=col)
    legend(x=1.02,y = 1,legend = c(paste0("AUROC = ",round(AUC,2)),
                                   paste0("95% CI ",round(ci95[1],2)," - ",round(ci95[3],2)),
                                   paste0("Sensitivity = ",round(sensitivity,2)),
                                   paste0("Specificity = ",round(specificity,2))),
           xpd=T,bty="n",cex=0.8,text.col = col)
  }
  return(list(AUC=AUC,perf_ROC=perf_ROC,CI95=ci95,Diagnosis=diagnosis,Accuracy=accuracy,Sensitivity=sensitivity,Specificity=specificity))
}

# validation cohort cross validation
VCCV <- function(RFmodel,newdata,newdataGroup,plot=FALSE,cvcol="#bf0808",RFROC,mdcol="#003366"){
  require("ROCR")
  require("pROC")
  RFmodel <- RFmodel$RFmodel
  
  mapID <- match(rownames(RFmodel$importance),rownames(newdata))
  otu.validataion <- newdata[mapID,]
  rownames(otu.validataion) <- rownames(RFmodel$importance)
  otu.validataion[is.na(otu.validataion)] <- 0
  
  CV.validataion <- predict(RFmodel,newdata = t(otu.validataion))
  CV.validataion.vote <- predict(RFmodel,newdata = t(otu.validataion),type = "vote")
  
  pred.CV=prediction(CV.validataion.vote[,2],as.numeric(newdataGroup)-1)
  perf_AUC.CV=performance(pred.CV,"auc") #Calculate the AUC value
  AUC.CV=perf_AUC.CV@y.values[[1]]
  perf_ROC.CV=performance(pred.CV,"tpr","fpr") #plot the actual ROC curve
  
  ci95.CV <- ci.auc(newdataGroup,CV.validataion.vote[,2])
  diagnosis.CV <- table(CV.validataion,newdataGroup)
  accuracy.CV <- sum(diag(diagnosis.CV))/sum(diagnosis.CV)
  cutff.CV <- which.min(abs(perf_ROC.CV@alpha.values[[1]] - accuracy.CV))
  sensitivity.CV <- perf_ROC.CV@y.values[[1]][cutff.CV]
  specificity.CV <- 1-perf_ROC.CV@x.values[[1]][cutff.CV]
  
  if(plot==TRUE){
    par(las=1,lwd=1.5,mar=c(5.1,4.1,4.1,10.1))
    plot(perf_ROC.CV,type="n",ylab="Sensitivity",xlab="1-Specificity")
    abline(a = 0,b = 1,lwd=1.5,lty=3,col="gray50")
    #model AUROC
    lines(RFROC$perf_ROC@x.values[[1]],RFROC$perf_ROC@y.values[[1]],col=mdcol)
    points(1-RFROC$Specificity,RFROC$Sensitivity,pch=17,col=mdcol)
    legend(x=1.02,y = 1,legend = c(paste0("AUROC = ",round(RFROC$AUC,2)),
                                   paste0("95% CI ",round(RFROC$CI95[1],2)," - ",round(RFROC$CI95[3],2)),
                                   paste0("Sensitivity = ",round(RFROC$Sensitivity,2)),
                                   paste0("Specificity = ",round(RFROC$Specificity,2))),
           xpd=T,bty="n",cex=0.8,text.col = mdcol,title = "Model")
    
    #validation AUROC
    lines(perf_ROC.CV@x.values[[1]],perf_ROC.CV@y.values[[1]],col=cvcol,lty=2)
    points(1-specificity.CV,sensitivity.CV,pch=2,col=cvcol)
    legend(x=1.02,y = 0.3,legend = c(paste0("AUROC = ",round(AUC.CV,2)),
                                   paste0("95% CI ",round(ci95.CV[1],2)," - ",round(ci95.CV[3],2)),
                                   paste0("Sensitivity = ",round(sensitivity.CV,2)),
                                   paste0("Specificity = ",round(specificity.CV,2))),
           xpd=T,bty="n",cex=0.8,text.col = cvcol,title = "Validation Dataset")
  }
  
  return(list(AUC=AUC.CV,perf_ROC=perf_ROC.CV,CI95=ci95.CV,Diagnosis=diagnosis.CV,Accuracy=accuracy.CV,Sensitivity=sensitivity.CV,Specificity=specificity.CV))
}
