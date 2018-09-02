library(ROCR)
###########################################################################
#################### Compute AUC when Y ~ 1
###########################################################################
# resp: response variable
# cv_dat: data for computing AUC of null model
nulAUC <- function(resp,cv_dat){
  f=as.formula(paste(resp,"~ 1"))
  res <- c()
  # Use Leave-one-out cross validation (LOOCV)
  for (i in 1:nrow(cv_dat)){
    test.dat <- cv_dat[i,]
    train.dat <- cv_dat[-i,]
    model = glm(f, data = train.dat, family = "binomial")
    pred_tr = predict(model, train.dat, type="response")
    trauc <- performance(prediction(pred_tr,train.dat[,match(resp,colnames(train.dat))]),"auc")@y.values[[1]][1]
    res <- c(res, trauc)
  }
  initrauc <- mean(res);initeauc <- 0.5
  nulres <-matrix(c(initrauc, initeauc), nrow = 1)
  return(nulres)
}

###########################################################################
#################### Forward step with LOOCV
###########################################################################
# cv_dat: data for forward step
# com: indices of variables that can be included in the existing model
# resp: response variable
# var: variables included in the existing model
forward_step <- function(cv_dat,com,resp,var){
  forward_ls <- NULL
  for (co in com){
    res <- c();pred <- c()
    f=as.formula(paste(resp,"~",paste(c(var,colnames(cv_dat)[co]),collapse = " + ")))
    for (i in 1:nrow(cv_dat)){
      test.dat <- cv_dat[i,]
      train.dat <- cv_dat[-i,]
      model = glm(f, data = train.dat, family = "binomial")
      pred_tr = predict(model, train.dat, type="response")
      trauc <- performance(prediction(pred_tr,train.dat[,match(resp,colnames(train.dat))]),"auc")@y.values[[1]][1]
      pred <- c(pred, predict(model, test.dat, type="response"))
      res <- c(res, trauc)
    }
    auc.value <- performance(prediction(pred,cv_dat[,match(resp,colnames(cv_dat))]),"auc")@y.values[[1]][1]
    tmp <-c(paste(c(var,colnames(train.dat)[co]),collapse = " + "),mean(res),auc.value)
    forward_ls <- rbind(forward_ls, tmp)
  }
  return(forward_ls)
}

###########################################################################
#################### Backward step with LOOCV
###########################################################################
# var_ls: indices of variables included in the backward model
# cv_dat: data for forward step
# resp: response variable
# backward_ls: NULL list
backward_step <- function(var_ls, cv_dat, resp, backward_ls){
  for (u in 1:nrow(var_ls)){
    f=as.formula(paste(resp,"~",paste(colnames(cv_dat)[var_ls[u,]],collapse = ' + ')))
    pred <- c();res <- c()
    for (i in 1:nrow(cv_dat)){
      test.dat <- cv_dat[i,]
      train.dat <- cv_dat[-i,]
      model = glm(f, data = train.dat, family = "binomial")
      pred_tr = predict(model, train.dat, type="response")
      trauc <- performance(prediction(pred_tr,train.dat[,match(resp,colnames(train.dat))]),"auc")@y.values[[1]][1]
      pred <- c(pred, predict(model, test.dat, type="response"))
      res <- c(res, trauc)
    }
    auc.value <- performance(prediction(pred,cv_dat[,match(resp,colnames(cv_dat))]),"auc")@y.values[[1]][1]
    tmp <-c(paste(colnames(cv_dat)[var_ls[u,]],collapse = ' + '),mean(res),auc.value)
    backward_ls <- rbind(backward_ls, tmp)
  }
  return(backward_ls)
}

###########################################################################
#################### test AUC-based stepwise selection
###########################################################################
# totvar: index of independent variables used for stepwise selection (in dataset)
# cv_dat: dataset for variable selection
# resp: response variable
# cutoff: cutoff to select the best model between updated model and existing model
# addres1: address to save the output file
testAUCstepwise <- function(totvar, cv_dat, resp, cutoff, addres1){
  imtres <- NULL;inires <-nulAUC(resp,cv_dat)
  # inires: initial AUC of null model
  var <- '';trauc <- inires[1];teauc <- inires[2]
  while (length(setdiff(strsplit(var[[1]], split=c(" "))[[1]],'+'))<=length(totvar)){
    com <- combn(totvar,1)
    ## get list of variables not included in the model if the number of variables > 1 in the model
    if(var!="") com <- com[-which(com %in% which(colnames(cv_dat)%in%setdiff(strsplit(var[[1]], split=c(" "))[[1]],'+')))]

    ## forward step
    forward_ls <- forward_step(cv_dat,com,resp,var)
    
    # select forward model with maximum test AUC
    forward.teauc1<-max(forward_ls[,3])
    forward.var1 <- forward_ls[which.max(forward_ls[,3]),1];forward.trauc1 <-forward_ls[which.max(forward_ls[,3]),2]
    forward.newstep <-matrix(c(forward.var1,forward.trauc1,forward.teauc1),nrow=1)
    
    # For forward step, select updated model only if updated AUC > existing AUC + cutoff
    if (forward.newstep[3] > (as.numeric(teauc) + as.numeric(cutoff))){
      imtres <- rbind(imtres, forward.newstep)
      colnames(forward.newstep)<-c('Variable','trainAUC','testAUC')
      eval(parse(text = paste("write.csv(forward.newstep,'./",addres1,"/intermediate_var_forward",nrow(imtres),".csv',row.names = F)",sep = "")))
      var <- forward.var1
      trauc <- forward.trauc1
      teauc <- forward.teauc1
    } else{
      mat<-matrix(imtres[nrow(imtres),],nrow=1)
      colnames(mat)<-c('Variable','trainAUC','testAUC')
      colnames(imtres)<-c('Variable','trainAUC','testAUC')
      eval(parse(text = paste("write.csv(imtres,'./",addres1,"/intermediate_var_stepwise.csv',row.names = F)",sep = "")))
      break
    }
    
    backward_ls <- NULL
    # var_list: variables included in the selected model
    var_list <- setdiff(strsplit(var[[1]], split=c(" "))[[1]],'+')[which(setdiff(strsplit(var[[1]], split=c(" "))[[1]],'+')!="")]
    # get indices of variables included in the backward model if the number of variables > 1
    if (length(var_list)>1){
      var_ls <- NULL
      for (k in 1:length(var_list)){
        var_ls <- rbind(var_ls, c(match(var_list[-k],colnames(cv_dat))))
      }
    } else{
      var_ls <- 1
    }
    ## compare the model with null model when the number of variables = 1 in the selected model
    if (!(is.matrix(var_ls))){
      f=as.formula(paste(resp,"~ 1"))
      res <- c()
      for (i in 1:nrow(cv_dat)){
        test.dat <- cv_dat[i,]
        train.dat <- cv_dat[-i,]
        model = glm(f, data = train.dat, family = "binomial")
        pred_tr = predict(model, train.dat, type="response")
        trauc <- performance(prediction(pred_tr,train.dat[,match(resp,colnames(train.dat))]),"auc")@y.values[[1]][1]
        res <- c(res, trauc)
      }
      auc.value <- 0.5
      tmp <-c(" ",mean(res),auc.value)
      backward_ls <- rbind(backward_ls, tmp)
    } else{
      ## backward step
      backward_ls <- backward_step(var_ls, cv_dat, resp, backward_ls)
    }
    # select backward model with maximum test AUC among multiple backward models
    backward.teauc1<-max(backward_ls[,3])
    backward.var1 <- backward_ls[which.max(backward_ls[,3]),1];backward.trauc1 <-backward_ls[which.max(backward_ls[,3]),2]
    backward.newstep <-matrix(c(backward.var1,backward.trauc1,backward.teauc1),nrow=1)
    
    # For backward step, select updated model only if updated AUC > existing AUC - cutoff
    if (backward.newstep[3] > (as.numeric(teauc)-cutoff)){
      imtres <- rbind(imtres, backward.newstep)
      colnames(backward.newstep)<-c('Variable','trainAUC','testAUC')
      eval(parse(text = paste("write.csv(backward.newstep,'./",addres1,"/intermediate_var_backward",nrow(imtres),".csv',row.names = F)",sep = "")))
      var <- backward.var1
      trauc <- backward.trauc1
      teauc <- backward.teauc1
    } 
  }
  return(mat)
}

###########################################################################
#################### main function to execute test AUC based stepwise
###########################################################################
main <- function(input, cutoff, exp_var, resp, addres){
  ## totvar: index of independent variables used for stepwise selection (in dataset)
  ifelse(is.null(exp_var),totvar <- seq(ncol(input))[-match(resp,colnames(input))],totvar <- seq(ncol(input))[-match(c(exp_var,resp),colnames(input))])
  addres1 <- paste0(addres,cutoff)
  total_res <-testAUCstepwise(totvar, input, resp, cutoff, addres1)
  eval(parse(text = paste("write.csv(total_res,'./",addres1,"/stepwise_result.csv',row.names = F)",sep = "")))
}

