# This file is licensed under GPLv2 (see LICENSE for further information)
# This file contains adapted versions of functions present in the mpm package, written by Luc Wouters
# See https://cran.r-project.org/web/packages/mpm/index.html for the original library
require(randomForest)

#do not add spacing to the following line
MAGISTA_dir="<<<MAGISTA_PATH>>>"


# The following functions are adapted from the package mpm
mpm.logtransf = function(data, logrepl = 1e-09){
    logrepl <- max(1e-09, logrepl)
    if (any(data <= 0, na.rm = TRUE)) {
        #warning(paste("Non-positive data replaced by", logrepl,
        #              "in computing mpm.logtransf"),
        #        call. = FALSE)
        data[data <= 0] <- logrepl
    }
    return(log(data))
}
mpm.pca_axes = function(training_data){
    ncol1 = ncol(training_data)
    complete_cases = complete.cases(training_data)
    training_data = training_data[complete_cases,]

    #remove variables which do not vary
    tmpData <- as.matrix(training_data)
    datavar = diag(var(tmpData))
    to_keep = (abs(datavar)>1e-09)
    training_data = training_data[,names(to_keep)]


    # remove variables which are perfectly correlated to eahc other
    cr = cor(training_data[,2:ncol(training_data)])
    cr[is.na(cr)]=0
    redundant = logical(ncol(training_data))
    for(currow in 1:(ncol(training_data)-2)){
        for(curcol in (currow+1):(ncol(training_data)-1)){
            if(abs(cr[currow,curcol])==1){
                redundant[curcol+1]=TRUE
            }
        }
    }
    training_data = training_data[,!redundant]

    if(ncol1!=ncol(training_data)){
        warning(paste("Number of columns reduced from",as.character(ncol1),"to",as.character(ncol(training_data))), call. = F)
    }

    mpm_result = mpm.2(training_data, center = "column", normal = "column")

    origdata_N = as.matrix(training_data[,mpm_result$col.names])
    origdata_log = mpm.logtransf(origdata_N,1e-09)
    colmeans = colMeans(origdata_log)
    origdata_center = origdata_log
    for(i in 1:ncol(origdata_N)){
        origdata_center[,i] = origdata_log[,i] - colmeans[i]
    }
    origdata_var = colSums(origdata_center^2)/nrow(origdata_center)
    normfactor = sqrt(origdata_var)
    origdata_norm = origdata_center
    for(i in 1:ncol(origdata_N)){
        origdata_norm[,i] = origdata_center[,i] / normfactor[i]
    }
    denom = sqrt(nrow(origdata_N))*sqrt(ncol(origdata_N))
    return(list(vt=mpm_result$SVD$vt,
                d =mpm_result$SVD$d,
                denom=denom,
                normfactor=normfactor,
                colmeans = colmeans,
                columns=mpm_result$col.names))
}
mpm.pca_predict = function(pca_axes, target_data){
    x_N = as.matrix(target_data[,pca_axes$columns])
    x_log = mpm.logtransf(x_N,1e-09)
    x_center = x_log
    for(i in 1:ncol(x_log)){
        x_center[,i] = x_log[,i] - pca_axes$colmeans[i]
    }
    x_norm = x_center
    for(i in 1:ncol(x_norm)){
        x_norm[,i] = x_center[,i] / pca_axes$normfactor[i]
    }
    actualD = diag(1.0/pca_axes$d)
    WData = x_norm / pca_axes$denom
    result = WData %*% t(pca_axes$vt) %*% actualD
    result =as.data.frame(result)
    for(f in 1:ncol(result)){
        colnames(result)[f] = paste("V",as.character(f),sep="")
    }

    return(result)
}
mpm.2 = function (data, logtrans = TRUE, logrepl = 1e-09,
                  center = c("double","row", "column", "global", "none"),
                  normal = c("global","row", "column", "none"),
                  closure = c("none", "row", "column", "global", "double"),
                  row.weight = c("constant", "mean", "median", "max", "logmean", "RW"),
                  col.weight = c("constant", "mean", "median", "max", "logmean", "CW"),
                  CW = NULL,
                  RW = NULL,
                  pos.row = NULL,
                  pos.column = NULL,
                  samples = c("cols","rows"))
{
    if (missing(data))
        stop("Argument \"data\" is missing, with no default")
    NData <- as.matrix(data[, -1])
    datavar = diag(var(NData))
    logrepl <- max(1e-09, logrepl)
    to_keep = (abs(datavar)>logrepl)
    if(length(to_keep) != ncol(NData)){
        warning("At least one column had a variance below 1e-9 and was removed")
    }
    NData = NData[,to_keep]
    if(is.null(pos.row)){
        pos.row = rep(FALSE, nrow(NData))
    }
    if(is.null(pos.column)){
        pos.column = rep(FALSE, ncol(NData))
    }
    if(is.null(RW)){
        RW = rep(1, nrow(data))
    }
    if(is.null(CW)){
        CW = rep(1, ncol(NData))
    }

    if (any(is.na(NData)) || !is.numeric(NData))
        stop("Data must be numeric without NA's")
    if (length(pos.row) != dim(NData)[1])
        stop("Length of pos.row argument not equal to number of rows in table")
    if (length(pos.column) != dim(NData)[2])
        stop(paste("Length of pos.column argument (",
                   as.character(length(pos.column)),
                   ") not equal to number of columns in table (",
                   as.character(dim(NData)[2])))
    if (length(RW) != nrow(NData))
        stop("Length of RW not equal to number of rows in table")
    if (length(CW) != ncol(NData))
        stop("Length of CW not equal to number of columns in table")
    center <- match.arg(center)
    normal <- match.arg(normal)
    closure <- match.arg(closure)
    row.weight <- match.arg(row.weight)
    col.weight <- match.arg(col.weight)
    samples <- match.arg(samples)
    Row.Names <- as.character(data[, 1])
    Col.names <- colnames(NData)
    if(is.null(Col.names)){
        Col.names <- as.character(1:ncol(NData))
    }

    RData <- NData
    RData[pos.row, ] <- NA
    RData[, pos.column] <- NA
    logtransf <- function(x, logrepl, comp) {
        if (any(x <= 0, na.rm = TRUE)) {
            warning(paste("Non-positive data replaced by", logrepl,
                          "in computing", comp, "in: spectralmap."),
                    call. = FALSE)
            x[x <= 0] <- logrepl
        }
        return(log(x))
    }
    LData <- if (logtrans)
        logtransf(NData, logrepl, comp = "logarithms")
    else NData
    RM <- rowMeans(NData[, !pos.column])
    CM <- colMeans(NData[!pos.row, ])
    Wn <- pmax(0, switch(row.weight,
                         constant = rep(1, length(RM)),
                         mean = apply(RData, 1, mean, na.rm = TRUE),
                         median = apply(RData, 1, median, na.rm = TRUE),
                         max = apply(RData, 1, max, na.rm = TRUE),
                         logmean = apply(logtransf(RData, logrepl, comp = "logmean weights"),1, mean, na.rm = TRUE),
                         RW = RW))
    Wp <- pmax(0, switch(col.weight,
                         constant = rep(1, length(CM)),
                         mean = apply(RData, 2, mean, na.rm = TRUE),
                         median = apply(RData, 2, median, na.rm = TRUE),
                         max = apply(RData, 2, max, na.rm = TRUE),
                         logmean = apply(logtransf(RData, logrepl, comp = "logmean weights"), 2, mean, na.rm = TRUE),
                         CW = CW))
    Wn[pos.row] <- 0
    Wp[pos.column] <- 0
    Wn <- Wn/sum(Wn)
    Wp <- Wp/sum(Wp)
    if (closure != "none" && any(LData < 0))
        warning("Closure operation with non-positive data")
    Tn <- rowSums(LData[, !pos.column], na.rm = TRUE)
    Tp <- colSums(LData[!pos.row, ], na.rm = TRUE)
    Tt <- sum(LData[!pos.row, !pos.column], na.rm = TRUE)
    ClData <- switch(closure,
                     none = LData,
                     row = sweep(LData, 1, Tn, "/"),
                     column = sweep(LData, 2, Tp, "/"),
                     global = LData/Tt,
                     double = Tt * sweep(sweep(LData, 1, Tn, "/"), 2, Tp, "/"))
    if (any(!is.finite(ClData)))
        stop("Division by 0 in closure operation")
    Mp <- colSums(sweep(ClData, 1, Wn, "*"))
    Mn <- rowSums(sweep(ClData, 2, Wp, "*"))
    Mg <- sum(Mp * Wp)
    CData <- switch(center,
                    double = Mg + sweep(sweep(ClData,2, Mp), 1, Mn),
                    row = sweep(ClData, 1, Mn),
                    column = sweep(ClData, 2, Mp),
                    global = ClData - Mg,
                    none = ClData)
    Vp <- colSums(sweep(CData^2, 1, Wn, "*"))
    Vn <- rowSums(sweep(CData^2, 2, Wp, "*"))
    Vg <- sum(Vp * Wp)
    SData <- switch(normal,
                    global = CData/sqrt(Vg),
                    row = sweep(CData, 1, sqrt(Vn), "/"),
                    column = sweep(CData, 2, sqrt(Vp),"/"),
                    none = CData)
    WData <- sweep(sweep(SData, 1, sqrt(Wn), "*"), 2, sqrt(Wp),"*")
    svd.res <- La.svd(WData)
    eigen <- svd.res$d^2
    contrib <- eigen/sum(eigen)
    r <- list(TData = SData, row.names = Row.Names, col.names = Col.names,
              closure = closure, center = center, normal = normal, logtrans = logtrans, logrepl = logrepl,
              row.weight = row.weight, col.weight = col.weight, eigen = eigen,
              contrib = contrib, Rm = RM, Cm = CM, Wn = Wn, Wp = Wp, Mp = Mp, Mn = Mn, Vn = Vn, Vp=Vp, CW=CW, RW=RW,Tn=Tn, Tp=Tp, Tt=Tt,
              SVD = svd.res, pos.column = pos.column, pos.row = pos.row, samples=samples,
              call = match.call(), CData=CData, LData=LData, NData=NData, WData=WData)
    class(r) <- "mpm"
    return(r)
}


#Utility functions
colnames_containing_pattern = function(data,pattern,fixed=F){
    result = (1:ncol(data))*0
    result = as.logical(result)
    for(i in 1:ncol(data)){
        result[i] = (length(grep(pattern,colnames(data)[i],fixed = fixed))>0)
    }
    return(colnames(data)[which(result)])
}
colnames_without_pattern = function(data,pattern,fixed=F){
    result = (1:ncol(data))*0
    result = as.logical(result)
    for(i in 1:ncol(data)){
        result[i] = (length(grep(pattern,colnames(data)[i],fixed = fixed))==0)
    }
    return(colnames(data)[which(result)])
}
numeric_cols = function(data){
    result = (1:ncol(data))*0
    result = as.logical(result)
    for(i in 1:ncol(data)){
        result[i] = is.numeric(data[,i])
    }
    return(colnames(data)[which(result)])
}
cv.RFpredict = function(formula, data, folds, subgroup=NULL, ...){
    input_order = sample(nrow(data))
    prediction = data.frame(value=numeric(length=nrow(data)),
                            times_predicted=numeric(length=nrow(data)),
                            fold_id=numeric(length=nrow(data)))
    if(is.null(subgroup)){
        rows_per_subset = ceiling(nrow(data)/folds)
        cat("|",file=stderr())
        for(i in 1:nfolds){
            cat(" ",file=stderr())
        }
        cat("|\r|",file=stderr())
        for(i in 1:folds){
            cat(".",file=stderr())
            selected_rows = logical(length=nrow(data))
            selected_rows[(rows_per_subset*(i-1)+1):min((rows_per_subset*(i)),nrow(data))] = TRUE
            test_subset = data[selected_rows,]
            train_subset = data[!selected_rows,]
            tmpmodel = randomForest(formula,train_subset, ...)
            tmpresult = predict(tmpmodel, test_subset)
            prediction$value[selected_rows] = tmpresult
            prediction$times_predicted[selected_rows] = prediction$times_predicted[selected_rows]+1
        }
    } else {
        usg=unique(subgroup)
        if(length(usg)>1){
            usg=sample(usg)
        }
        nsubgroups = length(usg)
        nsubgroup_splits = min(folds,nsubgroups)
        ignored_groups=max(c( (1.0*nsubgroups)/(1.0*nsubgroup_splits),1.0))
        x=0
        folds=1 #reduce folds to 1
        x = nsubgroup_splits*folds
        cat("|",file=stderr())
        for(i in 1:x){
            cat(" ",file=stderr())
        }
        cat("|\r|",file=stderr())
        rows_per_subset = ceiling(nrow(data)/folds)
        for(g in 1:nsubgroup_splits){
            #the mask ensures that nsubgroup_splits groups are completely during each round of cross-validation
            mask = logical(length=nrow(data))
            first_test_group = as.integer(floor((g-1)*ignored_groups)+1)
            last_test_group = as.integer(floor(g*ignored_groups))
            for(j in first_test_group:last_test_group){
                mask[subgroup==usg[j]]=TRUE
            }
            prediction$fold_id[mask]=g
            for(i in 1:folds){
                selected_rows = logical(length=nrow(data))
                if(folds==1){
                    selected_rows[1:nrow(data)] = FALSE
                } else {
                    selected_rows[(rows_per_subset*(i-1)+1):min((rows_per_subset*(i)),nrow(data))] = TRUE
                }
                selected_rows[mask] = TRUE
                test_subset = data[selected_rows,]
                train_subset = data[!selected_rows,]
                tmpmodel = randomForest(formula,train_subset, ...)
                tmpresult = predict(tmpmodel, test_subset)
                prediction$value[selected_rows] = prediction$value[selected_rows] + tmpresult
                prediction$times_predicted[selected_rows] = prediction$times_predicted[selected_rows]+1
                cat("#",file=stderr())
            }

        }
        prediction$value = prediction$value / prediction$times_predicted
    }
    cat("\n",file=stderr())
    prediction = data.frame(value = prediction[,c("value")], fold_id = as.factor(prediction[,c("fold_id")]))
    return(prediction)
}

MAGISTA_error_plot = function(data,x,y,color, ...){
    #textfamily="Arial"
    textsize=20
    data = data[complete.cases(data[,c(x,y)]),]
    g = mean(data[complete.cases(data),y])
    tmp_cor = cor(x=data[,x],y=data[,y],use = "complete.obs")
    delta = data[,x]-data[,y]
    SSE = sum(delta*delta, na.rm = T)
    AE = sqrt(SSE/nrow(data))
    R2 = table_rsquared(data[,y],data[,x],"model")
    ex = paste("atop({{R^2}[x%~%y]}==",as.character(floor(R2*1000)/1000),",",
               "{RMSE==",as.character(floor(AE*100)/100),"})",
               sep="")

    minx = 0
    maxx = 100
    miny = 0
    maxy = 100
    labelpos_x = minx+(maxx-minx)*0.7
    labelpos_y = minx+(maxx-minx)*0.15
    textdf = data.frame(x=c(labelpos_x),y=c(labelpos_y),t=c(ex))

    return(
            ggplot(data)+
                geom_point(aes_string(x=x,y=y, color=color))+geom_abline(intercept=0,slope=1) +
                geom_smooth(aes_string(x=x,y=y),method="lm", formula = y~x)+
                geom_text(data=textdf,aes_string(x="x",y="y",label="t"),parse = T,size = textsize*0.35) +
                xlim(c(0,100)) + ylim(c(0,100)) +
                theme_bw()
        )

}

# MAGISTA core functions
train_MAGISTA = function(data, predictor_groups, targets, variables_to_ignore=c(), random_seed=123456, error_plot_prefix=NULL){
    set.seed(random_seed)
    result = list()
    gnames = names(predictor_groups)
    allbin_group=NULL
    if("bin_group" %in% colnames(data)){
        allbin_group = data$bin_group
    }
    dataset_color=character(length=nrow(data))
    if("dataset" %in% colnames(data)){
        dataset_color=data$dataset
    }
    if(is.null(gnames)){
        gnames = paste("Subset",1:length(predictor_groups),sep=".")
    }
    for(cn in colnames(data)){
        if(length(grep("Rating",cn))>0){
            for(v in unique(data[,cn])){
                data[,paste(cn,as.character(v),sep=".")]=(data[,cn]==v)*1.0
            }
        }
        if(length(grep("bin.part",cn))>0){
            if("bin.size" %in% colnames(data)){
                if(abs(max(data[,cn],na.rm=T))>1){
                    data[,cn] = data[,cn]/data$bin.size
                }
            }
        }
    }
    vargroup_rating = numeric(length = length(predictor_groups))
    names(vargroup_rating) = gnames
    for(group in gnames){
        curcols = predictor_groups[[group]]
        vargroup_rating[group] = sum(complete.cases(data[,curcols]))
    }
    possible_ratings = sort(unique(vargroup_rating))
    res = list()
    error_plot=F
    if(!is.null(error_plot_prefix)){
        error_plot=T
        cat("error plot will be printed\n")
        error_plot_prefix=as.character(error_plot_prefix)
    }
    for(target in targets){
        curmodel = list()
        bin_rating = character(length=nrow(data))
        for(rating in possible_ratings){
            cat(paste("Generating submodel : ",target, " [",rating,"]\n",sep=""),file=stderr())
            tcols = names(vargroup_rating)[vargroup_rating>=rating]
            variables = unique(unlist(predictor_groups[tcols]))
            for(badvar  in variables_to_ignore){
                variables=variables[variables!=badvar]
            }
            indata = data[,variables]
            rows2keep = complete.cases(indata)
            indata = indata[rows2keep,]
            bin_rating[as.logical(rows2keep*(bin_rating==""))]=rating
            bin_group=NULL
            if(!is.null(allbin_group)){
                bin_group=allbin_group[rows2keep]
            }
            transformation = mpm.pca_axes(indata)
            train.fit = mpm.pca_predict(transformation, indata)
            indata = cbind(indata, train.fit)
            indata = cbind(data[rows2keep,target],indata)
            colnames(indata)[1] = target
            rfm = randomForest(as.formula(paste(target,"~",".")),data = indata, na.action = na.omit, ntree = 200, maxnodes = 60, importance=TRUE, localImp = TRUE)

            rf_model_cv = cv.RFpredict(as.formula(paste(target,"~",".")), data=indata, folds=10, na.action = na.omit, importance=T, subgroup=bin_group)
            rf_model_cv$actual = data[rows2keep, target]
            indata$predictionerror = rf_model_cv$actual - rf_model_cv$value

            if(error_plot){
                plot_data = data.frame(
                    actual = rf_model_cv$actual,
                    value = rf_model_cv$value,
                    cv.pass = rf_model_cv$fold_id
                )
                ggsave(paste(error_plot_prefix,"err_profile",target,rating,"pdf",sep='.'),
                       MAGISTA_error_plot(plot_data,x="actual",y="value",color="cv.pass")
                )
            }

            errpred_model = randomForest(as.formula(paste("predictionerror ~ . -",target)),data = indata, na.action = na.omit, ntree = 200, maxnodes = 60, importance=TRUE, localImp = TRUE)
            #rf_model_cv$error = predict(errpred_model,data=indata)

            #rf_correction = lm( actual ~ poly(value,3), rf_model_cv, na.action = na.omit)
            rf_correction = lm( actual ~ value, rf_model_cv, na.action = na.omit)
            #rf_correction = lm( actual ~ value * error, rf_model_cv, na.action = na.omit)

            submodel = list(name=paste("R",rating,sep=""),
                            rating=rating,
                            incols=variables,
                            resname=target,
                            pcaxes=transformation,
                            uncertainity_model=errpred_model,
                            correction=rf_correction,
                            model=rfm,
                            valuerange=c(min(data[rows2keep,target]),max(data[rows2keep,target])))
            curmodel[[submodel$name]]=submodel
        }
        res[[target]] = curmodel
    }
    return(res)
}

predict_MAGISTA = function(model, data){
    tvars = names(model)
    result = as.data.frame(matrix(nrow = nrow(data),ncol=length(model)+1))
    for(cn in colnames(data)){
        if(length(grep("Rating",cn))>0){
            for(v in unique(data[,cn])){
                data[,paste(cn,as.character(v),sep=".")]=(data[,cn]==v)*1.0
            }
        }
        if(length(grep("bin.part",cn))>0){
            if("bin.size" %in% colnames(data)){
                if(abs(max(data[,cn],na.rm=T))>1){
                    data[,cn] = data[,cn]/data$bin.size
                }
            }
        }
    }
    colnames(result)[1:length(model)] = names(model)
    colnames(result)[length(model)+1] = "rating"
    for(target in colnames(result)[1:length(model)]){
        subnames = names(model[[target]])
        ratings = numeric(length(subnames))
        names(ratings) = subnames
        for(cname in subnames){
            ratings[cname] = model[[target]][[cname]]$rating
        }
        ratings = sort(ratings)
        result[,target] = NA
        for(rating in ratings){
            torun = is.na(result[,target])
            if(sum(torun)>0){
                result$rating[torun] = rating
                indata = data[torun,]
                cname = names(ratings)[which(ratings==rating)]
                submodel = model[[target]][[cname]]
                missing_cols = !(model[[target]][[cname]]$pcaxes$columns %in% colnames(indata))
                if(any(missing_cols)){
                    stop(paste("missing columns:",paste(model[[target]][[cname]]$pcaxes$columns[missing_cols],collapse = ",")))
                }
                train.fit = mpm.pca_predict(submodel$pcaxes, indata)
                indata = cbind(indata, train.fit)

                missing_cols = !(rownames(model[[target]][[cname]]$model$importance) %in% colnames(indata))
                if(any(missing_cols)){
                    stop(paste("missing columns:",paste(model[[target]][[cname]]$uncertainity_model$columns[missing_cols],collapse = ",")))
                }

                missing_cols = !(rownames(model[[target]][[cname]]$model$importance) %in% colnames(indata))
                resvalues = predict(submodel$model, indata)
                if(cname %in% colnames(indata)){
                    indata$predictionerror = predict(submodel$uncertainity_model,indata)
                } else {
                    indata$predictionerror = NaN
                }
                df = data.frame(value = resvalues,error=indata$predictionerror)
                resvalues = predict(submodel$correction, df)
                resvalues[resvalues < submodel$valuerange[1]] = submodel$valuerange[1]
                resvalues[resvalues > submodel$valuerange[2]] = submodel$valuerange[2]
                result[torun,target] = resvalues
            }
        }
    }
    return(result)
}


#Main Script
MAGISTA_prepare_generic = function(modeltype, training_dataname,name){
    training_data = read.csv(training_dataname,header=T)
    targetvars = c("purity","completeness","F1","contamination")
    cols_desc_id = c("Bin.Id","dataset","method")
    cols_GC = colnames_containing_pattern(training_data,"GC.",fixed=T)
    cols_GC = cols_GC[cols_GC!="GC.bin.part"]
    cols_f4c = colnames_containing_pattern(training_data,"w5k_f4c.",fixed=T)
    cols_m3c = colnames_containing_pattern(training_data,"w10k_m3c.",fixed=T)
    cols_m4c = colnames_containing_pattern(training_data,"w10k_m4c.",fixed=T)
    cols_k4s = colnames_containing_pattern(training_data,"w50k_k4s.",fixed=T)
    if (modeltype=="MAGISTIC"){
        #include the depth to which checkM has reached as a usable value
        training_data$checkm.SpecificityRating=8-as.integer(factor(gsub("root.*","d",gsub("__.*","",training_data$checkm.Marker.lineage)),levels=c("s","g","f","o","c","p","k","d")))
        cols_checkm = colnames_containing_pattern(training_data,"checkm.",fixed=T)
        cols_checkm = cols_checkm[ sapply(training_data[1,cols_checkm],is.numeric)]
    }
    cols_dist = c(cols_GC,cols_f4c,cols_m3c, cols_m4c, cols_k4s)

    colsets= list(GC=cols_GC,
                  f4c=cols_f4c,
                  m3c=cols_m3c,
                  m4c=cols_m4c,
                  k4s=cols_k4s)

    for(f in names(colsets)){
        colsets[[f]]=numeric_cols(training_data[,colsets[[f]]])
        colsets[[f]]=c("bin.size",colsets[[f]])
        colsets[[f]]=c("Bin.Id",colsets[[f]])
    }
    if(modeltype=="MAGISTIC"){
        colsets$checkm = cols_checkm
        MAGISTIC = train_MAGISTA(training_data,colsets,targetvars,"Bin.Id")
    }
    if(modeltype=="MAGISTA"){
        MAGISTA = train_MAGISTA(training_data,colsets,targetvars,"Bin.Id")
    }
    save(list=c(modeltype),file=paste(MAGISTA_dir,"data/",modeltype,".",name,".RData",sep=""))
}
MAGISTA_prepare = function(training_dataname,name="default"){
    MAGISTA_prepare_generic("MAGISTA", training_dataname, name)
}
MAGISTIC_prepare = function(training_dataname,name="default"){
    MAGISTA_prepare_generic("MAGISTIC", training_dataname, name)
}
MAGISTA_testrun = function(){
    require(moments)
    require(dplyr)
    cat("Beginning Test Run\n", file=stdout())
    cat(paste("Temporarily changed working directory to: ", MAGISTA_dir,"data/\n",sep=""),file=stderr())
    setwd(paste(MAGISTA_dir,"/data/",sep=""))
    if(!file.exists("MAGISTA.default.RData")){
        MAGISTA_prepare("training_dataset.csv",name="default")
        cat("MAGISTA Model ready\n",file=stderr())
    }
    if(!file.exists("MAGISTIC.default.RData")){
        MAGISTIC_prepare("training_dataset.csv",name="default")
        cat("MAGISTIC Model ready\n",file=stderr())
    }
    load("MAGISTA.default.RData")
    load("MAGISTIC.default.RData")
    test_data = read.csv("test_dataset.csv",header=T)
    result_MAGISTA = predict_MAGISTA(MAGISTA, test_data)
    result_MAGISTIC = predict_MAGISTA(MAGISTIC, test_data)
    cat("Test run successful\n", file=stdout())
}

MAGISTA_run_generic = function(modeltype,infile,outfile,model){
    modelname=model
    cwd=getwd()
    model=paste(MAGISTA_dir,"data/",modeltype,".",modelname,".RData",sep="")
    if (!file.exists(model)){
        cat(paste("[WARN] Model '",model,"' not found\n",sep=""),file=stderr())
        cat("       A new model will be generated from the input data.\n",file=stderr())
        MAGISTA_prepare_generic(modeltype,infile, modelname)
    }
    load(file=model)
    inputdata = read.csv(infile,header=T)
    if (modeltype=="MAGISTA"){
        outputdata = predict_MAGISTA(MAGISTA, inputdata)
    }
    if (modeltype=="MAGISTIC"){
        predict_MAGISTA(MAGISTIC, inputdata)
    }
    idname=NULL
    for (potential_idname in c("Bin.Id","Bin.id","Id","ID","id","Bin","bin","name","Name","X")) {
        if ( potential_idname %in% colnames(inputdata)) {
            idname = potential_idname
        }
    }
    if (!is.null(idname)){
        tmp = data.frame(Bin.id=inputdata[,idname])
        colnames(tmp)[1] = idname
        outputdata = cbind(tmp,outputdata)
    } else {
        outputdata = cbind(data.frame(Bin.id=paste("Bin_",1:nrow(outputdata),sep="")),outputdata)
    }
    setwd(cwd)
    write.csv(outputdata,outfile,row.names = F,quote=F)
}

MAGISTA_run = function(infile,outfile=paste(gsub(".csv","",gsub("_input","",infile)),"_output.csv"),model="default"){
    cat("Run parameters\n", file=stderr())
    cat(" Model type: MAGISTA\n", file=stderr())
    cat(paste(" Input file name:",infile,"\n"), file=stderr())
    cat(paste(" Output file name:",outfile,"\n"), file=stderr())
    cat(paste(" Model name:",model,"\n"), file=stderr())
    MAGISTA_run_generic("MAGISTA",infile,outfile,model)
}
MAGISTIC_run = function(infile,outfile=paste(gsub(".csv","",gsub("_input","",infile)),"_output.csv"),model="default"){
    cat("Run parameters\n", file=stderr())
    cat(" Model type: MAGISTIC\n", file=stderr())
    cat(paste(" Input file name:",infile,"\n"), file=stderr())
    cat(paste(" Output file name:",outfile,"\n"), file=stderr())
    cat(paste(" Model name:",model,"\n"), file=stderr())
    MAGISTA_run_generic("MAGISTIC",infile,outfile,model)
}

args = commandArgs(trailingOnly=TRUE)

proper_input=T
cat("This is MAGISTA v0.2\n", file=stderr())
cat("Supplied arguments:\n  ", file=stderr())
cat(paste(args,collapse=","), file=stderr())
cat("\n", file=stderr())
if( length(args)==1 && args[1]=="FIRSTRUN") {
    MAGISTA_testrun()
} else if(length(args)>1) {
    if(!grep(pattern=".csv",x=args[2],fixed=T)){
        cat("[Error] Input (second argument) should be a .csv file\n", file=stdout())
        proper_input=F
    }
    if(args[1]!="MAGISTA" && args[1]!="MAGISTIC"){
        cat("[Error] Approach (first argument) must be 'MAGISTA' or 'MAGISTIC'\n", file=stdout())
        proper_input=F
    }
    if(proper_input){
        if(length(args)==2){
            if(args[1]=="MAGISTA"){
                MAGISTA_run(args[2])
            } else if (args[1]=="MAGISTIC") {
                MAGISTIC_run(args[2])
            }
        }
        if(length(args)==3){
            if(args[1]=="MAGISTA"){
                MAGISTA_run(args[2],args[3])
            } else if (args[1]=="MAGISTIC") {
                MAGISTIC_run(args[2],arg[3])
            }
        }
        if(length(args)==4){
            if(args[1]=="MAGISTA"){
                MAGISTA_run(args[2],args[3],args[4])
            } else if (args[1]=="MAGISTIC") {
                MAGISTIC_run(args[2],arg[3],args[4])
            }
        }
    }
} else {
    proper_input=F
    if(!interactive()){
        cat(paste("Only",length(args),"arguments supplied\n"), file=stdout())
        cat("Usage: Rscript MAGISTA.R <FIRSTRUN|MAGISTA|MAGISTIC> [<input file>] [<output file>] [<model name>]\n", file=stdout())
        cat("       Please note that <input file> and <output file> must be .csv files, and <model name> should not contain extensions.\n", file=stdout())
        cat("       If <model name> does not match any existing model names, an attempt will be made to create a new model in\n       ", file=stdout())
        cat(MAGISTA_dir,file=stdout())
        cat("\n       and based on the input data - in this case, <output file> will contain the prediction for the input data\n", file=stdout())
        quit(status=1)
    }
}


