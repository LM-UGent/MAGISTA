library("tidyverse")
library("reshape2")
library("ggplot2")
library("extrafont")
library("randomForest")
library("randomForestExplainer")
library("caret")
if(!exists("fonts_imported")){
    font_import() #run this only once
    windowsFonts() #run this only once
    fonts_imported=T
}

order.DataFrame = function(df,by=colnames(df)[1]){
    return(df[order(df[,by]),])
}

weighted_confusion_matrix= function(predicted, observed, case_weights=NULL){
    if(is.null(case_weights)){
        case_weights=1
    }
    upred = unique(predicted)
    uobs = unique(observed)
    result = matrix(nrow=length(upred), ncol=length(uobs), dimnames = list(predicted=upred,observed=uobs))

    for(up in upred){
        for(uo in uobs){
            result[up,uo] = sum( (predicted==up)*(observed==uo)*case_weights )
        }
    }
    return(result)
}
weighted_accuracy= function(predicted, observed, case_weights=NULL){
    if(is.null(case_weights)){
        case_weights=1
    }
    return(sum((predicted==observed)*case_weights)*1.0/(sum(case_weights)*1.0))
}

weighted_mean = function(x, weights = 1, na.rm=T){
    if(length(weights)==1){
        weights = numeric(length(x))*0+weights
    }
    if(na.rm){
        weights = weights[!is.na(x)]
        x = x[!is.na(x)]
    }
    totweight = sum(weights)
    return(sum(x*weights)/totweight)
}

weighted_var = function(x, weights = 1, na.rm=T) {
    if(length(weights)==1){
        weights = numeric(length(x))*0+weights
    }
    if(na.rm){
        weights = weights[!is.na(x)]
        x = x[!is.na(x)]
    }
    totweight = sum(weights)
    xmean = sum(x*weights)/totweight
    cs_x = (x - xmean)
    return(sum((cs_x)^2*weights)/totweight)
}

weighted_cov = function(x, y ,weights = 1, na.rm=T){
    if(length(weights)==1){
        weights = numeric(length(x))*0+weights
    }
    if(na.rm){
        weights = weights[!is.na(x)]
        x = x[!is.na(x)]
    }
    totweight = sum(weights)
    xmean = sum(x*weights)/totweight
    ymean = sum(y*weights)/totweight
    totweight = sum(weights)
    cx = (x - xmean)
    cy = (y - ymean)
    vx = sum((cx)^2*weights)
    vy = sum((cy)^2*weights)
    return( sum(cx*cy*weights)/totweight )
}


weighted_cor = function(x, y ,weights = 1, na.rm=T){
    return( weighted_cov(x,y,weights,na.rm)/(sqrt(weighted_var(x,weights,na.rm))*sqrt(weighted_var(y,weights,na.rm))) )
}
table_rsquared = function(predicted, actual, calcmode="lm", weights = 1){
    #observed values are the actual values of the bin
    #expected values are those predicted by the model
    cc = complete.cases(data.frame(predicted,actual))
    x = predicted[cc]
    y = actual[cc]
    if(length(weights)==1 && weights == "uniform"){
        #sort the observations
        new_order = order(y)
        y = y[new_order]
        x = x[new_order]
        # calculate weights in such a way that sum(weights) == length(data)
        # and weight(x_n) is proportional to the distance to its neighbors

        # calculate distances - weights will be proportial to them
        dy = y[2:length(y)] - y[1:(length(y)-1)]
        weights = numeric(length(y))
        weights[1:(length(y)-1)] =  weights[1:(length(y)-1)] + dy
        weights[2:length(y)] =  weights[2:length(y)-1] + dy
        weights[2:(length(y)-1)] = weights[2:(length(y)-1)]/2

        totweight = sum(weights)
        factor = (length(y)*1.0)/totweight
        weights = weights*factor
    }
    if(calcmode=="lm"){
        return(weighted_cor(predicted,actual,weights)^2)
    }
    if(calcmode=="model"){
        g = mean(y)
        delta = y-x
        R2 = 1 - sum((delta^2)*weights)/sum( (y-g)^2 * weights, na.rm = T)
        return(R2)
    }
    return(NULL)
}
table_rmse = function(expected, observed, weights = 1){
    if(length(weights)==1){
        weights = numeric(length(expected))*0+weights
    }
    cc = complete.cases(data.frame(expected,observed))
    weights = weights[cc]
    if(weights == "uniform"){
        #sort the observations
        new_order = order(y)
        y = y[new_order]
        x = x[new_order]
        # calculate weights in such a way that sum(weights) == length(data)
        # and weight(x_n) is proportional to the distance to its neighbors

        # calculate distances - weights will be proportial to them
        dy = y[2:length(y)] - y[1:(length(y)-1)]
        weights = numeric(length(y))
        weights[1:(length(y)-1)] =  weights[1:(length(y)-1)] + dy
        weights[2:length(y)] =  weights[2:length(y)-1] + dy
        weights[2:(length(y)-1)] = weights[2:(length(y)-1)]/2

        totweight = sum(weights)
        factor = (length(y)*1.0)/totweight
        weights = weights*factor
    }
    delta = (expected[cc]-observed[cc])^2
    return(sqrt(sum(delta*weights)/sum(weights)))
}
table_avgchisq = function(expected, observed, weights = 1){
    if(class(weights) != "character") {
        if(length(weights)==1){
            weights = numeric(length(expected))*0+weights
        }
    }
    cc = complete.cases(data.frame(expected,observed))
    weights = weights[cc]
    if(class(weights) != "character") {
        if(weights == "uniform"){
            #sort the observations
            new_order = order(y)
            y = y[new_order]
            x = x[new_order]
            # calculate weights in such a way that sum(weights) == length(data)
            # and weight(x_n) is proportional to the distance to its neighbors

            # calculate distances - weights will be proportial to them
            dy = y[2:length(y)] - y[1:(length(y)-1)]
            weights = numeric(length(y))
            weights[1:(length(y)-1)] =  weights[1:(length(y)-1)] + dy
            weights[2:length(y)] =  weights[2:length(y)-1] + dy
            weights[2:(length(y)-1)] = weights[2:(length(y)-1)]/2

            totweight = sum(weights)
            factor = (length(y)*1.0)/totweight
            weights = weights*factor
        }
    }
    delta = (expected[cc]-observed[cc])^2
    return(sum( (delta/expected[cc])*weights )/sum(weights))
}
table_accuracy = function(expected, observed){
    return(sum(expected==observed)/length(expected))
}
table_avgerr = function(expected, observed){
    return(sum(expected-observed)*1.0/length(expected))
}

maxidx <- function(arr) {
    return(which(arr == max(arr)))
}

shape = function(x){
    if(is.null(nrow(x)) || is.null(ncol(x))){
        return(length(x))
    }
    return(c(nrow(x),ncol(x)))
}

factor_mode = function(x, na.rm=F) {
    if(na.rm){
        x = x[!is.na(x)]
    }
    uniqx = unique(x)
    uniqx[which.max(tabulate(match(x, uniqx)))]
}

uniform_sample_indices = function(data_vals, alongcol, nbins){
    i=0
    rmcol=F
    if(length(alongcol)==nrow(data_vals)){
        data_vals = cbind(target=alongcol, data_vals)
        alongcol = "target"
        rmcol=T
    }
    if(!(alongcol %in% colnames(data_vals))){
        alongcol=1
    }
    data_vals = data_vals[!is.na(data_vals[,alongcol]),]
    minvalue = min(data_vals[,alongcol])
    maxvalue = max(data_vals[,alongcol])
    data_vals = data_vals[order(data_vals[,alongcol]),]
    mincount = nrow(data_vals)

    for(i in 1:nbins){
        curmin = (i-1)*(maxvalue-minvalue)/(1.0*nbins)
        curmax = (i)*(maxvalue-minvalue)/(1.0*nbins)
        curminrow = sum(data_vals[,alongcol]<=curmin)
        curmaxrow = sum(data_vals[,alongcol]<=curmax)
        if(curminrow<1){
            curminrow=1
        }
        if(curmaxrow>nrow(data_vals)){
            curmaxrow=nrow(data_vals)
        }
        if(mincount > curmaxrow-curminrow + 1 && curmaxrow-curminrow + 1 > 0){
            mincount = curmaxrow-curminrow + 1
        }
    }
    result = logical(nrow(data_vals))
    for(i in 1:nbins){
        curmin = (i-1)*(maxvalue-minvalue)/(1.0*nbins)
        curmax = (i)*(maxvalue-minvalue)/(1.0*nbins)
        curminrow = sum(data_vals[,alongcol]<=curmin)
        curmaxrow = sum(data_vals[,alongcol]<=curmax)
        if(curminrow<1){
            curminrow=1
        }
        if(curmaxrow>nrow(data_vals)){
            curmaxrow=nrow(data_vals)
        }
        result[sample(curminrow:curmaxrow,mincount)] = TRUE
    }

    return(result)
}

uniform_sampling = function(data_vals, alongcol, nbins){
    i=0
    rmcol=F
    if(length(alongcol)==nrow(data_vals)){
        data_vals = cbind(target=alongcol, data_vals)
        alongcol = "target"
        rmcol=T
    }
    if(!(alongcol %in% colnames(data_vals))){
        alongcol=1
    }
    data_vals = data_vals[!is.na(data_vals[,alongcol]),]
    minvalue = min(data_vals[,alongcol])
    maxvalue = max(data_vals[,alongcol])
    data_vals = data_vals[order(data_vals[,alongcol]),]
    mincount = nrow(data_vals)

    for(i in 1:nbins){
        curmin = (i-1)*(maxvalue-minvalue)/(1.0*nbins)
        curmax = (i)*(maxvalue-minvalue)/(1.0*nbins)
        curminrow = sum(data_vals[,alongcol]<=curmin)
        curmaxrow = sum(data_vals[,alongcol]<=curmax)
        if(curminrow<1){
            curminrow=1
        }
        if(curmaxrow>nrow(data_vals)){
            curmaxrow=nrow(data_vals)
        }
        if(mincount > curmaxrow-curminrow + 1 && curmaxrow-curminrow + 1 > 0){
            mincount = curmaxrow-curminrow + 1
        }
    }
    new_data = data.frame()
    for(i in 1:nbins){
        curmin = (i-1)*(maxvalue-minvalue)/(1.0*nbins)
        curmax = (i)*(maxvalue-minvalue)/(1.0*nbins)
        curminrow = sum(data_vals[,alongcol]<=curmin)
        curmaxrow = sum(data_vals[,alongcol]<=curmax)
        if(curminrow<1){
            curminrow=1
        }
        if(curmaxrow>nrow(data_vals)){
            curmaxrow=nrow(data_vals)
        }
        tmp_data = data_vals[sample(curminrow:curmaxrow,mincount),]
        new_data = bind_rows(new_data,tmp_data)
    }

    if(rmcol){
        new_data = new_data[,2:ncol(new_data)]
    }
    return(new_data)
}
uniform_category_sampling = function(data_vals, alongcol){
    i=0
    mincount = nrow(data_vals)
    uvals = unique(data_vals[,alongcol])
    for(i in uvals){
        count = sum(data_vals[,alongcol] == i)
        if(mincount > count){
            mincount = count
        }
    }
    new_data = data.frame()
    for(i in uvals){
        tmp_data = data_vals[sample(which(data_vals[,alongcol] == i),mincount),]
        new_data = bind_rows(new_data,tmp_data)
    }
    return(new_data)
}
pur2cat = function(x){
    if(x>95){
        return(1)
    }
    if(x>90){
        return(2)
    }
    if(x>85){
        return(3)
    }
    return(4)
}
comp2cat = function(x){
    if(x>90){
        return(1)
    }
    if(x>70){
        return(2)
    }
    if(x>50){
        return(3)
    }
    return(4)
}

numericcol = function(data_src){
    colfilter = sapply(data_src[1,], is.numeric)
    return(colfilter)
}
namepatterncol = function(original, pattern, inverse=F, use.regex=F){
    cols = logical(ncol(original))
    cols[grep(pattern=pattern,x=colnames(original),fixed=!use.regex)] = T
    if(inverse){
        cols = !cols
    }
    return(cols)
}
idofcol = function(data,colname){
    cnames = colnames(data)
    return(which(colname == cnames))
}

test_regression = function(original, predicted, color){
    datatoplot = data.frame(o=original,p=predicted,c=color)

}

regression_plot = function(data,x,y, textfamily="Arial", textsize=20, pointsize=1, transparency=NULL, labelx=NULL,labely=NULL, minv = NULL, maxv = NULL, showlm=FALSE, only.values=F, equalxy=T, weights=NULL, verbose = F, label_extreme = F, label_col = NULL, ...){
    if(label_extreme == T){
        if(is.null(label_col)){
            label_col = colnames(data)[1]
        }
        label_extreme = 2
    }
    if(verbose){
        message(paste("input data has",nrow(data),"rows and",ncol(data),"cols"))
        message("Colnames ares:")
        message(toString(colnames(data)))
    }
    data = data[complete.cases(data[,c(x,y)]),]
    if(verbose){
        message(paste("input data has",nrow(data),"remaining after removing incomplete cases"))
    }
    xydata=data.frame(actual=data[,x], predicted=data[,y])
    lmcolor="blue"
    if(!is.null(showlm) && !is.logical(showlm)){
        lmcolor=as.character(showlm)
        showlm=T
    }
    g = mean(data[complete.cases(data),y])
    delta = data[,x]-data[,y]
    if(is.null(weights)){
        SSE = sum(delta*delta, na.rm = T)
        AE = sqrt(SSE/nrow(data))
        R2 = table_rsquared(data[,y],data[,x],"model")
    } else {
        if(is.character(weights)){
            weights = data[,weights]
        }
        weights[weights<0]=0
        SSE = sum(delta*delta * weights, na.rm = T)
        AE = sqrt(SSE/sum(weights))
        R2 = table_rsquared(data[,y],data[,x],"model", weights = weights)
    }
    if(nrow(data)>=2){
        linmodel=lm(formula=as.formula("predicted~actual"),data=xydata,weights = weights)
        y_lm = predict(linmodel,xydata)
    } else {
        if(nrow(data)==0){
            linmodel=lm(formula = as.formula("predicted~actual"),data=data.frame(predicted=c(0,1),actual=c(0,1)))
        } else {
            linmodel=lm(formula = as.formula("predicted~actual"),data=data.frame(predicted=c(xydata$predicted,xydata$predicted),actual=c(xydata$actual,xydata$actual)))
        }
        y_lm = xydata$actual
    }
    delta_lm = xydata$actual - y_lm
    if(is.null(weights)){
        SSE_lm = sum(delta_lm*delta_lm, na.rm = T)
        AE_lm = sqrt(SSE/nrow(data))
        R2_lm = table_rsquared(data[,y],data[,x],"lm")
    } else {
        if(is.character(weights)){
            weights = data[,weights]
        }
        SSE_lm = sum(delta_lm*delta_lm * weights, na.rm = T)
        AE_lm = sqrt(SSE_lm/sum(weights))
        R2_lm = table_rsquared(data[,y],data[,x],"lm", weights = weights)
    }
    cor_lm = sqrt(R2_lm)

    #ex = paste("{R^2}[y%~%x]==",floor(R2*1000)/1000,sep="")

    if(only.values){
        return(data.frame(R2_xy=c(R2),R2_lm=c(R2_lm),cor_lm=c(cor_lm),RMSE_xy=c(AE),RMSE_lm=c(AE_lm)))
    }

    ex1 = paste("atop({{R^2}[y%~%x]==",as.character(floor(R2*1000)/1000),"},",
                "{{RMSE}[y%~%x]==",as.character(floor(AE*100)/100),"})",
                sep="")
    ex2 = paste("atop({{cor}[lm]==",as.character(floor(cor_lm*1000)/1000),"},",
                "{{RMSE}[lm]==",as.character(floor(AE_lm*100)/100),"})",
                sep="")

    minx = min(data[,x],na.rm=T)
    maxx = max(data[,x],na.rm=T)
    miny = min(data[,y],na.rm=T)
    maxy = max(data[,y],na.rm=T)
    if(!is.null(minv)){
        minx=minv
        miny=minv
    }
    if(!is.null(maxv)){
        maxx=maxv
        maxy=maxv
    }

    data = data[data[,x]>=minx,]
    data = data[data[,x]<=maxx,]
    data = data[data[,y]>=miny,]
    data = data[data[,y]<=maxy,]

    valuerange=c(min(minx,miny),max(maxx,maxy))
    if(equalxy){
        minx= valuerange[1]
        maxx= valuerange[2]
        miny= valuerange[1]
        maxy= valuerange[2]
    }
    ex = NULL
    if(is.null(labelx)){
        labelpos1_x = minx+(maxx-minx)*0.7
        labelpos2_x = minx+(maxx-minx)*0.3
    } else {
        labelpos1_x= labelx
        labelpos2_x= labelx
    }
    if(is.null(labely)){
        labelpos1_y = miny+(maxy-miny)*0.15
        labelpos2_y = miny+(maxy-miny)*0.85
    } else {
        labelpos1_y= labely
        labelpos2_y= labely
    }
    if(labelpos1_y==labelpos2_y && labelpos1_x==labelpos2_x){
        ex = paste("atop(",ex1,",",ex2,")",sep="")
        textdf = data.frame(x=c(labelpos_x),y=c(labelpos_y),t=c(ex1))
    }

    textdf_1 = data.frame(x=c(labelpos1_x),y=c(labelpos1_y),t=c(ex1))
    textdf_2 = data.frame(x=c(labelpos2_x),y=c(labelpos2_y),t=c(ex2))


    resplot = ggplot(data)
    if(is.null(transparency)){
        if(is.numeric(pointsize)){
            resplot = resplot + geom_point(aes_string(x=x,y=y, ...),size=pointsize)
        } else {
            resplot = resplot + geom_point(aes_string(x=x,y=y, ...))
        }
    } else {
        if(is.numeric(pointsize)){
            resplot = resplot + geom_point(aes_string(x=x,y=y, ...),size=pointsize,alpha=1-transparency)
        } else {
            resplot = resplot + geom_point(aes_string(x=x,y=y, ...),alpha=1-transparency)
        }
    }
    resplot = resplot + geom_abline(intercept=0,slope=1, color="black", size=1.2)
    resplot = resplot + geom_abline(intercept=linmodel$coefficients[1],slope=linmodel$coefficients[2], color=lmcolor, size=1.2)
    if(is.null(ex)){
        if(showlm){
            resplot = resplot + geom_text(data=textdf_1,aes(x=x,y=y,label=t),parse = T,size = textsize*0.35, family=textfamily, color="black" )
            resplot = resplot + geom_text(data=textdf_2,aes(x=x,y=y,label=t),parse = T,size = textsize*0.35, family=textfamily, color=lmcolor)
        } else {
            resplot = resplot + geom_text(data=textdf_1,aes(x=x,y=y,label=t),parse = T,size = textsize*0.35, family=textfamily )
        }
    } else {
        resplot = resplot + geom_text(data=textdf,aes(x=x,y=y,label=t),parse = T,size = textsize*0.35, family=textfamily )
    }
    if(is.numeric(label_extreme)){
        if(verbose){
            message(paste("adding",sum((abs(data[,x]-data[,y])>(label_extreme*AE))),"labels for col:",label_col))
        }
        resplot = resplot + geom_text(aes_string(x=x,y=y,label=label_col),data=data[abs(data[,x]-data[,y])>(label_extreme*AE),])
    }
    if(equalxy){
        resplot = resplot + xlim(valuerange) + ylim(valuerange)
    }
    return(resplot)
}


test_prediction_regress = function(prediction,relevant_data,class_table,classtype,cclass,test_set,graphpath,graph_prefix,test=T){
    tmp_cor = cor(x=relevant_data[test_set,1],y=prediction,use = "pairwise.complete.obs")
    delta = (relevant_data[test_set,1]-prediction)
    SSE = sum(delta*delta, na.rm = T)
    AE = sqrt(SSE/length(prediction))
    R2 = 1 - var(delta)/var(prediction)
    print(c(cclass,tmp_cor,AE))
    result = data.frame(desc = c(paste(as.character(classtype),as.character(cclass))), val = c(tmp_cor))

    labelpos_x = mean(relevant_data[test_set,1],na.rm=T)
    labelpos_y = mean(prediction[test_set],na.rm=T)
    if(is.null(class_table)){
        data_to_plot = data.frame(prediction=prediction,original=relevant_data[test_set,1],classval=as.character(classtype))
    } else {
        data_to_plot = data.frame(prediction=prediction,original=relevant_data[test_set,1],classval=class_table[test_set,classtype])
    }
    ex = paste("atop({R[x==y]}^2==",R2,",",
               "{R[lm]}^2==",as.character(tmp_cor*tmp_cor,use="pairwise.complete.obs"),")",
               sep="")
    textdf = data.frame(x=c(labelpos_x),y=c(labelpos_y),t=c(ex))

    axlims = seq(0,1,0.2)
    axmax = 1
    if(max(relevant_data[,1])>80){
        axlims = seq(0,100,20)
        axmax = 100
    } else {
        if(max(relevant_data[,1])>4){
            axlims = seq(0,5,1)
            axmax = 5
        }
    }

    p1 = ggplot(data_to_plot)+
        geom_point(aes(x=original,y=prediction,color=classval))+
        scale_x_continuous(breaks=axlims, limits = c(0,axmax)) +
        scale_y_continuous(breaks=axlims, limits = c(0,axmax)) +
        geom_smooth(aes(x=original,y=prediction),method="lm")+
        geom_abline(intercept=0,slope=1) +
        labs(caption=paste("Average error:",as.character(AE), "(based on sum of square errors)")) +
        geom_text(data=textdf,aes(x=x,y=y,label=t),parse = T) +
        ggtitle(paste(as.character(classtype),as.character(cclass))) +
        theme_minimal()

    if(test==T){
        plotname = paste(graphpath, graph_prefix,"_TEST_",colnames(relevant_data)[1],"_",as.character(classtype),"_",as.character(cclass),".png",sep="")
    } else {
        plotname = paste(graphpath, graph_prefix,"_",colnames(relevant_data)[1],"_",as.character(classtype),"_",as.character(cclass),".png",sep="")
    }
    ggsave(plot=p1,filename=plotname,dpi=300,width = 8, height = 5, units = "cm", scale=2)
    return(result)
}
test_prediction_classif =function(prediction,relevant_data,class_table,classtype,cclass,test_set,graphpath,graph_prefix){
    misclassifications=sum(relevant_data[test_set,1] != prediction, na.rm = T)
    print(c(cclass,misclassifications/length(prediction)))

    classes = unique(c(as.character(prediction),as.character(relevant_data[,1])))
    confusion_matrix = matrix(nrow = length(classes), ncol=length(classes))
    colnames(confusion_matrix) = classes
    rownames(confusion_matrix) = classes

    for(r in classes){
        for(c in classes){
            confusion_matrix[r,c] = sum((prediction==r)*(relevant_data[test_set,1]==c))
        }
    }

    data_to_plot = melt(confusion_matrix)
    colnames(data_to_plot)=c("prediction","original","count")
    p1 = ggplot(data_to_plot)+
        geom_raster(aes(x=original,y=prediction,fill=count))+
        geom_text(aes(x=original,y=prediction,label=count))+
        scale_fill_gradient(low="grey90", high="red") +
        labs(caption=paste("Misclassification rate:",as.character(misclassifications/length(prediction)))) +
        theme_minimal()
    plotname = paste(graphpath, graph_prefix,"_TEST_",colnames(relevant_data)[1],"_",as.character(classtype),"_",as.character(cclass),".png",sep="")
    ggsave(plot=p1,filename=plotname,dpi=300,width = 8, height = 5, units = "cm", scale=2)

    return(data.frame(desc = c(paste(as.character(classtype),as.character(cclass))),val=c(misclassifications/nrow(relevant_data))))
}


discretize = function(x,binw=NA,nbins=NA, minx=NA, maxx=NA){
    if(is.na(binw) && is.na(nbins)){
        error("please specify a number of bins or bin width, but not both")
    }
    if(!is.na(binw) && !is.na(nbins)){
        error("please specify a number of bins or bin width, but not both")
    }
    if(is.na(minx)){
        minx=min(x)
    }
    if(is.na(maxx)){
        maxx = max(x)
    }
    if(minx>=maxx){
        error("minimal values hsould be strictly below maximal value")
    }

    if(is.na(binw)){
        binw = (maxx-minx)/(nbins*1.0)
    } else{
        nbins = as.integer(ceiling((maxx-minx)/binw))
    }

    result = data.frame(bin=numeric(nbins),count=integer(nbins))
    for(i in 1:nbins){
        threshold = minx+(1.0*i)*binw
        result$bin[i] = threshold
        result$count[i] = sum(x<threshold)
        x = x[x>=threshold]
    }
    return(result)
}

get_FP_TP = function(df,threshold, col1=2, col2=3){
    part1 = df[df[,1]<=threshold,]
    part2 = df[df[,1]>threshold,]

    TP = sum(part1[,col1])
    FP = sum(part1[,col2])
    TN = sum(part2[,col2])
    FN = sum(part2[,col1])

    TPR = TP/(TP+FN)
    FPR = FP/(FP+TN)

    return(c(FPR,TPR))
}
get_MCC = function(df,threshold, col1=2, col2=3){
    part1 = df[df[,"bin"]<=threshold,]
    part2 = df[df[,"bin"]>threshold,]

    TP = sum(part1[,col1],na.rm = T)*1.0
    FP = sum(part1[,col2],na.rm = T)*1.0
    TN = sum(part2[,col2],na.rm = T)*1.0
    FN = sum(part2[,col1],na.rm = T)*1.0

    TPR = TP/(TP+FN)
    FPR = FP/(FP+TN)

    MCC = TP*TN - FP*FN
    MCC = MCC/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    return(MCC)
}
create_ROC = function(df,firstcol=2,secondcol=3){
    # df should consist of three columns: bin, good, bad (names do not matter, order does)
    # this function will iteratively go through all bin values and check wether using them a thresholds produces good results
    # note: it is assumed that "good" values are on the left and "bad" values are on the right
    df = as.data.frame(df)
    xy = data.frame(matrix(nrow=length(df[,1]), ncol = 3))
    colnames(xy) = c("bin","FPR","TPR")
    j = 0
    for(i in df[,1]){
        j=j+1
        xy[j,1] = i
        xy[j,2:3] = get_FP_TP(df,i, firstcol, secondcol)
    }
    return(xy)
}
create_MCC = function(df,firstcol=2,secondcol=3){
    # df should consist of three columns: bin, good, bad (names do not matter, order does)
    # this function will iteratively go through all bin values and check wether using them a thresholds produces good results
    # note: it is assumed that "good" values are on the left and "bad" values are on the right
    xy = data.frame(matrix(nrow=length(df[,1]), ncol = 3))
    colnames(xy) = c("bin","Cutoff","MCC")
    j = 0
    for(i in df[,1]){
        j=j+1
        xy[j,1] = i
        xy[j,2] = i
        xy[j,3] = get_MCC(df,i, firstcol, secondcol)
    }
    return(xy)
}

Area_Under_Curve = function(x,y){
    if(length(x)!=length(y)){
        error("Inputs do not have the same dimensions! (Area_Under_Curve)")
    }
    area=0
    for(i in 2:length(x)){
        area = area + (x[i]-x[i-1])*(y[i]+y[i-1])/2.0
    }
    return(area)
}

ROC_analysis = function(df,predictor1,predictor2,actual_class,
                        discretization1=50, discretization2=50,
                        bad_level = levels(df[,actual_class])[1],
                        good_level = levels(df[,actual_class])[2],
                        minp1=min(df[,predictor1]), minp2=min(df[,predictor2]),
                        maxp1=max(df[,predictor1]), maxp2=max(df[,predictor2])){

    step1 = (maxp1-minp1)/discretization1
    step2 = (maxp2-minp2)/discretization2

    discretization1 = discretization1-1
    discretization2 = discretization2-1

    bad_below1_below2 = matrix(nrow=discretization1,ncol=discretization2)
    bad_below1_above2= matrix(nrow=discretization1,ncol=discretization2)
    bad_above1_below2= matrix(nrow=discretization1,ncol=discretization2)
    bad_above1_above2= matrix(nrow=discretization1,ncol=discretization2)

    good_below1_below2 = matrix(nrow=discretization1,ncol=discretization2)
    good_below1_above2= matrix(nrow=discretization1,ncol=discretization2)
    good_above1_below2= matrix(nrow=discretization1,ncol=discretization2)
    good_above1_above2= matrix(nrow=discretization1,ncol=discretization2)


    for(i in 1:discretization1){
        threshold1 = minp1 + step1*i
        for(j in 1:discretization2){
            threshold2 = minp2 + step2*j

            bad_below1_below2[i,j] = sum(as.logical((df[,predictor1]<=threshold1)*(df[,predictor2]<=threshold2)*(df[,actual_class]==bad_level)))
            bad_below1_above2[i,j] = sum(as.logical((df[,predictor1]<=threshold1)*(df[,predictor2]>threshold2)*(df[,actual_class]==bad_level)))
            bad_above1_below2[i,j] = sum(as.logical((df[,predictor1]>threshold1)*(df[,predictor2]<=threshold2)*(df[,actual_class]==bad_level)))
            bad_above1_above2[i,j] = sum(as.logical((df[,predictor1]>threshold1)*(df[,predictor2]>threshold2)*(df[,actual_class]==bad_level)))

            good_below1_below2[i,j] = sum(as.logical((df[,predictor1]<=threshold1)*(df[,predictor2]<=threshold2)*(df[,actual_class]==good_level)))
            good_below1_above2[i,j] = sum(as.logical((df[,predictor1]<=threshold1)*(df[,predictor2]>threshold2)*(df[,actual_class]==good_level)))
            good_above1_below2[i,j] = sum(as.logical((df[,predictor1]>threshold1)*(df[,predictor2]<=threshold2)*(df[,actual_class]==good_level)))
            good_above1_above2[i,j] = sum(as.logical((df[,predictor1]>threshold1)*(df[,predictor2]>threshold2)*(df[,actual_class]==good_level)))
        }
    }


    # determine whether or not good classes correspond to low or high values of the variables
    high_is_good_1 = T
    high_is_good_2 = T
    if(sum((good_below1_above2+good_below1_below2)>(good_above1_above2+good_above1_below2))>(discretization1*discretization2/2)){
        high_is_good_1 = F
    }
    if(sum((good_below1_below2+good_above1_below2)>(good_above1_above2+good_below1_above2))>(discretization1*discretization2/2)){
        high_is_good_2 = F
    }
    if (high_is_good_1==T && high_is_good_2==T){
        TPs = good_above1_above2
        FPs = bad_above1_above2
        TNs = bad_below1_above2 + bad_above1_below2 + bad_below1_below2
        FNs = good_below1_above2 + good_above1_below2 + good_below1_below2
    } else if (high_is_good_1==T && high_is_good_2==F){
        TPs = good_above1_below2
        FPs = bad_above1_below2
        TNs = bad_below1_above2 + bad_above1_above2 + bad_below1_below2
        FNs = good_below1_above2 + good_above1_above2 + good_below1_below2
    } else if (high_is_good_1==F && high_is_good_2==T){
        TPs = good_below1_above2
        FPs = bad_below1_above2
        TNs = bad_above1_above2 + bad_above1_below2 + bad_below1_below2
        FNs = good_above1_above2 + good_above1_below2 + good_below1_below2
    } else if (high_is_good_1==F && high_is_good_2==F){
        TPs = good_below1_below2
        FPs = bad_below1_below2
        TNs = bad_below1_above2 + bad_above1_above2 + bad_above1_below2
        FNs = good_below1_above2 + good_above1_above2 + good_above1_below2
    }

    thresholds1 = (1:discretization1)*step1+minp1
    thresholds2 = (1:discretization2)*step2+minp2

    TPR = TPs/(TPs + FNs)
    FPR = FPs/(FPs + TNs)

    colnames(TPR) = as.character(thresholds2)
    rownames(TPR) = as.character(thresholds1)
    colnames(FPR) = as.character(thresholds2)
    rownames(FPR) = as.character(thresholds1)

    AUROC1 = numeric(discretization2)
    AUROC2 = numeric(discretization1)

    for(i in 1:discretization2){
        AUROC1[i] = Area_Under_Curve(x=FPR[,i],y=TPR[,i])
    }
    for(i in 1:discretization1){
        AUROC2[i] = Area_Under_Curve(x=FPR[i,],y=TPR[i,])
    }

    return(list(TPR=TPR, FPR=FPR, AUROC1=AUROC1, AUROC2=AUROC2, thresholds1=thresholds1, thresholds2=thresholds2))
}

heatmap_plot = function(data, rescale.rows=F, rescale.cols=F, ...){
    data = as.data.frame(data)
    for(column in colnames(data)){
        if(!is.numeric(data[,column])){
            data[,column] = as.integer(data[,column])
        }
    }
    if(rescale.rows && !rescale.cols){
        for(i in 1:nrow(data)){
            data[i,] = data[i,]/max(data[i,])
        }
    }else if(rescale.cols && !rescale.rows){
        for(i in 1:nrow(data)){
            data[i,] = data[i,]/max(data[i,])
        }
    } else {
        data = data/max(data)
    }
    if(!is.null(row.names(data))){
        data = cbind(row.names(data),data)
    }
    newdata = melt(data=data,id.vars=c(1))

    return(geom_tile(fztz=newdata,aes(x)))
}


all_equal = function(A,B){
    for(i in ncol(A)){
        for(j in nrow(A)){
            if(A[j,i] != B[j,i]){
                return(FALSE)
            }
        }
    }
    return(TRUE)
}

matrix_inverse = function(x){
    eye = diag(x=1,nrow=nrow(x),ncol=ncol(x))
    if(all_equal(t(x) %*% x, eye) == TRUE){
        return(t(x))
    }
    else{
        return(solve(x))
    }
}

mpm.logtransf = function(data, logrepl = 1e-09){
    logrepl <- max(1e-09, logrepl)
    if (any(data <= 0, na.rm = TRUE)) {
        warning(paste("Non-positive data replaced by", logrepl,
                      "in computing mpm.logtransf"),
                call. = FALSE)
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
    training_data = training_data[,to_keep]


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

mpm.sma_axes = function(training_data){
    ncol1 = ncol(training_data)
    complete_cases = complete.cases(training_data)
    training_data = training_data[complete_cases,]

    #remove variables which do not vary
    tmpData <- as.matrix(training_data)
    datavar = diag(var(tmpData))
    to_keep = (abs(datavar)>1e-09)
    training_data = training_data[,to_keep]


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

    mpm_result = mpm.2(training_data, row.weight = "mean", col.weight = "mean")

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

mpm.pca_fit = function(training_data, target_data){
    pca_axes = mpm.pca_axes(training_data)
    return(mpm.pca_predict(pca_axes,target_data))
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

mpm.2.logpca = function(data){
    return(mpm.2(data, center = "column", normal = "column", samples="rows"))
}
mpm.2.sma = function(data){
    return(mpm.2(data, row.weight = "mean", col.weight = "mean", samples="rows"))
}


predict_mpm2 = function(model, data){
    logtrans = model$logtrans
    logrepl = model$logrepl
    center = model$center
    normal = model$normal
    closure = model$closure
    row.weight = model$row.weight
    col.weight = model$col.weight
    CW = model$CW
    RW = model$RW
    pos.row = model$pos.row
    pos.column =  model$pos.column
    samples = model$samples

    if(samples == "cols"){
        data = data[model$row.names,]
        for(r in 1:nrow(data)){
            if(any(is.na(data[r,]))){
                stop(paste("Row",as.character(r),"contains NAs"))
            }
            if(!any(as.numeric(data[r,]) != data[r,])){
                data[r,] = as.numeric(data[r,])
            } else {
                stop(paste("Could not coerce row",as.character(r),rownames(data)[r],"to numeric"))
            }
        }
        data = as.matrix(data)
    }
    if(samples == "rows"){
        data = data[,model$col.names,pos.column]
        for(i in 1:ncol(data)){
            if(any(is.na(data[,i]))){
                warning(paste("Column",as.character(i),"contains NAs"))
            }
            if(!any(as.numeric(data[,i]) != data[,i])){
                data[,i] = as.numeric(data[,i])
            } else {
                tmp = data.frame(numeric= as.numeric(data[,i]), previous = data[,i])
                str(tmp)
                print(which(tmp$numeric != tmp$previous))
                stop(paste("Could not coerce column",as.character(i),colnames(data)[i],"to numeric"))
            }
        }
        data = as.matrix(data)
    }
    col.names = colnames(data)
    if(is.null(col.names)){
        col.names = as.character(1:ncol(data))
    }
    row.names = rownames(data)
    if(is.null(row.names)){
        row.names = as.character(1:nrow(data))
    }

    pos.row = rep(FALSE, nrow(data))
    pos.column = rep(FALSE, ncol(data))

    if (missing(data))
        stop("Argument \"data\" is missing, with no default")
    NData = data
    if(any(is.na(NData))){
        stop("Data must be without NA's")
    }
    if (!is.numeric(NData)){
        stop("Data must be numeric")
    }
    if (length(RW) != nrow(NData) && samples=="cols")
        stop("Length of RW not equal to number of rows in table")
    if (length(CW) != ncol(NData) && samples=="rows")
        stop("Length of CW not equal to number of columns in table")

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

    # All new samples have the same weight,  variables weights are defined by the model
    if(samples == "rows"){
        RM <- rowMeans(NData)
        CM = model$Cm
        Wn <- rep(1, length(RM))
        Wp = model$Wp

    }
    if(samples == "cols"){
        RM = model$RM
        CM <- colMeans(NData)
        Wn = model$Wn
        Wp <- rep(1, length(CM))
    }
    # Closure: make sure the sum of the data == 1
    # These factors cannot be computed for new samples but assuming variable distribution is the same
    # the variable factors can be re-used
    if (closure != "none" && any(LData < 0))
        warning("Closure operation with non-positive data")
    if(samples == "rows"){
        Tn <- rowSums(LData[, ], na.rm = TRUE)
        Tp <- model$Tp
    }
    if(samples == "cols"){
        Tn <- model$Tn
        Tp <- colSums(LData[, ], na.rm = TRUE)
    }
    Tt = model$Tt

    ClData <- switch(closure,
                     none = LData,
                     row = sweep(LData, 1, Tn, "/"),
                     column = sweep(LData, 2, Tp, "/"),
                     global = LData/Tt,
                     double = Tt * sweep(sweep(LData, 1, Tn, "/"), 2, Tp, "/"))
    if (any(!is.finite(ClData)))
        stop("Division by 0 in closure operation")

    # Centering : make sure the the (weighted) mean of data == 0
    if(samples == "rows"){
        Mp <- model$Mp
        Mn <- rowSums(sweep(ClData, 2, Wp, "*"))
    }
    if(samples == "cols"){
        Mp <- colSums(sweep(ClData, 1, Wn, "*"))
        Mn <- model$Mn
    }
    Mg <- sum(Mp * Wp)

    CData <- switch(center,
                    double = Mg + sweep(sweep(ClData,2, Mp), 1, Mn),
                    row = sweep(ClData, 1, Mn),
                    column = sweep(ClData, 2, Mp),
                    global = ClData - Mg,
                    none = ClData)

    if(samples == "rows"){
        Vp <- model$Vp
        Vn <- rowSums(sweep(CData^2, 2, Wp, "*"))
    }
    if(samples == "cols"){
        Vp <- colSums(sweep(CData^2, 1, Wn, "*"))
        Vn <- model$Vn
    }
    Vg <- sum(Vp * Wp)
    SData <- switch(normal,
                    global = CData/sqrt(Vg),
                    row = sweep(CData, 1, sqrt(Vn), "/"),
                    column = sweep(CData, 2, sqrt(Vp),"/"),
                    none = CData)
    WData <- sweep(sweep(SData, 1, sqrt(Wn), "*"), 2, sqrt(Wp),"*")

    actualD = diag(1.0/model$SVD$d)
    result = WData %*% t(model$SVD$vt) %*% actualD
    new_SVD = model$SVD
    new_SVD$u = result
    result = as.data.frame(result)
    if(samples=="rows"){
        for(f in 1:ncol(result)){
            colnames(result)[f] = paste("V",as.character(f),sep="")
        }
    }
    if(samples=="cols"){
        for(f in 1:ncol(result)){
            rownames(result)[f] = paste("V",as.character(f),sep="")
        }
    }
    r <- list(TData = SData, row.names = row.names, col.names = col.names,
              closure = closure, center = center, normal = normal, logtrans = logtrans, logrepl = logrepl,
              row.weight = row.weight, col.weight = col.weight, eigen = model$eigen,
              contrib = model$contrib, Rm = RM, Cm = CM, Wn = Wn, Wp = Wp, Mp = Mp, Mn = Mn, Vn = Vn,CW=CW, RW=RW,Tn=Tn, Tp=Tp, Tt=Tt,
              SVD = new_SVD, pos.column = pos.column, pos.row = pos.row, samples=samples,
              call = match.call(), CData=CData, LData=LData, NData=NData, WData=WData)
    class(r) <- "mpm"

    return(list(result = result, mpm = r))
}

sma.coords = function(x,scale="uvc"){
    dim = 1:length(x$SVD$d)
    rot = rep(-1, length(dim))
    Wn <- x$Wn
    Wp <- x$Wp
    Z <- x$TData
    d <- x$SVD$d
    U <- x$SVD$u
    V <- t(x$SVD$vt)
    fact.scale <- switch(scale, singul = rep(0.5, 2), eigen = rep(1, 2), uvc = c(1, 0), uvr = c(0, 1))
    S <- sweep(crossprod(t(as.matrix(sweep(Z, 2, sqrt(Wp), "*"))), V[, dim]), 2, d[dim]^(fact.scale[1] - 1), "*")
    S <- S * sqrt(sum(x$eigen))/sqrt(sum(x$eigen)^fact.scale[1])
    S[, d[dim] < 1e-06] <- 0
    L <- sweep(crossprod(as.matrix(sweep(Z, 1, sqrt(Wn), "*")), U[, dim]), 2, d[dim]^(fact.scale[2] - 1), "*")
    L <- L * sqrt(sum(x$eigen))/sqrt(sum(x$eigen)^fact.scale[2])
    L[, d[dim] < 1e-06] <- 0
    S <- S * matrix(rot, ncol = ncol(S), nrow = nrow(S), byrow = TRUE)
    L <- L * matrix(rot, ncol = ncol(L), nrow = nrow(L), byrow = TRUE)
    return(list(S=as.data.frame(S),L=as.data.frame(L)))
}

save_files = function(source_data,target_data,target_colnames, file_prefix, complete.cases.only=T){
    for(cname in target_colnames){
        outpt = cbind(source_data, target_data[,cname])
        colnames(outpt)[ncol(outpt)]=as.character(cname)
        if(complete.cases.only){
            outpt = outpt[complete.cases(outpt),]
        }
        write.csv(outpt,paste(file_prefix,"_",as.character(cname),".csv",sep=""),row.names = F)
    }
}

numeric_cols = function(data){
    result = (1:ncol(data))*0
    result = as.logical(result)
    for(i in 1:ncol(data)){
        result[i] = is.numeric(data[,i])
    }
    return(colnames(data)[which(result)])
}

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


cvpredict = function(data, folds, training_method, ...){
    result = numeric(length = nrow(data))
    first = 1
    for(i in 1:folds){
        cat(paste(as.character(i),"\r"))
        last=as.integer(round(i * nrow(data)/(1.0*folds)))
        test_set = first:last
        if(first==1){
            training_set = (last+1):nrow(data)
        } else if(last==nrow(data)) {
            training_set = 1:(first-1)
        } else {
            training_set = c(1:(first-1),(last+1):nrow(data))
        }
        training_data = data[training_set,]
        test_data = data[test_set,]
        model = training_method(data=training_data, ...)
        tmp = predict(model, test_data)
        if(is.factor(tmp) && !is.factor(result)){
            result = factor(rep(NA,nrow(data)),levels = levels(tmp))
        }
        result[test_set] = tmp
        first = last+1
    }
    return(result)
}


train_multimodel = function(data, variables_to_predict, global_predictors, submodel_predictors, failure_vals){
    result = list()
    model_names= c("")
    for(target_var in variables_to_predict){
        submodels = list()
        model_names = names(submodel_predictors)
        predicted = data.frame(matrix(ncol=length(model_names)+1,nrow=nrow(data), dimnames=list(NULL, c(target_var,model_names))))
        for(submodel in model_names){
            selected_vars = c(global_predictors,submodel_predictors[[submodel]])
            indata = data[,c(target_var,selected_vars)]
            curmodel = randomForest(formula(paste(target_var,"~",".")),data = indata, na.action = na.omit)
            submodels[[submodel]] = list(model = curmodel, vars=selected_vars)
            predicted[,submodel] = predict(curmodel, indata)
        }
        predicted[,target_var]=data[,target_var]
        predicted[is.na(predicted)]=failure_vals[target_var]
        final_model = randomForest(formula(paste(target_var,"~",".")),data = predicted, na.action = na.omit)
        result[[target_var]] = list(submodels = submodels, finalmodel = final_model, modelnames = model_names)
    }
    if(is.null(names(failure_vals))){
        names(failure_vals) = variables_to_predict
    }
    result = list(RFs=result,NArep=failure_vals,varnames=variables_to_predict)
    class(result) = "RFMM"
    return(result)
}

train_MAGISTA = function(data, predictor_groups, targets, variables_to_ignore=c()){
    result = list()
    gnames = names(predictor_groups)
    if(is.null(gnames)){
        gnames = paste("Subset",1:length(predictor_groups),sep=".")
    }

    vargroup_rating = numeric(length = length(predictor_groups))
    names(vargroup_rating) = gnames
    for(group in gnames){
        curcols = predictor_groups[[group]]
        vargroup_rating[group] = sum(complete.cases(data[,curcols]))
    }
    possible_ratings = sort(unique(vargroup_rating))
    res = list()
    for(target in targets){
        curmodel = list()
        for(rating in possible_ratings){
            cat(paste(target,":",rating,"\n"))
            tcols = names(vargroup_rating)[vargroup_rating>=rating]
            variables = unique(unlist(predictor_groups[tcols]))
            for(badvar  in variables_to_ignore){
                variables=variables[variables!=badvar]
            }
            indata = data[,variables]
            rows2keep = complete.cases(indata)
            indata = indata[rows2keep,]
            transformation = mpm.pca_axes(indata)
            train.fit = mpm.pca_predict(transformation, indata)
            indata = cbind(indata, train.fit)
            indata = cbind(data[rows2keep,target],indata)
            colnames(indata)[1] = target
            rfm = randomForest(as.formula(paste(target,"~",".")),data = indata, na.action = na.omit, ntree = 200, maxnodes = 60)

            rf_contamination_model_cv = cv.RFpredict( as.formula(paste(target,"~",".")), data=indata, folds=10, na.action = na.omit, importance=T)
            actual = data[rows2keep, target]
            rf_correction = lm( actual ~ value, rf_contamination_model_cv, na.action = na.omit)

            submodel = list(name=paste("R",rating,sep=""),
                            rating=rating,
                            incols=variables,
                            resname=target,
                            pcaxes=transformation,
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
            print(sum(torun))
            if(sum(torun)>0){
                result$rating[torun] = rating
                indata = data[torun,]
                cname = names(ratings)[which(ratings==rating)]
                submodel = model[[target]][[cname]]
                train.fit = mpm.pca_predict(submodel$pcaxes, indata)
                indata = cbind(indata, train.fit)
                resvalues = predict(submodel$model, indata)
                df = data.frame(value = resvalues)
                resvalues = predict(submodel$correction, df)
                resvalues[resvalues < submodel$valuerange[1]] = submodel$valuerange[1]
                resvalues[resvalues > submodel$valuerange[2]] = submodel$valuerange[2]
                result[torun,target] = resvalues
            }
        }
    }
    return(result)
}
predict_multimodel = function(data, multimodel){
    result = data.frame(matrix(ncol=length(multimodel$varnames),nrow=nrow(data), dimnames=list(NULL, multimodel$varnames)))
    for(targetvar in multimodel$varnames){
        model_names = multimodel$RFs[[targetvar]]$modelnames
        predicted = data.frame(matrix(ncol=length(model_names),nrow=nrow(data), dimnames=list(NULL, model_names)))
        for(submodel_name in model_names){
            submodel = multimodel$RFs[[targetvar]]$submodels[[submodel_name]]
            indata = data[,submodel$vars]
            predicted[,submodel_name] = predict(submodel$model, data)
        }
        predicted[is.na(predicted)]=multimodel$NArep[targetvar]
        result[,targetvar] = predict(multimodel$RFs[[targetvar]]$finalmodel, predicted)
    }
    return(result)
}
predict.RFMM = function(multimodel, data){
    return(predict_multimodel(data,multimodel))
}
group_by_prefix = function(names,separator){
    result = list()
    for(name in names){
        pos = str_locate(name,separator)[1]
        if(!is.na(pos)){
            curprefix = substr(name,1,pos-1)
            result[[curprefix]]=c(result[[curprefix]],name)
        } else{
            result[["__noprefix__"]] = c(result[["__noprefix__"]],name)
        }
    }
    return(result)
}

cv.RFpredict = function(formula, data, folds, ...){
    input_order = sample(nrow(data))
    prediction = data.frame(value=numeric(length=nrow(data)))
    rows_per_subset = ceiling(nrow(data)/folds)
    for(i in 1:folds){
        selected_rows = logical(length=nrow(data))
        selected_rows[(rows_per_subset*(i-1)+1):min((rows_per_subset*(i)),nrow(data))] = TRUE
        test_subset = data[selected_rows,]
        train_subset = data[!selected_rows,]
        tmpmodel = randomForest(formula,train_subset, ...)
        tmpresult = predict(tmpmodel, test_subset)
        print(c(nrow(test_subset),nrow(train_subset),length(tmpresult)))
        prediction$value[selected_rows] = tmpresult
    }
    return(prediction)
}

classifier_accuracy = function(prediction, reference, average_per_class=F){
    if(!average_per_class){
        return(sum(prediction==reference, na.rm = T)/(1.0*length(reference)))
    } else {
        curavg = 0
        uniqvals = unique(reference)
        nclass = length(uniqvals)
        for(i in uniqvals){
            tmpacc = sum( (prediction==reference) * (reference==i), na.rm = T)/(1.0*sum(reference==i,na.rm=T))
            curavg = curavg + tmpacc
        }
        return(curavg/(1.0*nclass))
    }
}

train.nnns = function(data, target_var, expl_vars, nns_amount=1, normalize=T, optimize_nn_count = T, nn_count=sqrt(nrow(data))){
    if(optimize_nn_count == T){
        optimum_found=F
        low_step = as.integer(floor(nn_count/2))
        high_step = as.integer(nn_count)
        nn_count = as.integer(floor(nn_count))
        prediction = cvpredict(data=data,10,train.nnns,
                               target_var=target_var,
                               expl_vars = expl_vars,
                               nns_amount=nns_amount,
                               normalize=normalize,
                               optimize_nn_count = F,
                               nn_count = nn_count)
        cur_accuracy = classifier_accuracy(prediction, data[,target_var], average_per_class = T)
        print(c(nn_count,cur_accuracy, high_step, low_step))
        while(!optimum_found){
            if(low_step >= nn_count){
                low_step = nn_count - 1
            }
            if(high_step + nn_count > nrow(data)){
                high_step = nrow(data) - nn_count
            }
            prediction = cvpredict(data=data,10,train.nnns,
                                   target_var=target_var,
                                   expl_vars = expl_vars,
                                   nns_amount=nns_amount,
                                   normalize=normalize,
                                   optimize_nn_count = F,
                                   nn_count = nn_count - low_step)
            low_accuracy = classifier_accuracy(prediction, data[,target_var], average_per_class = T)
            prediction = cvpredict(data=data,10,train.nnns,
                                   target_var=target_var,
                                   expl_vars = expl_vars,
                                   nns_amount=nns_amount,
                                   normalize=normalize,
                                   optimize_nn_count = F,
                                   nn_count = nn_count + high_step)
            high_accuracy = classifier_accuracy(prediction, data[,target_var], average_per_class = T)
            if(high_accuracy > cur_accuracy && high_accuracy > low_accuracy){
                nn_count = nn_count + high_step
                cur_accuracy = high_accuracy
            } else {
                if(low_accuracy > cur_accuracy){
                    nn_count = nn_count - low_step
                    cur_accuracy = low_accuracy
                    if(low_step >= nn_count){
                        low_step = as.integer(floor(nn_count/2))
                    }
                } else {
                    low_step = as.integer(floor(low_step/2))
                    high_step = as.integer(floor(high_step/2))
                }
            }
            if(high_step<=1 && low_step<=1){
                optimum_found = T
            } else {
                if(low_step == 0){
                    low_step = 1
                }
            }
            print(c(nn_count,cur_accuracy, high_step, low_step))
        }
    }

    if(length(target_var) == nrow(data)){
        data = cbind(data.frame(TARGET_VAR = target_var),data)
        target_var = "TARGET_VAR"
    }
    if(is.factor(data[,target_var])){
        sum_func = factor_mode
        factor_var = T
    } else{
        sum_func = median
        factor_var = F
    }
    normby = numeric(length(expl_vars)) + 1.0
    names(normby) = expl_vars
    if(normalize){
        x = 0
        for(colname in expl_vars){
            x = x + 1
            normby[x] = sqrt(var(data[,colname], use="complete.obs"))
            data[,colname] = data[,colname]/normby[x]
        }
    }
    data = data[,c(target_var,expl_vars)]
    x = as.matrix(dist(data[,expl_vars], method = "manhattan"))
    new_vals = data[,target_var]
    for(i in 1:nrow(data)){
        nns = order(x[i,])[(1:(nns_amount+1))]
        new_vals[i] = sum_func(data[nns,target_var],na.rm = T)
    }
    data[,target_var]=new_vals
    result = list(vars=expl_vars, ref=data, tvar=target_var, stdev = normby, sum_func=sum_func, nns_count = nn_count)
    class(result)="nnns"
    return(result)
}

correspondance_table = function(pred, ref){
    if(length(pred)!=length(ref)){
        error("Lengths do not match")
        return(NULL)
    }
    result = matrix(nrow=length(unique(pred)),ncol=length(unique(ref)))
    refvals = unique(ref)
    predvals = unique(pred)
    colnames(result) = paste("ref",as.character(refvals))
    rownames(result) = paste("pred",as.character(predvals))
    names(refvals) = colnames(result)
    names(predvals) = rownames(result)
    for(cname in colnames(result)){
        for(rname in rownames(result)){
            result[rname,cname] = sum( (ref==refvals[cname]) * (pred==predvals[rname]) )
        }
    }
    return(result)

}

predict.nnns = function(model,data, nns_count=model$nns_count){
    x = 0
    data = data[,model$vars]
    for(colname in model$vars){
        x = x + 1
        data[,colname] = data[,colname]/model$stdev[x]
    }
    tmp = rbind(data[,model$vars], model$ref[,model$vars])
    x = as.matrix(dist(tmp, method = "manhattan"))
    x = x[1:nrow(data),(nrow(data)+1):(nrow(data)+nrow(model$ref))]
    if(is.factor(model$ref[,model$tvar])){
        new_vals = factor(x=rep(NA,nrow(data)),levels = levels(model$ref[,model$tvar]))
    }else {
        new_vals = numeric(nrow(data))
    }
    for(i in 1:nrow(data)){
        nns = order(x[i,])[1:nns_count]
        new_vals[i] = model$sum_func(model$ref[nns,model$tvar])
    }
    return(new_vals)
}

train.BQM = function(data, target_variable, predictors, bestvalue=max(data[,target_variable]), worstvalue=min(data[,target_variable]), maxvars=10, subdivide=NULL, test_data=NULL, extravar=NULL, add_correction=T){
    complete_cases = complete.cases(data[,predictors])

    rf_in = data[,c(target_variable,predictors)]

    # linear model correction based on cross-validation - using the exact same output but ignoring this step
    cvres=NULL
    if(class(data[,target_variable])=="factor" || class(data[,target_variable])=="character"){
        fact_mode = T
        if(bestvalue=="max" && !("max" %in% levels(data[,target_variable]))){
            bestvalue = levels(data[,target_variable])[length(levels(data[,target_variable]))]
        }
        if(bestvalue=="min" && !("min" %in% levels(data[,target_variable]))){
            bestvalue = levels(data[,target_variable])[1]
        }
        if(worstvalue=="max" && !("max" %in% levels(data[,target_variable]))){
            worstvalue = levels(data[,target_variable])[length(levels(data[,target_variable]))]
        }
        if(worstvalue=="min" && !("min" %in% levels(data[,target_variable]))){
            worstvalue = levels(data[,target_variable])[1]
        }
        lowvalue = worstvalue
        highvalue = bestvalue
        narep = worstvalue
        linear_correction = NULL
        # extrapolate NA values assuming NA <=> "2x worse that the worst non-NA"
        predreps = numeric(length = length(predictors))
        names(predreps) = predictors
        for(pred in predictors){
            curworstavg = mean(data[data[,target_variable]==worstvalue,pred],na.rm=T)
            curbestavg = mean(data[data[,target_variable]==bestvalue,pred],na.rm=T)
            predreps[pred] = 2*curworstavg - curbestavg
            rf_in[is.na(rf_in[,pred]),pred] = predreps[pred]
            if(!is.null(test_data)){
                test_data[is.na(test_data[,pred]),pred] = predreps[pred]
            }
        }
    } else {
        fact_mode = F
        #parse input
        if(bestvalue=="max"){
            bestvalue = max(data[,target_variable], na.rm = T)
        }
        if(bestvalue=="min"){
            bestvalue = min(data[,target_variable], na.rm = T)
        }
        if(worstvalue=="max"){
            worstvalue = max(data[,target_variable], na.rm = T)
        }
        if(worstvalue=="min"){
            worstvalue = min(data[,target_variable], na.rm = T)
        }
        if(bestvalue < worstvalue){
            lowvalue = bestvalue
            highvalue = worstvalue
        }
        else{
            lowvalue = worstvalue
            highvalue = bestvalue
        }
        narep = worstvalue - (bestvalue-worstvalue)

        identity_lm = lm(actual~predicted,data.frame(actual=rf_in[complete_cases,target_variable],predicted=rf_in[complete_cases,target_variable]))
        if(add_correction){
            cvpred = data.frame(
                predicted = cvpredict(
                    data = rf_in[,],
                    folds=10,
                    training_method = train.BQM,
                    target_variable = target_variable,
                    predictors = predictors,
                    bestvalue=bestvalue,
                    worstvalue=worstvalue,
                    maxvars=maxvars,
                    subdivide=subdivide,
                    test_data=test_data,
                    extravar=extravar,
                    add_correction=F
                ),
                actual = rf_in[,target_variable])
            linear_correction = lm(predicted~actual,cvpred)
            cvres = (cvpred$predicted -  linear_correction$coefficients[1]) / linear_correction$coefficients[2]
            cat('A')
        } else{
            linear_correction = identity_lm
        }
        predreps = numeric(length = length(predictors))
        names(predreps) = predictors
        # extrapolate NA values assuming NA <=> "2x worse that the worst non-NA"
        for(pred in predictors){
            lmodel = lm(as.formula(paste(pred,"~",target_variable)), data=rf_in, na.action = na.omit)
            tmp = data.frame(x=worstvalue)
            colnames(tmp)=c(target_variable)
            predreps[pred] = predict(lmodel,tmp)[1]
            rf_in[is.na(rf_in[,pred]),pred] = predreps[pred]
            if(!is.null(test_data)){
                test_data[is.na(test_data[,pred]),pred] = predreps[pred]
            }
        }
    }
    rf_in[is.na(rf_in[,target_variable]),target_variable] = narep

    # add PCA transform data
    axes = mpm.pca_axes(rf_in[complete_cases,predictors])
    newvars = mpm.pca_predict(axes, rf_in[,predictors])
    rf_in = cbind(rf_in,newvars)

    # use the same PCA transform for test data, if present
    if(!is.null(test_data)){
        newvars2 = mpm.pca_predict(axes, test_data[,predictors])
        test_data = cbind(test_data, newvars2)
    }

    invars = predictors
    predictors = colnames(rf_in[2:ncol(rf_in)])

    pre_model = NULL
    #if(!is.null(pre)){
    #    tmpdata = data.frame(x=data[,pre])
    #    colnames(tmpdata)[1] = pre
    #    tmpdata = cbind(tmpdata,rf_in)
    #    pre_formula = as.formula(paste(pre,"~",paste(collapse = "+", predictors)))
    #    pre_model = randomForest(pre_formula, tmpdata)
    #    tmp = cv.RFpredict(pre_formula, tmpdata, 3)
    #    colnames(tmp)[1] = paste("predicted_",pre,sep="")
    #   rf_in =  cbind(rf_in,tmp)
    #    predictors = colnames(rf_in[2:ncol(rf_in)])
    #}
    data = rf_in

    if(maxvars<=0){
        maxvars = length(predictors)
    }
    formula_to_use = as.formula(paste(target_variable,"~",paste(collapse = " + ", predictors)))

    has_divisions=F


    preliminary_model = randomForest(formula_to_use, data=rf_in[complete_cases,], na.action = na.omit, importance=T)

    imp = measure_importance(preliminary_model)
    topvars = important_variables(imp, k = min(maxvars,length(predictors)))

    nareps = numeric(length = length(topvars))
    names(nareps) = topvars

    if(fact_mode){
        for(pred in topvars){
            curworstavg = mean(data[data[,target_variable]==worstvalue,pred],na.rm=T)
            curbestavg = mean(data[data[,target_variable]==bestvalue,pred],na.rm=T)
            predreps[pred] = 2*curworstavg - curbestavg
            rf_in[is.na(rf_in[,pred]),pred] = predreps[pred]
            if(!is.null(test_data)){
                test_data[is.na(test_data[,pred]),pred] = predreps[pred]
            }
        }
    } else {
        for(pred in topvars){
            lmodel = lm(as.formula(paste(pred,"~",target_variable)), data=rf_in, na.action = na.omit)
            tmp = data.frame(x=worstvalue)
            colnames(tmp)=c(target_variable)
            nareps[pred] = predict(lmodel,data)
            rf_in[is.na(rf_in[,pred]),pred] = nareps[pred]
        }
    }
    new_formula_to_use = as.formula(paste(target_variable,"~",paste(collapse = " + ", topvars)))
    new_model = randomForest(new_formula_to_use, data=rf_in[complete_cases,])
    backup_model = randomForest(new_formula_to_use, data=rf_in)

    submodels = NULL
    acceptable_levels=NULL
    if(is.numeric(subdivide) && subdivide>0 && !fact_mode){
        cutoffs_low=numeric(length=subdivide+1.)
        cutoffs_high=numeric(length=subdivide+1.)
        has_divisions=T
        subdivide = as.integer(subdivide)
        # generate the overlaps - overlaps should cover half of each adjacent fragment
        overlap_size = (highvalue - lowvalue)/(subdivide+1.0)
        division_size = overlap_size*2
        div_id = as.integer((data[,target_variable]-lowvalue)/overlap_size)
        div_id[div_id>=subdivide] = subdivide

        subdiv_id1 = div_id-1
        subdiv_id1[subdiv_id1<0]=0
        subdiv_id2 = div_id
        subdiv_id2[subdiv_id2==subdivide] = subdivide-1

        acceptable_levels = 0:subdivide
        subdiv_id = factor(div_id,levels=acceptable_levels)

        #subdiv_id2=factor(subdiv_id2, levels=acceptable_levels)
        #subdiv_id1=factor(subdiv_id1, levels=acceptable_levels)

        cdata = cbind(data,subdiv_id)

        cmodel = randomForest(as.formula(paste("subdiv_id ~",paste(collapse ="+", predictors))), cdata, na.action = na.omit)
        #predictor_rating = cmodel$importance[,1]
        #topvars = rownames(cmodel$importance)[order(predictor_rating, decreasing=T)][1:3]
        #cmodel = randomForest(as.formula(paste("subdiv_id ~",paste(collapse ="+", topvars))), cdata, na.action = na.omit)

        submodels = list(rep(NA,length(acceptable_levels-1)))
        cpred = as.integer(as.character(predict(cmodel,cdata)))
        for(i in 1:(length(acceptable_levels)-1)){
            subset = as.logical((cpred==i) + (cpred==i-1))
            tmpdata = cdata[subset,]
            #lowervals = as.logical(cpred<i-1)
            #highervals = as.logical(cpred>i)
            #tmpdata[highervals,target_variable] = highvalue #max(cdata[subset,target_variable])
            #tmpdata[lowervals,target_variable] = lowvalue #min(cdata[subset,target_variable])
            tmp = randomForest(as.formula(paste(target_variable,"~",paste(collapse ="+", predictors))), tmpdata, na.action = na.omit)
            submodels[[i]] = tmp
        }


        if(is.null(test_data)){
            test_data = cdata
        } else{
            tdiv_id = as.integer((test_data[,target_variable]-lowvalue)/overlap_size)
            tdiv_id[tdiv_id>=subdivide] = subdivide


            subdiv_id1 = tdiv_id-1
            subdiv_id1[subdiv_id1<0]=0
            subdiv_id2 = tdiv_id
            subdiv_id2[subdiv_id2==subdivide] = subdivide-1

            #subdiv_id2=factor(subdiv_id2, levels=acceptable_levels)
            #subdiv_id1=factor(subdiv_id1, levels=acceptable_levels)
            subdiv_id = factor(tdiv_id,levels=acceptable_levels)
            test_data = cbind(test_data,subdiv_id)
        }

        classpred = predict(cmodel, test_data)
        cpred = as.integer(as.character(classpred))
        result = numeric(length = nrow(test_data))
        result[]=0
        for(i in 0:(length(acceptable_levels))){
            to_parse = as.logical((cpred==i)+(cpred==i-1))
            tmodel = i
            if(tmodel==0){
                tmodel = tmodel + 1
            }
            if(tmodel==length(acceptable_levels)){
                tmodel = tmodel - 1
            }
            curmodel = submodels[[tmodel]]
            result[to_parse] = result[to_parse] + as.numeric(as.character(predict(curmodel, test_data[to_parse,])))/2.0
        }
        xy = data.frame(x=test_data[,target_variable], y=result,z=classpred)

        confusion_matrix = matrix(nrow = subdivide+1, ncol = subdivide+1)
        colnames(confusion_matrix) = as.character(acceptable_levels)
        rownames(confusion_matrix) = as.character(acceptable_levels)
        for(ref in acceptable_levels){
            TPFN = (subdiv_id==ref)
            for(pred in acceptable_levels){
                TPFP = (classpred==pred)
                confusion_matrix[as.character(pred),as.character(ref)] = sum(as.logical(TPFN*TPFP),na.rm = T)
            }
        }
    } else {
        if(!is.null(subdivide) && !is.numeric(subdivide)){
            cdata = data
            colnames(cdata)[subdivide] = "subdiv_id"
            cdata["subdiv_id",] = as.factor(cdata["subdiv_id",])
            cmodel = randomForest(as.formula(paste("subdiv_id ~",paste(collapse ="+", predictors))))
        } else {
            cdata = data
            cmodel = 1
        }

    }
    result = list(model = new_model,
                  backup_model = backup_model,
                  targetvar=target_variable,
                  invars=invars,
                  vars = topvars,
                  nareps=nareps,
                  mpmaxes=axes,
                  narep = narep,
                  cmodel = cmodel,
                  submodels = submodels,
                  subdivisions=acceptable_levels,
                  has_divisions = has_divisions,
                  pre_model = pre_model,
                  limits=c(lowvalue,highvalue),
                  linear_correction = linear_correction,
                  fact_mode = fact_mode,
                  cvres=cvres)
    class(result) = "BQM"
    return(result)
}

predict.BQM = function(model,data, full_output=F, simple=F) {
    # reorder the data so that it matches the model
    fact_mode = model$fact_mode
    rf_in = data[, model$invars]
    for(varname in model$invars){
        rf_in[is.na(rf_in[,varname]),varname] = model$nareps[varname]
    }
    complete_cases= complete.cases(rf_in)
    extravars = mpm.pca_predict(model$mpmaxes, rf_in[,])
    rf_in = cbind(rf_in, extravars)
    if(!is.null(model$pre_model)){
        tmp = data.frame(value=predict(model$pre_model, rf_in))
        colnames(tmp)[1] = paste("predicted_",model$pre,sep="")
        cbind(rf_in,tmp)
    }
    rf_in = rf_in[,model$vars]
    for(varname in model$vars){
        rf_in[is.na(rf_in[,varname]),varname] = model$nareps[varname]
    }
    rf_out = predict(model$backup_model, rf_in)
    rf_out[complete_cases] = predict(model$model, rf_in[complete_cases,])
    rf_out[is.na(rf_out)]=model$narep

    if(model$has_divisions && !simple &&!fact_mode){
        classpred = predict(model$cmodel, rf_in)
        cpred = as.integer(as.character(classpred))
        rf_out = numeric(length = nrow(rf_in))
        rf_out[]=0
        for(i in 0:(length(model$subdivisions))){
            to_parse = as.logical((cpred==i)+(cpred==i-1))
            tmodel = i
            if(tmodel==0){
                tmodel = tmodel + 1
            }
            if(tmodel==length(model$subdivision)){
                tmodel = tmodel - 1
            }
            curmodel = model$submodels[[tmodel]]
            rf_out[to_parse] = rf_out[to_parse] + as.numeric(as.character(predict(curmodel, rf_in[to_parse,])))/2.0
        }
        if(full_output){
            result = data.frame(pred=rf_out, binclass = cpred)
        } else {
            result = rf_out
        }
    } else{
        if(full_output) {
            result = data.frame(pred=rf_out, binclass = 0)
            if(!is.null(model$pre_model)){
                result = cbind(result,tmp)
            }
        } else{
            result = rf_out
        }
    }
    if(full_output){
        tmp = data.frame(predicted=result$pred)
        #result$pred = predict(model$linear_correction, tmp)
        if(!is_null(model$linear_correction)) {
            result$pred = (tmp$predicted -  model$linear_correction$coefficients[1]) / model$linear_correction$coefficients[2]
        }
        if(!fact_mode){
            result$pred[result$pred<model$limits[1]] = model$limits[1]
            result$pred[result$pred>model$limits[2]] = model$limits[2]
        }
    } else {
        tmp = data.frame(predicted=rf_out)
        #result = predict(model$linear_correction, tmp)
        if(!is_null(model$linear_correction)) {
            rf_out = (tmp$predicted -  model$linear_correction$coefficients[1]) / model$linear_correction$coefficients[2]
        }
        if(!fact_mode){
            rf_out[rf_out<model$limits[1]] = model$limits[1]
            rf_out[rf_out>model$limits[2]] = model$limits[2]
        }
        result = rf_out
    }
    return(result)
}

train.glmclass = function(targetvar,predictors, data, mpm="none", naclass=NULL){
    result = list()
    result$pre=NULL
    if(!is.null(naclass)){
        result$NAclass = naclass
    } else{
        result$NAclass = levels(actual)[min(as.integer(actual))]
    }
    if(length(predictors)==1){
        if(is.numeric(predictors)){
            predictors=colnames(data)[predictors]
        }
    } else {
        predictors = colnames(data[,predictors])
    }
    result$invars=predictors
    data = data[, c(targetvar,predictors)]
    data = data[complete.cases(data),]
    if(mpm!="none"){
        for(f in predictors){
            data[,f] = as.numeric(data[,f])
            if(any(is.na(data[,f]))){
                stop(paste("Column",f,"contains NA's"))
            }
        }
    }
    formula = as.formula(paste(targetvar,"~",paste(collapse="+",predictors)))
    if(mpm=="sma"){
        result$pre = mpm.2(data[,predictors], row.weight = "mean", col.weight = "mean", samples="rows")
        data_data = sma.coords(predict_mpm2(result$pre,data)$mpm)$S
        formula = as.formula(paste(targetvar,"~",paste(collapse="+",colnames(data_data))))
        data = cbind(data[,targetvar], data_data)
        colnames(data)[1]=targetvar
    }
    result$model = cva.glmnet(formula, data = data, family="multinomial")
    class(result) = "GLMQC"
    return(result)
}
predict.GLMQC = function(model, data){
    data = data[,model$invars]
    ccs = complete.cases(data)
    data = data[ccs,]
    if(!is.null(model$pre)){
        data = sma.coords(predict_mpm2(model$pre,data)$mpm)$S
    }
    resultsUtils <- predict(model$model,newdata=data, type="response", alpha=1)
    idx2 <- apply(resultsUtils, c(1), maxidx)
    prediction <- colnames(resultsUtils)[idx2]
    totcases=length(ccs)
    NAcases = totcases-sum(ccs)
    result = c(rep(NA,NAcases),prediction)
    result[ccs]=prediction
    result[!ccs]=model$NAclass
    return(result)
}
train.svmclass = function(targetvar,predictors, data, mpm="none", naclass=NULL){
    result = list()
    result$pre=NULL
    if(!is.null(naclass)){
        result$NAclass = naclass
    } else{
        result$NAclass = levels(data[,targetvar])[min(as.integer(data[,targetvar]))]
    }
    result$outlvls=levels(data[,targetvar])
    if(length(predictors)==1){
        if(is.numeric(predictors)){
            predictors=colnames(data)[predictors]
        }
    } else {
        predictors = colnames(data[,predictors])
    }
    result$invars=predictors
    data = data[, c(targetvar,predictors)]
    data = data[complete.cases(data),]
    if(mpm!="none"){
        for(f in predictors){
            data[,f] = as.numeric(data[,f])
            if(any(is.na(data[,f]))){
                stop(paste("Column",f,"contains NA's"))
            }
        }
    }
    formula = as.formula(paste(targetvar,"~",paste(collapse="+",predictors)))
    if(mpm=="sma"){
        result$pre = mpm.2(data[,predictors], row.weight = "mean", col.weight = "mean", samples="rows")
        data_data = sma.coords(predict_mpm2(result$pre,data)$mpm)$S
        formula = as.formula(paste(targetvar,"~",paste(collapse="+",colnames(data_data))))
        data = cbind(data[,targetvar], data_data)
        colnames(data)[1]=targetvar
    }
    tune = tune.svm(formula, data= data, gamma=10^(-6:-1), cost=10^(1:4))
    result$model = svm(formula, data = data, gamma = tune$best.parameters$gamma,cost=tune$best.parameters$cost)
    class(result) = "SVMQC"
    return(result)
}
predict.SVMQC = function(model, data){
    data = data[,model$invars]
    ccs = complete.cases(data)
    data = data[ccs,]
    if(!is.null(model$pre)){
        data = sma.coords(predict_mpm2(model$pre,data)$mpm)$S
    }
    prediction=predict(model$model,data)
    NAclass = which(model$outlvls==model$NAclass)
    totcases=length(ccs)
    NAcases = totcases-sum(ccs)
    result = c(rep(NA,NAcases),prediction)
    result[ccs]=prediction
    result[!ccs]=NAclass
    result=factor(result,levels=1:length(model$outlvls),labels=model$outlvls)
    return(result)
}
train.nbclass = function(targetvar,predictors, data, mpm="none", naclass=NULL){
    result = list()
    result$pre=NULL
    if(!is.null(naclass)){
        result$NAclass = naclass
    } else{
        result$NAclass = levels(actual)[min(as.integer(actual))]
    }
    if(length(predictors)==1){
        if(is.numeric(predictors)){
            predictors=colnames(data)[predictors]
        }
    } else {
        predictors = colnames(data[,predictors])
    }
    result$invars=predictors
    data = data[, c(targetvar,predictors)]
    data = data[complete.cases(data),]
    if(mpm!="none"){
        for(f in predictors){
            data[,f] = as.numeric(data[,f])
            if(any(is.na(data[,f]))){
                stop(paste("Column",f,"contains NA's"))
            }
        }
    }
    formula = as.formula(paste(targetvar,"~",paste(collapse="+",predictors)))
    if(mpm=="sma"){
        result$pre = mpm.2(data[,predictors], row.weight = "mean", col.weight = "mean", samples="rows")
        data_data = sma.coords(predict_mpm2(result$pre,data)$mpm)$S
        formula = as.formula(paste(targetvar,"~",paste(collapse="+",colnames(data_data))))
        data = cbind(data[,targetvar], data_data)
        colnames(data)[1]=targetvar
    }
    result$model = naiveBayes(formula, data=data)
    class(result) = "NBQC"
    return(result)
}
predict.NBQC = function(model, data){
    data = data[,model$invars]
    ccs = complete.cases(data)
    data = data[ccs,]
    if(!is.null(model$pre)){
        data = sma.coords(predict_mpm2(model$pre,data)$mpm)$S
    }
    prediction = predict(model$model,data)
    totcases=length(ccs)
    NAcases = totcases-sum(ccs)
    result = c(rep(NA,NAcases),prediction)
    result[ccs]=prediction
    result[!ccs]=model$NAclass
    return(result)
}
train.knnclass = function(targetvar,predictors, data, mpm="none", naclass=NULL){
    result = list()
    result$pre=NULL
    if(!is.null(naclass)){
        result$NAclass = naclass
    } else{
        result$NAclass = levels(actual)[min(as.integer(actual))]
    }
    if(length(predictors)==1){
        if(is.numeric(predictors)){
            predictors=colnames(data)[predictors]
        }
    } else {
        predictors = colnames(data[,predictors])
    }
    result$invars=predictors
    data = data[, c(targetvar,predictors)]
    data = data[complete.cases(data),]
    if(mpm!="none"){
        for(f in predictors){
            data[,f] = as.numeric(data[,f])
            if(any(is.na(data[,f]))){
                stop(paste("Column",f,"contains NA's"))
            }
        }
    }
    formula = as.formula(paste(targetvar,"~",paste(collapse="+",predictors)))
    if(mpm=="sma"){
        result$pre = mpm.2(data[,predictors], row.weight = "mean", col.weight = "mean", samples="rows")
        data_data = sma.coords(predict_mpm2(result$pre,data)$mpm)$S
        formula = as.formula(paste(targetvar,"~",paste(collapse="+",colnames(data_data))))
        data = cbind(data[,targetvar], data_data)
        colnames(data)[1]=targetvar
    }
    result$model = knn3(formula, data=data)
    class(result) = "KNNQC"
    return(result)
}
predict.KNNQC = function(model, data){
    data = data[,model$invars]
    ccs = complete.cases(data)
    data = data[ccs,]
    if(!is.null(model$pre)){
        data = sma.coords(predict_mpm2(model$pre,data)$mpm)$S
    }
    prediction = predict(model$model,data)
    top_pick = character(nrow(prediction))
    for(i in 1:nrow(prediction)){
        top_pick[i] = which(prediction[i,]==max(prediction[i,]))
    }
    top_pick = factor(top_pick,levels=1:ncol(prediction),labels=colnames(prediction))
    totcases=length(ccs)
    NAcases = totcases-sum(ccs)
    result = c(rep(NA,NAcases),top_pick)
    result[ccs]=top_pick
    result[!ccs]=model$NAclass
    return(result)
}
train.RFclass = function(targetvar,predictors, data, mpm="none", naclass=NULL){
    result = list()
    result$pre=NULL
    if(!is.null(naclass)){
        result$NAclass = naclass
    } else{
        result$NAclass = levels(actual)[min(as.integer(actual))]
    }
    if(length(predictors)==1){
        if(is.numeric(predictors)){
            predictors=colnames(data)[predictors]
        }
    } else {
        predictors = colnames(data[,predictors])
    }
    result$invars=predictors
    data = data[, c(targetvar,predictors)]
    data = data[complete.cases(data),]
    if(mpm!="none"){
        for(f in predictors){
            data[,f] = as.numeric(data[,f])
            if(any(is.na(data[,f]))){
                stop(paste("Column",f,"contains NA's"))
            }
        }
    }
    formula = as.formula(paste(targetvar,"~",paste(collapse="+",predictors)))
    if(mpm=="sma"){
        result$pre = mpm.2(data[,predictors], row.weight = "mean", col.weight = "mean", samples="rows")
        data_data = sma.coords(predict_mpm2(result$pre,data)$mpm)$S
        formula = as.formula(paste(targetvar,"~",paste(collapse="+",colnames(data_data))))
        data = cbind(data[,targetvar], data_data)
        colnames(data)[1]=targetvar
    }
    result$model = randomForest(formula, data=data)
    class(result) = "RFQC"
    return(result)
}
predict.RFQC = function(model, data){
    data = data[,model$invars]
    ccs = complete.cases(data)
    data = data[ccs,]
    if(!is.null(model$pre)){
        data = sma.coords(predict_mpm2(model$pre,data)$mpm)$S
    }
    prediction = predict(model$model,data)
    totcases=length(ccs)
    NAcases = totcases-sum(ccs)
    result = c(rep(NA,NAcases),prediction)
    result[ccs]=prediction
    result[!ccs]=model$NAclass
    return(result)
}




