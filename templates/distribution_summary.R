suppressWarnings(suppressMessages(library(moments)))
suppressWarnings(suppressMessages(library(dplyr)))

args = commandArgs(trailingOnly=TRUE)

if( length(args)==2) {
 write_summary = function(filename) {
  in_data = data.frame(distance=c(1.0))
  measurable = F
  result = tryCatch({
      in_data = read.csv(filename, header=T)
  }, error = function(e) {
      return(integer(0))
  })
  if(length(result)==0){
    cat("Input data could not be read\n")
  } else {
    if (nrow(in_data)<3){
      cat("Input data has less than 3 rows\n")
    } else {
      measurable = T
    }
  }
  if (measurable) {
    variable = in_data[,1]
    result = data.frame(
      dummy=c(""),
      name=c(filename),
      ncmp=c(nrow(in_data)),
      avg=c(mean(in_data[,1])),
      sdev=c(sd(in_data[,1])),
      skew=c(skewness(in_data[,1])),
      kurt=c(kurtosis(in_data[,1])),
      med=c(median(in_data[,1])),
      p025=c(quantile(in_data[,1],0.025)),
      p975=c(quantile(in_data[,1],0.975)),
      p050=c(quantile(in_data[,1],0.050)),
      p950=c(quantile(in_data[,1],0.950)),
      p100=c(quantile(in_data[,1],0.100)),
      p900=c(quantile(in_data[,1],0.900)),
      stringsAsFactors = FALSE
      )
  } else {
    result = data.frame(
      dummy=c(""),
      name=c(filename),
      ncmp=c(0),
      avg=c(NA),
      sdev=c(NA),
      skew=c(NA),
      kurt=c(NA),
      med=c(NA),
      p025=c(NA),
      p975=c(NA),
      p050=c(NA),
      p950=c(NA),
      p100=c(NA),
      p900=c(NA)
      )
  }
  return(result)
 }
 x=data.frame()
 for( line in readLines(args[1]) ){
  x = bind_rows(x, write_summary(line))
 }
 write.csv(x[2:ncol(x)],args[2], row.names=F)
} else {
 cat("Usage: Rscript <name of the script> <list of distribs> <output file name>\n")
}

