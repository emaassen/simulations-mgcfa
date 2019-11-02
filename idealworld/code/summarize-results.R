# if results.summarized is not saved correctly ----------------------------
# we need to load all individual files and summarize those.
# rm(list = ls())
# options(scipen=999)
# x = c("data.table")
# sapply(x,library,character.only=T)
# results.summarized <- c()
# 
# # Load all individual result files 
# filelist = list.files(pattern = "result_")
# datalist = lapply(filelist, fread)
# 
# for (i in seq_along(datalist)){
#   
#   # assign NA to string variables first, otherwise mean cannot be calculated
#   datalist[[i]] <- sapply(datalist[[i]], as.numeric)
#   
#   # take mean of all iterations per condition
#   results.summarized.new <- apply(datalist[[i]],2,mean,na.rm=T)
#   
#   # assign it to summarized results
#   results.summarized <- rbind(results.summarized,results.summarized.new)
#   
# }
# 
# # assign column names to results.summarized and sort by cond
# colnames(results.summarized) <- colnames(datalist[[1]]) 
# results.summarized <- results.summarized[order(results.summarized[,1]),]
# 
# # write summarized results to file
# write.table(results.summarized, file = "idealworld_summarized_results1.dat", row.names=F, col.names=T)
