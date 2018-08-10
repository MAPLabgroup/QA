# tabulate cni spike data
library(rjson)

setwd("~/host/Work/Circmaze/QA/")

spiketps<-function(file,thresh){
  
  require(rjson)
  qa<-fromJSON(file=file)
  #qa<-qa[[1]]
  zs<-qa$'timeseries zscore'
  tps<-lapply(zs,function(x)which(abs(x)>thresh))
  tps<-sort(unique(unlist(tps)))
  return(tps)
}

basedir<-'~/host/Work/Circmaze/QA/cm10'
thresh<-4

sessions<-list.files(basedir,full.names=TRUE)
spikes<-lapply(sessions,function(x)spiketps(x,thresh))
names(spikes)<-sessions
# this should be a list of your spikes for each scan, in the order that the files are specified in the sessions variable

setwd(basedir)
sink('test_spike_timepoints.txt') 
spikes
sink(file=NULL, append=FALSE)