# tabulate cni spike data
library(rjson)
library(ggplot2)

countspikes<-function(file,thresh){
  
  require(rjson)
  qa<-fromJSON(file=file)
  #qa<-qa[[1]]
  zs<-qa$'timeseries zscore'
  spikes<-sapply(zs,function(x)sum(abs(x)>thresh))
  tps<-sapply(zs,function(x)sum(x!=0))
  return(list(spikes=spikes,tps=tps))
  
}

tabulatesubs<-function(sessions,thresh){
  
  dat<-data.frame(subj=numeric(),scan=numeric(),nspikes=numeric())
  
  for (s in 1:length(sessions)){
    sess<-sessions[s]
    scans<-as.list(list.files(sess,full.names=TRUE))
    scans<-scans[grepl('json',scans)]
    spikedat<-lapply(scans,function(x){tmp<-countspikes(x,thresh)}) # for now collapse across slices
    spikecount<-sapply(spikedat,function(x)sum(x$spikes))
    tpcount<-sapply(spikedat,function(x)max(x$tps))
    d<-cbind(s,scan=1:length(spikecount),spikecount,tpcount)
    dat<-rbind(dat,d)
  }
  dat$thresh<-thresh
  return(dat)
}

basedir<-'~/host/Work/Circmaze/QA/'
sessions<-list.files(basedir,full.names=TRUE)
sessions<-sessions[grepl('/cm',sessions)] # keep this constraint
sessions<-sessions[1:15] # remove this once files have .json extension
thresh6<-tabulatesubs(sessions,6)
thresh5<-tabulatesubs(sessions,5)
thresh4<-tabulatesubs(sessions,4)
thresh3<-tabulatesubs(sessions,3)

# plot scan x sub x thresh
master<-rbind(thresh3,thresh4,thresh5,thresh6)
master$spikeprop<-master$spikecount/master$tpcount

g<-ggplot(master,aes(x=scan,y=spikeprop)) + geom_bar(stat='identity') + facet_grid(s~thresh,labeller=label_both)
print(g)
ggsave('~/host/Work/Circmaze/QA/qa_scan_sub_thresh.pdf')

# plot sub x thresh
agg<-aggregate(cbind(spikecount,tpcount)~s+thresh,master,FUN=sum)
agg$spikeprop<-agg$spikecount/agg$tpcount

g<-ggplot(agg,aes(x=factor(thresh),y=spikeprop)) + geom_boxplot() + geom_jitter(aes(color=factor(s)),size=3)#geom_jitter is variant on geom_point
print(g)
ggsave('~/host/Work/Circmaze/QA/qa_sub_thresh.pdf')