
################################################################################
#####    calculate local PCAs and summary statistics for each window      ######
################################################################################


library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
samplelist_file<-args[1]
regionslist_file<-args[2]
genotypes_file<-args[3]
prefix=args[4]

samplelist<-read.table(samplelist_file,header=F,col.names=c("sample"))
regionslist<-read.table(regionslist_file,header=F,col.names=c("chrom","windowstart","windowend","region"))
genotypes<-read.table(genotypes_file,header=F,col.names=c("region","sample","genotype","freq_group"))

PCAfit_stats<-data.frame()
for (r in 1:nrow(regionslist)) {
  regionline<-regionslist[r,]
  regionstring=paste0(regionline$chrom,".",regionline$windowstart,".",regionline$windowend)
  windowname=regionline$region
  PCAfile=paste0(prefix,".",regionstring,".cov")
  PCAeigen<-eigen(read.table(PCAfile,header = F))
  PCAeigen_vecs<-bind_cols(samplelist,
                           as.data.frame(PCAeigen$vectors[,1:2]) %>% rename_with(function(x){str_replace(x,"V","PC")})) %>%
    left_join(genotypes %>% filter(region==windowname)) %>%
    mutate(genotype_value=genotype,genotype=factor(genotype))
  PCAeigen_vals<-round(PCAeigen$values[1:2]/sum(PCAeigen$values)*100,2)
  ###plot PCA
  PCAplot<-PCAeigen_vecs %>% ggplot(aes(x=PC1,y=PC2,color=genotype,shape=freq_group)) + 
    geom_point() + theme_bw() + 
    xlab(paste0("PC1: ",PCAeigen_vals[1],"%")) + 
    ylab(paste0("PC2: ",PCAeigen_vals[2],"%")) + 
    ggtitle(regionstring)
  ggsave(paste0(PCAfile,".plot.pdf"),PCAplot,device = "pdf",width = 3,height=3,units = "in",dpi=600)
  PCAeigen_vecs<-PCAeigen_vecs %>% filter(!is.na(genotype))
  
  ###test 1: linear model fit
  model<-lm(PC1~genotype_value,data=PCAeigen_vecs)
  r.squared<-summary(model)$r.squared
  
  ###test 2: heterozygote fit
  mean_g0<-PCAeigen_vecs %>% filter(genotype==0) %>% summarize(mean=mean(PC1)) %>% pull(mean)
  mean_g2<-PCAeigen_vecs %>% filter(genotype==2) %>% summarize(mean=mean(PC1)) %>% pull(mean)
  center=(mean_g0+mean_g2)/2
  hets_scaled<-PCAeigen_vecs %>% filter(genotype==1) %>% 
    mutate(scaled=(PC1-center)/(mean_g2-center)) %>%
    summarize(mean_scaled=mean(scaled),
              var_scaled=var(scaled)) %>%
    mutate(zscore=mean_scaled/(sqrt(var_scaled)))
  
  ###test 3: fit of introgressed vs. non-introgressed samples
  gut_nons<-PCAeigen_vecs %>% filter(genotype==0 & freq_group=="guttatus") %>%
    mutate(deviation=(PC1-mean_g0)^2)
  nas_nons<-PCAeigen_vecs %>% filter(genotype==2 & freq_group=="nasutus") %>%
    mutate(deviation=(PC1-mean_g2)^2)
  intro.g0s<-PCAeigen_vecs %>% filter(genotype==0 & freq_group=="nasutus") %>%
    mutate(deviation=(PC1-mean_g0)^2)
  intro.g1s<-PCAeigen_vecs %>% filter(genotype==1 & freq_group %in% c("guttatus","nasutus")) %>% 
    mutate(deviation=(PC1-center)^2)
  intro.g2s<-PCAeigen_vecs %>% filter(genotype==2 & freq_group=="guttatus") %>%
    mutate(deviation=(PC1-mean_g2)^2)
  all_nons<-bind_rows(gut_nons,nas_nons) %>% summarize(sum.sq=sum(deviation),n=n()) %>% 
    mutate(variance=sum.sq/(n-1))
  all_intros<-bind_rows(intro.g0s,intro.g1s,intro.g2s) %>% summarize(sum.sq=sum(deviation),n=n()) %>% 
    mutate(variance=sum.sq/(n-1))
  
  PCAfit_stats<-PCAfit_stats %>% bind_rows(data.frame(region=regionstring,modelfit.r.squared=r.squared,
                                                      het.scaled.mean=hets_scaled$mean_scaled,
                                                      het.scaled.zscore=hets_scaled$zscore,
                                                      non.introgressed.variance=all_nons$variance,
                                                      introgressed.variance=all_intros$variance,
                                                      variance_factor=all_intros$variance/all_nons$variance))
}
write.table(PCAfit_stats,"PCA_fit_stats.txt",quote = F,row.names=F,col.names=T,sep="\t")
