
################################################################################
#####    Obtain windowed recombination rates from                         ######
#####   list of crossover events                                          ######
################################################################################

###Crossover data from https://doi.org/10.1371/journal.pgen.1011072.s010

thisfolder<-"./" ###path to local folder


recombinations<-read_excel(file.path(thisfolder,"JKK_recombination_events_IMlines.xlsx"),skip=2)
AIMs_counts<-read.table(file.path(thisfolder,"AIMs_July2025_panel33.AIMs_counts.txt"),
                        col.names=c("chr","pos","a1g","a2g","a1n","a2n"))
fai<-read.table(file.path(thisfolder,"Mguttatusvar_IM767_887_v2.0.fa.fai"),
                col.names=c("chr","chrsize","V3","V4","V5")) %>% select(chr,chrsize) %>%
  mutate(last_window_end=ceiling(chrsize/1000000)*1000000,
         lastwindow=paste0(chr,"_",format(last_window_end,scientific=F)),
         lastwindow_size=chrsize+1000000-last_window_end) %>% head(14)
                
get_windowsize<-function(windowID,base=1000000) {
  chromID=str_split_1(windowID,"_")[2]
  position=as.integer(str_split_1(windowID,"_")[3])
  chrom=paste0("Chr_",chromID)
  isEnd=windowID %in% fai$lastwindow
  if (isEnd) {
    return(fai %>% filter(chr==chrom) %>% pull(lastwindow_size))
  } else {
    return(base)
  }
}

recombinations_sizes<-recombinations %>% rename(left=`Marker left`,right=`Marker right`) %>% 
  select(Chromosome,left,right) %>%
  mutate(bracketsize=right-left,log10_bracketsize=log10(bracketsize),
         center=(right+left)/2) 

recombinations_sizes %>% 
  ggplot(aes(x=log10_bracketsize)) +
  geom_histogram() + theme_bw()
recombinations_sizes %>% 
  ggplot(aes(x=center)) + 
  facet_wrap(~Chromosome,ncol=2) + 
  geom_histogram(boundary=0,binwidth=1000000) + 
  theme_bw()

recombinations_windows<-recombinations_sizes %>% 
  mutate(bin_end=ceiling(center/1000000)*1000000) %>%
  rename("chr"="Chromosome") %>% 
  #mutate(chr=str_remove(Chromosome,"_")) %>%   ###if using fixed chr names
  group_by(chr,bin_end) %>% summarize(n_events=n()) %>%
  ungroup() %>%
  mutate(windowID=paste0(chr,"_",format(bin_end,scientific=F)),
         isEnd=(windowID %in% fai$lastwindow),
         windowsize=sapply(windowID,get_windowsize),
         events_per_bp=n_events/windowsize,
         M_per_bp=events_per_bp/(1373*2))

recombinations_windows_50kb<-recombinations_sizes %>% 
  mutate(bin_end=ceiling(center/50000)*50000) %>%
  rename("chr"="Chromosome") %>% 
  group_by(chr,bin_end) %>% summarize(n_events=n()) %>%
  ungroup() %>%
  mutate(windowID=paste0(chr,"_",format(bin_end,scientific=F)),
         isEnd=(windowID %in% fai$lastwindow),
         windowsize=sapply(windowID,get_windowsize,base=50000),
         events_per_bp=n_events/windowsize,
         M_per_bp=events_per_bp/(1373*2))

#recombinations_windows %>% ggplot(aes(x=bin_end,y=M_per_bp)) + 
#  facet_wrap(~Chromosome,nrow=7) + geom_line() + theme_bw()
#recombinations_windows_50kb %>% ggplot(aes(x=bin_end,y=M_per_bp)) + 
#  facet_wrap(~Chromosome,nrow=7) + geom_line() + theme_bw()

get_rate<-function(windowID) {
  return(recombinations_windows[recombinations_windows$windowID==windowID,]$M_per_bp)
}
windowrates<-recombinations_windows %>% ungroup() %>% select(windowID,M_per_bp)

Aims_counts_rates<-AIMs_counts %>% 
  mutate(previous=ifelse(!is.na(lag(chr)) & chr==lag(chr),lag(pos),0),
         size=pos-previous,
         center=(pos+previous)/2,
         window_end=ceiling(center/1000000)*1000000,
         windowID=paste0(chr,"_",format(window_end,scientific=F))) %>%
  left_join(windowrates) %>% 
  mutate(adj_rate=ifelse(is.na(M_per_bp),0.5/(1000000*1373*2),M_per_bp),
         M_distance=adj_rate*size)
         

Aims_counts_rates %>% select(chr,pos,a1g,a2g,a1n,a2n,M_distance) %>% 
  mutate(M_distance=format(M_distance,scientific=F,nsmall=10,width=12)) %>%
  write.table(file=file.path(thisfolder,"AIMs_July2025_panel33.AIMs_counts.realcMs.txt"),
              col.names=F,row.names=F,quote=F,sep="\t")

#Aims_counts_rates %>% group_by(chr,window_end) %>% summarize(n_AIMs=n()) %>%
#  ggplot(aes(x=window_end,y=n_AIMs)) + facet_wrap(~chr,ncol=2) + geom_col() + theme_bw()



