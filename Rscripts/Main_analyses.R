################################################################################
#####    Main analysis script for core manuscript analyses                ######
################################################################################

thisfolder<-"./" ###input/output folder


###libraries
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(car)
library(rcartocolor)

###Load in windowed ancestry data
Chrlist<-paste0("Chr",str_pad(c(1:14),2,pad="0"))
windowed_ancestry_allsamples<-read.table(file.path(thisfolder,"ancestry_IM767v2_withpanel/ancestry-probs_allchrs.tsv_rec.windowed.txt_transposed"),header=T)
windowed_ancestry_allsamples_long<-windowed_ancestry_allsamples %>% 
  pivot_longer(cols=-c(chrom.windowstart.windowend,nsites),
               names_to="sample",values_to="genotype") %>%
  separate(chrom.windowstart.windowend,into=c("chrom","windowstart","windowend"),sep=":",remove=F) %>%
  mutate(windowstart=as.integer(windowstart),windowend=as.integer(windowend),
         chrom=factor(chrom,levels=Chrlist))

###Load in global ancestry summaries
ancestry_globals<-read.table(file.path(thisfolder,"ancestry_IM767v2_withpanel/ancestry-probs_allchrs.tsv_rec.samplesummary.txt"),header=T) %>%
  mutate(sampleID=str_remove(sample,"_read_1.fastq")) %>%
  filter(sampleID!="KK047") ##filter one bad sample by name
filtered_globals<-ancestry_globals %>% filter(counted_sites/(counted_sites+missing_sites)>=0.25)

###ancestry_windowed_summaries
windowed_globals<-windowed_ancestry_allsamples_long %>% group_by(sample) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called_windows=n_called/(n_called+n_missing),
         nasutus_proportion_windows=genosum/(2*n_called),
         het_prop_windows=n_1/n_called) %>%
  mutate(sampleID=str_remove(sample,"_read_1.fastq")) %>%
  mutate(sampleID=str_replace_all(sampleID,fixed("."),"_")) %>%
  filter(sampleID!="KK047") ##filter one bad sample by name

###Load in angsd data
angsd_cov<-read.table(file.path(thisfolder,"angsd/genolike.25percent_calledsites.IM767_v2.IM767panel.highqualSNPs.thinned1kb.allChrs.cov"))
angsd_eigen<-eigen(angsd_cov)
angsd_bamkey<-read.table(file.path(thisfolder,"angsd/bamlist_25percent_calledsites.txt"),col.names="bampath") %>%
  separate(bampath,into=c("V1","V2","V3","V4","V5","bamfile"),sep="/") %>% 
  mutate(sampleID=str_remove(bamfile,".IM767v2.frdms.bam")) %>% 
  mutate(sampleID=str_replace_all(sampleID,"-","_")) %>%
  select(sampleID)
angsd_PCA<-bind_cols(angsd_bamkey,angsd_eigen$vectors[,1:10]) %>% filter(sampleID!="KK047")
colnames(angsd_PCA)<-c("sampleID",paste0("PC",c(1:10)))

###remember to check/exclude sookensis

###Load in ancestry transitions
transitions_key<-read.table(file.path(thisfolder,"ancestry_IM767v2_withpanel/ancestry_withpanel_samplelist.txt"),header=F,col.names="sample") %>%
  mutate(sampleID=str_remove(sample,"_read_1.fastq")) %>% 
  mutate(sampleID=str_replace_all(sampleID,"-","_")) %>%
  select(sampleID) %>% mutate(IDnum=c(1:nrow(.)))
transitions<-read.table(file.path(thisfolder,"ancestry_IM767v2_withpanel/ancestry-probs_allchrs.tsv_rec.txt_ancestrytransitions_allchrs"),
                        col.names=c("chrom","lowerbound","upperbound","IDnum")) %>%
  left_join(transitions_key,by="IDnum")
transitions_summary<-transitions %>% group_by(sampleID,chrom) %>% summarize(n=n()) %>%
  mutate(Chr_all=sum(n,na.rm=T)) %>%
  pivot_wider(names_from=chrom,values_from=n)
transitions_summary[is.na(transitions_summary)]<-0
transitions_summary_filtered<-transitions_summary %>% filter(sampleID %in% filtered_globals$sampleID) %>%
  rename_with(function(x){paste0("Transitions_",x)},.cols=-sampleID)

###Load in JKK recombination data and gene density info
recombinations<-read_csv(file.path(thisfolder,"JKK_recombination_events_IMlines.csv"),skip=2)
fai<-read.table(file.path(thisfolder,"Mguttatusvar_IM767_887_v2.0.fixed.fa.fai"),
                col.names=c("chr","chrsize","V3","V4","V5")) %>% select(chr,chrsize) %>%
  mutate(last_window_end=ceiling(chrsize/1000000)*1000000,
         lastwindow=paste0(chr,"_",format(last_window_end,scientific=F)),
         lastwindow_size=chrsize+1000000-last_window_end)

get_windowsize<-function(windowID,base=1000000) {
  chrom=str_split_1(windowID,"_")[1]
  position=as.integer(str_split_1(windowID,"_")[2])
  isEnd=windowID %in% fai$lastwindow
  if (isEnd) {
    return(fai[fai$chr==chrom,]$lastwindow_size)
  } else {
    return(base)
  }
}

recombinations_sizes<-recombinations %>% rename(left=`Marker left`,right=`Marker right`) %>% 
  select(Chromosome,left,right) %>%
  mutate(bracketsize=right-left,log10_bracketsize=log10(bracketsize),
         center=(right+left)/2) 

recombinations_windows<-recombinations_sizes %>% 
  mutate(bin_end=ceiling(center/1000000)*1000000) %>%
  group_by(Chromosome,bin_end) %>% summarize(n_events=n()) %>%
  mutate(chr=str_remove(Chromosome,"_"),
         windowID=paste0(chr,"_",format(bin_end,scientific=F)),
         isEnd=(windowID %in% fai$lastwindow),
         windowsize=sapply(windowID,get_windowsize),
         events_per_bp=n_events/windowsize,
         M_per_bp=events_per_bp/(1373*2))

recombinations_windows_50kb<-recombinations_sizes %>% 
  mutate(bin_end=ceiling(center/50000)*50000) %>%
  group_by(Chromosome,bin_end) %>% summarize(n_events=n()) %>%
  mutate(chr=str_remove(Chromosome,"_"),
         windowID=paste0(chr,"_",format(bin_end,scientific=F)),
         isEnd=(windowID %in% fai$lastwindow),
         windowsize=sapply(windowID,get_windowsize,base=50000),
         events_per_bp=n_events/windowsize,
         M_per_bp=events_per_bp/(1373*2))

genecounts<-read.table(file.path(thisfolder,"Mguttatusvar_IM767_887_v2.1.gene.gff3"),header=F) %>%
  filter(V3=="gene") %>% mutate(chr=str_remove(V1,"_"),windowend=ceiling(V4/50000)*50000) %>% select(chr,windowend) %>%
  group_by(chr,windowend) %>% summarize(ngenes=n())


###combine sample-level data and categorize
allopatric_DPRarea_pops=c("BFR","COP","GCH","MOC","RCF","RHI")
highelev_pops=c("HHT","HHR","TUO")

sampleset_info<-windowed_globals %>% 
  left_join(angsd_PCA,by="sampleID") %>% left_join(transitions_summary_filtered,by="sampleID") %>%
  mutate(group=case_when(
    str_ends(sampleID,"_sampled") & nasutus_proportion_windows > 0.5 ~ "panel_NAS",
    str_ends(sampleID,"_sampled") & str_starts(sampleID,"CAC") ~ "panel_CAC",
    str_ends(sampleID,"_sampled") ~ "panel_GUT",
    str_starts(sampleID,"CAC") | str_starts(sampleID,"Census") | str_starts(sampleID,"KK") | str_starts(sampleID,"Mom") ~ "CAC",
    str_starts(sampleID,"LM") ~ "LM",
    str_sub(sampleID,1,3) %in% allopatric_DPRarea_pops ~ "Southern_allopatric",
    str_sub(sampleID,1,3) %in% highelev_pops ~ "Southern_highelev",
    .default = "Southern_sympatric"
  ))

sampleset_info %>% group_by(group) %>% summarize(n())   

F1_like<-sampleset_info %>% filter(het_prop_windows>0.75,frac_called_windows>=0.25)
NAS<-sampleset_info %>% filter(nasutus_proportion_windows>0.85,frac_called_windows>=0.25)

enumerate_windows<-windowed_ancestry_allsamples %>% select(chrom.windowstart.windowend) %>% mutate(window_number=c(1:nrow(.))) 
chrom_blocks<-enumerate_windows %>%
  separate(chrom.windowstart.windowend,into=c("chrom","windowstart","windowend"),sep=":",remove=F) %>%
  mutate(blockend=(chrom!=lead(chrom))) %>% filter(blockend)
  
windowed_ancestry_allsamples_long_enumerated<-windowed_ancestry_allsamples_long %>% 
  left_join(enumerate_windows,by="chrom.windowstart.windowend")

colors_sp<-c("#993929","#CC892F","#FFDA35") #G,A,N

paint_ancestry<-function(data,chromends,write_sampleIDs=T) {
  paintplot<-data %>%
    mutate(genotype=factor(genotype,levels=c(0,1,2))) %>%
    ggplot(aes(x=window_number,y=sample,fill=genotype)) + 
    geom_tile() + 
    geom_vline(data=chromends,mapping=aes(xintercept=window_number),color="blue") + 
    scale_fill_manual(values=colors_sp) + 
    theme_bw() + 
    scale_x_continuous(expand = c(0,0))
  if (write_sampleIDs) {
    return(paintplot)
  } else {
    return(paintplot + scale_y_discrete(labels=NULL))
  }
}

###filtering sites with bad calls in NAS individuals
nas_site_summary<-windowed_ancestry_allsamples_long %>% filter(sample %in% NAS$sample) %>% 
  group_by(chrom,windowstart,windowend,chrom.windowstart.windowend) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called=n_called/(n_called+n_missing),
         nasutus_freq=genosum/(2*n_called),
         het_freq=n_1/n_called)

#mask_sites_nas<-nas_site_summary %>% filter(frac_called<0.5 || nasutus_freq<0.25) %>% pull(chrom.windowstart.windowend)
#keep_sites_nas<-nas_site_summary %>% filter(frac_called>=0.5 && nasutus_freq>=0.25) %>% pull(chrom.windowstart.windowend)

sookensis_exclusions<-F1_like %>% filter(group=="CAC",str_detect(sample,"2022")) %>% pull(sample)
full_site_summary<-windowed_ancestry_allsamples_long %>% filter(!(sample %in% sookensis_exclusions)) %>% 
  group_by(chrom,windowstart,windowend,chrom.windowstart.windowend) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called=n_called/(n_called+n_missing),
         nasutus_freq=genosum/(2*n_called),
         het_freq=n_1/n_called)

mask_sites_full<-full_site_summary %>% filter(frac_called<0.5) %>% pull(chrom.windowstart.windowend)
keep_sites_full<-full_site_summary %>% filter(frac_called>=0.5) %>% pull(chrom.windowstart.windowend)

###make final dataset
windowed_sitefilter<-windowed_ancestry_allsamples_long %>% filter(chrom.windowstart.windowend %in% keep_sites_full)

###get individual summaries from site-filtered data
individuals_sitefiltered<-windowed_sitefilter %>% group_by(sample) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called_windows=n_called/(n_called+n_missing),
         nasutus_proportion_windows=genosum/(2*n_called),
         het_prop_windows=n_1/n_called) %>%
  left_join(select(sampleset_info,sample,sampleID,group))

###check sookensis-looking samples
prefilter_list<-individuals_sitefiltered %>% 
  filter(sample!="KK047_read_1.fastq", ###two data runs don't match
         #sample!="DPR_plot2I_BL_read_1.fastq", ###clusters with northern samples
         !(sample %in% sookensis_exclusions),
         !(str_detect(group,"panel")))
final_individual_list<-individuals_sitefiltered %>% 
  filter(sample!="KK047_read_1.fastq", ###two data runs don't match
         sample!="DPR_plot2I_BL_read_1.fastq", ###clusters with northern samples
         !(sample %in% sookensis_exclusions),
         !(str_detect(group,"panel")),
         frac_called_windows>=0.5) %>% 
  left_join(angsd_PCA,by="sampleID") %>% left_join(transitions_summary_filtered,by="sampleID") %>%
  mutate(freq_group=ifelse(nasutus_proportion_windows<0.15,"guttatus",
                           ifelse(nasutus_proportion_windows>0.85,"nasutus",
                                  "hybrid")))


###final windows
final_window_set<-windowed_sitefilter %>% filter(sample %in% final_individual_list$sample)
enumerate_final_windows<-final_window_set %>% select(chrom.windowstart.windowend) %>% unique() %>% mutate(window_number=c(1:nrow(.)))
chrom_blocks_final_windows<-enumerate_final_windows %>%
  separate(chrom.windowstart.windowend,into=c("chrom","windowstart","windowend"),sep=":",remove=F) %>%
  mutate(blockend=(chrom!=lead(chrom))) %>% filter(blockend)
final_window_set_enumerated<-final_window_set %>% left_join(enumerate_final_windows) %>% left_join(final_individual_list %>% select(sample,group))

group_paints<-list()
samplelists<-list()
for (i in unique(final_window_set_enumerated$group)) {
  samplelists[[i]]<-final_individual_list %>% filter(group==i) %>% 
    arrange(nasutus_proportion_windows) %>% 
    mutate(sample_ordered=factor(sample,levels=sample)) %>% select(sample_ordered)
  subdata<-final_window_set_enumerated %>% filter(group==i) %>% mutate(sample=factor(sample,levels=samplelists[[i]]$sample_ordered))
  group_paints[[i]]<-paint_ancestry(subdata,chrom_blocks_final_windows,write_sampleIDs = F)
}

###define key groups
key_groups<-final_individual_list %>% group_by(group,freq_group) %>% summarize(n=n()) %>%
  mutate(keygroup=paste0(group,"_",freq_group,"_",n))
focal_nine<-key_groups %>% filter(n>=10)

###get_keygroup_frequencies
all_groups_window_freqs<-data.frame()
for (g in key_groups$keygroup) {
  thispop=key_groups %>% filter(keygroup==g) %>% pull(group)
  thisfreq=key_groups %>% filter(keygroup==g) %>% pull(freq_group)
  these_individuals<-final_individual_list %>% filter(group==thispop,freq_group==thisfreq) %>% pull(sample)
  window_freqs<-final_window_set_enumerated %>% filter(sample %in% these_individuals) %>% 
    group_by(chrom.windowstart.windowend,chrom,windowstart,windowend,window_number) %>%
    summarize(n_called=sum(!is.na(genotype)),
              n_missing=sum(is.na(genotype)),
              n_0=sum(genotype==0,na.rm=T),
              n_1=sum(genotype==1,na.rm=T),
              n_2=sum(genotype==2,na.rm=T),
              genosum=sum(genotype,na.rm=T)) %>%
    mutate(frac_called_samples=n_called/(n_called+n_missing),
           nasutus_frequency=genosum/(2*n_called),
           het_frequency=n_1/n_called,
           keygroup=g)
  all_groups_window_freqs<-bind_rows(all_groups_window_freqs,window_freqs)
}
panel_GUT_list<-sampleset_info %>% filter(group=="panel_GUT",nasutus_proportion_windows<0.01)
panel_GUT_list_includehybrids<-sampleset_info %>% filter(group=="panel_GUT")
panel_GUT_freqs<-windowed_ancestry_allsamples_long %>% 
  filter(chrom.windowstart.windowend %in% final_window_set$chrom.windowstart.windowend,
         sample %in% panel_GUT_list$sample) %>% 
  left_join(enumerate_final_windows) %>% 
  group_by(chrom.windowstart.windowend,chrom,windowstart,windowend,window_number) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called_samples=n_called/(n_called+n_missing),
         nasutus_frequency=genosum/(2*n_called),
         het_frequency=n_1/n_called,
         keygroup="panel_GUT")

key_groups_freqorder<-all_groups_window_freqs %>% filter(keygroup %in% focal_nine$keygroup) %>% 
  group_by(keygroup) %>% summarize(mean_freq=mean(nasutus_frequency,na.rm=T)) %>%
  arrange(mean_freq) %>% mutate(keygroup=factor(keygroup,levels=keygroup)) 

###only 'truly' allopatric pops
true_allos<-final_individual_list %>% filter(group=="Southern_allopatric") %>% mutate(Pop=str_sub(sampleID,1,3)) %>%
  filter(Pop %in% c("COP","GCH","RHI"))
pseudo_allos<-final_individual_list %>% filter(group=="Southern_allopatric") %>% mutate(Pop=str_sub(sampleID,1,3)) %>%
  filter(Pop %in% c("MOC","RCF","BFR"))
true_allos_freqs<-windowed_ancestry_allsamples_long %>% 
  filter(chrom.windowstart.windowend %in% final_window_set$chrom.windowstart.windowend,
         sample %in% true_allos$sample) %>% 
  left_join(enumerate_final_windows) %>% 
  group_by(chrom.windowstart.windowend,chrom,windowstart,windowend,window_number) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called_samples=n_called/(n_called+n_missing),
         nasutus_frequency=genosum/(2*n_called),
         het_frequency=n_1/n_called,
         keygroup="Southern_true_allopatric")
pseudo_allos_freqs<-windowed_ancestry_allsamples_long %>% 
  filter(chrom.windowstart.windowend %in% final_window_set$chrom.windowstart.windowend,
         sample %in% pseudo_allos$sample) %>% 
  left_join(enumerate_final_windows) %>% 
  group_by(chrom.windowstart.windowend,chrom,windowstart,windowend,window_number) %>%
  summarize(n_called=sum(!is.na(genotype)),
            n_missing=sum(is.na(genotype)),
            n_0=sum(genotype==0,na.rm=T),
            n_1=sum(genotype==1,na.rm=T),
            n_2=sum(genotype==2,na.rm=T),
            genosum=sum(genotype,na.rm=T)) %>%
  mutate(frac_called_samples=n_called/(n_called+n_missing),
         nasutus_frequency=genosum/(2*n_called),
         het_frequency=n_1/n_called,
         keygroup="Southern_pseudo_allopatric")

###integrate with recombination rates
rr_ready_50kb<-recombinations_windows_50kb %>% ungroup() %>% select(chr,bin_end,M_per_bp)
rr_ready_1Mb<-recombinations_windows %>% ungroup() %>% select(chr,bin_end,M_per_bp)

with_rr<-all_groups_window_freqs %>%
  mutate(windowend_1Mb=ceiling(windowend/1000000)*1000000) %>%
  left_join(rr_ready_50kb,by=c("chrom"="chr","windowend"="bin_end"),relationship="many-to-one") %>%
  rename("M_per_bp_50kb"="M_per_bp") %>% 
  left_join(rr_ready_1Mb,by=c("chrom"="chr","windowend_1Mb"="bin_end"),relationship="many-to-one") %>%
  rename("M_per_bp_1Mb"="M_per_bp") %>%
  left_join(genecounts,by=c("chrom"="chr","windowend"="windowend"),relationship="many-to-one") %>%
  mutate(M_per_bp_50kb=ifelse(is.na(M_per_bp_50kb),0,M_per_bp_50kb),
         M_per_bp_1Mb=ifelse(is.na(M_per_bp_1Mb),0,M_per_bp_1Mb),
         ngenes=ifelse(is.na(ngenes),0,ngenes)) %>%
  left_join(select(fai,chr,chrsize),by=c("chrom"="chr")) %>%
  mutate(windowmiddle=(windowend+windowstart)/2,
         dist_to_chr_end=ifelse(windowmiddle>(chrsize/2),chrsize-windowmiddle,windowmiddle),
         dist_to_chr_middle=ifelse(windowmiddle>(chrsize/2),windowmiddle-(chrsize/2),(chrsize/2)-windowmiddle),
         relative_chr_pos=dist_to_chr_middle/(chrsize/2))

rr_quantile_cutoffs<-quantile(rr_ready_1Mb$M_per_bp,probs=c(1:10)/10)
my_ecdf<-ecdf(rr_ready_1Mb$M_per_bp)
rr_quantile_cutoffs_50kb<-quantile(rr_ready_50kb$M_per_bp,probs=c(1:10)/10)
my_ecdf_50kb<-ecdf(rr_ready_50kb$M_per_bp)

##get correlations with recombination rate and each other
correlation_matrix<-with_rr %>% filter(keygroup %in% focal_nine$keygroup) %>% select(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos,keygroup,nasutus_frequency) %>%
  pivot_wider(id_cols=c(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos),names_from=keygroup,values_from=nasutus_frequency) %>%
  filter(!is.na(`Southern_highelev_hybrid_16`)) %>%  ###removes one window where one group has 0 called samples
  arrange(chrom,windowend)
correlation_matrix_withmissingness<-with_rr %>% filter(keygroup %in% focal_nine$keygroup) %>% select(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos,keygroup,nasutus_frequency,frac_called_samples) %>%
  rename(frequency=nasutus_frequency,missingness=frac_called_samples) %>%
  pivot_wider(id_cols=c(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos),names_from=keygroup,values_from=c(frequency,missingness),names_sep = "_") %>%
  filter(!is.na(`frequency_Southern_highelev_hybrid_16`)) %>%  ###removes one window where one group has 0 called samples
  arrange(chrom,windowend)
correlation_matrix_with_trueallos<-true_allos_freqs %>% ungroup() %>%
  select(chrom,windowend,nasutus_frequency) %>%
  rename("Southern_true_allopatric"="nasutus_frequency") %>% 
  right_join(
    pseudo_allos_freqs %>% ungroup() %>%
      select(chrom,windowend,nasutus_frequency) %>%
      rename("Southern_pseudo_allopatric"="nasutus_frequency")
  ) %>%
  right_join(correlation_matrix)
correlation_matrix_with_trueallos_missingness<-true_allos_freqs %>% ungroup() %>%
  select(chrom,windowend,nasutus_frequency,frac_called_samples) %>%
  rename("frequency_Southern_true_allopatric"="nasutus_frequency","missingness_Southern_true_allopatric"="frac_called_samples") %>%
  right_join(
    pseudo_allos_freqs %>% ungroup() %>%
      select(chrom,windowend,nasutus_frequency,frac_called_samples) %>%
      rename("frequency_Southern_pseudo_allopatric"="nasutus_frequency","missingness_Southern_pseudo_allopatric"="frac_called_samples")
  ) %>%
  right_join(correlation_matrix_withmissingness)

correlation_matrix_withsquares<-correlation_matrix %>% mutate(relative_chr_pos_squared=relative_chr_pos^2,
                                                              ngenes_squared=ngenes,
                                                              M_per_bp_1Mb_squared=M_per_bp_1Mb^2,
                                                              M_per_bp_50kb_squared=M_per_bp_50kb^2) %>%
  select(relative_chr_pos,ngenes,M_per_bp_1Mb,M_per_bp_50kb,
         relative_chr_pos_squared,ngenes_squared,M_per_bp_1Mb_squared,M_per_bp_50kb_squared,
         CAC_guttatus_137,CAC_hybrid_183,CAC_nasutus_58,LM_hybrid_56,
         Southern_allopatric_guttatus_61,Southern_sympatric_guttatus_174,Southern_sympatric_nasutus_75,
         Southern_highelev_guttatus_12,Southern_highelev_hybrid_16)

correlation_matrix_long<-with_rr %>% filter(keygroup %in% focal_nine$keygroup) %>% select(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos,keygroup,nasutus_frequency) %>%
  filter(!is.na(nasutus_frequency)) %>%
  arrange(chrom,windowend)

rescale_z<-function(x) {
  return((x-mean(x))/sqrt(var(x)))
}

correlation_matrix_long_rescaled<-correlation_matrix_long %>% mutate(rr_1Mb_rescaled=rescale_z(M_per_bp_1Mb),
                                            rr_50kb_rescaled=rescale_z(M_per_bp_50kb),
                                            ngenes_rescaled=rescale_z(ngenes))
each_rr_models<-list()
each_rr_newmodels<-list()
each_rr_fullinteract_models<-list()
each_missingness_models<-list()
correlation_matrix_with_residuals<-correlation_matrix_with_trueallos %>% select(-c(chrom,windowend,M_per_bp_50kb:relative_chr_pos))
correlation_matrix_with_residuals_missingness<-correlation_matrix_with_trueallos_missingness %>% select(-c(chrom,windowend,M_per_bp_50kb:relative_chr_pos),-starts_with("missingness_")) %>%
  rename_with(function(x){str_remove(x,"frequency_")},.cols=starts_with("frequency"))

for (i in c(focal_nine$keygroup,"Southern_pseudo_allopatric","Southern_true_allopatric")) {
  baseformula_string="relative_chr_pos + ngenes + M_per_bp_1Mb + M_per_bp_50kb"
  newmodel_string="I(relative_chr_pos^2) + ngenes + M_per_bp_1Mb + M_per_bp_50kb"
  fullinteract_string="relative_chr_pos * ngenes * M_per_bp_1Mb * M_per_bp_50kb *
                       I(relative_chr_pos^2) * I(ngenes^2) * I(M_per_bp_1Mb^2) * I(M_per_bp_50kb^2)"
  withmissingness_string="I(relative_chr_pos^2) + ngenes + M_per_bp_1Mb + M_per_bp_50kb + missingness_"
  basemodel<-lm(as.formula(paste0(i," ~ ",baseformula_string)),data=correlation_matrix_with_trueallos)
  newmodel<-lm(as.formula(paste0(i," ~ ",newmodel_string)),data=correlation_matrix_with_trueallos)
  fullinteractmodel<-lm(as.formula(paste0(i," ~ ",fullinteract_string)),data=correlation_matrix_with_trueallos)
  withmissingness_model<-lm(as.formula(paste0("frequency_",i," ~ ",withmissingness_string,i)),data=correlation_matrix_with_trueallos_missingness)
  correlation_matrix_with_residuals<-correlation_matrix_with_residuals %>% 
    add_column(!!paste0(i,".residuals"):=residuals(basemodel)) %>%
    add_column(!!paste0(i,".residuals_newmodel"):=residuals(newmodel))
  correlation_matrix_with_residuals_missingness<-correlation_matrix_with_residuals_missingness %>%
    add_column(!!paste0(i,".residuals"):=residuals(basemodel)) %>%
    add_column(!!paste0(i,".residuals_newmodel"):=residuals(newmodel)) %>%
    add_column(!!paste0(i,".residuals_missingness"):=residuals(withmissingness_model))
  each_rr_models[[i]]<-basemodel
  each_rr_newmodels[[i]]<-newmodel
  each_missingness_models[[i]]<-withmissingness_model
}

key_groups_customorder<-key_groups_freqorder %>% mutate(keygroup=factor(keygroup,levels=c(
  "LM_hybrid_56","CAC_hybrid_183","CAC_guttatus_137",
  "Southern_sympatric_guttatus_174","Southern_allopatric_guttatus_61",
  "Southern_highelev_hybrid_16","Southern_highelev_guttatus_12",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
))) %>% arrange(keygroup)
key_groups_customorder2<-key_groups_freqorder %>% mutate(keygroup=factor(keygroup,levels=c(
  "CAC_guttatus_137","CAC_hybrid_183","LM_hybrid_56",
  "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174",
  "Southern_highelev_guttatus_12","Southern_highelev_hybrid_16",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
))) %>% arrange(keygroup)
key_groups_customorder_withtrueallos<-c(
  "LM_hybrid_56","CAC_hybrid_183","CAC_guttatus_137",
  "Southern_sympatric_guttatus_174","Southern_allopatric_guttatus_61",
  "Southern_pseudo_allopatric","Southern_true_allopatric",
  "Southern_highelev_hybrid_16","Southern_highelev_guttatus_12",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
)
key_groups_customorder2_withtrueallos<-c(
  "CAC_guttatus_137","CAC_hybrid_183","LM_hybrid_56",
  "Southern_true_allopatric","Southern_pseudo_allopatric",
  "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174",
  "Southern_highelev_guttatus_12","Southern_highelev_hybrid_16",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
)
all_correlations_withcorrections<-as.data.frame(cor(correlation_matrix_with_residuals)) %>% 
  rownames_to_column(var="group1") %>%
  pivot_longer(cols=-group1,names_to="compare_to",values_to="r") %>%
  separate(compare_to,into=c("group2","comparison"),sep="\\.",fill = "right") %>%
  mutate(comparison=ifelse(is.na(comparison),"uncorrected",comparison)) %>%
  pivot_wider(names_from=comparison,values_from=r) %>%
  mutate(group1=factor(group1,levels=key_groups_customorder_withtrueallos),
         group2=factor(group2,levels=key_groups_customorder_withtrueallos),
         #whichval=ifelse(as.integer(group1)>as.integer(group2),residuals,uncorrected),
         #whichval2=ifelse(as.integer(group1)>as.integer(group2),residuals_newmodel,residuals_newmodel),
         uncorrected_prettyval=format(round(uncorrected,3),digits=3),
         residuals_newmodel_prettyval=format(round(residuals_newmodel,3),digits=3)) %>%
  filter(!is.na(group1),!is.na(group2))
all_correlations_withcorrections_missingness<-as.data.frame(cor(correlation_matrix_with_residuals_missingness)) %>% 
  rownames_to_column(var="group1") %>%
  pivot_longer(cols=-group1,names_to="compare_to",values_to="r") %>%
  separate(compare_to,into=c("group2","comparison"),sep="\\.",fill = "right") %>%
  mutate(comparison=ifelse(is.na(comparison),"uncorrected",comparison)) %>%
  pivot_wider(names_from=comparison,values_from=r) %>%
  mutate(group1=factor(group1,levels=key_groups_customorder_withtrueallos),
         group2=factor(group2,levels=key_groups_customorder_withtrueallos),
         #whichval=ifelse(as.integer(group1)>as.integer(group2),residuals_newmodel,uncorrected),
         #whichval2=ifelse(as.integer(group1)>as.integer(group2),residuals_missingness,residuals_missingness),
         residuals_newmodel_prettyval=format(round(residuals_newmodel,3),digits=3),
         residuals_missingness_prettyval=format(round(residuals_missingness,3),digits=3)) %>%
  filter(!is.na(group1),!is.na(group2))

###predict ancestry based on genome-structure data alone
ancestry_predictions<-correlation_matrix_with_trueallos %>% select(chrom,windowend)
ancestry_predictions_newmodel<-correlation_matrix_with_trueallos %>% select(chrom,windowend)
ancestry_predictions_missingness<-correlation_matrix_with_trueallos %>% select(chrom,windowend)
for (i in c(focal_nine$keygroup,"Southern_true_allopatric","Southern_pseudo_allopatric")) {
  ancestry_predictions<-ancestry_predictions %>% add_column(!!i:=predict(each_rr_models[[i]]))
  ancestry_predictions_newmodel<-ancestry_predictions_newmodel %>% add_column(!!i:=predict(each_rr_newmodels[[i]]))
  ancestry_predictions_missingness<-ancestry_predictions_missingness %>% add_column(!!i:=predict(each_missingness_models[[i]]))
}
all_groups_window_freqs_with_trueallos<-all_groups_window_freqs %>% bind_rows(true_allos_freqs) %>% bind_rows(pseudo_allos_freqs)
ancestry_predictions_ready<-ancestry_predictions %>% pivot_longer(-c(chrom,windowend),names_to="keygroup",values_to="prediction") %>%
  left_join(ancestry_predictions_newmodel %>% pivot_longer(-c(chrom,windowend),names_to="keygroup",values_to="prediction_newmodel")) %>%
  left_join(ancestry_predictions_missingness %>% pivot_longer(-c(chrom,windowend),names_to="keygroup",values_to="prediction_missingness")) %>%
  left_join(all_groups_window_freqs_with_trueallos %>% filter(keygroup %in% c(focal_nine$keygroup,"Southern_true_allopatric","Southern_pseudo_allopatric")),by = c("chrom","windowend","keygroup")) %>%
  pivot_longer(c(nasutus_frequency,prediction,prediction_newmodel,prediction_missingness),names_to="statistic",values_to="nasutus_frequency")

prediction_deltas_allmodels<-ancestry_predictions_ready %>%
  pivot_wider(names_from = "statistic",values_from="nasutus_frequency") %>%
  mutate(delta_predict=nasutus_frequency-prediction_newmodel,
         delta_predict_missingness=nasutus_frequency-prediction_missingness)

###tests of window outlier overlap -- add new allopatric groups
focal_outliers_95_uncorrected<-list()
focal_outliers_95_corrected<-list()
focal_outliers_95_missingness<-list()
permuted_outliers_95<-list()
for (i in c(key_groups_customorder2_withtrueallos)) {
  subdata<-prediction_deltas_allmodels %>% filter(keygroup==i,!is.na(delta_predict)) %>% select(chrom,windowend,chrom.windowstart.windowend,delta_predict,delta_predict_missingness,nasutus_frequency)
  meanvalue<-mean(subdata$nasutus_frequency)
  permuted<-data.frame()
  if (meanvalue>0.5) {
    cutoff_uncorrected=quantile(subdata$nasutus_frequency,0.05)
    cutoff_corrected=quantile(subdata$delta_predict,0.05)
    cutoff_corrected_missingness=quantile(subdata$delta_predict_missingness,0.05)
    focal_outliers_95_uncorrected[[i]]<-subdata %>% filter(nasutus_frequency<=cutoff_uncorrected)
    focal_outliers_95_corrected[[i]]<-subdata %>% filter(delta_predict<=cutoff_corrected)
    focal_outliers_95_missingness[[i]]<-subdata %>% filter(delta_predict_missingness<=cutoff_corrected_missingness)
    
  } else {
    cutoff_uncorrected=quantile(subdata$nasutus_frequency,0.95)
    cutoff_corrected=quantile(subdata$delta_predict,0.95)
    cutoff_corrected_missingness=quantile(subdata$delta_predict_missingness,0.95)
    focal_outliers_95_uncorrected[[i]]<-subdata %>% filter(nasutus_frequency>=cutoff_uncorrected)
    focal_outliers_95_corrected[[i]]<-subdata %>% filter(delta_predict>=cutoff_corrected)
    focal_outliers_95_missingness[[i]]<-subdata %>% filter(delta_predict_missingness>=cutoff_corrected_missingness)
  }
  for (p in 1:1000) {
    one_permute<-sample_n(subdata,nrow(focal_outliers_95_corrected[[i]])) %>% mutate(permutation=p)
    permuted<-permuted %>% bind_rows(one_permute)
  }
  permuted_outliers_95[[i]]<-permuted
}
count_overlaps<-data.frame()
overlap_permuter<-function(p,i,j,corrected=T,missingness=T) {
  if (missingness) {
    return(inner_join(focal_outliers_95_missingness[[i]],permuted_outliers_95[[j]] %>% filter(permutation==p),by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow())
  } else if (corrected) {
    return(inner_join(focal_outliers_95_corrected[[i]],permuted_outliers_95[[j]] %>% filter(permutation==p),by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow())
  } else {
    return(inner_join(focal_outliers_95_uncorrected[[i]],permuted_outliers_95[[j]] %>% filter(permutation==p),by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow())
  }
}
for (i in key_groups_customorder2_withtrueallos) {
  for (j in key_groups_customorder2_withtrueallos) {
    N_i<-prediction_deltas_allmodels %>% filter(keygroup==i,!is.na(delta_predict)) %>% nrow()
    N_j<-prediction_deltas_allmodels %>% filter(keygroup==j,!is.na(delta_predict)) %>% nrow()
    
    n_outliers_uncorrected_i<-focal_outliers_95_uncorrected[[i]] %>% nrow()
    n_outliers_corrected_i<-focal_outliers_95_corrected[[i]] %>% nrow()
    n_outliers_missingness_i<-focal_outliers_95_missingness[[i]] %>% nrow()
    n_outliers_uncorrected_j<-focal_outliers_95_uncorrected[[j]] %>% nrow()
    n_outliers_corrected_j<-focal_outliers_95_corrected[[j]] %>% nrow()
    n_outliers_missingness_j<-focal_outliers_95_missingness[[j]] %>% nrow()
    
    observed_overlap_uncorrected<-inner_join(focal_outliers_95_uncorrected[[i]],focal_outliers_95_uncorrected[[j]],by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow()
    observed_overlap_corrected<-inner_join(focal_outliers_95_corrected[[i]],focal_outliers_95_corrected[[j]],by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow()
    observed_overlap_missingness<-inner_join(focal_outliers_95_missingness[[i]],focal_outliers_95_missingness[[j]],by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow()
    
    expected_overlap_uncorrected<-n_outliers_uncorrected_i*n_outliers_uncorrected_j/(N_j)
    expected_overlap_corrected<-n_outliers_corrected_i*n_outliers_corrected_j/(N_j)
    expected_overlap_missingness<-n_outliers_missingness_i*n_outliers_missingness_j/(N_j)
    #compare observed i to permuted j
    permute_overlaps<-sapply(c(1:1000),overlap_permuter,i=i,j=j)
    mean_permute_overlaps<-mean(permute_overlaps)
    stderr_permute_overlaps<-sqrt(var(permute_overlaps)/1000)
    permuted_ecdf<-ecdf(permute_overlaps)
    pval_permutes_uncorrected<-1-permuted_ecdf(observed_overlap_uncorrected)
    pval_permutes_corrected<-1-permuted_ecdf(observed_overlap_corrected)
    pval_permutes_missingness<-1-permuted_ecdf(observed_overlap_missingness)
    n_permutes_above_observed_uncorrected<-sum(permute_overlaps>=observed_overlap_uncorrected)
    n_permutes_above_observed_corrected<-sum(permute_overlaps>=observed_overlap_corrected)
    n_permutes_above_observed_missingness<-sum(permute_overlaps>=observed_overlap_missingness)
    count_overlaps_line<-data.frame(i,j,N_i,N_j,n_outliers_uncorrected_i,n_outliers_corrected_i,n_outliers_missingness_i,n_outliers_uncorrected_j,n_outliers_corrected_j,n_outliers_missingness_j,
                                    observed_overlap_uncorrected,observed_overlap_corrected,observed_overlap_missingness,expected_overlap_uncorrected,expected_overlap_corrected,expected_overlap_missingness,
                                    mean_permute_overlaps,stderr_permute_overlaps,
                                    pval_permutes_uncorrected,pval_permutes_corrected,pval_permutes_missingness,
                                    n_permutes_above_observed_uncorrected,n_permutes_above_observed_corrected,n_permutes_above_observed_missingness)
    count_overlaps<-count_overlaps %>% bind_rows(count_overlaps_line)
  }
}

###Mantel regions corrected for IM767ref
Mantel_regions<-read_excel(file.path(thisfolder,"Mantel2023_regions.xlsx"),sheet="Mantel_convert_IM767")
Mantel_regions_withGOIs<-read_excel(file.path(thisfolder,"Mantel2023_regions.xlsx"),sheet="Mantel_convert_IM767_withGOIs") %>% 
  filter(is.na(Exclude_size))
Mantel_regions_GOIs<-read_excel(file.path(thisfolder,"Mantel2023_regions.xlsx"),sheet="Mantel_convert_IM767_withGOIs") %>% 
  filter(Group %in% c("Hybrid_lethality","Photoperiod")) %>%
  mutate(chrom=paste0("Chr",str_pad(LG,2,pad="0")))

Mantel_regions_nolarge<-Mantel_regions %>% filter(is.na(Exclude_size))
Mantel_regions_windows<-data.frame()
for (i in 1:nrow(Mantel_regions)) {
  line=Mantel_regions[i,]
  chr=paste0("Chr",str_pad(line$LG,2,pad="0"))
  first=line$first_windowend_IM767
  last=line$last_windowend_IM767
  for (w in c((first/50000):(last/50000))*50000){
    windowline=line %>% select(Name,Group,Predicted_effect) %>% 
      mutate(chrom=chr,region_size=last-first+50000,region_start=first,region_end=last,windowend=w)
    Mantel_regions_windows<-bind_rows(Mantel_regions_windows,windowline)
  }
}
Mantel_regions_windows_nolarge<-Mantel_regions_windows %>% filter(Name %in% Mantel_regions_nolarge$Name)

Mantel_regions_windows_withGOIs<-data.frame()
for (i in 1:nrow(Mantel_regions_withGOIs)) {
  line=Mantel_regions_withGOIs[i,]
  chr=paste0("Chr",str_pad(line$LG,2,pad="0"))
  first=line$first_windowend_IM767
  last=line$last_windowend_IM767
  for (w in c((first/50000):(last/50000))*50000){
    windowline=line %>% select(Name,Group,Predicted_effect) %>% 
      mutate(chrom=chr,region_size=last-first+50000,region_start=first,region_end=last,windowend=w)
    Mantel_regions_windows_withGOIs<-bind_rows(Mantel_regions_windows_withGOIs,windowline)
  }
}

Mantel_regions_withGOIs_freqs_withpredictions_newmodel<-Mantel_regions_windows_withGOIs %>% filter(Name!="pTAC14") %>% left_join(prediction_deltas_allmodels,by=c("chrom","windowend"),relationship="one-to-many") %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos))  ###remove pTAC14 because the window is already counted once
Mantel_regions_nolarge_freqs_withpredictions_newmodel<-Mantel_regions_windows_nolarge %>% left_join(prediction_deltas_allmodels,by=c("chrom","windowend"),relationship="one-to-many") %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos))


###do outliers occur in windows with QTL more often than expected?

outliers_in_windows_summary<-data.frame()
for (i in key_groups_customorder2_withtrueallos) {
  n_total_windows<-ancestry_predictions_ready %>% filter(keygroup==i,statistic=="nasutus_frequency") %>% nrow()
  n_QTLwindows<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>% filter(keygroup==i) %>% nrow()
  n_outlier_windows_uncorrected<-focal_outliers_95_uncorrected[[i]] %>% nrow()
  n_outlier_windows_corrected<-focal_outliers_95_corrected[[i]] %>% nrow()
  n_outlier_windows_missingness<-focal_outliers_95_missingness[[i]] %>% nrow()
  n_outliers_in_QTL_uncorrected<-focal_outliers_95_uncorrected[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_corrected<-focal_outliers_95_corrected[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_missingness<-focal_outliers_95_missingness[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  dataline<-data.frame(keygroup=i,n_total_windows=n_total_windows,n_QTLwindows=n_QTLwindows,
                       n_outlier_windows_uncorrected=n_outlier_windows_uncorrected,n_outlier_windows_corrected=n_outlier_windows_corrected,n_outlier_windows_missingness=n_outlier_windows_missingness,
                       n_outliers_in_QTL_uncorrected=n_outliers_in_QTL_uncorrected,n_outliers_in_QTL_corrected=n_outliers_in_QTL_corrected,n_outliers_in_QTL_missingness=n_outliers_in_QTL_missingness)
  outliers_in_windows_summary<-bind_rows(outliers_in_windows_summary,dataline)
}
do_chisq<-function(observed,total,expectedproportion,returnval="statistic") {
  test=chisq.test(c(observed,total-observed),p=c(expectedproportion,1-expectedproportion))
  if (returnval=="statistic") {
    return(test$statistic)
  } else if (returnval=="p.value") {
    return(test$p.value)
  }
}
testfunc<-function(row,column) {
  return(column[row])
}
outliers_in_windows_summary<-outliers_in_windows_summary %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos),
         expected_uncorrected=n_outlier_windows_uncorrected*n_QTLwindows/n_total_windows,
         expected_corrected=n_outlier_windows_corrected*n_QTLwindows/n_total_windows,
         expected_missingness=n_outlier_windows_missingness*n_QTLwindows/n_total_windows,
         expected_fraction_uncorrected=expected_uncorrected/n_outlier_windows_uncorrected,
         expected_fraction_corrected=expected_corrected/n_outlier_windows_corrected,
         expected_fraction_missingness=expected_missingness/n_outlier_windows_corrected,
         observed_fraction_uncorrected=n_outliers_in_QTL_uncorrected/n_outlier_windows_uncorrected,
         observed_fraction_corrected=n_outliers_in_QTL_corrected/n_outlier_windows_corrected,
         observed_fraction_missingness=n_outliers_in_QTL_missingness/n_outlier_windows_corrected) %>%
  rowwise() %>% mutate(
    chisq_uncorrected=do_chisq(observed=n_outliers_in_QTL_uncorrected,total=n_outlier_windows_uncorrected,expectedproportion=expected_fraction_uncorrected),
    chisq_corrected=do_chisq(observed=n_outliers_in_QTL_corrected,total=n_outlier_windows_corrected,expectedproportion=expected_fraction_corrected),
    chisq_missingness=do_chisq(observed=n_outliers_in_QTL_missingness,total=n_outlier_windows_missingness,expectedproportion=expected_fraction_missingness),
    chisq_pval_uncorrected=do_chisq(observed=n_outliers_in_QTL_uncorrected,total=n_outlier_windows_uncorrected,expectedproportion=expected_fraction_uncorrected,returnval="p.value"),
    chisq_pval_corrected=do_chisq(observed=n_outliers_in_QTL_corrected,total=n_outlier_windows_corrected,expectedproportion=expected_fraction_corrected,returnval="p.value"),
    chisq_pval_missingness=do_chisq(observed=n_outliers_in_QTL_missingness,total=n_outlier_windows_missingness,expectedproportion=expected_fraction_missingness,returnval="p.value"),
    chisq_pval_bonf_uncorrected=chisq_pval_uncorrected*9,
    chisq_pval_bonf_corrected=chisq_pval_corrected*9,
    chisq_pval_bonf_missingness=chisq_pval_missingness*9,
    sig_uncorrected=ifelse(chisq_pval_bonf_uncorrected<0.001,"***",ifelse(chisq_pval_bonf_uncorrected<0.01,"**",ifelse(chisq_pval_bonf_uncorrected<0.05,"*",ifelse(chisq_pval_bonf_uncorrected<0.1,".","")))),
    sig_corrected=ifelse(chisq_pval_bonf_corrected<0.001,"***",ifelse(chisq_pval_bonf_corrected<0.01,"**",ifelse(chisq_pval_bonf_corrected<0.05,"*",ifelse(chisq_pval_bonf_corrected<0.1,".","")))),
    sig_missingness=ifelse(chisq_pval_bonf_missingness<0.001,"***",ifelse(chisq_pval_bonf_missingness<0.01,"**",ifelse(chisq_pval_bonf_missingness<0.05,"*",ifelse(chisq_pval_bonf_missingness<0.1,".",""))))
  )
outliers_in_windows_summary_longer<-bind_rows(
  outliers_in_windows_summary %>% select(keygroup,n_total_windows,n_QTLwindows,ends_with("_uncorrected")) %>%
    rename_with(function(x){str_remove(x,"_uncorrected")}) %>% 
    mutate("dataset"="uncorrected"),
  outliers_in_windows_summary %>% select(keygroup,n_total_windows,n_QTLwindows,ends_with("_corrected")) %>%
    rename_with(function(x){str_remove(x,"_corrected")}) %>% 
    mutate("dataset"="corrected"),
  outliers_in_windows_summary %>% select(keygroup,n_total_windows,n_QTLwindows,ends_with("_missingness")) %>%
    rename_with(function(x){str_remove(x,"_missingness")}) %>% 
    mutate("dataset"="missingness")) %>% arrange(desc(dataset),keygroup)

###outliers in windows: use only N+ vs. N- QTL
Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>% filter(Predicted_effect=="N+")
Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>% filter(Predicted_effect=="N-")
outliers_in_windows_summary_split<-data.frame()
for (i in key_groups_customorder2_withtrueallos) {
  n_total_windows<-ancestry_predictions_ready %>% filter(keygroup==i,statistic=="nasutus_frequency") %>% nrow()
  n_QTLwindows_plus<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus %>% filter(keygroup==i) %>% nrow()
  n_QTLwindows_minus<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus %>% filter(keygroup==i) %>% nrow()
  n_outlier_windows_uncorrected<-focal_outliers_95_uncorrected[[i]] %>% nrow()
  n_outlier_windows_corrected<-focal_outliers_95_corrected[[i]] %>% nrow()
  n_outlier_windows_missingness<-focal_outliers_95_missingness[[i]] %>% nrow()
  n_outliers_in_QTL_uncorrected_plus<-focal_outliers_95_uncorrected[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_corrected_plus<-focal_outliers_95_corrected[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_missingness_plus<-focal_outliers_95_missingness[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_uncorrected_minus<-focal_outliers_95_uncorrected[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_corrected_minus<-focal_outliers_95_corrected[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_missingness_minus<-focal_outliers_95_missingness[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  dataline<-data.frame(keygroup=i,n_total_windows=n_total_windows,n_QTLwindows_plus=n_QTLwindows_plus,n_QTLwindows_minus=n_QTLwindows_minus,
                       n_outlier_windows_uncorrected=n_outlier_windows_uncorrected,n_outlier_windows_corrected=n_outlier_windows_corrected,n_outlier_windows_missingness=n_outlier_windows_missingness,
                       n_outliers_in_QTL_uncorrected_minus=n_outliers_in_QTL_uncorrected_minus,n_outliers_in_QTL_corrected_minus=n_outliers_in_QTL_corrected_minus,n_outliers_in_QTL_missingness_minus=n_outliers_in_QTL_missingness_minus,
                       n_outliers_in_QTL_uncorrected_plus=n_outliers_in_QTL_uncorrected_plus,n_outliers_in_QTL_corrected_plus=n_outliers_in_QTL_corrected_plus,n_outliers_in_QTL_missingness_plus=n_outliers_in_QTL_missingness_plus)
  outliers_in_windows_summary_split<-bind_rows(outliers_in_windows_summary_split,dataline)
}

outliers_in_windows_summary_split<-outliers_in_windows_summary_split %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos),
         expected_uncorrected_plus=n_outlier_windows_uncorrected*n_QTLwindows_plus/n_total_windows,
         expected_corrected_plus=n_outlier_windows_corrected*n_QTLwindows_plus/n_total_windows,
         expected_missingness_plus=n_outlier_windows_missingness*n_QTLwindows_plus/n_total_windows,
         expected_fraction_uncorrected_plus=expected_uncorrected_plus/n_outlier_windows_uncorrected,
         expected_fraction_corrected_plus=expected_corrected_plus/n_outlier_windows_corrected,
         expected_fraction_missingness_plus=expected_missingness_plus/n_outlier_windows_corrected,
         observed_fraction_uncorrected_plus=n_outliers_in_QTL_uncorrected_plus/n_outlier_windows_uncorrected,
         observed_fraction_corrected_plus=n_outliers_in_QTL_corrected_plus/n_outlier_windows_corrected,
         observed_fraction_missingness_plus=n_outliers_in_QTL_missingness_plus/n_outlier_windows_corrected,
         
         expected_uncorrected_minus=n_outlier_windows_uncorrected*n_QTLwindows_minus/n_total_windows,
         expected_corrected_minus=n_outlier_windows_corrected*n_QTLwindows_minus/n_total_windows,
         expected_missingness_minus=n_outlier_windows_missingness*n_QTLwindows_minus/n_total_windows,
         expected_fraction_uncorrected_minus=expected_uncorrected_minus/n_outlier_windows_uncorrected,
         expected_fraction_corrected_minus=expected_corrected_minus/n_outlier_windows_corrected,
         expected_fraction_missingness_minus=expected_missingness_minus/n_outlier_windows_corrected,
         observed_fraction_uncorrected_minus=n_outliers_in_QTL_uncorrected_minus/n_outlier_windows_uncorrected,
         observed_fraction_corrected_minus=n_outliers_in_QTL_corrected_minus/n_outlier_windows_corrected,
         observed_fraction_missingness_minus=n_outliers_in_QTL_missingness_minus/n_outlier_windows_corrected) %>%
  rowwise() %>% mutate(
    chisq_uncorrected_plus=do_chisq(observed=n_outliers_in_QTL_uncorrected_plus,total=n_outlier_windows_uncorrected,expectedproportion=expected_fraction_uncorrected_plus),
    chisq_corrected_plus=do_chisq(observed=n_outliers_in_QTL_corrected_plus,total=n_outlier_windows_corrected,expectedproportion=expected_fraction_corrected_plus),
    chisq_missingness_plus=do_chisq(observed=n_outliers_in_QTL_missingness_plus,total=n_outlier_windows_missingness,expectedproportion=expected_fraction_missingness_plus),
    chisq_uncorrected_minus=do_chisq(observed=n_outliers_in_QTL_uncorrected_minus,total=n_outlier_windows_uncorrected,expectedproportion=expected_fraction_uncorrected_minus),
    chisq_corrected_minus=do_chisq(observed=n_outliers_in_QTL_corrected_minus,total=n_outlier_windows_corrected,expectedproportion=expected_fraction_corrected_minus),
    chisq_missingness_minus=do_chisq(observed=n_outliers_in_QTL_missingness_minus,total=n_outlier_windows_missingness,expectedproportion=expected_fraction_missingness_minus),
    chisq_pval_uncorrected_plus=do_chisq(observed=n_outliers_in_QTL_uncorrected_plus,total=n_outlier_windows_uncorrected,expectedproportion=expected_fraction_uncorrected_plus,returnval="p.value"),
    chisq_pval_corrected_plus=do_chisq(observed=n_outliers_in_QTL_corrected_plus,total=n_outlier_windows_corrected,expectedproportion=expected_fraction_corrected_plus,returnval="p.value"),
    chisq_pval_missingness_plus=do_chisq(observed=n_outliers_in_QTL_missingness_plus,total=n_outlier_windows_missingness,expectedproportion=expected_fraction_missingness_plus,returnval="p.value"),
    chisq_pval_uncorrected_minus=do_chisq(observed=n_outliers_in_QTL_uncorrected_minus,total=n_outlier_windows_uncorrected,expectedproportion=expected_fraction_uncorrected_minus,returnval="p.value"),
    chisq_pval_corrected_minus=do_chisq(observed=n_outliers_in_QTL_corrected_minus,total=n_outlier_windows_corrected,expectedproportion=expected_fraction_corrected_minus,returnval="p.value"),
    chisq_pval_missingness_minus=do_chisq(observed=n_outliers_in_QTL_missingness_minus,total=n_outlier_windows_missingness,expectedproportion=expected_fraction_missingness_minus,returnval="p.value"),
    chisq_pval_bonf_uncorrected_plus=chisq_pval_uncorrected_plus*9,
    chisq_pval_bonf_corrected_plus=chisq_pval_corrected_plus*9,
    chisq_pval_bonf_missingness_plus=chisq_pval_missingness_plus*9,
    chisq_pval_bonf_uncorrected_minus=chisq_pval_uncorrected_minus*9,
    chisq_pval_bonf_corrected_minus=chisq_pval_corrected_minus*9,
    chisq_pval_bonf_missingness_minus=chisq_pval_missingness_minus*9,
    sig_uncorrected_plus=ifelse(chisq_pval_bonf_uncorrected_plus<0.001,"***",ifelse(chisq_pval_bonf_uncorrected_plus<0.01,"**",ifelse(chisq_pval_bonf_uncorrected_plus<0.05,"*",ifelse(chisq_pval_bonf_uncorrected_plus<0.1,".","")))),
    sig_corrected_plus=ifelse(chisq_pval_bonf_corrected_plus<0.001,"***",ifelse(chisq_pval_bonf_corrected_plus<0.01,"**",ifelse(chisq_pval_bonf_corrected_plus<0.05,"*",ifelse(chisq_pval_bonf_corrected_plus<0.1,".","")))),
    sig_missingness_plus=ifelse(chisq_pval_bonf_missingness_plus<0.001,"***",ifelse(chisq_pval_bonf_missingness_plus<0.01,"**",ifelse(chisq_pval_bonf_missingness_plus<0.05,"*",ifelse(chisq_pval_bonf_missingness_plus<0.1,".","")))),
    sig_uncorrected_minus=ifelse(chisq_pval_bonf_uncorrected_minus<0.001,"***",ifelse(chisq_pval_bonf_uncorrected_minus<0.01,"**",ifelse(chisq_pval_bonf_uncorrected_minus<0.05,"*",ifelse(chisq_pval_bonf_uncorrected_minus<0.1,".","")))),
    sig_corrected_minus=ifelse(chisq_pval_bonf_corrected_minus<0.001,"***",ifelse(chisq_pval_bonf_corrected_minus<0.01,"**",ifelse(chisq_pval_bonf_corrected_minus<0.05,"*",ifelse(chisq_pval_bonf_corrected_minus<0.1,".","")))),
    sig_missingness_minus=ifelse(chisq_pval_bonf_missingness_minus<0.001,"***",ifelse(chisq_pval_bonf_missingness_minus<0.01,"**",ifelse(chisq_pval_bonf_missingness_minus<0.05,"*",ifelse(chisq_pval_bonf_missingness_minus<0.1,".",""))))
  )
outliers_in_windows_summary_split_longer<-bind_rows(
  outliers_in_windows_summary_split %>% select(keygroup,n_total_windows,n_QTLwindows_plus,ends_with("_uncorrected_plus")) %>%
    rename_with(function(x){str_remove(x,"_uncorrected_plus")}) %>% rename("n_QTLwindows"="n_QTLwindows_plus") %>%
    mutate("dataset"="uncorrected","QTLset"="plus"),
  outliers_in_windows_summary_split %>% select(keygroup,n_total_windows,n_QTLwindows_minus,ends_with("_uncorrected_minus")) %>%
    rename_with(function(x){str_remove(x,"_uncorrected_minus")}) %>% rename("n_QTLwindows"="n_QTLwindows_minus") %>%
    mutate("dataset"="uncorrected","QTLset"="minus"),
  outliers_in_windows_summary_split %>% select(keygroup,n_total_windows,n_QTLwindows_plus,ends_with("_corrected_plus")) %>%
    rename_with(function(x){str_remove(x,"_corrected_plus")}) %>% rename("n_QTLwindows"="n_QTLwindows_plus") %>%
    mutate("dataset"="corrected","QTLset"="plus"),
  outliers_in_windows_summary_split %>% select(keygroup,n_total_windows,n_QTLwindows_minus,ends_with("_corrected_minus")) %>%
    rename_with(function(x){str_remove(x,"_corrected_minus")}) %>% rename("n_QTLwindows"="n_QTLwindows_minus") %>%
    mutate("dataset"="corrected","QTLset"="minus"),
  outliers_in_windows_summary_split %>% select(keygroup,n_total_windows,n_QTLwindows_plus,ends_with("_missingness_plus")) %>%
    rename_with(function(x){str_remove(x,"_missingness_plus")}) %>% rename("n_QTLwindows"="n_QTLwindows_plus") %>%
    mutate("dataset"="missingness","QTLset"="plus"),
  outliers_in_windows_summary_split %>% select(keygroup,n_total_windows,n_QTLwindows_minus,ends_with("_missingness_minus")) %>%
    rename_with(function(x){str_remove(x,"_missingness_minus")}) %>% rename("n_QTLwindows"="n_QTLwindows_minus") %>%
    mutate("dataset"="missingness","QTLset"="minus")) %>% 
  arrange(QTLset,desc(dataset),keygroup)


###t-tests on individual windows
ttests_Nplus<-list()
ttests_Nplus_withGOIs<-list()
ttests_Nplus_missingness<-list()
ttests_Nplus_withGOIs_missingness<-list()
ttests_Nminus<-list()
ttests_Nminus_withGOIs<-list()
ttests_Nminus_missingness<-list()
ttests_Nminus_withGOIs_missingness<-list()
ttest_results<-data.frame()
ttest_results_withGOIs<-data.frame()
ttest_results_missingness<-data.frame()
ttest_results_withGOIs_missingness<-data.frame()
for (i in c(focal_nine$keygroup,"Southern_true_allopatric","Southern_pseudo_allopatric")) {
  ttests_Nplus[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_withGOIs[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_withGOIs[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_missingness[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_missingness[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_withGOIs_missingness[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_withGOIs_missingness[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttest_results<-bind_rows(ttest_results,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus[[i]]$statistic,
    df=ttests_Nplus[[i]]$parameter,
    p.value=ttests_Nplus[[i]]$p.value,
    mean=ttests_Nplus[[i]]$estimate,
    stderr=ttests_Nplus[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus[[i]]$statistic,
    df=ttests_Nminus[[i]]$parameter,
    p.value=ttests_Nminus[[i]]$p.value,
    mean=ttests_Nminus[[i]]$estimate,
    stderr=ttests_Nminus[[i]]$stderr
  )) 
  ttest_results_withGOIs<-bind_rows(ttest_results_withGOIs,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_withGOIs[[i]]$statistic,
    df=ttests_Nplus_withGOIs[[i]]$parameter,
    p.value=ttests_Nplus_withGOIs[[i]]$p.value,
    mean=ttests_Nplus_withGOIs[[i]]$estimate,
    stderr=ttests_Nplus_withGOIs[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_withGOIs[[i]]$statistic,
    df=ttests_Nminus_withGOIs[[i]]$parameter,
    p.value=ttests_Nminus_withGOIs[[i]]$p.value,
    mean=ttests_Nminus_withGOIs[[i]]$estimate,
    stderr=ttests_Nminus_withGOIs[[i]]$stderr
  )) 
  ttest_results_missingness<-bind_rows(ttest_results_missingness,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_missingness[[i]]$statistic,
    df=ttests_Nplus_missingness[[i]]$parameter,
    p.value=ttests_Nplus_missingness[[i]]$p.value,
    mean=ttests_Nplus_missingness[[i]]$estimate,
    stderr=ttests_Nplus_missingness[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_missingness[[i]]$statistic,
    df=ttests_Nminus_missingness[[i]]$parameter,
    p.value=ttests_Nminus_missingness[[i]]$p.value,
    mean=ttests_Nminus_missingness[[i]]$estimate,
    stderr=ttests_Nminus_missingness[[i]]$stderr
  )) 
  ttest_results_withGOIs_missingness<-bind_rows(ttest_results_withGOIs_missingness,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_withGOIs_missingness[[i]]$statistic,
    df=ttests_Nplus_withGOIs_missingness[[i]]$parameter,
    p.value=ttests_Nplus_withGOIs_missingness[[i]]$p.value,
    mean=ttests_Nplus_withGOIs_missingness[[i]]$estimate,
    stderr=ttests_Nplus_withGOIs_missingness[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_withGOIs_missingness[[i]]$statistic,
    df=ttests_Nminus_withGOIs_missingness[[i]]$parameter,
    p.value=ttests_Nminus_withGOIs_missingness[[i]]$p.value,
    mean=ttests_Nminus_withGOIs_missingness[[i]]$estimate,
    stderr=ttests_Nminus_withGOIs_missingness[[i]]$stderr
  )) 
}
label_significance<-function(p) {
  return(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","n.s.")))))
}

ttest_results<-ttest_results %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_withGOIs<-ttest_results_withGOIs %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_missingness<-ttest_results_missingness %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_withGOIs_missingness<-ttest_results_withGOIs_missingness %>% mutate(p.bonf.corrected=p.value*18,
                                                          significance_uncorrected=label_significance(p.value),
                                                          significance_corrected=label_significance(p.bonf.corrected))



###t-tests on region means
ttests_Nplus_regions<-list()
ttests_Nplus_regions_withGOIs<-list()
ttests_Nplus_regions_missingness<-list()
ttests_Nplus_regions_withGOIs_missingness<-list()
ttests_Nminus_regions<-list()
ttests_Nminus_regions_withGOIs<-list()
ttests_Nminus_regions_missingness<-list()
ttests_Nminus_regions_withGOIs_missingness<-list()
ttest_results_regions<-data.frame()
ttest_results_regions_withGOIs<-data.frame()
ttest_results_regions_missingness<-data.frame()
ttest_results_regions_withGOIs_missingness<-data.frame()

for (i in c(focal_nine$keygroup,"Southern_true_allopatric","Southern_pseudo_allopatric")) {
  ttests_Nplus_regions[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_regions[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_regions_withGOIs[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_regions_withGOIs[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_regions_missingness[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% 
    group_by(Name) %>% summarize(delta_predict_missingness=mean(delta_predict_missingness)) %>% 
    pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_regions_missingness[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% 
    group_by(Name) %>% summarize(delta_predict_missingness=mean(delta_predict_missingness)) %>% 
    pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_regions_withGOIs_missingness[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% 
    group_by(Name) %>% summarize(delta_predict_missingness=mean(delta_predict_missingness)) %>% 
    pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_regions_withGOIs_missingness[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% 
    group_by(Name) %>% summarize(delta_predict_missingness=mean(delta_predict_missingness)) %>% 
    pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttest_results_regions<-bind_rows(ttest_results_regions,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_regions[[i]]$statistic,
    df=ttests_Nplus_regions[[i]]$parameter,
    p.value=ttests_Nplus_regions[[i]]$p.value,
    mean=ttests_Nplus_regions[[i]]$estimate,
    stderr=ttests_Nplus_regions[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_regions[[i]]$statistic,
    df=ttests_Nminus_regions[[i]]$parameter,
    p.value=ttests_Nminus_regions[[i]]$p.value,
    mean=ttests_Nminus_regions[[i]]$estimate,
    stderr=ttests_Nminus_regions[[i]]$stderr
  )) 
  ttest_results_regions_withGOIs<-bind_rows(ttest_results_regions_withGOIs,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_regions_withGOIs[[i]]$statistic,
    df=ttests_Nplus_regions_withGOIs[[i]]$parameter,
    p.value=ttests_Nplus_regions_withGOIs[[i]]$p.value,
    mean=ttests_Nplus_regions_withGOIs[[i]]$estimate,
    stderr=ttests_Nplus_regions_withGOIs[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_regions_withGOIs[[i]]$statistic,
    df=ttests_Nminus_regions_withGOIs[[i]]$parameter,
    p.value=ttests_Nminus_regions_withGOIs[[i]]$p.value,
    mean=ttests_Nminus_regions_withGOIs[[i]]$estimate,
    stderr=ttests_Nminus_regions_withGOIs[[i]]$stderr
  )) 
  ttest_results_regions_missingness<-bind_rows(ttest_results_regions_missingness,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_regions_missingness[[i]]$statistic,
    df=ttests_Nplus_regions_missingness[[i]]$parameter,
    p.value=ttests_Nplus_regions_missingness[[i]]$p.value,
    mean=ttests_Nplus_regions_missingness[[i]]$estimate,
    stderr=ttests_Nplus_regions_missingness[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_regions_missingness[[i]]$statistic,
    df=ttests_Nminus_regions_missingness[[i]]$parameter,
    p.value=ttests_Nminus_regions_missingness[[i]]$p.value,
    mean=ttests_Nminus_regions_missingness[[i]]$estimate,
    stderr=ttests_Nminus_regions_missingness[[i]]$stderr
  )) 
  ttest_results_regions_withGOIs_missingness<-bind_rows(ttest_results_regions_withGOIs_missingness,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_regions_withGOIs_missingness[[i]]$statistic,
    df=ttests_Nplus_regions_withGOIs_missingness[[i]]$parameter,
    p.value=ttests_Nplus_regions_withGOIs_missingness[[i]]$p.value,
    mean=ttests_Nplus_regions_withGOIs_missingness[[i]]$estimate,
    stderr=ttests_Nplus_regions_withGOIs_missingness[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_regions_withGOIs_missingness[[i]]$statistic,
    df=ttests_Nminus_regions_withGOIs_missingness[[i]]$parameter,
    p.value=ttests_Nminus_regions_withGOIs_missingness[[i]]$p.value,
    mean=ttests_Nminus_regions_withGOIs_missingness[[i]]$estimate,
    stderr=ttests_Nminus_regions_withGOIs_missingness[[i]]$stderr
  )) 
}
ttest_results_regions<-ttest_results_regions %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_regions_withGOIs<-ttest_results_regions_withGOIs %>% mutate(p.bonf.corrected=p.value*18,
                                                        significance_uncorrected=label_significance(p.value),
                                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_regions_missingness<-ttest_results_regions_missingness %>% mutate(p.bonf.corrected=p.value*18,
                                                        significance_uncorrected=label_significance(p.value),
                                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_regions_withGOIs_missingness<-ttest_results_regions_withGOIs_missingness %>% mutate(p.bonf.corrected=p.value*18,
                                                                          significance_uncorrected=label_significance(p.value),
                                                                          significance_corrected=label_significance(p.bonf.corrected))

outlier_stretches=data.frame()
for (i in key_groups_customorder2_withtrueallos) {
  x=focal_outliers_95_corrected[[i]] %>% 
    mutate(isbreak=(is.na(lag(chrom)) | chrom!=lag(chrom) | windowend!=lag(windowend)+50000)) %>% 
    group_by(chrom) %>% summarize(nstretches=sum(isbreak),nwindows=n()) %>% 
    summarize(nchroms=n(),nstretches=sum(nstretches,na.rm=T),nwindows=sum(nwindows,na.rm=T)) %>% 
    mutate(keygroup=i)
  outlier_stretches=bind_rows(outlier_stretches,x)
}

###create groupings based on type of comparisons
comparison_types<-count_overlaps %>% select(i,j) %>% unique() %>% 
  mutate(
    region_i=ifelse(str_starts(i,"Southern"),"Southern","Northern"),
    region_j=ifelse(str_starts(j,"Southern"),"Southern","Northern"),
    cohort_i=ifelse(str_detect(i,"guttatus"),"guttatus",
                    ifelse(str_detect(i,"allopatric"),"guttatus",
                           ifelse(str_detect(i,"nasutus"),"nasutus","hybrid"))),
    cohort_j=ifelse(str_detect(j,"guttatus"),"guttatus",
                    ifelse(str_detect(j,"allopatric"),"guttatus",
                           ifelse(str_detect(j,"nasutus"),"nasutus","hybrid"))),
    group_i=ifelse(str_detect(i,"allopatric"),"Southern_allopatric",
                   ifelse(str_detect(i,"sympatric"),"Southern_sympatric",
                          ifelse(str_detect(i,"highelev"),"Southern_highelev",
                                 ifelse(str_detect(i,"CAC"),"CAC","LM")))),
    group_j=ifelse(str_detect(j,"allopatric"),"Southern_allopatric",
                   ifelse(str_detect(j,"sympatric"),"Southern_sympatric",
                          ifelse(str_detect(j,"highelev"),"Southern_highelev",
                                 ifelse(str_detect(j,"CAC"),"CAC","LM"))))
  ) %>%
  mutate(comparison_type=ifelse(i==j,"self",
                                ifelse(group_i==group_j,"Within-group Across-cohort",
                                       ifelse(region_i==region_j,
                                              ifelse(cohort_i==cohort_j,"Within-region Within-cohort","Within-region Across-cohort"),
                                              ifelse(cohort_i==cohort_j,"Across-region Within-cohort","Across-region Across-cohort")))),
         region=ifelse(region_i==region_j,region_i,"Across"))

###Try an anova approach for all pairwise comparisons
###matrix should be correlation_matrix_with_trueallos_missingness
###tests effect of features+j on i
anova_pairwise<-function(i,j,matrix) {
  base_string="I(relative_chr_pos^2) + ngenes + M_per_bp_1Mb + M_per_bp_50kb + missingness_"
  full_model_string=paste0("frequency_",i," ~ ",base_string,i," + frequency_",j)
  model<-lm(as.formula(full_model_string),data=matrix)
  anova_results<-as.data.frame(anova(model)) %>% rownames_to_column("predictor") %>%
    mutate(var_explained=`Sum Sq`/sum(`Sum Sq`)) %>% rename("pval"="Pr(>F)") %>% select(predictor,var_explained,pval) %>% 
    filter(predictor!="Residuals") %>%
    mutate(predictor_group=ifelse(str_starts(predictor,"missingness"),"missingness",
                                  ifelse(str_starts(predictor,"frequency"),"testgroup","genome_structure"))) %>%
    group_by(predictor_group) %>% summarize(var_explained=sum(var_explained),pval=min(pval)) %>%
    pivot_wider(names_from=predictor_group,values_from=c(var_explained,pval)) %>%
    mutate(i=i,j=j)
  return(anova_results)
}
anova_results_all<-data.frame()
for (i in key_groups_customorder2_withtrueallos) {
  for (j in key_groups_customorder2_withtrueallos) {
    if (i!=j) {
      anova_results_all<-bind_rows(anova_results_all,anova_pairwise(i,j,correlation_matrix_with_trueallos_missingness))
    }
  }
}
anova_results_main<-anova_results_all %>% select(-ends_with("genome_structure"),-ends_with("missingness")) %>% select(i,j,var_explained_testgroup,pval_testgroup)
anova_results_pre1<-anova_results_all %>% select(-ends_with("testgroup"),-ends_with("missingness"),-j) %>% group_by(i) %>%
  summarize(pval_testgroup=max(pval_genome_structure),
            var_explained_testgroup=max(var_explained_genome_structure)) %>%
  mutate(j="genome_structure") %>% select(i,j,var_explained_testgroup,pval_testgroup)
anova_results_pre2<-anova_results_all %>% select(-ends_with("testgroup"),-ends_with("genome_structure"),-j) %>% group_by(i) %>%
  summarize(pval_testgroup=max(pval_missingness),
            var_explained_testgroup=max(var_explained_missingness)) %>%
  mutate(j="missingness") %>% select(i,j,var_explained_testgroup,pval_testgroup)
anova_results_ready<-bind_rows(anova_results_main,anova_results_pre1,anova_results_pre2) %>%
  mutate(pval_bonferroni=pval_testgroup*110,
         significance=label_significance(pval_bonferroni),
         significance=ifelse(significance=="n.s.","",significance),
         prettyval=paste0(significance,"\n",format(round(var_explained_testgroup*100,1),digits=1)),
         i=factor(i,levels=key_groups_customorder2_withtrueallos),
         j=factor(j,levels=c("genome_structure","missingness",key_groups_customorder2_withtrueallos)))

####raw coverage info
rawcoverage<-read.table(file.path(thisfolder,"Metadata/ancestry_fastqs_stats.txt"),header=T) %>%
  mutate(sample=str_remove(file,".gz")) %>% right_join(sampleset_info,by = "sample") %>% mutate(num_seqs=as.numeric(str_remove_all(num_seqs,",")))
rawcoverage %>% filter(group %in% c("Southern_allopatric","Southern_highelev","Southern_sympatric")) %>% summarize(mean(num_seqs),median(num_seqs),min(num_seqs),max(num_seqs))
rawcoverage %>% filter(group %in% c("CAC","LM")) %>% filter(!is.na(num_seqs)) %>% summarize(mean(num_seqs),median(num_seqs),min(num_seqs),max(num_seqs))

  