
################################################################################
#####    Code to run supplemental analysis on reduced set of windows     ######
################################################################################

###read summary local PCA stats from all CAC
PCA_fit_stats_allwindows_CAC<-read.table(file.path(thisfolder,"CAC.PCA_fit_stats.txt"),header=T)
PCA_fit_stats_allwindows_CAC_frequencies<-all_groups_window_freqs %>% filter(str_starts(keygroup,"CAC")) %>%
  mutate(region=paste0(chrom,".",windowstart,".",windowend)) %>%
  left_join(PCA_fit_stats_allwindows_CAC,by="region",relationship="many-to-one")
outliers_cutoffs<-PCA_fit_stats_allwindows_CAC_frequencies %>% group_by(keygroup) %>% summarize(cutoff95=quantile(nasutus_frequency,0.95),
                                                                                                cutoff05=quantile(nasutus_frequency,0.05)) %>%
  right_join(PCA_fit_stats_allwindows_CAC_frequencies,by="keygroup") %>%
  mutate(outlier_type=ifelse(nasutus_frequency>cutoff95,"top 5%",ifelse(nasutus_frequency<cutoff05,"bottom 5%","other")))

###read summary local PCA stats from all Southern
PCA_fit_stats_allwindows_southern<-read.table(file.path(thisfolder,"southern_core.PCA_fit_stats.txt"),header=T)
PCA_fit_stats_allwindows_southern_frequencies<-all_groups_window_freqs %>% filter(keygroup %in% c("Southern_sympatric_guttatus_174","Southern_allopatric_guttatus_61","Southern_sympatric_nasutus_75")) %>%
  mutate(region=paste0(chrom,".",windowstart,".",windowend)) %>%
  left_join(PCA_fit_stats_allwindows_southern,by="region",relationship="many-to-one")
outliers_cutoffs_southern<-PCA_fit_stats_allwindows_southern_frequencies %>% group_by(keygroup) %>% summarize(cutoff95=quantile(nasutus_frequency,0.95),
                                                                                                              cutoff05=quantile(nasutus_frequency,0.05)) %>%
  right_join(PCA_fit_stats_allwindows_southern_frequencies,by="keygroup") %>%
  mutate(outlier_type=ifelse(nasutus_frequency>cutoff95,"top 5%",ifelse(nasutus_frequency<cutoff05,"bottom 5%","other")))
shared_outliers<-outliers_cutoffs_southern %>% filter(outlier_type=="top 5%") %>% group_by(region) %>% 
  summarize(ngroups=n()) %>% filter(ngroups==2) %>% pull(region)
with_shared_outliers<-outliers_cutoffs_southern %>% mutate(shared=(region %in% shared_outliers))

###get list of excluded and good windows
###total windows to start: 3139
###final reduced window set: 2069
exclude_windows_CAC<-PCA_fit_stats_allwindows_CAC %>% filter(modelfit.r.squared<0.9 | het.scaled.zscore<(-2) | het.scaled.zscore>2)
exclude_windows_southern<-PCA_fit_stats_allwindows_southern %>% filter(modelfit.r.squared<0.9 | het.scaled.zscore<(-2) | het.scaled.zscore>2)
all_exclude<-c(exclude_windows_CAC$region, exclude_windows_southern$region) %>% sort() %>% unique()
keep_windows<-PCA_fit_stats_allwindows_CAC %>% filter(!(region %in% all_exclude)) %>% pull(region)

nasutus_hetz<-all_groups_window_freqs %>% filter(keygroup %in% c("CAC_nasutus_58","Southern_sympatric_nasutus_75")) %>% 
  mutate(region=paste0(chrom,".",windowstart,".",windowend)) %>% group_by(region) %>%
  summarize(n_hets=sum(n_1,na.rm=T))
nasutus_hetz %>% ggplot(aes(x=n_hets)) + geom_histogram(binwidth=1) + theme_bw()
nasutus_hetz %>% filter(region %in% keep_windows) %>% ggplot(aes(x=n_hets)) + geom_histogram(binwidth=1) + theme_bw()
nasutus_highhetz<-nasutus_hetz %>% filter(n_hets>2) %>% pull(region)

keep_windows_nohighhets<-keep_windows[!(keep_windows %in% nasutus_highhetz)]

reduced_dataset_window_freqs<-all_groups_window_freqs %>% mutate(region=paste0(chrom,".",windowstart,".",windowend)) %>%
  filter(region %in% keep_windows_nohighhets)
reduced_dataset_windowlist<-reduced_dataset_window_freqs %>% pull(chrom.windowstart.windowend) %>% unique()


###make final reduced dataset
windowed_sitefilter_reduced<-windowed_ancestry_allsamples_long %>%
  filter(chrom.windowstart.windowend %in% reduced_dataset_windowlist)

###get individual summaries from site-filtered data
individuals_sitefiltered_reduced<-windowed_sitefilter_reduced %>% group_by(sample) %>%
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
final_individual_list_reduced<-individuals_sitefiltered_reduced %>% 
  filter(sample!="KK047_read_1.fastq", ###two data runs don't match
         sample!="DPR_plot2I_BL_read_1.fastq", ###clusters with northern samples
         !(sample %in% sookensis_exclusions),
         !(str_detect(group,"panel")),
         frac_called_windows>=0.5) %>% 
  left_join(angsd_PCA,by="sampleID") %>% left_join(transitions_summary_filtered,by="sampleID") %>%
  mutate(freq_group=ifelse(nasutus_proportion_windows<0.15,"guttatus",
                           ifelse(nasutus_proportion_windows>0.85,"nasutus",
                                  "hybrid")))
final_individual_list_reduced %>% group_by(group,freq_group) %>% summarize(n=n()) %>% pivot_wider(names_from=freq_group,values_from=n)

panelinfo_reduced<-individuals_sitefiltered_reduced %>% filter(str_detect(group,"panel")) %>%
  mutate(Line=str_remove(sampleID,".sampled")) %>% select(Line,nasutus_proportion_windows,group) %>% arrange(group,nasutus_proportion_windows,Line)
write.table(panelinfo_reduced,file.path(thisfolder,"Metadata/gutnas_panel_info_reduced.txt"),sep="\t",row.names=F,col.names=T,quote=F)
southern_allopatric_info_reduced<-individuals_sitefiltered_reduced %>% filter(group=="Southern_allopatric") %>% 
  mutate(Pop=str_sub(sampleID,1,3)) %>% arrange(Pop,nasutus_proportion_windows)
southern_allopatric_info_reduced %>% group_by(Pop) %>% summarize(n_total=n(),n_above_1=sum(nasutus_proportion_windows>=0.01,na.rm=T))

###final windows
final_window_set_reduced<-windowed_sitefilter_reduced %>% filter(sample %in% final_individual_list_reduced$sample)
enumerate_final_windows_reduced<-final_window_set_reduced %>% select(chrom.windowstart.windowend) %>% unique() %>% mutate(window_number=c(1:nrow(.)))
chrom_blocks_final_windows_reduced<-enumerate_final_windows_reduced %>%
  separate(chrom.windowstart.windowend,into=c("chrom","windowstart","windowend"),sep=":",remove=F) %>%
  mutate(blockend=(chrom!=lead(chrom))) %>% filter(blockend)
final_window_set_enumerated_reduced<-final_window_set_reduced %>% left_join(enumerate_final_windows_reduced) %>% left_join(final_individual_list_reduced %>% select(sample,group))

group_paints_reduced<-list()
samplelists_reduced<-list()
for (i in unique(final_window_set_enumerated_reduced$group)) {
  samplelists_reduced[[i]]<-final_individual_list_reduced %>% filter(group==i) %>% 
    arrange(nasutus_proportion_windows) %>% 
    mutate(sample_ordered=factor(sample,levels=sample)) %>% select(sample_ordered)
  subdata<-final_window_set_enumerated_reduced %>% filter(group==i) %>% mutate(sample=factor(sample,levels=samplelists[[i]]$sample_ordered))
  group_paints_reduced[[i]]<-paint_ancestry(subdata,chrom_blocks_final_windows_reduced,write_sampleIDs = F)
}


###define key groups
key_groups_reduced<-final_individual_list_reduced %>% group_by(group,freq_group) %>% summarize(n=n()) %>%
  mutate(keygroup=paste0(group,"_",freq_group,"_",n))
focal_nine_reduced<-key_groups_reduced %>% filter(n>=10)

###get_keygroup_frequencies
all_groups_window_freqs_reduced<-data.frame()
for (g in key_groups_reduced$keygroup) {
  thispop=key_groups_reduced %>% filter(keygroup==g) %>% pull(group)
  thisfreq=key_groups_reduced %>% filter(keygroup==g) %>% pull(freq_group)
  these_individuals<-final_individual_list_reduced %>% filter(group==thispop,freq_group==thisfreq) %>% pull(sample)
  window_freqs<-final_window_set_enumerated_reduced %>% filter(sample %in% these_individuals) %>% 
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
  all_groups_window_freqs_reduced<-bind_rows(all_groups_window_freqs_reduced,window_freqs)
}
panel_GUT_freqs_reduced<-windowed_ancestry_allsamples_long %>% 
  filter(chrom.windowstart.windowend %in% final_window_set_reduced$chrom.windowstart.windowend,
         sample %in% panel_GUT_list$sample) %>% 
  left_join(enumerate_final_windows_reduced) %>% 
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

key_groups_freqorder_reduced<-all_groups_window_freqs_reduced %>% filter(keygroup %in% focal_nine_reduced$keygroup) %>% 
  group_by(keygroup) %>% summarize(mean_freq=mean(nasutus_frequency,na.rm=T)) %>%
  arrange(mean_freq) %>% mutate(keygroup=factor(keygroup,levels=keygroup))

###only 'truly' allopatric pops
true_allos_reduced<-final_individual_list_reduced %>% filter(group=="Southern_allopatric") %>% mutate(Pop=str_sub(sampleID,1,3)) %>%
  filter(Pop %in% c("COP","GCH","RHI"))
pseudo_allos_reduced<-final_individual_list_reduced %>% filter(group=="Southern_allopatric") %>% mutate(Pop=str_sub(sampleID,1,3)) %>%
  filter(Pop %in% c("MOC","RCF","BFR"))
true_allos_freqs_reduced<-windowed_ancestry_allsamples_long %>% 
  filter(chrom.windowstart.windowend %in% final_window_set_reduced$chrom.windowstart.windowend,
         sample %in% true_allos_reduced$sample) %>% 
  left_join(enumerate_final_windows_reduced) %>% 
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
pseudo_allos_freqs_reduced<-windowed_ancestry_allsamples_long %>% 
  filter(chrom.windowstart.windowend %in% final_window_set_reduced$chrom.windowstart.windowend,
         sample %in% pseudo_allos_reduced$sample) %>% 
  left_join(enumerate_final_windows_reduced) %>% 
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

with_rr_reduced<-all_groups_window_freqs_reduced %>%
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

correlation_matrix_reduced<-with_rr_reduced %>% filter(keygroup %in% focal_nine_reduced$keygroup) %>% select(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos,keygroup,nasutus_frequency) %>%
  pivot_wider(id_cols=c(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos),names_from=keygroup,values_from=nasutus_frequency) %>%
  arrange(chrom,windowend)
correlation_matrix_withmissingness_reduced<-with_rr_reduced %>% filter(keygroup %in% focal_nine_reduced$keygroup) %>% select(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos,keygroup,nasutus_frequency,frac_called_samples) %>%
  rename(frequency=nasutus_frequency,missingness=frac_called_samples) %>%
  pivot_wider(id_cols=c(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos),names_from=keygroup,values_from=c(frequency,missingness),names_sep = "_") %>%
  arrange(chrom,windowend)
correlation_matrix_with_trueallos_reduced<-true_allos_freqs_reduced %>% ungroup() %>%
  select(chrom,windowend,nasutus_frequency) %>%
  rename("Southern_true_allopatric"="nasutus_frequency") %>% 
  right_join(
    pseudo_allos_freqs_reduced %>% ungroup() %>%
      select(chrom,windowend,nasutus_frequency) %>%
      rename("Southern_pseudo_allopatric"="nasutus_frequency")
  ) %>%
  right_join(correlation_matrix_reduced)
correlation_matrix_with_trueallos_missingness_reduced<-true_allos_freqs_reduced %>% ungroup() %>%
  select(chrom,windowend,nasutus_frequency,frac_called_samples) %>%
  rename("frequency_Southern_true_allopatric"="nasutus_frequency","missingness_Southern_true_allopatric"="frac_called_samples") %>%
  right_join(
    pseudo_allos_freqs_reduced %>% ungroup() %>%
      select(chrom,windowend,nasutus_frequency,frac_called_samples) %>%
      rename("frequency_Southern_pseudo_allopatric"="nasutus_frequency","missingness_Southern_pseudo_allopatric"="frac_called_samples")
  ) %>%
  right_join(correlation_matrix_withmissingness_reduced)

correlation_matrix_withsquares_reduced<-correlation_matrix_reduced %>% mutate(relative_chr_pos_squared=relative_chr_pos^2,
                                                              ngenes_squared=ngenes,
                                                              M_per_bp_1Mb_squared=M_per_bp_1Mb^2,
                                                              M_per_bp_50kb_squared=M_per_bp_50kb^2) %>%
  select(relative_chr_pos,ngenes,M_per_bp_1Mb,M_per_bp_50kb,
         relative_chr_pos_squared,ngenes_squared,M_per_bp_1Mb_squared,M_per_bp_50kb_squared,
         starts_with("CAC"),starts_with("LM"),starts_with("Southern"))

correlation_matrix_long_reduced<-with_rr_reduced %>% filter(keygroup %in% focal_nine_reduced$keygroup) %>% select(chrom,windowend,M_per_bp_50kb,M_per_bp_1Mb,ngenes,relative_chr_pos,keygroup,nasutus_frequency) %>%
  filter(!is.na(nasutus_frequency)) %>%
  arrange(chrom,windowend)

correlation_matrix_long_rescaled_reduced<-correlation_matrix_long_reduced %>% mutate(rr_1Mb_rescaled=rescale_z(M_per_bp_1Mb),
                                            rr_50kb_rescaled=rescale_z(M_per_bp_50kb),
                                            ngenes_rescaled=rescale_z(ngenes))
fullmodel_norescale_reduced<-lmer(nasutus_frequency~relative_chr_pos+ngenes+M_per_bp_1Mb+M_per_bp_50kb+(1|keygroup),data=correlation_matrix_long_reduced)
anova(fullmodel_norescale_reduced,type="II")
fullmodel_reduced<-lmer(nasutus_frequency~relative_chr_pos+ngenes_rescaled+rr_1Mb_rescaled+rr_50kb_rescaled+(1|keygroup),data=correlation_matrix_long_rescaled_reduced)
anova(fullmodel_reduced,type="I")
anova(fullmodel_reduced,type="III")

###get residuals for full interaction with squares model
each_rr_models_reduced<-list()
each_rr_newmodels_reduced<-list()
each_rr_fullinteract_models_reduced<-list()
each_missingness_models_reduced<-list()
correlation_matrix_with_residuals_reduced<-correlation_matrix_with_trueallos_reduced %>% select(-c(chrom,windowend,M_per_bp_50kb:relative_chr_pos))
correlation_matrix_with_residuals_missingness_reduced<-correlation_matrix_with_trueallos_missingness_reduced %>% select(-c(chrom,windowend,M_per_bp_50kb:relative_chr_pos),-starts_with("missingness_")) %>%
  rename_with(function(x){str_remove(x,"frequency_")},.cols=starts_with("frequency"))

for (i in c(focal_nine_reduced$keygroup,"Southern_pseudo_allopatric","Southern_true_allopatric")) {
  baseformula_string="relative_chr_pos + ngenes + M_per_bp_1Mb + M_per_bp_50kb"
  newmodel_string="I(relative_chr_pos^2) + ngenes + M_per_bp_1Mb + M_per_bp_50kb"
  fullinteract_string="relative_chr_pos * ngenes * M_per_bp_1Mb * M_per_bp_50kb *
                       I(relative_chr_pos^2) * I(ngenes^2) * I(M_per_bp_1Mb^2) * I(M_per_bp_50kb^2)"
  withmissingness_string="I(relative_chr_pos^2) + ngenes + M_per_bp_1Mb + M_per_bp_50kb + missingness_"
  basemodel<-lm(as.formula(paste0(i," ~ ",baseformula_string)),data=correlation_matrix_with_trueallos_reduced)
  newmodel<-lm(as.formula(paste0(i," ~ ",newmodel_string)),data=correlation_matrix_with_trueallos_reduced)
  fullinteractmodel<-lm(as.formula(paste0(i," ~ ",fullinteract_string)),data=correlation_matrix_with_trueallos_reduced)
  withmissingness_model<-lm(as.formula(paste0("frequency_",i," ~ ",withmissingness_string,i)),data=correlation_matrix_with_trueallos_missingness_reduced)
  correlation_matrix_with_residuals_reduced<-correlation_matrix_with_residuals_reduced %>% 
    add_column(!!paste0(i,".residuals"):=residuals(basemodel)) %>%
    add_column(!!paste0(i,".residuals_newmodel"):=residuals(newmodel))
  correlation_matrix_with_residuals_missingness_reduced<-correlation_matrix_with_residuals_missingness_reduced %>%
    add_column(!!paste0(i,".residuals"):=residuals(basemodel)) %>%
    add_column(!!paste0(i,".residuals_newmodel"):=residuals(newmodel)) %>%
    add_column(!!paste0(i,".residuals_missingness"):=residuals(withmissingness_model))
  each_rr_models_reduced[[i]]<-basemodel
  each_rr_newmodels_reduced[[i]]<-newmodel
  each_missingness_models_reduced[[i]]<-withmissingness_model
}

key_groups_customorder_reduced<-key_groups_freqorder_reduced %>% mutate(keygroup=factor(keygroup,levels=c(
  "LM_hybrid_56","CAC_hybrid_164","CAC_guttatus_156",
  "Southern_sympatric_guttatus_174","Southern_allopatric_guttatus_61",
  "Southern_highelev_hybrid_15","Southern_highelev_guttatus_12",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
))) %>% arrange(keygroup)
key_groups_customorder2_reduced<-key_groups_freqorder_reduced %>% mutate(keygroup=factor(keygroup,levels=c(
  "CAC_guttatus_156","CAC_hybrid_164","LM_hybrid_56",
  "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174",
  "Southern_highelev_guttatus_12","Southern_highelev_hybrid_15",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
))) %>% arrange(keygroup)
key_groups_customorder_withtrueallos_reduced<-c(
  "LM_hybrid_56","CAC_hybrid_164","CAC_guttatus_156",
  "Southern_sympatric_guttatus_174","Southern_allopatric_guttatus_61","Southern_true_allopatric",
  "Southern_highelev_hybrid_15","Southern_highelev_guttatus_12",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
)
key_groups_customorder2_withtrueallos_reduced<-c(
  "CAC_guttatus_156","CAC_hybrid_164","LM_hybrid_56",
  "Southern_true_allopatric","Southern_pseudo_allopatric",
  "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174",
  "Southern_highelev_guttatus_12","Southern_highelev_hybrid_15",
  "CAC_nasutus_58",
  "Southern_sympatric_nasutus_75"
)


all_correlations_withcorrections_reduced<-as.data.frame(cor(correlation_matrix_with_residuals_reduced)) %>% 
  rownames_to_column(var="group1") %>%
  pivot_longer(cols=-group1,names_to="compare_to",values_to="r") %>%
  separate(compare_to,into=c("group2","comparison"),sep="\\.",fill = "right") %>%
  mutate(comparison=ifelse(is.na(comparison),"uncorrected",comparison)) %>%
  pivot_wider(names_from=comparison,values_from=r) %>%
  mutate(group1=factor(group1,levels=key_groups_customorder_withtrueallos_reduced),
         group2=factor(group2,levels=key_groups_customorder_withtrueallos_reduced),
         #whichval=ifelse(as.integer(group1)>as.integer(group2),residuals,uncorrected),
         #whichval2=ifelse(as.integer(group1)>as.integer(group2),residuals_newmodel,residuals_newmodel),
         uncorrected_prettyval=format(round(uncorrected,2),digits=2),
         residuals_newmodel_prettyval=format(round(residuals_newmodel,2),digits=2)) %>%
  filter(!is.na(group1),!is.na(group2))
all_correlations_withcorrections_missingness_reduced<-as.data.frame(cor(correlation_matrix_with_residuals_missingness_reduced)) %>% 
  rownames_to_column(var="group1") %>%
  pivot_longer(cols=-group1,names_to="compare_to",values_to="r") %>%
  separate(compare_to,into=c("group2","comparison"),sep="\\.",fill = "right") %>%
  mutate(comparison=ifelse(is.na(comparison),"uncorrected",comparison)) %>%
  pivot_wider(names_from=comparison,values_from=r) %>%
  mutate(group1=factor(group1,levels=key_groups_customorder_withtrueallos_reduced),
         group2=factor(group2,levels=key_groups_customorder_withtrueallos_reduced),
         #whichval=ifelse(as.integer(group1)>as.integer(group2),residuals_newmodel,uncorrected),
         #whichval2=ifelse(as.integer(group1)>as.integer(group2),residuals_missingness,residuals_missingness),
         residuals_newmodel_prettyval=format(round(residuals_newmodel,2),digits=2),
         residuals_missingness_prettyval=format(round(residuals_missingness,2),digits=2)) %>%
  filter(!is.na(group1),!is.na(group2))

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
anova_results_all_reduced<-data.frame()
for (i in key_groups_customorder2_withtrueallos_reduced) {
  for (j in key_groups_customorder2_withtrueallos_reduced) {
    if (i!=j) {
      anova_results_all_reduced<-bind_rows(anova_results_all_reduced,anova_pairwise(i,j,correlation_matrix_with_trueallos_missingness_reduced))
    }
  }
}
anova_results_main_reduced<-anova_results_all_reduced %>% select(-ends_with("genome_structure"),-ends_with("missingness")) %>% select(i,j,var_explained_testgroup,pval_testgroup)
anova_results_pre1_reduced<-anova_results_all_reduced %>% select(-ends_with("testgroup"),-ends_with("missingness"),-j) %>% group_by(i) %>%
  summarize(pval_testgroup=max(pval_genome_structure),
            var_explained_testgroup=max(var_explained_genome_structure)) %>%
  mutate(j="genome_structure") %>% select(i,j,var_explained_testgroup,pval_testgroup)
anova_results_pre2_reduced<-anova_results_all_reduced %>% select(-ends_with("testgroup"),-ends_with("genome_structure"),-j) %>% group_by(i) %>%
  summarize(pval_testgroup=max(pval_missingness),
            var_explained_testgroup=max(var_explained_missingness)) %>%
  mutate(j="missingness") %>% select(i,j,var_explained_testgroup,pval_testgroup)
anova_results_ready_reduced<-bind_rows(anova_results_main_reduced,anova_results_pre1_reduced,anova_results_pre2_reduced) %>%
  mutate(pval_bonferroni=pval_testgroup*110,
         significance=label_significance(pval_bonferroni),
         significance=ifelse(significance=="n.s.","",significance),
         prettyval=paste0(significance,"\n",format(round(var_explained_testgroup*100,1),digits=1)),
         i=factor(i,levels=key_groups_customorder2_withtrueallos_reduced),
         j=factor(j,levels=c("genome_structure","missingness",key_groups_customorder2_withtrueallos_reduced)))


###predict ancestry based on genome-structure data alone
ancestry_predictions_reduced<-correlation_matrix_with_trueallos_reduced %>% select(chrom,windowend)
ancestry_predictions_newmodel_reduced<-correlation_matrix_with_trueallos_reduced %>% select(chrom,windowend)
ancestry_predictions_missingness_reduced<-correlation_matrix_with_trueallos_reduced %>% select(chrom,windowend)
for (i in c(focal_nine_reduced$keygroup,"Southern_true_allopatric","Southern_pseudo_allopatric")) {
  ancestry_predictions_reduced<-ancestry_predictions_reduced %>% add_column(!!i:=predict(each_rr_models_reduced[[i]]))
  ancestry_predictions_newmodel_reduced<-ancestry_predictions_newmodel_reduced %>% add_column(!!i:=predict(each_rr_newmodels_reduced[[i]]))
  ancestry_predictions_missingness_reduced<-ancestry_predictions_missingness_reduced %>% add_column(!!i:=predict(each_missingness_models_reduced[[i]]))
}
all_groups_window_freqs_with_trueallos_reduced<-all_groups_window_freqs_reduced %>% bind_rows(true_allos_freqs_reduced) %>% bind_rows(pseudo_allos_freqs_reduced)
ancestry_predictions_ready_reduced<-ancestry_predictions_reduced %>% pivot_longer(-c(chrom,windowend),names_to="keygroup",values_to="prediction") %>%
  left_join(ancestry_predictions_newmodel_reduced %>% pivot_longer(-c(chrom,windowend),names_to="keygroup",values_to="prediction_newmodel")) %>%
  left_join(ancestry_predictions_missingness_reduced %>% pivot_longer(-c(chrom,windowend),names_to="keygroup",values_to="prediction_missingness")) %>%
  left_join(all_groups_window_freqs_with_trueallos_reduced %>% filter(keygroup %in% c(focal_nine_reduced$keygroup,"Southern_true_allopatric","Southern_pseudo_allopatric")),by = c("chrom","windowend","keygroup")) %>%
  pivot_longer(c(nasutus_frequency,prediction,prediction_newmodel,prediction_missingness),names_to="statistic",values_to="nasutus_frequency")

prediction_deltas_allmodels_reduced<-ancestry_predictions_ready_reduced %>%
  pivot_wider(names_from = "statistic",values_from="nasutus_frequency") %>%
  mutate(delta_predict=nasutus_frequency-prediction_newmodel,
         delta_predict_missingness=nasutus_frequency-prediction_missingness)

###tests of window outlier overlap 
focal_outliers_95_uncorrected_reduced<-list()
focal_outliers_95_corrected_reduced<-list()
focal_outliers_95_missingness_reduced<-list()
permuted_outliers_95_reduced<-list()
for (i in c(key_groups_customorder2_withtrueallos_reduced)) {
  subdata<-prediction_deltas_allmodels_reduced %>% filter(keygroup==i,!is.na(delta_predict)) %>% select(chrom,windowend,chrom.windowstart.windowend,delta_predict,delta_predict_missingness,nasutus_frequency)
  meanvalue<-mean(subdata$nasutus_frequency)
  permuted<-data.frame()
  if (meanvalue>0.5) {
    cutoff_uncorrected=quantile(subdata$nasutus_frequency,0.05)
    cutoff_corrected=quantile(subdata$delta_predict,0.05)
    cutoff_corrected_missingness=quantile(subdata$delta_predict_missingness,0.05)
    focal_outliers_95_uncorrected_reduced[[i]]<-subdata %>% filter(nasutus_frequency<=cutoff_uncorrected)
    focal_outliers_95_corrected_reduced[[i]]<-subdata %>% filter(delta_predict<=cutoff_corrected)
    focal_outliers_95_missingness_reduced[[i]]<-subdata %>% filter(delta_predict_missingness<=cutoff_corrected_missingness)
    
  } else {
    cutoff_uncorrected=quantile(subdata$nasutus_frequency,0.95)
    cutoff_corrected=quantile(subdata$delta_predict,0.95)
    cutoff_corrected_missingness=quantile(subdata$delta_predict_missingness,0.95)
    focal_outliers_95_uncorrected_reduced[[i]]<-subdata %>% filter(nasutus_frequency>=cutoff_uncorrected)
    focal_outliers_95_corrected_reduced[[i]]<-subdata %>% filter(delta_predict>=cutoff_corrected)
    focal_outliers_95_missingness_reduced[[i]]<-subdata %>% filter(delta_predict_missingness>=cutoff_corrected_missingness)
  }
  for (p in 1:1000) {
    one_permute<-sample_n(subdata,nrow(focal_outliers_95_corrected_reduced[[i]])) %>% mutate(permutation=p)
    permuted<-permuted %>% bind_rows(one_permute)
  }
  permuted_outliers_95_reduced[[i]]<-permuted
}
count_overlaps_reduced<-data.frame()
overlap_permuter_reduced<-function(p,i,j,corrected=T,missingness=T) {
  if (missingness) {
    return(inner_join(focal_outliers_95_missingness_reduced[[i]],permuted_outliers_95_reduced[[j]] %>% filter(permutation==p),by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow())
  } else if (corrected) {
    return(inner_join(focal_outliers_95_corrected_reduced[[i]],permuted_outliers_95_reduced[[j]] %>% filter(permutation==p),by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow())
  } else {
    return(inner_join(focal_outliers_95_uncorrected_reduced[[i]],permuted_outliers_95_reduced[[j]] %>% filter(permutation==p),by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow())
  }
}
for (i in key_groups_customorder2_withtrueallos_reduced) {
  for (j in key_groups_customorder2_withtrueallos_reduced) {
    #print(paste0(i," vs. ",j))
    N_i<-prediction_deltas_allmodels_reduced %>% filter(keygroup==i,!is.na(delta_predict)) %>% nrow()
    N_j<-prediction_deltas_allmodels_reduced %>% filter(keygroup==j,!is.na(delta_predict)) %>% nrow()
    
    n_outliers_uncorrected_i<-focal_outliers_95_uncorrected_reduced[[i]] %>% nrow()
    n_outliers_corrected_i<-focal_outliers_95_corrected_reduced[[i]] %>% nrow()
    n_outliers_missingness_i<-focal_outliers_95_missingness_reduced[[i]] %>% nrow()
    n_outliers_uncorrected_j<-focal_outliers_95_uncorrected_reduced[[j]] %>% nrow()
    n_outliers_corrected_j<-focal_outliers_95_corrected_reduced[[j]] %>% nrow()
    n_outliers_missingness_j<-focal_outliers_95_missingness_reduced[[j]] %>% nrow()
    
    observed_overlap_uncorrected<-inner_join(focal_outliers_95_uncorrected_reduced[[i]],focal_outliers_95_uncorrected_reduced[[j]],by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow()
    observed_overlap_corrected<-inner_join(focal_outliers_95_corrected_reduced[[i]],focal_outliers_95_corrected_reduced[[j]],by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow()
    observed_overlap_missingness<-inner_join(focal_outliers_95_missingness_reduced[[i]],focal_outliers_95_missingness_reduced[[j]],by=c("chrom","windowend","chrom.windowstart.windowend")) %>% nrow()
    
    expected_overlap_uncorrected<-n_outliers_uncorrected_i*n_outliers_uncorrected_j/(N_j)
    expected_overlap_corrected<-n_outliers_corrected_i*n_outliers_corrected_j/(N_j)
    expected_overlap_missingness<-n_outliers_missingness_i*n_outliers_missingness_j/(N_j)
    #compare observed i to permuted j
    permute_overlaps<-sapply(c(1:1000),overlap_permuter_reduced,i=i,j=j)
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
    count_overlaps_reduced<-count_overlaps_reduced %>% bind_rows(count_overlaps_line)
  }
}

write.table(count_overlaps_reduced,file=file.path(thisfolder,"Outlier_deltapredict_overlap_comparisons_reduced.txt"),quote = F,sep = "\t",row.names=F)


###Mantel regions corrected for IM767ref
Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced<-Mantel_regions_windows_withGOIs %>% filter(Name!="pTAC14") %>% left_join(prediction_deltas_allmodels_reduced,by=c("chrom","windowend"),relationship="one-to-many") %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos_reduced))  ###remove pTAC14 because the window is already counted once
Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced<-Mantel_regions_windows_nolarge %>% left_join(prediction_deltas_allmodels_reduced,by=c("chrom","windowend"),relationship="one-to-many") %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos_reduced))


###do outliers occur in windows with QTL more often than expected?

outliers_in_windows_summary_reduced<-data.frame()
for (i in key_groups_customorder2_withtrueallos_reduced) {
  n_total_windows<-ancestry_predictions_ready_reduced %>% filter(keygroup==i) %>% nrow()
  n_QTLwindows<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>% filter(keygroup==i) %>% nrow()
  n_outlier_windows_uncorrected<-focal_outliers_95_uncorrected_reduced[[i]] %>% nrow()
  n_outlier_windows_corrected<-focal_outliers_95_corrected_reduced[[i]] %>% nrow()
  n_outliers_in_QTL_uncorrected<-focal_outliers_95_uncorrected_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_corrected<-focal_outliers_95_corrected_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  dataline<-data.frame(keygroup=i,n_total_windows=n_total_windows,n_QTLwindows=n_QTLwindows,
                       n_outlier_windows_uncorrected=n_outlier_windows_uncorrected,n_outlier_windows_corrected=n_outlier_windows_corrected,
                       n_outliers_in_QTL_uncorrected=n_outliers_in_QTL_uncorrected,n_outliers_in_QTL_corrected=n_outliers_in_QTL_corrected)
  outliers_in_windows_summary_reduced<-bind_rows(outliers_in_windows_summary_reduced,dataline)
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
Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus_reduced<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>% filter(Predicted_effect=="N+")
Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus_reduced<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>% filter(Predicted_effect=="N-")
outliers_in_windows_summary_split_reduced<-data.frame()
for (i in key_groups_customorder2_withtrueallos_reduced) {
  n_total_windows<-ancestry_predictions_ready_reduced %>% filter(keygroup==i,statistic=="nasutus_frequency") %>% nrow()
  n_QTLwindows_plus<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus_reduced %>% filter(keygroup==i) %>% nrow()
  n_QTLwindows_minus<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus_reduced %>% filter(keygroup==i) %>% nrow()
  n_outlier_windows_uncorrected<-focal_outliers_95_uncorrected_reduced[[i]] %>% nrow()
  n_outlier_windows_corrected<-focal_outliers_95_corrected_reduced[[i]] %>% nrow()
  n_outlier_windows_missingness<-focal_outliers_95_missingness_reduced[[i]] %>% nrow()
  n_outliers_in_QTL_uncorrected_plus<-focal_outliers_95_uncorrected_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_corrected_plus<-focal_outliers_95_corrected_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_missingness_plus<-focal_outliers_95_missingness_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nplus_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_uncorrected_minus<-focal_outliers_95_uncorrected_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_corrected_minus<-focal_outliers_95_corrected_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  n_outliers_in_QTL_missingness_minus<-focal_outliers_95_missingness_reduced[[i]] %>% inner_join(Mantel_regions_withGOIs_freqs_withpredictions_newmodel_Nminus_reduced %>% filter(keygroup==i),by=c("chrom","windowend")) %>% nrow()
  dataline<-data.frame(keygroup=i,n_total_windows=n_total_windows,n_QTLwindows_plus=n_QTLwindows_plus,n_QTLwindows_minus=n_QTLwindows_minus,
                       n_outlier_windows_uncorrected=n_outlier_windows_uncorrected,n_outlier_windows_corrected=n_outlier_windows_corrected,n_outlier_windows_missingness=n_outlier_windows_missingness,
                       n_outliers_in_QTL_uncorrected_minus=n_outliers_in_QTL_uncorrected_minus,n_outliers_in_QTL_corrected_minus=n_outliers_in_QTL_corrected_minus,n_outliers_in_QTL_missingness_minus=n_outliers_in_QTL_missingness_minus,
                       n_outliers_in_QTL_uncorrected_plus=n_outliers_in_QTL_uncorrected_plus,n_outliers_in_QTL_corrected_plus=n_outliers_in_QTL_corrected_plus,n_outliers_in_QTL_missingness_plus=n_outliers_in_QTL_missingness_plus)
  outliers_in_windows_summary_split_reduced<-bind_rows(outliers_in_windows_summary_split_reduced,dataline)
}

outliers_in_windows_summary_split_reduced<-outliers_in_windows_summary_split_reduced %>%
  mutate(keygroup=factor(keygroup,levels=key_groups_customorder2_withtrueallos_reduced),
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
outliers_in_windows_summary_split_longer_reduced<-bind_rows(
  outliers_in_windows_summary_split_reduced %>% select(keygroup,n_total_windows,n_QTLwindows_plus,ends_with("_uncorrected_plus")) %>%
    rename_with(function(x){str_remove(x,"_uncorrected_plus")}) %>% rename("n_QTLwindows"="n_QTLwindows_plus") %>%
    mutate("dataset"="uncorrected","QTLset"="plus"),
  outliers_in_windows_summary_split_reduced %>% select(keygroup,n_total_windows,n_QTLwindows_minus,ends_with("_uncorrected_minus")) %>%
    rename_with(function(x){str_remove(x,"_uncorrected_minus")}) %>% rename("n_QTLwindows"="n_QTLwindows_minus") %>%
    mutate("dataset"="uncorrected","QTLset"="minus"),
  outliers_in_windows_summary_split_reduced %>% select(keygroup,n_total_windows,n_QTLwindows_plus,ends_with("_corrected_plus")) %>%
    rename_with(function(x){str_remove(x,"_corrected_plus")}) %>% rename("n_QTLwindows"="n_QTLwindows_plus") %>%
    mutate("dataset"="corrected","QTLset"="plus"),
  outliers_in_windows_summary_split_reduced %>% select(keygroup,n_total_windows,n_QTLwindows_minus,ends_with("_corrected_minus")) %>%
    rename_with(function(x){str_remove(x,"_corrected_minus")}) %>% rename("n_QTLwindows"="n_QTLwindows_minus") %>%
    mutate("dataset"="corrected","QTLset"="minus"),
  outliers_in_windows_summary_split_reduced %>% select(keygroup,n_total_windows,n_QTLwindows_plus,ends_with("_missingness_plus")) %>%
    rename_with(function(x){str_remove(x,"_missingness_plus")}) %>% rename("n_QTLwindows"="n_QTLwindows_plus") %>%
    mutate("dataset"="missingness","QTLset"="plus"),
  outliers_in_windows_summary_split_reduced %>% select(keygroup,n_total_windows,n_QTLwindows_minus,ends_with("_missingness_minus")) %>%
    rename_with(function(x){str_remove(x,"_missingness_minus")}) %>% rename("n_QTLwindows"="n_QTLwindows_minus") %>%
    mutate("dataset"="missingness","QTLset"="minus")) %>% 
  arrange(QTLset,desc(dataset),keygroup)

###t-tests on individual windows
ttests_Nplus_reduced<-list()
ttests_Nplus_withGOIs_reduced<-list()
ttests_Nminus_reduced<-list()
ttests_Nminus_withGOIs_reduced<-list()
ttest_results_reduced<-data.frame()
ttest_results_withGOIs_reduced<-data.frame()
for (i in focal_nine_reduced$keygroup) {
  ttests_Nplus_reduced[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_reduced[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_withGOIs_reduced[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_withGOIs_reduced[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttest_results_reduced<-bind_rows(ttest_results_reduced,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_reduced[[i]]$statistic,
    df=ttests_Nplus_reduced[[i]]$parameter,
    p.value=ttests_Nplus_reduced[[i]]$p.value,
    mean=ttests_Nplus_reduced[[i]]$estimate,
    stderr=ttests_Nplus_reduced[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_reduced[[i]]$statistic,
    df=ttests_Nminus_reduced[[i]]$parameter,
    p.value=ttests_Nminus_reduced[[i]]$p.value,
    mean=ttests_Nminus_reduced[[i]]$estimate,
    stderr=ttests_Nminus_reduced[[i]]$stderr
  )) 
  ttest_results_withGOIs_reduced<-bind_rows(ttest_results_withGOIs_reduced,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_withGOIs_reduced[[i]]$statistic,
    df=ttests_Nplus_withGOIs_reduced[[i]]$parameter,
    p.value=ttests_Nplus_withGOIs_reduced[[i]]$p.value,
    mean=ttests_Nplus_withGOIs_reduced[[i]]$estimate,
    stderr=ttests_Nplus_withGOIs_reduced[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_withGOIs_reduced[[i]]$statistic,
    df=ttests_Nminus_withGOIs_reduced[[i]]$parameter,
    p.value=ttests_Nminus_withGOIs_reduced[[i]]$p.value,
    mean=ttests_Nminus_withGOIs_reduced[[i]]$estimate,
    stderr=ttests_Nminus_withGOIs_reduced[[i]]$stderr
  )) 
}
label_significance<-function(p) {
  return(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","n.s.")))))
}

ttest_results_reduced<-ttest_results_reduced %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_withGOIs_reduced<-ttest_results_withGOIs_reduced %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))



###t-tests on region means
ttests_Nplus_regions_reduced<-list()
ttests_Nplus_regions_withGOIs_reduced<-list()
ttests_Nminus_regions_reduced<-list()
ttests_Nminus_regions_withGOIs_reduced<-list()
ttest_results_regions_reduced<-data.frame()
ttest_results_regions_withGOIs_reduced<-data.frame()
for (i in focal_nine_reduced$keygroup) {
  ttests_Nplus_regions_reduced[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_regions_reduced[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttests_Nplus_regions_withGOIs_reduced[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_regions_withGOIs_reduced[[i]]<-Mantel_regions_withGOIs_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% 
    group_by(Name) %>% summarize(delta_predict=mean(delta_predict)) %>% 
    pull(delta_predict) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttest_results_regions_reduced<-bind_rows(ttest_results_regions_reduced,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_regions_reduced[[i]]$statistic,
    df=ttests_Nplus_regions_reduced[[i]]$parameter,
    p.value=ttests_Nplus_regions_reduced[[i]]$p.value,
    mean=ttests_Nplus_regions_reduced[[i]]$estimate,
    stderr=ttests_Nplus_regions_reduced[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_regions_reduced[[i]]$statistic,
    df=ttests_Nminus_regions_reduced[[i]]$parameter,
    p.value=ttests_Nminus_regions_reduced[[i]]$p.value,
    mean=ttests_Nminus_regions_reduced[[i]]$estimate,
    stderr=ttests_Nminus_regions_reduced[[i]]$stderr
  )) 
  ttest_results_regions_withGOIs_reduced<-bind_rows(ttest_results_regions_withGOIs_reduced,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_regions_withGOIs_reduced[[i]]$statistic,
    df=ttests_Nplus_regions_withGOIs_reduced[[i]]$parameter,
    p.value=ttests_Nplus_regions_withGOIs_reduced[[i]]$p.value,
    mean=ttests_Nplus_regions_withGOIs_reduced[[i]]$estimate,
    stderr=ttests_Nplus_regions_withGOIs_reduced[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_regions_withGOIs_reduced[[i]]$statistic,
    df=ttests_Nminus_regions_withGOIs_reduced[[i]]$parameter,
    p.value=ttests_Nminus_regions_withGOIs_reduced[[i]]$p.value,
    mean=ttests_Nminus_regions_withGOIs_reduced[[i]]$estimate,
    stderr=ttests_Nminus_regions_withGOIs_reduced[[i]]$stderr
  )) 
}
ttest_results_regions_reduced<-ttest_results_regions_reduced %>% mutate(p.bonf.corrected=p.value*18,
                                        significance_uncorrected=label_significance(p.value),
                                        significance_corrected=label_significance(p.bonf.corrected))
ttest_results_regions_withGOIs_reduced<-ttest_results_regions_withGOIs_reduced %>% mutate(p.bonf.corrected=p.value*18,
                                                        significance_uncorrected=label_significance(p.value),
                                                        significance_corrected=label_significance(p.bonf.corrected))

outlier_stretches_reduced=data.frame()
for (i in key_groups_customorder_reduced$keygroup) {
  x=focal_outliers_95_corrected_reduced[[i]] %>% 
    mutate(isbreak=(is.na(lag(chrom)) | chrom!=lag(chrom) | windowend!=lag(windowend)+50000)) %>% 
    group_by(chrom) %>% summarize(nstretches=sum(isbreak)) %>% 
    summarize(nchroms=n(),nstretches=sum(nstretches)) %>% 
    mutate(keygroup=i)
  outlier_stretches_reduced=bind_rows(outlier_stretches_reduced,x)
}

