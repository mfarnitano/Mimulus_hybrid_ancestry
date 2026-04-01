
################################################################################
#####    Code to produce supplemental figures and tables                  ######
################################################################################

thisfolder<-"./"

####################
#### TABLE S1
####################
table_S1_nofilt<-prefilter_list %>% mutate(freq_group=ifelse(nasutus_proportion_windows<0.15,"guttatus",
                                                             ifelse(nasutus_proportion_windows>0.85,"nasutus",
                                                                    "hybrid"))) %>% 
  group_by(group,freq_group) %>% summarize(n=n()) %>% pivot_wider(names_from=freq_group,values_from=n)
write.table(table_S1_nofilt,file.path(thisfolder,"Metadata/initial_individuals_summary.txt"),quote = F,row.names = F,col.names=T,sep="\t")

table_S1_filt<-final_individual_list %>% group_by(group,freq_group) %>% summarize(n=n()) %>% pivot_wider(names_from=freq_group,values_from=n)
write.table(table_S1_filt,file.path(thisfolder,"Metadata/final_individuals_summary.txt"),quote = F,row.names = F,col.names=T,sep="\t")

####################
#### TABLE S2
####################
panelinfo<-individuals_sitefiltered %>% filter(str_detect(group,"panel")) %>%
  mutate(Line=str_remove(sampleID,".sampled")) %>% select(Line,nasutus_proportion_windows,group) %>% arrange(group,nasutus_proportion_windows,Line)
write.table(panelinfo,file.path(thisfolder,"Metadata/gutnas_panel_info.txt"),sep="\t",row.names=F,col.names=T,quote=F)
southern_allopatric_info<-individuals_sitefiltered %>% filter(group=="Southern_allopatric") %>% 
  mutate(Pop=str_sub(sampleID,1,3)) %>% arrange(Pop,nasutus_proportion_windows)
southern_allopatric_info %>% group_by(Pop) %>% summarize(n_total=n(),n_above_1=sum(nasutus_proportion_windows>=0.01,na.rm=T))


################################################################################
#####    Figure S1                                                        ######
################################################################################


###distribution of rsq_PC1 and Zhet values
PCA_fit_stats_allwindows_CAC<-read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch3_genomescans/Outputs_July2025/CAC.PCA_fit_stats.txt",header=T)
PCA_fit_stats_allwindows_CAC_frequencies<-all_groups_window_freqs %>% filter(str_starts(keygroup,"CAC")) %>%
  mutate(region=paste0(chrom,".",windowstart,".",windowend)) %>%
  left_join(PCA_fit_stats_allwindows_CAC,by="region",relationship="many-to-one")
outliers_cutoffs_CAC<-PCA_fit_stats_allwindows_CAC_frequencies %>% group_by(keygroup) %>% summarize(cutoff95=quantile(nasutus_frequency,0.95),
                                                                                                    cutoff05=quantile(nasutus_frequency,0.05)) %>%
  right_join(PCA_fit_stats_allwindows_CAC_frequencies,by="keygroup") %>%
  mutate(outlier_type=ifelse(nasutus_frequency>cutoff95,"top 5%",ifelse(nasutus_frequency<cutoff05,"bottom 5%","other")))
PCA_fit_stats_allwindows_southern<-read.table("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch3_genomescans/Outputs_July2025/southern_core.PCA_fit_stats.txt",header=T)
PCA_fit_stats_allwindows_southern_frequencies<-all_groups_window_freqs %>% filter(keygroup %in% c("Southern_sympatric_guttatus_174","Southern_allopatric_guttatus_61","Southern_sympatric_nasutus_75")) %>%
  mutate(region=paste0(chrom,".",windowstart,".",windowend)) %>%
  left_join(PCA_fit_stats_allwindows_southern,by="region",relationship="many-to-one")
outliers_cutoffs_southern<-PCA_fit_stats_allwindows_southern_frequencies %>% group_by(keygroup) %>% summarize(cutoff95=quantile(nasutus_frequency,0.95),
                                                                                                              cutoff05=quantile(nasutus_frequency,0.05)) %>%
  right_join(PCA_fit_stats_allwindows_southern_frequencies,by="keygroup") %>%
  mutate(outlier_type=ifelse(nasutus_frequency>cutoff95,"top 5%",ifelse(nasutus_frequency<cutoff05,"bottom 5%","other")))
rsq_vs_ancestry_both<-outliers_cutoffs_CAC %>% bind_rows(outliers_cutoffs_southern) %>% 
  ggplot(aes(x=nasutus_frequency,y=modelfit.r.squared)) + 
  facet_wrap(~keygroup,labeller = as_labeller(c("CAC_guttatus_137" = "N-CAC-sym guttatus",
                                                "CAC_hybrid_183" = "N-CAC-sym hybrid",
                                                "CAC_nasutus_58" = "N-CAC-sym nasutus",
                                                "Southern_allopatric_guttatus_61" = "S-Foothills-allo guttatus",
                                                "Southern_sympatric_guttatus_174" = "S-Foothills-sym guttatus",
                                                "Southern_sympatric_nasutus_75" = "S-Foothills-sym nasutus"))) + 
  geom_point() + theme_bw() + 
  geom_hline(yintercept=0.9,linetype="dotted",color="red") + 
  ylab("Rsquared (PC1 v. ancestry)") + xlab("Ancestry frequency")
rsq_histogram<-outliers_cutoffs_CAC %>% bind_rows(outliers_cutoffs_southern) %>% 
  filter(keygroup %in% c("CAC_guttatus_137","Southern_allopatric_guttatus_61")) %>%
  ggplot(aes(y=modelfit.r.squared)) + 
  facet_wrap(~keygroup,nrow=2,labeller = as_labeller(c("CAC_guttatus_137" = "Northern", "Southern_allopatric_guttatus_61" = "Southern"))) + 
  geom_histogram(binwidth=0.05,boundary=0.9) + theme_bw() + 
  geom_hline(yintercept=0.9,linetype="dotted",color="red") + 
  ylab("Rsquared (PC1 v. ancestry)") + xlab("Window count")
hetscore_vs_ancestry_both<-outliers_cutoffs_CAC %>% bind_rows(outliers_cutoffs_southern) %>% 
  ggplot(aes(x=nasutus_frequency,y=het.scaled.zscore)) + 
  facet_wrap(~keygroup,labeller = as_labeller(c("CAC_guttatus_137" = "N-CAC-sym guttatus",
                                                "CAC_hybrid_183" = "N-CAC-sym hybrid",
                                                "CAC_nasutus_58" = "N-CAC-sym nasutus",
                                                "Southern_allopatric_guttatus_61" = "S-Foothills-allo guttatus",
                                                "Southern_sympatric_guttatus_174" = "S-Foothills-sym guttatus",
                                                "Southern_sympatric_nasutus_75" = "S-Foothills-sym nasutus"))) + 
  geom_point() + theme_bw() + 
  geom_hline(yintercept=2,linetype="dotted",color="red") + 
  geom_hline(yintercept=-2,linetype="dotted",color="red") + 
  ylab("Z_het") + xlab("Ancestry frequency")
hetscore_histogram<-outliers_cutoffs_CAC %>% bind_rows(outliers_cutoffs_southern) %>% 
  filter(keygroup %in% c("CAC_guttatus_137","Southern_allopatric_guttatus_61")) %>%
  ggplot(aes(y=het.scaled.zscore)) + 
  facet_wrap(~keygroup,nrow=2,labeller = as_labeller(c("CAC_guttatus_137" = "Northern", "Southern_allopatric_guttatus_61" = "Southern"))) +
  geom_histogram(binwidth=1,boundary=0) + theme_bw() + 
  geom_hline(yintercept=2,linetype="dotted",color="red") + 
  geom_hline(yintercept=-2,linetype="dotted",color="red") + 
  ylab("Z_het") + xlab("Window count")

PCinvestigation_summary<-plot_grid(rsq_vs_ancestry_both,rsq_histogram,hetscore_vs_ancestry_both,hetscore_histogram,nrow=2,rel_widths = c(3,1))
ggsave(file.path(thisfolder,"Paper_figs/LocalPCA_investigations_summary.pdf"),PCinvestigation_summary,device="pdf",width=12,height=8,units="in",dpi=600)

################################################################################
#####    Table S3                                                        ######
################################################################################


###correlations with recombination excluding centromeric regions

###do it all with higher-resolution recombination rate
genes_per_cM<-correlation_matrix %>% mutate(genes_percM=ngenes/(M_per_bp_50kb*5000000)) %>% filter(!is.na(genes_percM),!is.infinite((genes_percM)))
rsq_genes_per_cM<-round(summary(lm(CAC_hybrid_183~genes_percM,data=genes_per_cM))$r.squared,3)
raw_genespercM<-ggplot(genes_per_cM,aes(x=genes_percM,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(color="blue") + geom_smooth(method="lm",color="red") + 
  annotate(geom="text",label=paste0("R^2=",rsq_genes_per_cM),color="red",x=200,y=0.8) + 
  theme_bw() + xlab("gene density per cM") + ylab("CAC hybrid ancestry frequency") + 
  ggtitle("All windows")
genes_per_cM_newfilter<-correlation_matrix %>% filter(M_per_bp_50kb>=1e-8) %>% 
  mutate(genes_percM=ngenes/(M_per_bp_50kb*5000000)) %>% filter(!is.na(genes_percM),!is.infinite((genes_percM)))
rsq_genes_per_cM_newfilter<-round(summary(lm(CAC_hybrid_183~genes_percM,data=genes_per_cM_newfilter))$r.squared,3)
newfilter_genespercM<-ggplot(genes_per_cM_newfilter,aes(x=genes_percM,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(color="blue") + geom_smooth(method="lm",color="red") + 
  annotate(geom="text",label=paste0("R^2=",rsq_genes_per_cM_newfilter),color="red",x=150,y=0.8) + 
  theme_bw() + xlab("gene density per cM") + ylab("CAC hybrid ancestry frequency") + 
  ggtitle("Exclude <1cM/Mb")
rsq_genes_per_cM_newfilter_nooutliers<-round(summary(lm(CAC_hybrid_183~genes_percM,data=genes_per_cM_newfilter %>% filter(genes_percM<=100)))$r.squared,3)
newfilter_genespercM_nooutliers<-genes_per_cM_newfilter %>% filter(genes_percM<=100) %>% 
  ggplot(aes(x=genes_percM,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(color="blue") + geom_smooth(method="lm",color="red") + 
  annotate(geom="text",label=paste0("R^2=",rsq_genes_per_cM_newfilter_nooutliers),color="red",x=80,y=0.8) + 
  theme_bw() + xlab("gene density per cM") + ylab("CAC hybrid ancestry frequency") + 
  ggtitle("Exclude <1cM/Mb, exclude density>100")

genedensity_plots<-plot_grid(raw_genespercM,newfilter_genespercM,newfilter_genespercM_nooutliers,nrow=3)
ggsave(file.path(thisfolder,"genedensity_correlations.png"),genedensity_plots,device = "png",width=4,height=12,units = "in",dpi=600)

nrow(correlation_matrix) #3139
nrow(correlation_matrix %>% filter(M_per_bp_50kb==0)) ###585
nrow(correlation_matrix %>% filter(M_per_bp_50kb<1e-8)) ###777
nrow(correlation_matrix %>% filter(relative_chr_pos<0.25)) ###281
nrow(correlation_matrix %>% filter(M_per_bp_50kb<1e-8,relative_chr_pos<0.25)) ###177
correlation_matrix %>% filter(M_per_bp_50kb>=1e-8,relative_chr_pos<0.10) %>% group_by(chrom) %>% summarize(n())
updated_filter_long<-correlation_matrix %>% filter(M_per_bp_50kb>=1e-8,relative_chr_pos>=0.25) %>% 
  mutate(genes_percM=ngenes/(M_per_bp_50kb*5000000)) %>% 
  pivot_longer(cols=CAC_guttatus_137:Southern_sympatric_nasutus_75,names_to = "keygroup",values_to="nasutus_ancestry") %>%
  pivot_longer(cols=c(relative_chr_pos,ngenes,M_per_bp_1Mb,M_per_bp_50kb,genes_percM),names_to="statistic",values_to="value")
unfiltered_long<-correlation_matrix %>%
  mutate(genes_percM=ngenes/(M_per_bp_50kb*5000000)) %>%
  pivot_longer(cols=CAC_guttatus_137:Southern_sympatric_nasutus_75,names_to = "keygroup",values_to="nasutus_ancestry") %>%
  pivot_longer(cols=c(relative_chr_pos,ngenes,M_per_bp_1Mb,M_per_bp_50kb,genes_percM),names_to="statistic",values_to="value")
all_rsqs_updated<-data.frame()
for (k in unique(updated_filter_long$keygroup)) {
  for (s in unique(updated_filter_long$statistic)) {
    subdata<-updated_filter_long %>% filter(keygroup==k,statistic==s,!is.infinite(value))
    rsq<-round(summary(lm(nasutus_ancestry~value,data=subdata))$adj.r.squared,3)
    pval<-round(summary(lm(nasutus_ancestry~value,data=subdata))$coefficients[2,4],5)
    all_rsqs_updated<-all_rsqs_updated %>% bind_rows(data.frame(keygroup=k,statistic=s,rsq=rsq,pval=pval,withfilter=T))
    subdata<-unfiltered_long %>% filter(keygroup==k,statistic==s,!is.infinite(value))
    rsq<-round(summary(lm(nasutus_ancestry~value,data=subdata))$adj.r.squared,3)
    pval<-round(summary(lm(nasutus_ancestry~value,data=subdata))$coefficients[2,4],5)
    all_rsqs_updated<-all_rsqs_updated %>% bind_rows(data.frame(keygroup=k,statistic=s,rsq=rsq,pval=pval,withfilter=F))
  }
  
}
all_rsqs_updated_wide<-all_rsqs_updated %>% mutate(bonf=pval*9,sig=ifelse(bonf<0.001,"***",
                                                                          ifelse(bonf<0.01,"**",
                                                                                 ifelse(bonf<0.05,"*",
                                                                                        ifelse(bonf<0.1,".","")))),
                                                   rsq_sig=paste0(rsq,sig)) %>% 
  select(keygroup,statistic,withfilter,rsq_sig) %>%
  pivot_wider(names_from="statistic",values_from="rsq_sig") %>%
  arrange(withfilter,keygroup)
write.table(all_rsqs_updated_wide,file.path(thisfolder,"Paper_tables/genomic_feature_correlations_before_and_after.txt"),row.names=F,col.names=T,sep = "\t",quote=F)

bothfilter<-bind_rows(updated_filter_long %>% mutate(withfilter="Exclude pericentromic"),
                      unfiltered_long %>% mutate(withfilter="All data")) %>%
  mutate(statistic=factor(case_match(statistic,
                                     "relative_chr_pos" ~ "Relative position",
                                     "ngenes" ~ "genes per 50kb",
                                     "M_per_bp_1Mb" ~ "recombination rate (1Mb)",
                                     "M_per_bp_50kb" ~ "recombination rate (50kb)",
                                     "genes_percM" ~ "genes per cM")))
ggplot(bothfilter,aes(x=value,y=nasutus_ancestry,color=keygroup)) + 
  geom_smooth(method="lm") + 
  facet_grid(rows=vars(withfilter),cols=vars(statistic),scales="free")

################################################################################
#####    Figure S2                                                        ######
################################################################################

group_histograms_facetted_ready_reduced<-final_individual_list_reduced %>% 
  mutate(title=case_match(group,
                          "CAC" ~ "Northern - CAC - sympatric",
                          "LM" ~ "Northern - LM - sympatric",
                          "Southern_allopatric" ~ "Southern - Foothills - allopatric",
                          "Southern_sympatric" ~ "Southern - Foothills - sympatric",
                          "Southern_highelev" ~ "Southern - Montane - sympatric"))
group_histograms_facetted_reduced<-group_histograms_facetted_ready_reduced %>%
  ggplot(aes(x=nasutus_proportion_windows)) + theme_bw() + 
  facet_wrap(~title,ncol=1,scales="free_y") + 
  geom_histogram(data=subset(group_histograms_facetted_ready_reduced,freq_group=="guttatus"),binwidth=0.05,boundary=0,fill=colors_sp[1]) + 
  geom_histogram(data=subset(group_histograms_facetted_ready_reduced,freq_group=="hybrid"),binwidth=0.05,boundary=0,fill=colors_sp[2]) + 
  geom_histogram(data=subset(group_histograms_facetted_ready_reduced,freq_group=="nasutus"),binwidth=0.05,boundary=0,fill=colors_sp[3]) + 
  scale_x_continuous(limits=c(0,1),expand=c(0,0),name="Hybrid index") + 
  scale_y_continuous(name="Number of samples") + 
  theme(text=element_text(size=12))

ggsave(file.path(thisfolder,"Paper_figs/FigS2_histograms_reduced.pdf"),group_histograms_facetted_reduced,device="pdf",width=5,height=15,units="in",dpi=600)


####Table of correlations
write.table(cor(correlation_matrix_withsquares_reduced),file = file.path(thisfolder,"correlations_table_reduced.txt"),quote=F,sep="\t")
write.table(model_test_variances_reduced %>% filter(model=="newmodel"),file=file.path(thisfolder,"model_variances_AICs_table_reduced.txt"),quote=F,sep="\t")

################################################################################
#####    Figure S3                                                       ######
################################################################################


###Figure S3: (analog to Fig 3) frequency histograms
colorkey=colors_groups[c(1,1,2,3,3,3,4,5,5,1,4)]
neworder_reduced<-c("CAC_guttatus_156","CAC_hybrid_164","CAC_nasutus_58",
                    "LM_hybrid_56",
                    "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174","Southern_sympatric_nasutus_75",
                    "Southern_highelev_guttatus_12","Southern_highelev_hybrid_15")
hist_grouplabels_reduced<-c(
  "N-CAC-sym guttatus, n=156",
  "N-CAC-sym hybrid, n=164",
  "N-CAC-sym nasutus, n=58",
  "N-LM-sym hybrid, n=56",
  "S-Foothills-allo guttatus, n=61",
  "S-Foothills-sym guttatus, n=174",
  "S-Foothills-sym nasutus, n=75",
  "S-Montane-sym guttatus, n=12",
  "S-Montane-sym hybrid, n=15"
)
freq_hists_all_ready_reduced<-all_groups_window_freqs_reduced %>% ungroup() %>% filter(!is.na(nasutus_frequency)) %>%
  mutate(grouplabel=factor(keygroup,levels=neworder_reduced)) %>% filter(!is.na(grouplabel)) %>%
  mutate(grouplabel=dplyr::recode(grouplabel, !!!setNames(hist_grouplabels_reduced,neworder_reduced))) %>%
  mutate(freq_bin=ceiling(nasutus_frequency*200)/200) %>% group_by(grouplabel,freq_bin) %>%
  summarize(n_windows=n())
maxmins_reduced<-freq_hists_all_ready_reduced %>% group_by(grouplabel) %>% summarize(max=max(freq_bin),min=min(freq_bin)) %>%
  pivot_longer(cols=c(max,min),names_to = "which",values_to="arrow_x") %>% 
  mutate(xadjust=ifelse(which=="max",arrow_x+0.02,arrow_x-0.02)) %>%
  left_join(freq_hists_all_ready_reduced %>% group_by(grouplabel) %>% summarize(maxy=max(n_windows)),relationship="many-to-one")
freq_hists_all_reduced<-freq_hists_all_ready_reduced %>% ggplot(aes(x=freq_bin,y=n_windows,fill=grouplabel)) + 
  facet_wrap(~grouplabel,scale="free_y",ncol=1) + 
  geom_col() + 
  theme_bw() + ylab("Window count") + 
  scale_fill_discrete(palette=colors_groups[c(1,1,1,2,3,4,4,5,5)],guide="none") + 
  geom_segment(data=maxmins_reduced,aes(x=xadjust,xend=arrow_x,y=maxy/2,yend=0),arrow=arrow(type="closed",length = unit(0.05, "in")),linetype="dashed") + 
  scale_x_continuous(name="Ancestry frequency",limits=c(-0.05,1.05))

###Fig S3 black-and-white version
freq_hists_all_bw_reduced<-freq_hists_all_ready_reduced %>% ggplot(aes(x=freq_bin,y=n_windows)) + 
  facet_wrap(~grouplabel,scale="free_y",ncol=1) + 
  geom_rect(data=maxmins_reduced %>% filter(which=="min"),aes(xmin=0,xmax=arrow_x,ymin=0,ymax=maxy,x=NULL,y=NULL),fill="lightgrey") + 
  geom_rect(data=maxmins_reduced %>% filter(which=="max"),aes(xmin=arrow_x,xmax=1,ymin=0,ymax=maxy,x=NULL,y=NULL),fill="lightgrey") + 
  geom_col(fill="black") + 
  theme_bw() + ylab("Window count") + 
  scale_x_continuous(name="Ancestry frequency",limits=c(-0.05,1.05)) + 
  theme(text=element_text(size=12))
ggsave(file.path(thisfolder,"Paper_figs/histograms_allele_frequencies_bygroup_facetted_reduced_bw.pdf"),freq_hists_all_bw_reduced,device="pdf",width = 4,height=8,units="in",dpi = 600)  


################################################################################
#####    Figure S4                                                        ######
################################################################################


###Figure S4: Analysis of variance (compare to Figure 3B) for reduced dataset

rect_colors<-brewer.pal(5,name="Set1")
c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
boxes=data.frame(xmin=c(0.6,0.6,0.6,0.6,2.6,3.6,3.6,3.6,5.6,5.6,7.6,7.6),
                 xmax=c(2.4,2.4,3.4,7.4,3.4,7.4,5.4,5.4,7.4,7.4,9.4,9.4),
                 ymin=c(0.6,2.6,3.6,7.6,0.6,0.6,3.6,5.6,3.6,5.6,0.6,7.6),
                 ymax=c(2.4,3.4,7.4,9.4,2.4,3.4,5.4,7.4,5.4,7.4,7.4,9.4),
                 colorgroup=factor(c(1,2,4,5,2,4,2,3,3,1,5,4)))
boxes2<-boxes %>% mutate(xmin=xmin-0.05,xmax=xmax+0.05,ymin=ymin-0.05,ymax=ymax+0.05) 
boxes_shift<-boxes %>% mutate(xmin=xmin+2,xmax=xmax+2)
boxes2_shift<-boxes2 %>% mutate(xmin=xmin+2,xmax=xmax+2)

Figure3_combined_ready_reduced<-all_correlations_withcorrections_reduced %>% 
  mutate(i=factor(group1,levels=key_groups_customorder2_withtrueallos_reduced),
         j=factor(group2,levels=c("genome_structure","missingness",key_groups_customorder2_withtrueallos_reduced)),
         var_explained_raw=uncorrected^2,
         var_explained_prettyval_raw=format(round(var_explained_raw*100,1),digits=1)) %>%
  select(i,j,var_explained_raw,var_explained_prettyval_raw) %>%
  full_join(anova_results_ready_reduced,by=c("i","j")) %>%
  filter(!(i %in% c("Southern_true_allopatric","Southern_pseudo_allopatric")),
         !(j %in% c("Southern_true_allopatric","Southern_pseudo_allopatric")),
         !(as.character(i)==as.character(j))) %>%
  mutate(whichval=ifelse(as.integer(i)>(as.integer(j)-2),var_explained_testgroup,var_explained_raw),
         whichprettyval=ifelse(as.integer(i)>(as.integer(j)-2),prettyval,var_explained_prettyval_raw))
tile_varexplained_merged_reduced<-Figure3_combined_ready_reduced %>% 
  ggplot(aes(x=j,y=i,label=whichprettyval,fill=whichval)) + 
  scale_fill_gradient(low="white",high="darkgreen",guide="none") + 
  geom_tile() + geom_text() + theme_bw() + 
  geom_vline(xintercept=2.5,linetype="solid") + 
  geom_vline(xintercept=c(5.5,7.5,9.5),linetype="dashed") + 
  geom_hline(yintercept=c(3.5,5.5,7.5),linetype="dashed") + 
  geom_abline(slope=1,intercept=-2,linetype="solid",linewidth=2) + 
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x=element_text(angle=90))
tile_varexplained_merged_nonames_reduced<-tile_varexplained_merged_reduced + 
  scale_x_discrete(labels=NULL) + scale_y_discrete(labels=NULL) + 
  ylab(NULL) + xlab(NULL)
tile_varexplained_merged_nonames_box_reduced<-tile_varexplained_merged_nonames_reduced + 
  geom_rect(data=boxes2_shift,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=colorgroup,fill=NULL,label=NULL),fill="transparent",linewidth=1,linetype="dashed") + 
  scale_color_manual(values=rect_colors,name=NULL,
                     labels=c("Within-site","Nearby-site","Within-region","Between-region","Between-species")) + 
  coord_equal(ratio=1)
tile_varexplained_merged_nonames_box_reduced
ggsave(file.path(thisfolder,"Paper_figs/tile_varexplained_anovas_plusraw_nonames_boxes_reduced.pdf"),tile_varexplained_merged_nonames_box_reduced,device="pdf",width = 7,height=5,units="in",dpi = 600)  


################################################################################
#####    Figure S5                                                        ######
################################################################################


### Figure S5: outlier overlaps
###outlier overlap comparing raw to missingness (RAW now in bottom right triangle!)
###tile of outlier overlaps with missingness correction
tile_overlaps_missingness_boxes_reduced<-count_overlaps_reduced %>% mutate(i=factor(i,levels=key_groups_customorder2_reduced$keygroup),
                                                                           j=factor(j,levels=key_groups_customorder2_reduced$keygroup)) %>%
  filter(!is.na(i),!is.na(j)) %>%
  mutate(observed_overlap=ifelse(as.integer(j)>=as.integer(i),observed_overlap_missingness,observed_overlap_uncorrected),
         pval_permutes=ifelse(as.integer(j)>=as.integer(i),pval_permutes_missingness,pval_permutes_uncorrected),
         pval_bonf_corrected=pval_permutes*110,significance=ifelse(i=="CAC_nasutus_58","",
                                                                   ifelse(pval_bonf_corrected<0.001,"***",
                                                                          ifelse(pval_bonf_corrected<0.01,"**",
                                                                                 ifelse(pval_bonf_corrected<0.05,"*",
                                                                                        ifelse(pval_bonf_corrected<0.1,".",""))))),
         tile_label=ifelse(as.integer(i)==as.integer(j),
                           paste0(observed_overlap_missingness,"\n",observed_overlap_uncorrected),
                           paste0(observed_overlap,"\n",significance))) %>%
  ggplot(aes(x=i,y=j,fill=observed_overlap,label=tile_label)) + 
  geom_tile() + geom_text() + 
  geom_rect(data=boxes2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=colorgroup,x=NULL,y=NULL,fill=NULL,label=NULL),fill="transparent",linewidth=1,linetype="dashed") + 
  scale_color_manual(values=rect_colors,name=NULL,
                     labels=c("Within-site","Nearby-site","Within-region","Between-region","Between-species")) + 
  scale_fill_gradient(low="white",high="darkgreen",guide="none") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) + 
  ylab(NULL) + xlab(NULL) + 
  geom_hline(yintercept=c(3.5,5.5,7.5)) + geom_vline(xintercept=c(3.5,5.5,7.5)) + 
  coord_equal(ratio=1)
tile_overlaps_missingness_nonames_boxes_reduced<-tile_overlaps_missingness_boxes_reduced + 
  scale_x_discrete(labels=NULL) + scale_y_discrete(labels=NULL)
ggsave(file.path(thisfolder,"Paper_figs/deltapredict_outliers_overlap_matrix_missingness_boxes_reduced.pdf"),tile_overlaps_missingness_boxes_reduced,device="pdf",width = 8,height=7,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/deltapredict_outliers_overlap_matrix_missingness_nonames_boxes_reduced.pdf"),tile_overlaps_missingness_nonames_boxes_reduced,device="pdf",width = 6,height=5,units="in",dpi = 600)  


################################################################################
#####    Figure S6                                                        ######
################################################################################


###Figure S6: outliers vs. QTLs -- violin plots (compare to main Fig 5)

ttests_Nplus_missingness_reduced<-list()
ttests_Nminus_missingness_reduced<-list()
ttest_results_missingness_reduced<-data.frame()
for (i in focal_nine_reduced$keygroup) {
  ttests_Nplus_missingness_reduced[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N+") %>% pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "greater")
  ttests_Nminus_missingness_reduced[[i]]<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>%
    filter(keygroup==i,Predicted_effect=="N-") %>% pull(delta_predict_missingness) %>%
    t.test(conf.level=0.95,alternative = "less")
  ttest_results_missingness_reduced<-bind_rows(ttest_results_missingness_reduced,data.frame(
    keygroup=i,Predicted_effect="N+",
    t=ttests_Nplus_missingness_reduced[[i]]$statistic,
    df=ttests_Nplus_missingness_reduced[[i]]$parameter,
    p.value=ttests_Nplus_missingness_reduced[[i]]$p.value,
    mean=ttests_Nplus_missingness_reduced[[i]]$estimate,
    stderr=ttests_Nplus_missingness_reduced[[i]]$stderr
  )) %>% bind_rows(data.frame(
    keygroup=i,Predicted_effect="N-",
    t=ttests_Nminus_missingness_reduced[[i]]$statistic,
    df=ttests_Nminus_missingness_reduced[[i]]$parameter,
    p.value=ttests_Nminus_missingness_reduced[[i]]$p.value,
    mean=ttests_Nminus_missingness_reduced[[i]]$estimate,
    stderr=ttests_Nminus_missingness_reduced[[i]]$stderr
  )) 
}
label_significance<-function(p) {
  return(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","")))))
}

ttest_results_missingness_reduced<-ttest_results_missingness_reduced %>% mutate(p.bonf=p.value*18,
                                                                                significance_uncorrected=label_significance(p.value),
                                                                                significance_corrected=label_significance(p.bonf))

write.table(ttest_results_reduced %>% arrange(p.value),file.path(thisfolder,"Mantel_deltapredict_ttests_windows_reduced.txt"),quote = F,sep = "\t",row.names=F)
write.table(ttest_results_withGOIs_reduced %>% arrange(p.value),file=file.path(thisfolder,"Mantel_deltapredict_ttests_windows_withGOIs_reduced.txt"),quote = F,sep = "\t",row.names=F)
write.table(ttest_results_regions_reduced %>% arrange(p.value),file=file.path(thisfolder,"Mantel_deltapredict_ttests_regions.txt_reduced"),quote = F,sep = "\t",row.names=F)
write.table(ttest_results_regions_withGOIs_reduced %>% arrange(p.value),file=file.path(thisfolder,"Mantel_deltapredict_ttests_regions_withGOIs_reduced.txt"),quote = F,sep = "\t",row.names=F)

Mantel_regions_deltapredict_windows_missingness_reduced<-Mantel_regions_nolarge_freqs_withpredictions_newmodel_reduced %>% 
  filter(keygroup %in% key_groups_customorder2_reduced$keygroup) %>%
  filter(!is.na(keygroup)) %>%
  ggplot(aes(x=keygroup,y=delta_predict_missingness)) + 
  facet_wrap(~Predicted_effect,nrow=1) + 
  geom_violin(scale="width",quantile.linetype="solid") + #geom_point(position = position_dodge(width=0.75)) + 
  stat_summary(fun=mean,geom="point",position = position_dodge(width=0.9),shape=4) + 
  geom_text(data=ttest_results_missingness_reduced,aes(x=keygroup,y=0.5,label=significance_corrected,color=NULL),color="black") + 
  #scale_color_manual(values=colors_groups[c(1,1,2,3,4,5,5,1,4)],guide=NULL) + 
  scale_x_discrete(name=NULL) + 
  geom_hline(yintercept=0,linetype="dashed") + ylab(NULL) + 
  theme_bw() + theme(axis.text.x=element_text(angle=90))
ggsave(file.path(thisfolder,"Paper_figs/violin_outliers_vs_QTL_missingness_reduced.pdf"),Mantel_regions_deltapredict_windows_missingness_reduced,device="pdf",width = 7,height=5,units="in",dpi = 600)  

Mantel_regions_deltapredict_windows_missingness_nonames_reduced<-Mantel_regions_deltapredict_windows_missingness_reduced + 
  scale_x_discrete(labels=NULL,name=NULL) + 
  geom_hline(yintercept=0,linetype="dashed") + ylab(NULL) + 
  theme_bw()
ggsave(file.path(thisfolder,"Paper_figs/violin_outliers_vs_QTL_missingness_reduced_nonames.pdf"),Mantel_regions_deltapredict_windows_missingness_nonames_reduced,device="pdf",width = 7,height=5,units="in",dpi = 600)  


################################################################################
#####    Table S4                                                         ######
################################################################################

###Breakdown of Southern allopatric populations

southern_pops_breakdown<-final_individual_list %>% filter(group %in% c("Southern_allopatric","Southern_sympatric","Southern_highelev"),
                                                          freq_group=="guttatus") %>%
  mutate(population=str_sub(sample,1,3)) %>% 
  mutate(population=ifelse(population=="Tn5","DPR",population)) %>% 
  group_by(group,population) %>%
  summarize(n_samples=n(),
            mean_proportion=mean(nasutus_proportion_windows),
            samples_above_1percent=sum(nasutus_proportion_windows>=0.01)) %>%
  arrange(desc(mean_proportion))

write.table(southern_pops_breakdown,file.path(thisfolder,"Paper_tables/southern_populations_breakdown.txt"),row.names = F,col.names=T,sep="\t",quote=F)

################################################################################
#####    Figure S7                                                        ######
################################################################################


###relationship between allopatric pops ancestry frequency and distance to sympatry
library(geosphere)
dist_ready<-DPR_mapdata %>% select(Code,Latitude,Longitude,cohort) %>% filter(cohort %in% c("sympatric","allopatric"))
get_distance<-function(lat1,lat2,long1,long2) {
  return(distm(c(long1,lat1),c(long2,lat2),fun = distHaversine))
}
dist_pairwise<-cross_join(dist_ready %>% rename(Code_1=Code,Latitude_1=Latitude,Longitude_1=Longitude,cohort_1=cohort),
                         dist_ready %>% rename(Code_2=Code,Latitude_2=Latitude,Longitude_2=Longitude,cohort_2=cohort)) %>%
  filter(Code_1!=Code_2) %>% rowwise() %>% mutate(distance=get_distance(Latitude_1,Latitude_2,Longitude_1,Longitude_2)) %>% 
  ungroup() %>% 
  select(Code_1,Code_2,cohort_1,cohort_2,distance) %>% filter(cohort_1=="allopatric",cohort_2=="sympatric") %>% 
  group_by(Code_1) %>% summarize(min_distance=min(distance)) %>% rename(population = Code_1) %>%
  inner_join(southern_pops_breakdown,by="population")

summary(lm(mean_proportion~min_distance,data=dist_pairwise)) ###rsquared=0.7714, adj.R2=0.7257, F=16.87 on 1 and 5 df, p=0.009287
allopatric_distance_to_sympatry<-dist_pairwise %>% ggplot(aes(x=min_distance,y=mean_proportion)) + 
  geom_point() + geom_smooth(method="lm") + theme_bw() + 
  ylab("Mean hybrid index") + xlab("Minimum distance to sympatric locale")
write.table(dist_pairwise %>% select(population,min_distance),file.path(thisfolder,"Paper_tables/allopatric_distance_to_sympatry.txt"),quote=F,row.names=F,col.names=T,sep="\t")

ggsave(file.path(thisfolder,"Paper_figs/allopatric_distance_to_sympatry.pdf"),allopatric_distance_to_sympatry,device="pdf",width = 4,height=3,units="in",dpi = 600)

################################################################################
#####    Figure S8                                                        ######
################################################################################


###correlation matrices with 'true-allopatric' split
boxes=data.frame(xmin=c(0.6,0.6,0.6,0.6,2.6,3.6,3.6,3.6,5.6,5.6,7.6,7.6),
                 xmax=c(2.4,2.4,3.4,7.4,3.4,7.4,5.4,5.4,7.4,7.4,9.4,9.4),
                 ymin=c(0.6,2.6,3.6,7.6,0.6,0.6,3.6,5.6,3.6,5.6,0.6,7.6),
                 ymax=c(2.4,3.4,7.4,9.4,2.4,3.4,5.4,7.4,5.4,7.4,7.4,9.4),
                 colorgroup=factor(c(1,2,4,5,2,4,2,3,3,1,5,4)))
boxes2<-boxes %>% mutate(xmin=xmin-0.05,xmax=xmax+0.05,ymin=ymin-0.05,ymax=ymax+0.05) 
boxes_shift<-boxes %>% mutate(xmin=xmin+2,xmax=xmax+2)
boxes_withtrueallos<-data.frame(xmin=c(0.6,0.6,0.6,0.6,2.6,3.6,3.6,3.6,7.6,7.6,9.6,9.6),
                                xmax=c(2.4,2.4,3.4,9.4,3.4,9.4,7.4,7.4,9.4,9.4,11.4,11.4),
                                ymin=c(0.6,2.6,3.6,9.6,0.6,0.6,3.6,7.6,3.6,7.6,0.6,9.6),
                                ymax=c(2.4,3.4,9.4,11.4,2.4,3.4,7.4,9.4,7.4,9.4,9.4,11.4),
                                colorgroup=factor(c(1,2,4,5,2,4,2,3,3,1,5,4)))
boxes_withtrueallos2<-boxes_withtrueallos %>% mutate(xmin=xmin-0.05,xmax=xmax+0.05,ymin=ymin-0.05,ymax=ymax+0.05) 
boxes_withtrueallos2_shift<-boxes_withtrueallos2 %>% mutate(xmin=xmin+2,xmax=xmax+2)

FigureS8_combined_ready<-all_correlations_withcorrections %>% 
  mutate(i=factor(group1,levels=key_groups_customorder2_withtrueallos),
         j=factor(group2,levels=c("genome_structure","missingness",key_groups_customorder2_withtrueallos)),
         var_explained_raw=uncorrected^2,
         var_explained_prettyval_raw=format(round(var_explained_raw*100,1),digits=1)) %>%
  select(i,j,var_explained_raw,var_explained_prettyval_raw) %>%
  full_join(anova_results_ready,by=c("i","j")) %>%
  filter(!(as.character(i)==as.character(j))) %>%
  mutate(whichval=ifelse(as.integer(i)>(as.integer(j)-2),var_explained_testgroup,var_explained_raw),
         whichprettyval=ifelse(as.integer(i)>(as.integer(j)-2),prettyval,var_explained_prettyval_raw))
tile_varexplained_merged_withtrueallos<-FigureS8_combined_ready %>% 
  ggplot(aes(x=j,y=i,label=whichprettyval,fill=whichval)) + 
  scale_fill_gradient(low="white",high="darkgreen",guide="none") + 
  geom_tile() + geom_text() + theme_bw() + 
  geom_vline(xintercept=2.5,linetype="solid") + 
  geom_vline(xintercept=c(5.5,7.5,9.5,11.5),linetype="dashed") + 
  geom_hline(yintercept=c(3.5,5.5,7.5,9.5),linetype="dashed") + 
  geom_abline(slope=1,intercept=-2,linetype="solid",linewidth=2) + 
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x=element_text(angle=90))
tile_varexplained_merged_withtrueallos_nonames<-tile_varexplained_merged_withtrueallos + 
  scale_x_discrete(labels=NULL) + scale_y_discrete(labels=NULL) + 
  ylab(NULL) + xlab(NULL)
tile_varexplained_merged_withtrueallos_nonames_box<-tile_varexplained_merged_withtrueallos_nonames + 
  geom_rect(data=boxes_withtrueallos2_shift,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=colorgroup,fill=NULL,label=NULL),fill="transparent",linewidth=1,linetype="dashed") + 
  scale_color_manual(values=rect_colors,name=NULL,
                     labels=c("Within-site","Nearby-site","Within-region","Between-region","Between-species")) + 
  coord_equal(ratio=1)
tile_varexplained_merged_withtrueallos_nonames_box
ggsave(file.path(thisfolder,"Paper_figs/tile_varexplained_anovas_plusraw_withtrueallos.pdf"),tile_varexplained_merged_withtrueallos,device="pdf",width = 8,height=6,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/tile_varexplained_anovas_plusraw_nonames_boxes_withtrueallos.pdf"),tile_varexplained_merged_withtrueallos_nonames_box,device="pdf",width = 8,height=6,units="in",dpi = 600)  

################################################################################
#####    Figure S9                                                        ######
################################################################################

###outlier overlap tiles with true-allos
tile_overlaps_missingness_boxes_withtrueallos<-count_overlaps %>% mutate(i=factor(i,levels=key_groups_customorder2_withtrueallos),
                                                           j=factor(j,levels=key_groups_customorder2_withtrueallos)) %>%
  filter(!is.na(i),!is.na(j)) %>%
  mutate(observed_overlap=ifelse(as.integer(j)>=as.integer(i),observed_overlap_missingness,observed_overlap_uncorrected),
         pval_permutes=ifelse(as.integer(j)>=as.integer(i),pval_permutes_missingness,pval_permutes_uncorrected),
         pval_bonf_corrected=pval_permutes*110,significance=ifelse(pval_bonf_corrected<0.001,"***",
                                                                   ifelse(pval_bonf_corrected<0.01,"**",
                                                                          ifelse(pval_bonf_corrected<0.05,"*",
                                                                                 ifelse(pval_bonf_corrected<0.1,".","")))),
         tile_label=ifelse(as.integer(i)==as.integer(j),
                           paste0(observed_overlap_missingness,"\n",observed_overlap_uncorrected),
                           paste0(observed_overlap,"\n",significance))) %>%
  ggplot(aes(x=i,y=j,fill=observed_overlap,label=tile_label)) + 
  geom_tile() + geom_text() + 
  geom_rect(data=boxes_withtrueallos2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=colorgroup,x=NULL,y=NULL,fill=NULL,label=NULL),fill="transparent",linewidth=1,linetype="dashed") + 
  scale_color_manual(values=rect_colors,name=NULL,
                     labels=c("Within-site","Nearby-site","Within-region","Between-region","Between-species")) + 
  scale_fill_gradient(low="white",high="darkgreen",guide="none") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) + 
  ylab(NULL) + xlab(NULL) + 
  geom_hline(yintercept=c(3.5,5.5,7.5)) + geom_vline(xintercept=c(3.5,5.5,7.5)) + 
  coord_equal(ratio=1)
tile_overlaps_missingness_withtrueallos_nonames_boxes<-tile_overlaps_missingness_boxes_withtrueallos + 
  scale_x_discrete(labels=NULL) + scale_y_discrete(labels=NULL)
ggsave(file.path(thisfolder,"Paper_figs/deltapredict_outliers_overlap_matrix_missingness_boxes_withtrueallos.pdf"),tile_overlaps_missingness_boxes_withtrueallos,device="pdf",width = 8,height=7,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/deltapredict_outliers_overlap_matrix_missingness_nonames_boxes_withtrueallos.pdf"),tile_overlaps_missingness_withtrueallos_nonames_boxes,device="pdf",width = 7,height=6,units="in",dpi = 600)  

################################################################################
#####    Table S5                                                         ######
################################################################################

###Table of percentiles for raw ancestry frequencies and frequency residuals at GOI windows
quantile_table<-data.frame()
for (i in neworder) {
  subdata_freqs<-all_groups_window_freqs %>% filter(keygroup==i) %>% arrange(nasutus_frequency) %>% select(keygroup,chrom,windowend,nasutus_frequency)
  subdata_deltas<-prediction_deltas_allmodels %>% filter(keygroup==i) %>% arrange(delta_predict_missingness) %>% select(keygroup,chrom,windowend,delta_predict_missingness)
  ecdf_freqs<-ecdf(subdata_freqs$nasutus_frequency)
  ecdf_deltas<-ecdf(subdata_deltas$delta_predict_missingness)
  subdata_quantiles<-full_join(subdata_freqs,subdata_deltas,by = c("keygroup","chrom","windowend")) %>%
    mutate(quantile_freqs=ecdf_freqs(nasutus_frequency),quantile_deltas=ecdf_deltas(delta_predict_missingness))
  get_GOI_quantiles<-Mantel_regions_GOIs %>% select(Name,chrom,last_windowend_IM767) %>% rename("windowend"="last_windowend_IM767") %>%
    left_join(subdata_quantiles)
  get_surrounding_quantiles_means<-choose_windows_limited %>% select(chrom,windowend) %>% unique() %>%
    left_join(subdata_quantiles) %>% group_by(chrom) %>% summarize(neighborhoodmean_nasutus_frequency=mean(nasutus_frequency),
                                                                   neighborhoodmean_delta_predict_missingness=mean(delta_predict_missingness),
                                                                   neighborhoodmean_quantile_freqs=mean(quantile_freqs),
                                                                   neighborhoodmean_quantile_deltas=mean(quantile_deltas))
  table_ready<-left_join(get_GOI_quantiles,get_surrounding_quantiles_means,by="chrom")
  quantile_table<-bind_rows(quantile_table,table_ready)
}
quantile_table_ready<-quantile_table %>% mutate(grouplabel=factor(keygroup,levels=neworder)) %>% filter(!is.na(grouplabel)) %>%
  mutate(grouplabel=dplyr::recode(grouplabel, !!!setNames(GOI_grouplabels,neworder)),
         window_midpoint=windowend-25000) 

write.table(quantile_table_ready,file=file.path(thisfolder,"Paper_tables/GOI_quantiles.txt"),quote = F,sep="\t",col.names=T,row.names=F)


################################################################################
#####    Figures S10-S13                                                ######
################################################################################


###Supp figure with remaining genome scans not in Figure 6
supp_scans_ready<-ancestry_predictions_ready %>% mutate(keygroup=factor(keygroup,levels=key_groups_customorder2$keygroup),
                                      cohort=factor(ifelse(str_detect(keygroup,"nasutus"),"M. nasutus",
                                                           ifelse(str_detect(keygroup,"guttatus"),"M. guttatus","Hybrid")),
                                                    levels=c("M. nasutus","Hybrid","M. guttatus")),
                                      window_midpoint=windowend-25000) %>% 
  filter(statistic %in% c("prediction_missingness","nasutus_frequency")) %>%
  mutate(grouplabel=factor(keygroup,levels=neworder)) %>% filter(!is.na(grouplabel)) %>%
  mutate(grouplabel=dplyr::recode(grouplabel, !!!setNames(GOI_grouplabels,neworder)))  %>% 
  select(chrom,windowend,keygroup,statistic,nasutus_frequency,cohort,window_midpoint,grouplabel) %>%
  pivot_wider(names_from = "statistic",values_from="nasutus_frequency") %>%
  left_join(delta_quantiles %>% select(chrom,windowend,keygroup,quantile_freqs,quantile_deltas),by=c("chrom","windowend","keygroup"))
scan1_2<-supp_scans_ready %>% filter(chrom %in% c("Chr01","Chr02")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scan3_4<-supp_scans_ready %>% filter(chrom %in% c("Chr03","Chr04")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scan5_6<-supp_scans_ready %>% filter(chrom %in% c("Chr05","Chr06")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scan7_8<-supp_scans_ready %>% filter(chrom %in% c("Chr07","Chr08")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scan9_10<-supp_scans_ready %>% filter(chrom %in% c("Chr09","Chr10")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scan11_12<-supp_scans_ready %>% filter(chrom %in% c("Chr11","Chr12")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scan13_14<-supp_scans_ready %>% filter(chrom %in% c("Chr13","Chr14")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free") + 
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
scans_A<-plot_grid(scan1_2,scan3_4,ncol=1,align="lr",axis="hv")
scans_B<-plot_grid(scan5_6,scan7_8,ncol=1,align="lr",axis="hv")
scans_C<-plot_grid(scan9_10,scan11_12,ncol=1,align="lr",axis="hv")
scans_D<-plot_grid(scan13_14,ncol=1,align="lr",axis="hv")

ggsave(file.path(thisfolder,"Paper_figs/suppfig_scans_chr1_4.pdf"),scans_A,device="pdf",width = 10,height=10,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/suppfig_scans_chr5_8.pdf"),scans_B,device="pdf",width = 10,height=10,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/suppfig_scans_chr9_12.pdf"),scans_C,device="pdf",width = 10,height=10,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/suppfig_scans_chr13_14.pdf"),scans_D,device="pdf",width = 10,height=5,units="in",dpi = 600)  

