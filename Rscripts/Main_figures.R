
################################################################################
#####    Code to produce main figures and tables                          ######
################################################################################

library(tidyverse)
library(readxl)
library(geodata)
library(tidyterra)
library(sf)
library(ggrepel)
library(ggspatial)
library(patchwork)
library(cowplot)
library(ggtext)
library(RColorBrewer)

colors_sp<-c("#993929","#CC892F","#FFDA35") #G,A,N
colors_regions<-c("#005FA3","#78290F")
colors_groups<-c("#005FA3","#27B7CE","#78290F","#FF8B1F","#CAB9A5")
thisfolder<-"./"


################################################################################
#####    Figure 1                                                         ######
################################################################################


###broad map
DPR_mapdata<-read_excel(file.path(thisfolder,"Metadata/Pop_codes_DPR_area_sequencing.xlsx")) %>%
  select(Code,Latitude,Longitude,Elevation,Species_sequenced) %>% 
  mutate(cohort=ifelse(Elevation>1000,"high-elevation sympatric",ifelse(Species_sequenced=="sym","sympatric","allopatric"))) %>%
  mutate(cohort=factor(cohort,levels=c("allopatric","sympatric","high-elevation sympatric")))

area_locations<-data.frame(Code=c("Northern region","Southern region"),lat=c(45.71133219,37.82925), long=c(-121.3636758,-120.3386))


CAC_site_locations<-read_excel(file.path(thisfolder,"Metadata/CAC-site-locations.xlsx")) %>%
  filter(Site %in% c("CAC","LM")) %>% rename("lat"="latitude","long"="longitude")

us <- gadm(country="USA",level=1,path = file.path(thisfolder,"gadm/USA_outline"))
canada <- gadm(country="CAN",level=1,path = file.path(thisfolder,"gadm/CAN_outline"))
mexico <- gadm(country="MEX",level=1,path = file.path(thisfolder,"gadm/MEX_outline"))

us.states<-us[us$NAME_1 %in% c("Oregon","Washington","California"),]
us.states.plus<-us[us$NAME_1 %in% c("California","Oregon","Washington","Idaho","Nevada","Arizona","Montana","Utah","Wyoming"),]
can.states<-canada[canada$NAME_1 %in% c("British Columbia","Alberta"),]
mex.states<-mexico[mexico$NAME_1 %in% c("Baja California","Sonora"),]

###GBIF download citation: GBIF.org (05 September 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.sj7c59
gutnas<-read_xlsx(file.path(thisfolder,"Metadata/GBIF_gutnas_northamerica.xlsx")) %>% 
  select(gbifID,species,countryCode,stateProvince,decimalLongitude,decimalLatitude)
gutnassort<-arrange(gutnas,species)

###Broad map
westcoast<-ggplot(us.states.plus) + geom_spatvector(fill=NA) +
  geom_spatvector(data=can.states,fill=NA) + geom_spatvector(data=mex.states,fill=NA) + 
  theme_bw() + coord_sf(xlim=c(-128,-110),ylim=c(30,50)) + 
  geom_point(data=gutnassort,aes(x=decimalLongitude,y=decimalLatitude,color=species),group=NA,alpha=0.1) + 
  scale_color_manual(values=colors_sp[c(1,3)],name="GBIF\noccurrence\nrecords",labels=c(expression(italic("M. guttatus")),expression(italic("M. nasutus")))) + 
  #scale_color_manual(values=colors_sp[c(1,3)],guide="none") + 
  geom_rect(xmin=-121.42,xmax=-121.34,ymin=45.67,ymax=45.73,color="black",fill=NA) + 
  geom_rect(xmin=-120.5,xmax=-119.7,ymin=37.5,ymax=38.1,color="black",fill=NA) + 
  geom_label_repel(data=area_locations,aes(x=long,y=lat,group=group,label=Code,fill=Code,color=NA),group=NA,color="white",nudge_y = 1) + 
  scale_fill_manual(values=colors_regions,guide="none") + 
  ylab("Latitude") + xlab("Longitude") + theme(axis.text.x=element_text(angle=45,hjust=1,size=12),text=element_text(size=12),
                                               legend.position="bottom",legend.direction="vertical",legend.justification="top") + 
  theme(legend.background=element_blank()) + 
  guides(color=guide_legend(override.aes=list(alpha=1)))

###Northern region map
northern_region<-ggplot(us.states.plus) + geom_spatvector(fill=NA) +
  theme_bw() + coord_sf(xlim=c(-121.42,-121.34),ylim=c(45.67,45.73)) + 
  annotation_scale(bar_cols="black") +
  #geom_point(data=CAC_site_locations,aes(x=long,y=lat,group=group,fill=Site),group=NA,shape=21,size=3) + 
  #scale_fill_manual(values=colors_groups[1:2],name="Northern region",labels=c("CC-sym","LM-sym")) + 
  geom_point(data=CAC_site_locations,aes(x=long,y=lat,group=group,shape=Site),group=NA,fill="transparent",size=2) + 
  scale_shape_manual(values=c(23,24),name=NULL,labels=c("CAC-sym","LM-sym")) + 
  ylab("Latitude") + xlab("Longitude") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        text=element_text(size=12),
        legend.background = element_blank(),
        panel.grid = element_blank()
        )

southern_region<-ggplot(us.states.plus) + geom_spatvector(fill=NA) +
  theme_bw() + coord_sf(xlim=c(-120.5,-119.7),ylim=c(37.5,38.1)) + 
  annotation_scale(bar_cols="black") +
  #geom_point(data=DPR_mapdata,aes(x=Longitude,y=Latitude,group=group,fill=cohort),group=NA,shape=21,size=3) + 
  #scale_fill_manual(values=colors_groups[3:5],name="Southern region",labels=c("FH-allo","FH-sym","MO-sym")) + 
  geom_point(data=DPR_mapdata,aes(x=Longitude,y=Latitude,group=group,shape=cohort),group=NA,fill="transparent",size=2) + 
  scale_shape_manual(values=c(21,22,25),name=NULL,labels=c("Foothills-allo","Foothills-sym","Montane-sym")) + 
  ylab("Latitude") + xlab("Longitude") +  
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        text=element_text(size=12),
        legend.background = element_blank(),
        panel.grid = element_blank()
        )

group_histograms<-list()
for (i in unique(final_individual_list$group)) {
  group_subdata<-final_individual_list %>% filter(group==i)
  title=ifelse(i=="Southern_highelev","Southern Montane",
               ifelse(str_detect(i,"Southern_"),str_replace(i,"Southern_","Southern Foothills "),
                      str_replace(i,"_"," ")))
  group_histograms[[i]]<-ggplot(group_subdata,aes(x=nasutus_proportion_windows)) + theme_bw() + 
    geom_histogram(data=subset(group_subdata,freq_group=="guttatus"),binwidth=0.05,boundary=0,fill=colors_sp[1]) + 
    geom_histogram(data=subset(group_subdata,freq_group=="hybrid"),binwidth=0.05,boundary=0,fill=colors_sp[2]) + 
    geom_histogram(data=subset(group_subdata,freq_group=="nasutus"),binwidth=0.05,boundary=0,fill=colors_sp[3]) + 
    scale_x_continuous(limits=c(0,1),expand=c(0,0),name="Hybrid index") + 
    scale_y_continuous(name="Number of samples") + 
    ggtitle(title) + theme(title=element_text(hjust = 0.5,size=12),text=element_text(size=12))
}
group_histograms_facetted_ready<-final_individual_list %>% 
  mutate(title=case_match(group,
                          "CAC" ~ "Northern - CAC - sympatric",
                          "LM" ~ "Northern - LM - sympatric",
                          "Southern_allopatric" ~ "Southern - Foothills - allopatric",
                          "Southern_sympatric" ~ "Southern - Foothills - sympatric",
                          "Southern_highelev" ~ "Southern - Montane - sympatric"))
group_histograms_facetted<-group_histograms_facetted_ready %>%
  ggplot(aes(x=nasutus_proportion_windows)) + theme_bw() + 
  facet_wrap(~title,ncol=1,scales="free_y") + 
  geom_histogram(data=subset(group_histograms_facetted_ready,freq_group=="guttatus"),binwidth=0.05,boundary=0,fill=colors_sp[1]) + 
  geom_histogram(data=subset(group_histograms_facetted_ready,freq_group=="hybrid"),binwidth=0.05,boundary=0,fill=colors_sp[2]) + 
  geom_histogram(data=subset(group_histograms_facetted_ready,freq_group=="nasutus"),binwidth=0.05,boundary=0,fill=colors_sp[3]) + 
  scale_x_continuous(limits=c(0,1),expand=c(0,0),name="Hybrid index") + 
  scale_y_continuous(name="Number of samples") + 
  theme(text=element_text(size=12))

percent_variance<-round(angsd_eigen$values*100/sum(angsd_eigen$values),2)[1:2]
PCA<-final_individual_list %>% 
  mutate(title=factor(case_match(group,
                                 "CAC" ~ "N-CAC-sym",
                                 "LM" ~ "N-LM-sym",
                                 "Southern_allopatric" ~ "S-Foothills-allo",
                                 "Southern_sympatric" ~ "S-Foothills-sym",
                                 "Southern_highelev" ~ "S-Montane-sym"),
                      levels=c("N-CAC-sym","N-LM-sym","S-Foothills-allo","S-Foothills-sym","S-Montane-sym"))) %>%
  ggplot(aes(x=PC1,y=PC2,fill=freq_group,shape=title)) + 
  geom_point() + 
  xlab(paste0("PC1: ",percent_variance[1],"% of variance")) + 
  ylab(paste0("PC2: ",percent_variance[2],"% of variance")) + 
  scale_fill_manual(values=colors_sp,name=NULL) + theme_bw() + 
  scale_shape_manual(values=c(23,24,21,22,25),name=NULL) + 
  theme(text=element_text(size=12),
        legend.background = element_blank(),
        panel.grid = element_blank()) + 
  guides(
    shape = guide_legend(
      label.theme = element_text(face = "italic")),
    fill = guide_legend(override.aes = list(shape=21,size=3)))

full_fig1_plotlist<-list(westcoast,northern_region,southern_region,PCA,group_histograms_facetted)
full_fig1_plotlist_nolegends<-list(westcoast + theme(legend.position="none"),
                                   northern_region + theme(legend.position="none"),
                                   southern_region + theme(legend.position="none"),
                                   PCA + theme(legend.position="none"),
                                   group_histograms_facetted + theme(legend.position="none"))
design<-"
ABE
ACE
DDE
DDE
"
full_fig1<-full_fig1_plotlist %>% wrap_plots(design = design,widths = c(1,NA,1)) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face="bold"),
                                            panel.background = element_rect(fill = "transparent", colour = NA),  
                                            plot.background = element_rect(fill = "transparent", colour = NA))
full_fig1_nolegends<-full_fig1_plotlist_nolegends %>% wrap_plots(design = design,widths = c(1,NA,1)) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face="bold"),
                                            panel.background = element_rect(fill = "transparent", colour = NA),  
                                            plot.background = element_rect(fill = "transparent", colour = NA))

ggsave(file.path(thisfolder,"Paper_figs/Fig1_threemaps_PCA_plushistograms_new3.png"),full_fig1,device="png",width=10,height=10,units="in",dpi=600)
ggsave(file.path(thisfolder,"Paper_figs/Fig1_threemaps_PCA_plushistograms_nolegends.png"),full_fig1_nolegends,device="png",width=10,height=10,units="in",dpi=600)

################################################################################
#####    Table 1                                                         ######
################################################################################

####Table of correlations
write.table(cor(correlation_matrix_withsquares),file = file.path(thisfolder,"Paper_figs/correlations_table.txt"),quote=F,sep="\t")
write.table(model_test_variances %>% filter(model=="newmodel"),file=file.path(thisfolder,"Paper_figs/model_variances_AICs_table.txt"),quote=F,sep="\t")

###correlations including gene density/cM
correlation_matrix_withdensity<-correlation_matrix %>% 
  mutate(relative_chr_pos_sq=relative_chr_pos^2,genes_percM=ngenes/(M_per_bp_1Mb*5000000)) %>% 
  filter(!is.na(genes_percM),!is.infinite((genes_percM)))
r_values<-cor(correlation_matrix_withdensity %>% select(-chrom,-windowend)) %>% as.data.frame() %>%
  select(relative_chr_pos,relative_chr_pos_sq,ngenes,M_per_bp_1Mb,M_per_bp_50kb,genes_percM) %>% 
  rownames_to_column("x") %>%
  mutate(x=factor(x,levels=c("relative_chr_pos","relative_chr_pos_sq","ngenes","M_per_bp_1Mb","M_per_bp_50kb","genes_percM",
                             "CAC_guttatus_137","CAC_hybrid_183","CAC_nasutus_58","LM_hybrid_56",
                             "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174",
                             "Southern_sympatric_nasutus_75","Southern_highelev_guttatus_12","Southern_highelev_hybrid_16"
  ))) %>%
  arrange(x)
View(r_values)
write.table(r_values,file=file.path(thisfolder,"Table1_with_genedensity.txt"),quote=F,row.names=F,col.names=T,sep="\t")              

################################################################################
#####    Figure 2                                                         ######
################################################################################

###Histograms of admixed allele counts per window, for each group
allele_count_histograms<-list()
allele_frequency_histograms<-list()
allele_frequency_histograms_main<-list()
colorkey=colors_groups[c(1,1,2,3,3,3,4,5,5,1,4)]
neworder<-c("CAC_guttatus_137","CAC_hybrid_183","CAC_nasutus_58",
            "LM_hybrid_56",
            "Southern_allopatric_guttatus_61","Southern_sympatric_guttatus_174","Southern_sympatric_nasutus_75",
            "Southern_highelev_guttatus_12","Southern_highelev_hybrid_16")
hist_grouplabels<-c(
  "N-CAC-sym guttatus, n=137",
  "N-CAC-sym hybrid, n=183",
  "N-CAC-sym nasutus, n=75",
  "N-LM-sym hybrid, n=56",
  "S-Foothills-allo guttatus, n=61",
  "S-Foothills-sym guttatus, n=174",
  "S-Foothills-sym nasutus, n=75",
  "S-Montane-sym guttatus, n=12",
  "S-Montane-sym hybrid, n=16"
)
freq_hists_all_ready<-all_groups_window_freqs %>% ungroup() %>% filter(!is.na(nasutus_frequency)) %>%
  mutate(grouplabel=factor(keygroup,levels=neworder)) %>% filter(!is.na(grouplabel)) %>%
  mutate(grouplabel=dplyr::recode(grouplabel, !!!setNames(hist_grouplabels,neworder))) %>%
  mutate(freq_bin=ceiling(nasutus_frequency*200)/200) %>% group_by(grouplabel,freq_bin) %>%
  summarize(n_windows=n())
maxmins<-freq_hists_all_ready %>% group_by(grouplabel) %>% summarize(max=max(freq_bin),min=min(freq_bin)) %>%
  pivot_longer(cols=c(max,min),names_to = "which",values_to="arrow_x") %>% 
  mutate(xadjust=ifelse(which=="max",arrow_x+0.02,arrow_x-0.02)) %>%
  left_join(freq_hists_all_ready %>% group_by(grouplabel) %>% summarize(maxy=max(n_windows)),relationship="many-to-one")
freq_hists_all<-freq_hists_all_ready %>% ggplot(aes(x=freq_bin,y=n_windows,fill=grouplabel)) + 
  facet_wrap(~grouplabel,scale="free_y",ncol=1) + 
  geom_col() + 
  theme_bw() + ylab("Window count") + 
  scale_fill_discrete(palette=colors_groups[c(1,1,1,2,3,4,4,5,5)],guide="none") + 
  geom_segment(data=maxmins,aes(x=xadjust,xend=arrow_x,y=maxy/2,yend=0),arrow=arrow(type="closed",length = unit(0.05, "in")),linetype="dashed") + 
  scale_x_continuous(name="Ancestry frequency",limits=c(-0.05,1.05)) + 
  theme(text=element_text(size=12))

###Fig 2A black-and-white version
freq_hists_all_bw<-freq_hists_all_ready %>% ggplot(aes(x=freq_bin,y=n_windows)) + 
  facet_wrap(~grouplabel,scale="free_y",ncol=1) + 
  geom_rect(data=maxmins %>% filter(which=="min"),aes(xmin=0,xmax=arrow_x,ymin=0,ymax=maxy,x=NULL,y=NULL),fill="lightgrey") + 
  geom_rect(data=maxmins %>% filter(which=="max"),aes(xmin=arrow_x,xmax=1,ymin=0,ymax=maxy,x=NULL,y=NULL),fill="lightgrey") + 
  geom_col(fill="black") + 
  theme_bw() + ylab("Window count") + 
  scale_x_continuous(name="Ancestry frequency",limits=c(-0.05,1.05)) + 
  theme(text=element_text(size=12))


###Figure 2B: correlations for CAC
rsq_position<-round(summary(lm(CAC_hybrid_183~relative_chr_pos,data=correlation_matrix))$adj.r.squared,3)
rsq_position_sq<-round(summary(lm(CAC_hybrid_183~I(relative_chr_pos^2),data=correlation_matrix))$adj.r.squared,3)
raw_position<-ggplot(correlation_matrix,aes(x=relative_chr_pos,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(method="lm",color="red") + geom_smooth(method="lm",formula="y ~ x + I(x^2)",color="blue") + 
  theme_bw() + xlab("relative chromosomal position") + ylab("N-CAC-sym hybrid\nancestry frequency") + 
  theme(text=element_text(size=12))
rsq_recomb<-round(summary(lm(CAC_hybrid_183~M_per_bp_1Mb,data=correlation_matrix))$adj.r.squared,3)
rsq_recomb_50kb<-round(summary(lm(CAC_hybrid_183~M_per_bp_50kb,data=correlation_matrix))$adj.r.squared,3)
raw_recomb<-ggplot(correlation_matrix,aes(x=M_per_bp_1Mb,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(method="lm",color="red") + 
  theme_bw() + xlab("Recombination rate (M/bp)") + ylab("N-CAC-sym hybrid\nancestry frequency") + 
  theme(text=element_text(size=12))
rsq_genes<-round(summary(lm(CAC_hybrid_183~ngenes,data=correlation_matrix))$adj.r.squared,3)
raw_genes<-ggplot(correlation_matrix,aes(x=ngenes,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(method="lm",color="red") + 
  theme_bw() + xlab("genes per 50kb") + ylab("N-CAC-sym hybrid\nancestry frequency") + 
  theme(text=element_text(size=12))

genes_per_cM<-correlation_matrix %>% mutate(genes_percM=ngenes/(M_per_bp_50kb*5000000)) %>% filter(!is.na(genes_percM),!is.infinite((genes_percM)))
rsq_genes_per_cM<-round(summary(lm(CAC_hybrid_183~genes_percM,data=genes_per_cM))$r.squared,3)
raw_genespercM<-ggplot(genes_per_cM,aes(x=genes_percM,y=CAC_hybrid_183)) + geom_point() + 
  geom_smooth(method="lm",color="red") + 
  theme_bw() + xlab("gene density per cM") + ylab("N-CAC-sym hybrid\nancestry frequency") + 
  theme(text=element_text(size=12))

combine<-plot_grid(raw_position,raw_recomb,raw_genes,raw_genespercM,ncol=1)

Figure2_twocols_bw<-plot_grid(freq_hists_all_bw,combine,ncol=2)
ggsave(file.path(thisfolder,"Paper_figs/Figure2_twoparts_bw.pdf"),Figure2_twocols_bw,device="pdf",width=6,height=8,units="in",dpi=600)

################################################################################
#####    Figure 3                                                         ######
################################################################################

write.table(anova_results_all,file=file.path(thisfolder,"Paper_tables/anova_pairwise_results.txt"),quote=F,sep="\t",row.names=F)

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

###Combine uncorrected and corrected into opposing diagonals
Figure3_combined_ready<-all_correlations_withcorrections %>% 
  mutate(i=factor(group1,levels=key_groups_customorder2_withtrueallos),
         j=factor(group2,levels=c("genome_structure","missingness",key_groups_customorder2_withtrueallos)),
         var_explained_raw=uncorrected^2,
         var_explained_prettyval_raw=format(round(var_explained_raw*100,1),digits=1)) %>%
  select(i,j,var_explained_raw,var_explained_prettyval_raw) %>%
  full_join(anova_results_ready,by=c("i","j")) %>%
  filter(!(i %in% c("Southern_true_allopatric","Southern_pseudo_allopatric")),
         !(j %in% c("Southern_true_allopatric","Southern_pseudo_allopatric")),
         !(as.character(i)==as.character(j))) %>%
  mutate(whichval=ifelse(as.integer(i)>(as.integer(j)-2),var_explained_testgroup,var_explained_raw),
         whichprettyval=ifelse(as.integer(i)>(as.integer(j)-2),prettyval,var_explained_prettyval_raw))
tile_varexplained_merged<-Figure3_combined_ready %>% 
  ggplot(aes(x=j,y=i,label=whichprettyval,fill=whichval)) + 
  scale_fill_gradient(low="white",high="darkgreen",guide="none") + 
  geom_tile() + geom_text() + theme_bw() + 
  geom_vline(xintercept=2.5,linetype="solid") + 
  geom_vline(xintercept=c(5.5,7.5,9.5),linetype="dashed") + 
  geom_hline(yintercept=c(3.5,5.5,7.5),linetype="dashed") + 
  geom_abline(slope=1,intercept=-2,linetype="solid",linewidth=2) + 
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x=element_text(angle=90))
tile_varexplained_merged_nonames<-tile_varexplained_merged + 
  scale_x_discrete(labels=NULL) + scale_y_discrete(labels=NULL) + 
  ylab(NULL) + xlab(NULL)
tile_varexplained_merged_nonames_box<-tile_varexplained_merged_nonames + 
  geom_rect(data=boxes2_shift,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=colorgroup,fill=NULL,label=NULL),fill="transparent",linewidth=1,linetype="dashed") + 
  scale_color_manual(values=rect_colors,name=NULL,
                     labels=c("Within-site","Nearby-site","Within-region","Between-region","Between-species")) + 
  coord_equal(ratio=1)
tile_varexplained_merged_nonames_box
ggsave(file.path(thisfolder,"Paper_figs/tile_varexplained_anovas_plusraw_nonames_boxes.pdf"),tile_varexplained_merged_nonames_box,device="pdf",width = 7,height=5,units="in",dpi = 600)  

################################################################################
#####    Figure 4                                                         ######
################################################################################

###Figure 4A: Outlier overlap comparing raw to missingness (RAW now in bottom right triangle!)

tile_overlaps_missingness_boxes<-count_overlaps %>% mutate(i=factor(i,levels=key_groups_customorder2$keygroup),
                                                           j=factor(j,levels=key_groups_customorder2$keygroup)) %>%
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
  geom_rect(data=boxes2,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=colorgroup,x=NULL,y=NULL,fill=NULL,label=NULL),fill="transparent",linewidth=1,linetype="dashed") + 
  scale_color_manual(values=rect_colors,name=NULL,
                     labels=c("Within-site","Nearby-site","Within-region","Between-region","Between-species")) + 
  scale_fill_gradient(low="white",high="darkgreen",guide="none") + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90)) + 
  ylab(NULL) + xlab(NULL) + 
  geom_hline(yintercept=c(3.5,5.5,7.5)) + geom_vline(xintercept=c(3.5,5.5,7.5)) + 
  coord_equal(ratio=1)
tile_overlaps_missingness_nonames_boxes<-tile_overlaps_missingness_boxes + 
  scale_x_discrete(labels=NULL) + scale_y_discrete(labels=NULL)
ggsave(file.path(thisfolder,"Paper_figs/deltapredict_outliers_overlap_matrix_missingness_boxes.pdf"),tile_overlaps_missingness_boxes,device="pdf",width = 8,height=7,units="in",dpi = 600)  
ggsave(file.path(thisfolder,"Paper_figs/deltapredict_outliers_overlap_matrix_missingness_nonames_boxes.pdf"),tile_overlaps_missingness_nonames_boxes,device="pdf",width = 6,height=5,units="in",dpi = 600)  

###Figure 4b-c: effects of model correction on overlap
###effect of corrections on outlier overlap
outlier_change_residuals<-count_overlaps %>% 
  mutate(i=factor(i,levels=key_groups_customorder2$keygroup),
         j=factor(j,levels=key_groups_customorder2$keygroup)) %>%
  filter(!(i==j),!is.na(i),!is.na(j)) %>%
  mutate(direction=factor(sign(observed_overlap_corrected-observed_overlap_uncorrected))) %>%
  ggplot(aes(x=observed_overlap_uncorrected,y=observed_overlap_uncorrected,yend=observed_overlap_corrected)) +
  geom_segment(color="red") + geom_point(color="black") + geom_point(aes(y=observed_overlap_corrected),color="red") + 
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey") + 
  geom_hline(yintercept=7.85,linetype="dotted",color="grey") + 
  geom_vline(xintercept=7.85,linetype="dotted",color="grey") + 
  xlab("overlap of uncorrected outliers") + 
  ylab("overlap of corrected residual outliers") + 
  theme_bw()
outlier_change_violin<-count_overlaps %>% 
  mutate(i=factor(i,levels=key_groups_customorder2$keygroup),
         j=factor(j,levels=key_groups_customorder2$keygroup)) %>%
  filter(!(i==j),!is.na(i),!is.na(j)) %>%
  mutate(change=observed_overlap_corrected-observed_overlap_uncorrected) %>%
  ggplot(aes(x=1,y=change)) + 
  geom_hline(yintercept=0,linetype="dotted",color="grey") + 
  geom_violin(color="red",fill="transparent") + geom_jitter(color="red",height=0,width=0.1) + theme_bw() + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  ylim(-30,30)
outlier_change_missingness<-count_overlaps %>% 
  mutate(i=factor(i,levels=key_groups_customorder2$keygroup),
         j=factor(j,levels=key_groups_customorder2$keygroup)) %>%
  filter(!(i==j),!is.na(i),!is.na(j)) %>%
  mutate(direction=factor(sign(observed_overlap_missingness-observed_overlap_corrected))) %>%
  ggplot(aes(x=observed_overlap_uncorrected,y=observed_overlap_corrected,yend=observed_overlap_missingness)) +
  geom_segment(color="blue") + geom_point(color="red") + geom_point(aes(y=observed_overlap_missingness),color="blue") + 
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey") + 
  geom_hline(yintercept=7.85,linetype="dotted",color="grey") + 
  geom_vline(xintercept=7.85,linetype="dotted",color="grey") + 
  xlab("overlap of uncorrected outliers") + 
  ylab("after missingness correction") +
  theme_bw()
outlier_change_both<-count_overlaps %>% 
  mutate(i=factor(i,levels=key_groups_customorder2$keygroup),
         j=factor(j,levels=key_groups_customorder2$keygroup)) %>%
  filter(!(i==j),!is.na(i),!is.na(j)) %>%
  mutate(direction=factor(sign(observed_overlap_missingness-observed_overlap_corrected))) %>%
  ggplot(aes(x=observed_overlap_uncorrected,y=observed_overlap_corrected,yend=observed_overlap_missingness)) +
  geom_point(color="black",aes(y=observed_overlap_uncorrected)) + 
  geom_segment(color="red",linetype="dashed",aes(y=observed_overlap_uncorrected,yend=observed_overlap_corrected)) +
  geom_segment(color="blue",linetype="dashed") + geom_point(color="red") + geom_point(aes(y=observed_overlap_missingness),color="blue") + 
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey") + 
  geom_hline(yintercept=7.85,linetype="dotted",color="grey") + 
  geom_vline(xintercept=7.85,linetype="dotted",color="grey") + 
  xlab("overlap of uncorrected outliers") + 
  ylab("after missingness correction") +
  theme_bw()

outlier_change_violin_missingness<-count_overlaps %>% 
  mutate(i=factor(i,levels=key_groups_customorder2$keygroup),
         j=factor(j,levels=key_groups_customorder2$keygroup)) %>%
  filter(!(i==j),!is.na(i),!is.na(j)) %>%
  mutate(change=observed_overlap_missingness-observed_overlap_corrected) %>%
  ggplot(aes(x=1,y=change)) + 
  geom_hline(yintercept=0,linetype="dotted",color="grey") + 
  geom_violin(color="blue",fill="transparent") + geom_jitter(color="blue",height=0,width=0.1) + theme_bw() + 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  ylim(-30,30)
outlier_change_four<-plot_grid(outlier_change_residuals,outlier_change_violin,
                               outlier_change_missingness,outlier_change_violin_missingness,
                               ncol=2,rel_widths=c(3,1),axis="tb",align="h")

outlier_change_three<-plot_grid(outlier_change_both+ylab(NULL)+xlab(NULL),
                                outlier_change_violin,outlier_change_violin_missingness,
                                ncol=3,rel_widths=c(3,1,1),axis="tb",align="h")
ggsave(file.path(thisfolder,"Paper_figs/outlier_model_effects_combined.pdf"),outlier_change_three,device="pdf",width=8,height=4,units="in",dpi=600)


################################################################################
#####    Table 2                                                         ######
################################################################################

###do outliers occur in windows with QTL more often than expected?
write.table(outliers_in_windows_summary_longer %>% filter(!(keygroup %in% c("Southern_true_allopatric","Southern_pseudo_allopatric"))),file=file.path(thisfolder,"Paper_tables/outliers_in_windows_summary_longer.txt"),quote = F,sep = "\t",row.names=F)
write.table(outliers_in_windows_summary_split_longer %>% filter(!(keygroup %in% c("Southern_true_allopatric","Southern_pseudo_allopatric"))),file=file.path(thisfolder,"Paper_tables/outliers_in_windows_summary_split_longer.txt"),quote = F,sep = "\t",row.names=F)

################################################################################
#####    Figure 5                                                         ######
################################################################################


###Figure 5: outliers vs. QTLs -- violin plots (now post-missingness)

write.table(ttest_results %>% arrange(p.value),file=file.path(thisfolder,"Paper_figs/Mantel_deltapredict_ttests_windows_newmodel.txt"),quote = F,sep = "\t",row.names=F)
write.table(ttest_results_withGOIs %>% arrange(p.value),file=file.path(thisfolder,"Paper_figs/Mantel_deltapredict_ttests_windows_newmodel_withGOIs.txt"),quote = F,sep = "\t",row.names=F)
write.table(ttest_results_missingness %>% arrange(p.value),file=file.path(thisfolder,"Paper_figs/Mantel_deltapredict_ttests_windows_newmodel_missingness.txt"),quote = F,sep = "\t",row.names=F)
write.table(ttest_results_withGOIs_missingness %>% arrange(p.value),file=file.path(thisfolder,"Paper_figs/Mantel_deltapredict_ttests_windows_newmodel_withGOIs_missingness.txt"),quote = F,sep = "\t",row.names=F)

ttest_results_missingness_ready<-ttest_results_missingness %>% filter(keygroup %in% key_groups_customorder2$keygroup) %>%
  mutate(significance_corrected=ifelse(significance_corrected=="n.s.","",significance_corrected))

Mantel_regions_deltapredict_windows_missingness_nonames_bw<-Mantel_regions_nolarge_freqs_withpredictions_newmodel %>% 
  filter(keygroup %in% key_groups_customorder2$keygroup) %>%
  ggplot(aes(x=keygroup,y=delta_predict_missingness)) + 
  facet_wrap(~Predicted_effect,nrow=1) + 
  geom_violin(scale="width",quantile.linetype="solid") + 
  stat_summary(fun=mean,geom="point",position = position_dodge(width=0.9),shape=4) + 
  geom_text(data=ttest_results_missingness_ready,aes(x=keygroup,y=0.5,label=significance_corrected,color=NULL),color="black") + 
  scale_x_discrete(labels=NULL,name=NULL) + 
  geom_hline(yintercept=0,linetype="dashed") + ylab(NULL) + 
  theme_bw()
ggsave(file.path(thisfolder,"Paper_figs/violin_outliers_vs_QTL_missingness_nonames_bw.pdf"),Mantel_regions_deltapredict_windows_missingness_nonames_bw,device="pdf",width = 7,height=5,units="in",dpi = 600)  

################################################################################
#####    Figure 6                                                         ######
################################################################################

###Figure 6: ancestry at 4 GOIs and surrounding 40+1+40 windows (2Mb on either side)

choose_windows_full<-data.frame()
for (i in Mantel_regions_GOIs$Name) {
  line=Mantel_regions_GOIs %>% filter(Name==i)
  if (i=="pTAC13") {
    subdata<-ancestry_predictions_ready %>% filter(chrom==line$chrom)
  } else {
    subdata<-ancestry_predictions_ready %>% filter(chrom==line$chrom)
  }
  choose_windows_full<-choose_windows_full %>% bind_rows(subdata)
}
choose_windows_full<-choose_windows_full %>% mutate(keygroup=factor(keygroup,levels=key_groups_customorder2$keygroup),
                                          cohort=factor(ifelse(str_detect(keygroup,"nasutus"),"M. nasutus",
                                                               ifelse(str_detect(keygroup,"guttatus"),"M. guttatus","Hybrid")),
                                                        levels=c("M. nasutus","Hybrid","M. guttatus")),
                                          window_midpoint=windowend-25000)

###Ancestry at 4 GOIs and surrounding chrom, facet all groups
GOI_grouplabels<-c(
  "N-CAC-sym\nguttatus, n=137",
  "N-CAC-sym\nhybrid, n=183",
  "N-CAC-sym\nnasutus, n=58",
  "N-LM-sym\nhybrid, n=56",
  "S-Foothills-allo\nguttatus, n=61",
  "S-Foothills-sym\nguttatus, n=174",
  "S-Foothills-sym\nnasutus, n=75",
  "S-Montane-sym\nguttatus, n=12",
  "S-Montane-sym\nhybrid, n=16"
)
new_Fig6_GOIs_ready<-choose_windows_full %>% filter(statistic %in% c("prediction_missingness","nasutus_frequency")) %>%
  mutate(grouplabel=factor(keygroup,levels=neworder)) %>% filter(!is.na(grouplabel)) %>%
  mutate(grouplabel=dplyr::recode(grouplabel, !!!setNames(GOI_grouplabels,neworder))) 
new_Fig6_GOIs_ready_nonas<-new_Fig6_GOIs_ready %>% filter(!(grouplabel %in% c("N-CAC-sym\nnasutus, n=58","S-Foothills-sym\nnasutus, n=75"))) 

delta_quantiles<-data.frame()
for (i in neworder) {
  subdata_freqs<-all_groups_window_freqs %>% filter(keygroup==i) %>% arrange(nasutus_frequency) %>% select(keygroup,chrom,windowend,nasutus_frequency)
  subdata_deltas<-prediction_deltas_allmodels %>% filter(keygroup==i) %>% arrange(delta_predict_missingness) %>% select(keygroup,chrom,windowend,delta_predict_missingness)
  ecdf_freqs<-ecdf(subdata_freqs$nasutus_frequency)
  ecdf_deltas<-ecdf(subdata_deltas$delta_predict_missingness)
  subdata_quantiles<-full_join(subdata_freqs,subdata_deltas,by = c("keygroup","chrom","windowend")) %>%
    mutate(quantile_freqs=ecdf_freqs(nasutus_frequency),quantile_deltas=ecdf_deltas(delta_predict_missingness)) %>% 
    arrange(chrom,windowend)
  delta_quantiles<-bind_rows(delta_quantiles,subdata_quantiles)
}
new_Fig6_GOIs_ready_nonas_withquantiles<-new_Fig6_GOIs_ready_nonas %>% 
  select(chrom,windowend,keygroup,statistic,nasutus_frequency,cohort,window_midpoint,grouplabel) %>%
  pivot_wider(names_from = "statistic",values_from="nasutus_frequency") %>%
  left_join(delta_quantiles %>% select(chrom,windowend,keygroup,quantile_freqs,quantile_deltas),by=c("chrom","windowend","keygroup"))

split1<-new_Fig6_GOIs_ready_nonas_withquantiles %>% filter(chrom %in% c("Chr07","Chr08")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free_y") + 
  geom_vline(data=Mantel_regions_GOIs %>% filter(chrom %in% c("Chr07","Chr08")),aes(xintercept=first_windowend_IM767-25000,x=NULL,y=NULL,color=NULL,linetype=NULL),color="gold",linetype="solid",linewidth=1,alpha=0.8) +
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
ggsave(file.path(thisfolder,"Paper_figs/new_fig6_GOIs_nonas_withquantiles_zoom_lines_split1.pdf"),split1,device="pdf",width = 10,height=5,units="in",dpi = 600)  

split2<-new_Fig6_GOIs_ready_nonas_withquantiles %>% filter(chrom %in% c("Chr13","Chr14")) %>%
  ggplot(aes(x=window_midpoint,y=nasutus_frequency,color=quantile_deltas)) + theme_bw() + 
  facet_grid(grouplabel~chrom,scales="free_y") + 
  geom_vline(data=Mantel_regions_GOIs %>% filter(chrom %in% c("Chr13","Chr14")),aes(xintercept=first_windowend_IM767-25000,x=NULL,y=NULL,color=NULL,linetype=NULL),color="gold",linetype="solid",linewidth=1,alpha=0.8) +
  geom_line(aes(y=prediction_missingness,color=NULL),color="black") + 
  geom_line(aes(color=NULL),color="grey") + 
  geom_point(size=0.5) + 
  scale_color_gradient2(low="red",mid="grey",high="blue",midpoint=0.5,name="residual\nquantile",limits=c(0,1)) + 
  scale_x_continuous(labels=function(x){x/1000000},name="Window position (Mb)") + 
  ylab("Ancestry frequency") + 
  theme(
    strip.text.y = element_text(angle = 0)
  )
ggsave(file.path(thisfolder,"Paper_figs/new_fig6_GOIs_nonas_withquantiles_zoom_lines_split2.pdf"),split2,device="pdf",width = 10,height=5,units="in",dpi = 600)  




