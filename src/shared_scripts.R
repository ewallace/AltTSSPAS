# shared_scripts.R

library(tidyverse)
library(MASS)
select <- dplyr::select
library(broom)
library(cowplot)
library(pheatmap)
library(qvalue)
theme_set(theme_cowplot() + theme(strip.background=element_blank()))

## figure helpers
mybreaks = (1:9 * 10^rep(-1:5,each=9))
mylabels = mybreaks
mylabels[ c( c(3,4,6,7,8,9),c(3,4,6,7,8,9)+9,
             c(3,4,6,7,8,9)+18,c(3,4,6,7,8,9)+27,
             c(3,4,6,7,8,9)+36,c(3,4,6,7,8,9)+45,c(3,4,6,7,8,9)+54 ) ] <- ""

scale_y_log10nicer <- function(name = waiver(), 
                               breaks=mybreaks, 
                               labels=mylabels,...) {
    scale_y_log10(name=name,breaks=breaks,labels=labels,...)
}

## plotting functions ofr alternate counts

plot_altcountC <- function(df,Gname) {
    ## Plot alternate counts by Condition
    ggplot(data=filter(df,Gene_name==Gname),
           aes(x=Condition,y=Count,
               shape=Rep,colour=factor(Site))) +
        geom_point(position=position_dodge(width = 0.4)) + 
        scale_colour_brewer("Site",palette = "Set1") +
        coord_cartesian(expand=FALSE,clip="off") + 
        expand_limits(y=0,x=0.5) +
        labs(title=Gname) + 
        theme(axis.text.x=element_text(angle=90,vjust=0.5),
              axis.title.x=element_blank())
}

plot_altcountS <- function(df,Gname,samplenames=JEC21_TSS_samples_norep) {
    ## Plot alternate counts by Sample (Condition x Strain)
    ggplot(data=filter(df,Gene_name==Gname) %>%
               unite(Sample,Strain,Condition) %>%
               mutate(Sample=factor(Sample,
                                    levels=samplenames)),
           aes(x=Sample,y=Count,
               shape=Rep,colour=factor(Site))) +
        geom_point(position=position_dodge(width = 0.4)) + 
        scale_colour_brewer("Site",palette = "Set1") +
        coord_cartesian(expand=FALSE,clip="off") + 
        expand_limits(y=0,x=0.5) +
        labs(title=Gname) + 
        theme(axis.text.x=element_text(angle=90,vjust=0.5),
              axis.title.x=element_blank())
}

summarize_counts_bygene <- function(countdf,samplenames) {
    ## Summarize the counts and locations by gene across multiple sites
    countdf %>%
        mutate(Tot=rowSums(select(.,samplenames))) %>%
        group_by(Gene_name) %>%
        summarize(nSites=n_distinct(Site),
                  nLocations=n_distinct(location),
                  Tot=sum(Tot),
                  cdsAlt=("cdsE" %in% location | "cdsI" %in% location),
                  utr5Altx=sum(str_count(location,"utr5") ) >= 2,
                  utr3Altx=sum(str_count(location,"utr3") ) >= 2
        )
}



read_mitofates <- function(mffile,              
                           mitofates_colnames = 
                               c("Gene_name", "Prob_preseq", "Pred_preseq")
) {
    # Read mitofates output
    read_tsv(mffile,comment="#",skip=1,
             col_names=mitofates_colnames,
             col_types = "ccc") %>%
        mutate(Pred_preseq = (Pred_preseq == "Possessing mitochondrial presequence"))
}
