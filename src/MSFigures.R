#!/usr/bin/env Rscript

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

# Figure 2AB
{
  fig2a_df <- readxl::read_excel('data/SupplementaryTables.xlsx', sheet = 'Table S3') |>
    as.data.frame() |>
    filter(type == 'plasma')
  
  comparison1 <- list(c("pretreatment", "primary treatment"), 
                      c("pretreatment", "follow-up"),
                      c("primary treatment", "relapse"),
                      c("relapse", "follow-up") )
  
  fig2a <- fig2a_df %>% 
    mutate(label = ifelse(TP53_VAF > .3, eoc_code, '')) %>%
    ggplot(aes(x = factor(treatmentPhase, levels = c('pretreatment', 'primary treatment','relapse',  'follow-up')), 
               y = TP53_VAF)) + #ggboxplot(., x='stage', y = 'TP53_VAF') +
    geom_boxplot(aes(fill=treatmentPhase), alpha=.5, show.legend = F, 
                 outlier.shape = NA, width=.7) + 
    geom_jitter() + 
    ggrepel::geom_text_repel(aes(label=label), vjust=.1)+
    scale_fill_brewer(palette="BuGn") +
    labs(x='', y='TP53 VAF') +
    theme_pubr() +
    theme(title = element_text(size = 15), 
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.position = c(.85, .6)) +
    stat_compare_means(aes(label=..p.adj..), comparisons=comparison1) + 
    stat_compare_means(label.y = .8, size=5)
  
  fragment.df <- read.delim('data/DNAFragmentData.tsv')
  
  comparison2 <- list(c("pretreatment", "follow-up"),
                      c("primary treatment", "relapse"),
                      c("relapse", "follow-up") )
  
  fig2b <- fragment.df %>% 
    mutate(label = ifelse(median < 150 & TP53_VAF > .4, patientID, '')) %>%
    ggplot(aes(x = factor(treatmentPhase, levels = c('pretreatment', 'primary treatment','relapse',  'follow-up')), 
               y = median)) + 
    geom_boxplot(aes(fill=treatmentPhase), alpha=.5, 
                 show.legend = F, outlier.alpha = 0,
                 width = .7) + 
    geom_jitter(aes(size=TP53_VAF), alpha=.6) +
    ggrepel::geom_text_repel(aes(label=label), vjust=.1)+
    scale_fill_brewer(palette="PuBu") +
    labs(x='', y='Median DNA fragment size', size='ctDNA fraction') +
    theme_pubr() +
    theme(title = element_text(size = 15), 
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.position = c(.25, .2)) +
    stat_compare_means(aes(label=..p.adj..), comparisons=comparison2) + 
    stat_compare_means(label.y = 181, size=5)
  
  Fig2 <- ggarrange(fig2a, fig2b, nrow=1, ncol=2)
  
  #svg('results/figures/Fig2.svg', width = 15, height = 7)
  #Fig2
  #dev.off()
}

# Figure 2CD
{
  fig2c_df <- readxl::read_excel('data/SupplementaryTables.xlsx', sheet = 'Table S4') |>
    as.data.frame()
  
  fig2c <- fig2c_df %>% 
    mutate(gr = ifelse(ctDNA_TP53 <= .05, 'TP53 VAF: 1-5%', 
                       ifelse(ctDNA_TP53 <= .1, 'TP53 VAF: 5-10%', 'TP53 VAF: >10%')), 
           significant=ifelse(p_val <= .05, 'y', 'n')) %>% 
    ggplot(aes(x = factor(gr, levels = c('TP53 VAF: 1-5%', 'TP53 VAF: 5-10%', 'TP53 VAF: >10%')), 
               y = Kendall_correlation)) +
    geom_boxplot(aes(fill = gr), show.legend = F) + 
    scale_fill_brewer(palette = 'Set2') +
    geom_jitter(size=2, color='black',alpha=.5)+
    ylim(0,1) + labs(x = '', y ='Kendall Correlation') + 
    theme(title = element_text(size = 15), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.position = c(.8,.2))
  
  svg('results/figures/fig2c.svg', width = 5, height = 5)
  fig2c; dev.off()
  
  
  burden <- read.delim('data/CNASignalBurden.tsv')
  
  fig2d <- burden %>% 
    ggpubr::ggscatter(., x = 'TP53_VAF', y = 'FCS', cor.coef = T, cor.coef.size = 10, 
                      add = 'reg.line', size = .6, color = 'sample_type') + 
    labs(x = 'TP53 VAF', y = 'Focal Copy Number Score', color = 'Type') + 
    theme(title = element_text(size = 15), 
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size=10), 
          legend.title = element_text(size=15),
          legend.position = c(.75, .2))
  
  Fig2 <- ggarrange(fig2a, fig2b, fig2c, fig2d, nrow = 2, ncol = 2) 
            labels = 'AUTO', font.label = list(size=30))
  #svg('results/figures/Fig2.svg', width = 15, height = 15)
  #Fig2
  #dev.off()
}

# Figure 3
#-- based on gistic analysis
{
  sampleData <- readxl::read_excel('data/SupplementaryTables.xlsx', sheet = 'Table S3') |>
    as.data.frame() |>
    mutate(sampleID =  paste(patient, sampleID, sep = '_'))
  
  tis_sel <- (sampleData %>%
    filter(type == 'tissue', TP53_VAF > .0, treatmentPhase=='pretreatment'))$sampleID
  
  prim_tis_gistic <- read.delim('data/gistic_pretreatment_tissue.tsv') %>%
    filter(Gene.Symbol %in% c('MECOM', 'MYC', 'KRAS', 'CCNE1')) %>%
    dplyr::select(-Gene.ID) %>%
    tibble::column_to_rownames(var='Cytoband') %>%
    tidyr::gather(key=sample, value='CN_gistic', 2:ncol(.)) %>%
    mutate(CN = 2^CN_gistic*2,
           CNA = ifelse(CN_gistic >= .5, 'amp',
                    ifelse(CN_gistic < -.5, 'del', 'nor')))
  
  prim_tis_gistic_CNA <- prim_tis_gistic %>%
    filter(sample %in% tis_sel) %>%
    mutate(patient = sub('_.*', '', sample)) %>%
    dplyr::select(patient, gene = Gene.Symbol, CNA) %>%
    dplyr::group_by(gene, CNA) %>%
    summarise(N = n()) %>%
    as.data.frame()
  
  fig3a <- prim_tis_gistic_CNA %>%
    dplyr::group_by(gene) %>%
    mutate(pct=N/sum(N)) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    ggplot(aes(x=factor(gene, 
                        levels = c('MECOM', 'MYC', 'KRAS', 'CCNE1'),
                        labels = c('MECOM', 'MYC', 'KRAS', 'CCNE1')),
               y=N, fill=factor(CNA,
                                levels = c('nor', 'del', 'amp'),
                                labels = c('normal', 'loss', 'gain')))) +
    geom_bar(position = 'fill', stat='identity') +
    geom_text(aes(label = N, y = pct), size=8,
              show.legend = F, position = position_stack(vjust =0.5)) +
    scale_fill_manual(values = c("#dddddd", "#72a0e5", "#ff7f7f"))  +
    labs(x='', y = 'Proportion', fill='Status') +
    theme_pubr() +
    theme(legend.position = 'bottom',
          title = element_text(size=20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          legend.text=element_text(size=15))
  #svg('results/figures/fig3a.svg', height = 12, width = 8)
  #fig3a
  #dev.off()
  
}  
  
  
{
  # rescale for gene visualization
  segment_res <- read.delim('data/segment_result.tsv') |>
    left_join(oseq_samples |>
                dplyr::select(sample = sampleID, patientID = eoc_code)) |>
    mutate(sample = paste(patientID, sub('.*_', '', sample), sep = '_')) 
  
  # gene
  ensembl_genes_v96 <- read.delim('data/ensembl_genes_v96.csv') %>%
    filter(type == 'protein_coding') %>%
    dplyr::select(chr, start, end, gene=Gene, gene_id = ID, band) %>%
    filter(gene %in% c('MECOM', 'MYC', 'KRAS', 'CCNE1'))
  
  gene_gr <- with(ensembl_genes_v96, GRanges(chr, IRanges(start, end)))
  
  rescale_segment_gr <- with(rescale_segment, GRanges(chr, IRanges(startpos, endpos)))
  overlaps <- findOverlaps(rescale_segment_gr, gene_gr)
  
  df <- data.frame(
    rescale_segment[from(overlaps), ] %>%
      dplyr::select(-chr, -startpos, -endpos),
    ensembl_genes_v96[to(overlaps), ]
  )  
  
  rescale_segment <- segment_res %>% 
    mutate(sign = ifelse(logR > 0, 'pos', ifelse(logR < 0, 'neg', 'zero'))) %>% 
    dplyr::group_by(sample, sign) %>% 
    mutate(q3 = quantile(logR, .75),
           q1 = quantile(logR, .25)) %>% 
    dplyr::ungroup() %>%
    as.data.frame() %>%
    mutate(logR_scale = ifelse(logR >= 0, 
                               logR/q3, logR/(-q1)))
  
  
  samples4fig4 <- sampleData %>%
    filter(treatmentPhase %in% c('pretreatment', 'relapse'),
           TP53_VAF >= 0.05,
           !grepl('_r', sampleID)) %>%
    dplyr::group_by(patient, treatmentPhase, type) %>%
    dplyr::slice_max(TP53_VAF) %>%
    dplyr::group_by(patient) %>%
    mutate(N = n()) %>%
    dplyr::ungroup() %>%
    as.data.frame()  %>%
    filter(N>1)
  
  compare_df <- df %>%
    inner_join(samples4fig4 |>
                 dplyr::rename(sample=sampleID)) %>%
    dplyr::group_by( gene, type, patientID, TP53_VAF, treatmentPhase, N) %>% 
    summarise(logR=mean(logR_scale)) %>% 
    as.data.frame() %>% 
    mutate(var = paste0(type, '_', treatmentPhase)) %>%
    select(-type, -TP53_VAF, -treatmentPhase) %>% 
    tidyr::spread(., key=var, value=logR) %>%
    tidyr::gather(., phase, logR, 4:6)
  
  fig3b <- compare_df %>%
    ggplot(aes(factor(phase, levels = c('tissue_pretreatment', 'plasma_pretreatment', 'plasma_relapse'),
                      labels = c('Pretreatment\nTissue', 'Pretreatment\nPlasma', 'Relapse\nPlasma')), 
               factor(patientID,
                      levels = c('EOC912', 'EOC740', 'EOC736', 'EOC587', 'EOC482', 
                                 'EOC204', 'EOC415', 'EOC1067', 'EOC295', 'EOC198', 
                                 'EOC183', 'EOC165', 'EOC568', 'EOC1005', 'EOC1099',
                                 'EOC998', 'EOC989', 'EOC49','EOC426', 'EOC429',
                                 'EOC1030', 'EOC172', 'EOC677', 'EOC742', 'EOC3', 
                                 'EOC412', 'EOC321', 'EOC1120')), 
               fill=logR)) + 
    geom_tile() + 
    facet_wrap( ~ factor(gene, levels = c('MECOM', 'MYC', 'KRAS', 'CCNE1')), ncol=4) + 
    scale_fill_gradientn( colours = c("#3399ff", "white", "#ffc3a0", "#ff4040", "#800000"), 
                          na.value = 'grey', values = scales::rescale(c(-1.4, -.5, 0, 1, 3, 13))) + 
    labs(x='', y = 'patient', fill = 'scaled logR') +
    theme(title = element_text(size = 25),
          legend.position='bottom',
          legend.text = element_text(size=10, color='black'),
          strip.text.x = element_text(size=20, color='black'),
          axis.text.y = element_text(size = 15, color='black'),
          axis.text.x = element_text(size = 15, color='black'))
  
  #svg('results/figures/main/fig3b.svg', height = 14, width = 24)
  #fig3b
  #dev.off()
  

}
