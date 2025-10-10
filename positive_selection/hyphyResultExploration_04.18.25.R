
# BiocManager::install("ggtree")

library(tidyverse)
library(ggtree)
library(viridis)

setwd('/Volumes/WorkingBWP/WorkingBWP/bearComparativeGenomics_Spring2025/hyphy')

allsig <- read_tsv('absrelResults_04.26.2025/allsig_absrel_ap26.txt',col_names = c('all')) %>% 
  mutate(node = str_split_fixed(all,' ',3)[,2] %>% str_remove_all(',')) %>% 
  mutate(ogid = str_split_fixed(all,' ',3)[,1] %>% str_remove_all('all_absrel_logs/|.absrel.all.out.txt-[*]'))

allsig %>% group_by(node) %>% tally() %>% arrange(-n)

length(unique(allsig$ogid))

allsig.summary <- allsig %>% group_by(node) %>% tally() %>% arrange(-n)

bear_tree <- read.tree('custom_BL_rooted.unroot.nodelab.nwk')

bear_tree <- left_join(bear_tree,allsig.summary,by=c('label'='node'))
root_bear_tree <- treeio::root(bear_tree, outgroup = "panda", edgelabel = TRUE)

ggtree(bear_tree,aes(color=n),size=2) +
  ggrepel::geom_label_repel(aes(x=branch, label=n), fill='white',color='black',nudge_y = 0.3,nudge_x = -0.0003) + 
  geom_tiplab(size=5,color='black',offset = 0.0001) +
  labs(color='# of genes')+
  scale_color_viridis()

p0a <- ggtree(bear_tree,aes(color=n),size=1.5,layout='roundrect',show.legend=F) +
  geom_nodelab(size=3,color='black',nudge_x = 0.0008) +
  geom_tiplab(size=4,color='black',offset = 0.0001) +
  ggrepel::geom_label_repel(aes(x=branch, label=n), fill='white',color='black',nudge_y = 0.3,nudge_x = -0.0003) +
  labs(color='# of genes')+
  colorspace::scale_color_continuous_sequential('Rocket',rev=F,end = .8) 
p0a

ggtree(bear_tree) +
  geom_nodelab(size=5,color='black') +
  scale_color_viridis()

ggtree(bear_tree) + geom_text(aes(label=node))
flip(tree_view = p0a,node1=8,node2=7)

ggtree(bear_tree,aes(color=n),size=1.5,layout='roundrect',show.legend=F) +
  geom_nodelab(size=3,color='black',nudge_x = 0.0008) +
  geom_tiplab(size=4,color='black',offset = 0.0001) +
  # ggrepel::geom_label_repel(aes(x=branch, label=n), fill='white',color='black',nudge_y = 0.3,nudge_x = -0.0003) +
  labs(color='# of genes')+
  colorspace::scale_color_continuous_sequential('Rocket',rev=F,end = .8) 

# GO analysis -------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(AnnotationHub)

  
# Load brown bear database for GO analysis
hub <- AnnotationHub()
# query(hub, "Ursus arctos")
bear_annot <- hub[["AH118004"]]


# Prep background and gene ID info ----------------------------------------

singlecopy <- read_tsv('orthofinder_oldBrownAllNew_apr22/Orthogroups_SingleCopyOrthologues.txt',col_names = 'ogid') 

brown_og <- read_tsv('orthofinder_oldBrownAllNew_apr22/brown.existing.longestIso.prots.tsv') %>%
  dplyr::select(ogid=Orthogroup,id = brown.existing.longestIso.prots) %>%
  filter(ogid %in% singlecopy$ogid) %>%
  mutate(id = str_remove_all(id,'rna-')) %>% 
  unique() 
  
allsig.brownID <- allsig %>% left_join(brown_og)

# write_tsv(allsig.brownID,'allAbsrelSig_OrthoRefSeq.txt')

unique(allsig.brownID$node)

# Make named list of significant genes
all.sig <- list('brown'=as.list(allsig.brownID[which(allsig.brownID$node=='brown'),'id'])$id,
                'black'=as.list(allsig.brownID[which(allsig.brownID$node=='black'),'id'])$id,
                'polar'=as.list(allsig.brownID[which(allsig.brownID$node=='polar'),'id'])$id,
                'panda'=as.list(allsig.brownID[which(allsig.brownID$node=='panda'),'id'])$id,
                'asiatic'=as.list(allsig.brownID[which(allsig.brownID$node=='asiatic'),'id'])$id,
                'spectacled'=as.list(allsig.brownID[which(allsig.brownID$node=='spectacled'),'id'])$id,
                'sloth'=as.list(allsig.brownID[which(allsig.brownID$node=='sloth'),'id'])$id,
                'sun'=as.list(allsig.brownID[which(allsig.brownID$node=='sun'),'id'])$id,
                'Node8'=as.list(allsig.brownID[which(allsig.brownID$node=='Node8'),'id'])$id,
                'Node2'=as.list(allsig.brownID[which(allsig.brownID$node=='Node2'),'id'])$id,
                'Node10'=as.list(allsig.brownID[which(allsig.brownID$node=='Node10'),'id'])$id,
                'Node5'=as.list(allsig.brownID[which(allsig.brownID$node=='Node5'),'id'])$id,
                'Node3'=as.list(allsig.brownID[which(allsig.brownID$node=='Node3'),'id'])$id
)


# Turn off universe intersection. See: https://github.com/YuLab-SMU/clusterProfiler/issues/636
options(enrichment_force_universe = TRUE)


# Use compare cluster to run analysis on all at once
allbear_compGO.bp <- compareCluster(geneClusters = all.sig, 
                                 fun='enrichGO',
                                 universe = brown_og$id,
                                 OrgDb = bear_annot,
                                 ont = 'BP',
                                 keyType='REFSEQ',
                                 minGSSize=4,
                                 pAdjustMethod = 'BH',
                                 pvalueCutoff = 0.05)  %>% 
  simplify(cutoff=0.7, by="p.adjust", select_fun=min)

go.bp <- dotplot(allbear_compGO.bp,showCategory=10,) +
  colorspace::scale_fill_continuous_sequential('dark mint',rev=F) +
  xlab('') +
  ggtitle('Biological Process')+
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) + 
  theme_linedraw(base_size = 8)

allbear_compGO.mf <- compareCluster(geneClusters = all.sig, 
                                    fun='enrichGO',
                                    universe = brown_og$id,
                                    OrgDb = bear_annot,
                                    ont = 'MF',
                                    keyType='REFSEQ',
                                    minGSSize=4,
                                    pAdjustMethod = 'BH',
                                    pvalueCutoff = 0.05)  %>% 
  simplify(cutoff=0.7, by="p.adjust", select_fun=min)

go.mf <- dotplot(allbear_compGO.mf,showCategory=10,) +
  colorspace::scale_fill_continuous_sequential('dark mint',rev=F) +
  xlab('') +
  ggtitle('Molecular Function')+
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) + 
  theme_linedraw(base_size = 8)

allbear_compGO.cc <- compareCluster(geneClusters = all.sig, 
                                    fun='enrichGO',
                                    universe = brown_og$id,
                                    OrgDb = bear_annot,
                                    ont = 'CC',
                                    keyType='REFSEQ',
                                    minGSSize=4,
                                    pAdjustMethod = 'BH',
                                    pvalueCutoff = 0.05)  %>% 
  simplify(cutoff=0.7, by="p.adjust", select_fun=min)

go.cc <- dotplot(allbear_compGO.cc,showCategory=10,) +
  colorspace::scale_fill_continuous_sequential('dark mint',rev=F) +
  xlab('') +
  ggtitle('Cellular Component')+
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) + 
  theme_linedraw(base_size = 8)

go.bp + go.mf + go.cc

# allbear_compKEGG <- compareCluster(geneClusters = all.sig,
#                                    fun='enrichKEGG',
#                                    universe = brown_og$id,
#                                    organism     = 'uah',
#                                    keyType = 'REFSEQ',
#                                    pvalueCutoff = 0.05)

allbear_compKEG


# enrichGO(gene=all.sig$brown,
#          universe =  brown_og$id,
#          OrgDb = bear_annot,
#          ont = 'BP',
#          keyType = 'REFSEQ',
#          minGSSize = 2,
#          pAdjustMethod = 'BH',
#          pvalueCutoff = 0.05,
#          readable = T)

#


temp.bp <- as.data.frame(allbear_compGO.bp) %>% mutate(ontology_db = 'Biological Process')
temp.cc <- as.data.frame(allbear_compGO.cc) %>% mutate(ontology_db = 'Cellular Component')
temp.mf <- as.data.frame(allbear_compGO.mf) %>% mutate(ontology_db = 'Molecular Function')

allbear_compGO_allResTable <- temp.bp %>% bind_rows(temp.cc,temp.mf) %>% mutate(GeneRatio = paste(' ',GeneRatio,sep=''))

# write_csv(allbear_compGO_allResTable,'supp_tables/SuppTableS7_aBSREL_GOResults_05.26.25.csv')




# Read in candidate genes -------------------------------------------------
library(patchwork)

cands <- read_tsv('../candidate_genes/candidateGenes_bearCompGen_04.29.25.tsv')

cands %>% group_by(term) %>% tally()

# Filter and tally to those present in single copy ortholog set
singlcopy_symbols <- bitr(brown_og$id,fromType = 'REFSEQ',toType = 'SYMBOL',OrgDb = bear_annot)
sc_cands <- cands %>% filter(symbol %in% singlcopy_symbols$SYMBOL)

sc_cands %>% group_by(term) %>% tally()

#### Absrel Results
all.absrel.brownID <- enframe(all.sig) %>% unnest(cols = value)
absrel_temp.symbol <-  bitr(all.absrel.brownID$value,fromType = 'REFSEQ',toType = 'SYMBOL',OrgDb = bear_annot)

all.absrel.brownID.symbol <- all.absrel.brownID %>% left_join(absrel_temp.symbol,by = c('value'='REFSEQ'))

# write_tsv(all.absrel.brownID.symbol,'allAbsrelSig_RefSeqSymbol.txt')

# Testing enrichment analysis of candidate genes
test_panda <- all.absrel.brownID.symbol %>% filter(name=='sloth')
test_cand <- sc_cands %>% dplyr::select(term,geneID=symbol)
test <- enricher(test_panda$SYMBOL, 
         TERM2GENE = test_cand,
         universe = singlcopy_symbols$SYMBOL,
         pvalueCutoff = 0.05)
barplot(test)
dotplot(test)

all.absrel.sig <- list('brown'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='brown'),'SYMBOL'])$SYMBOL,
                'black'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='black'),'SYMBOL'])$SYMBOL,
                'polar'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='polar'),'SYMBOL'])$SYMBOL,
                'panda'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='panda'),'SYMBOL'])$SYMBOL,
                'asiatic'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='asiatic'),'SYMBOL'])$SYMBOL,
                'spectacled'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='spectacled'),'SYMBOL'])$SYMBOL,
                'sloth'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='sloth'),'SYMBOL'])$SYMBOL,
                'sun'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='sun'),'SYMBOL'])$SYMBOL,
                'Node8'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='Node8'),'SYMBOL'])$SYMBOL,
                'Node2'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='Node2'),'SYMBOL'])$SYMBOL,
                'Node10'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='Node10'),'SYMBOL'])$SYMBOL,
                'Node5'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='Node5'),'SYMBOL'])$SYMBOL,
                'Node3'=as.list(all.absrel.brownID.symbol[which(all.absrel.brownID.symbol$name=='Node3'),'SYMBOL'])$SYMBOL
)



all_candGO <- compareCluster(geneClusters = all.absrel.sig, 
               fun='enricher',
               universe = singlcopy_symbols$SYMBOL,
               TERM2GENE = test_cand,
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05)

temp <- as.data.frame(all_candGO)

dotplot(all_candGO)

# Simpler plotting

all.absrel.brownID.symbol.cand <- all.absrel.brownID.symbol %>% 
  filter(SYMBOL %in% cands$symbol) %>% left_join(cands,by=c('SYMBOL'='symbol')) %>% 
  janitor::clean_names() %>% 
  mutate(term = factor(term,levels=rev(c('lipid metabolic process','insulin receptor signaling pathway','circadian rhythm','digestion','sensory perception of taste','pigmentation',
                                         'skeletal system development','cranial skeletal system development','adipose tissue development','muscle organ development')))) %>% 
  mutate(group = case_when(
    str_detect(term,'lipid|insulin|circadian') ~ 'Metabolism and Hibernation',
    str_detect(term,'digestion|taste') ~ 'Diet',
    str_detect(term,'pigment') ~ 'Coloration',
    str_detect(term,'skeletal|adipose|muscle') ~ 'Anatomical structure and development',
  )) %>% 
  mutate(name=factor(name,levels=c('spectacled','panda','asiatic','Node3','sloth','Node5','sun','Node2','black','Node8','brown','Node10','polar')))

# write_tsv(all.absrel.brownID.symbol.cand,'aBSREL_CandidateGeneResults_05.01.25.tsv')


ggplot(all.absrel.brownID.symbol.cand,aes(y=term,x=symbol,color=group)) +
  geom_point(size=4) +
  # facet_wrap(~name,ncol=1)+ # when commented out, shows genes significant in at least one or more tests
  labs(x='Gene Symbol',y='GO Term',color='Candidate Gene Category')+
  ggtitle('Candidate genes under positive selection (BUSTED)') +
  scale_color_manual(values = c('Metabolism and Hibernation'='darkblue','Diet'='seagreen','Coloration'='black','Anatomical structure and development'='goldenrod'))+
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),plot.title = element_text(face='bold'),plot.title.position = 'plot')


ggplot(all.absrel.brownID.symbol.cand,aes(y=name,x=symbol,color=group)) +
  geom_point(size=4) +
  # facet_wrap(~group,scales = 'free_y')+ # when commented out, shows genes significant in at least one or more tests
  facet_grid(cols=vars(term),scales = 'free_x',space = 'free_x')+
  labs(y='Gene Symbol',x='Node',color='Candidate Gene Category')+
  ggtitle('Candidate genes under positive selection (aBSREL)') +
  scale_color_manual(values = c('Metabolism and Hibernation'='darkblue','Diet'='seagreen','Coloration'='black','Anatomical structure and development'='goldenrod'))+
  theme_linedraw() + theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),plot.title = element_text(face='bold'),plot.title.position = 'plot')

# Gene count summary
p0b <- all.absrel.brownID.symbol.cand %>% 
  group_by(name,term,group) %>% 
  tally() %>% 
  ggplot(aes(x=term,y=name,size=n,color=group)) +
  geom_point() +
  geom_text(aes(label=n),color='white',size=3)+
  scale_size_continuous(range = c(1,12),limits = c(0,25))+
  labs(x='GO Term',y='Node',color='Candidate Gene Category',size='# of Genes')+
  scale_color_manual(values = c('Metabolism and Hibernation'='darkblue','Diet'='seagreen','Coloration'='black','Anatomical structure and development'='goldenrod'))+
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
p0b

p0a + p0b + plot_layout(widths = c(2,1))


### Full plotting

p1 <- ggplot(subset(all.absrel.brownID.symbol.cand,str_detect(group,'structure|Color|Diet')),aes(y=name,x=symbol,color=group)) +
  geom_point(size=4,show.legend = T) +
  # facet_wrap(~group,scales = 'free_y')+ # when commented out, shows genes significant in at least one or more tests
  facet_grid(cols=vars(term),scales = 'free_x',space = 'free_x')+
  labs(y='Gene Symbol',x='',color='Candidate Gene Category')+
  ggtitle('Candidate genes under positive selection (aBSREL)') +
  scale_color_manual(values = c('Metabolism and Hibernation'='darkblue','Diet'='seagreen','Coloration'='black','Anatomical structure and development'='goldenrod'))+
  theme_linedraw() + theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),plot.title = element_text(face='bold'),plot.title.position = 'plot')

p2 <- ggplot(subset(all.absrel.brownID.symbol.cand,str_detect(group,'Metabolism')),aes(y=name,x=symbol,color=group)) +
  geom_point(size=4,show.legend = T) +
  # facet_wrap(~group,scales = 'free_y')+ # when commented out, shows genes significant in at least one or more tests
  facet_grid(cols=vars(term),scales = 'free_x',space = 'free_x')+
  labs(y='Gene Symbol',x='Node',color='Candidate Gene Category')+
  # ggtitle('Candidate genes under positive selection (aBSREL)') +
  scale_color_manual(values = c('Metabolism and Hibernation'='darkblue','Diet'='seagreen','Coloration'='black','Anatomical structure and development'='goldenrod'))+
  theme_linedraw() + theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),plot.title = element_text(face='bold'),plot.title.position = 'plot')

p1 + p2 + plot_layout(ncol = 1,guides = 'collect')

