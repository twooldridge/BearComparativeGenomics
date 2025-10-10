
library(tidyverse)
library(Biostrings)

# Read in list of single copy orthogroups
singles <- read_tsv('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/Orthogroups_SingleCopyOrthologues.txt',col_names = 'id')

# Read in all orthologs for existing brown bear annotation
brown_ortho <- read_tsv('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/brown.existing.longestIso.prots.tsv')

# Filter to retain ortholog info for only single-copy orthogroups
brown_ortho_filtered <- brown_ortho %>% 
  filter(Orthogroup %in% singles$id) %>% 
  mutate(Species = str_split_fixed(Species,'[.]',2)[,1])

# Make a table for each species with orthogroup, species, brown bear gene ID, egap gene ID
brown.ogs <- brown_ortho_filtered %>% filter(Species=='brown')  
black.ogs <- brown_ortho_filtered %>% filter(Species=='black')  
polar.ogs <- brown_ortho_filtered %>% filter(Species=='polar')  
panda.ogs <- brown_ortho_filtered %>% filter(Species=='panda')  
sun.ogs <- brown_ortho_filtered %>% filter(Species=='sun')  
spectacled.ogs <- brown_ortho_filtered %>% filter(Species=='spec')  %>% mutate(Species='spectacled')
sloth.ogs <- brown_ortho_filtered %>% filter(Species=='sloth')  
asiatic.ogs <- brown_ortho_filtered %>% filter(Species=='asiatic')  



# Read in longest isoform seqs --------------------------------------------
# For each species, read in longest isoform CDS fasta and filter/modify to retain single copy orthologs only

# Read in cds fasta
brown.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/brown.egapx.longestIso.cds.fna')  

# Clean up egap IDs
names(brown.iso) <- str_split_fixed(names(brown.iso),' ',2)[,1]
names(brown.iso) <- str_replace_all(names(brown.iso),':','_')

# Make table with egap IDs and Orthogroup IDs
brown.ids <- data.frame(row.names = brown.ogs$Orthologs,val=brown.ogs$Orthogroup)

# Filter to retain only single copy orthogroups
brown.iso.filtered <- brown.iso[names(brown.iso) %in% brown.ogs$Orthologs]

# Save as new table and prepend orthogroup ID with species id
brown.iso.filtered.og <- brown.iso.filtered
names(brown.iso.filtered.og) <- paste('brown_',brown.ids[names(brown.iso.filtered.og),],sep='') # rename with species orthogroup IDs
head(names(brown.iso.filtered.og))

# Repeat the process for each other species - probably will want to make this into a loop!
# black.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/black.egapx.longestIso.cds.fna')  
# names(black.iso) <- str_split_fixed(names(black.iso),' ',2)[,1]
# names(black.iso) <- str_replace_all(names(black.iso),':','_')
# black.ids <- data.frame(row.names = black.ogs$Orthologs,val=black.ogs$Orthogroup)
# black.iso.filtered <- black.iso[names(black.iso) %in% black.ogs$Orthologs]
# black.iso.filtered.og <- black.iso.filtered
# names(black.iso.filtered.og) <- paste('black_',black.ids[names(black.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(black.iso.filtered.og))
# 
# polar.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/polar.egapx.longestIso.cds.fna')  
# names(polar.iso) <- str_split_fixed(names(polar.iso),' ',2)[,1]
# names(polar.iso) <- str_replace_all(names(polar.iso),':','_')
# polar.ids <- data.frame(row.names = polar.ogs$Orthologs,val=polar.ogs$Orthogroup)
# polar.iso.filtered <- polar.iso[names(polar.iso) %in% polar.ogs$Orthologs]
# polar.iso.filtered.og <- polar.iso.filtered
# names(polar.iso.filtered.og) <- paste('polar_',polar.ids[names(polar.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(polar.iso.filtered.og))
# 
# panda.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/panda.egapx.longestIso.cds.fna')  
# names(panda.iso) <- str_split_fixed(names(panda.iso),' ',2)[,1]
# names(panda.iso) <- str_replace_all(names(panda.iso),':','_')
# panda.ids <- data.frame(row.names = panda.ogs$Orthologs,val=panda.ogs$Orthogroup)
# panda.iso.filtered <- panda.iso[names(panda.iso) %in% panda.ogs$Orthologs]
# panda.iso.filtered.og <- panda.iso.filtered
# names(panda.iso.filtered.og) <- paste('panda_',panda.ids[names(panda.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(panda.iso.filtered.og))
# 
# sun.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/sun.egapx.longestIso.cds.fna')  
# names(sun.iso) <- str_split_fixed(names(sun.iso),' ',2)[,1]
# names(sun.iso) <- str_replace_all(names(sun.iso),':','_')
# sun.ids <- data.frame(row.names = sun.ogs$Orthologs,val=sun.ogs$Orthogroup)
# sun.iso.filtered <- sun.iso[names(sun.iso) %in% sun.ogs$Orthologs]
# sun.iso.filtered.og <- sun.iso.filtered
# names(sun.iso.filtered.og) <- paste('sun_',sun.ids[names(sun.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(sun.iso.filtered.og))
# 
# spectacled.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/spec.egapx.longestIso.cds.fna')  
# names(spectacled.iso) <- str_split_fixed(names(spectacled.iso),' ',2)[,1]
# names(spectacled.iso) <- str_replace_all(names(spectacled.iso),':','_')
# spectacled.ids <- data.frame(row.names = spectacled.ogs$Orthologs,val=spectacled.ogs$Orthogroup)
# spectacled.iso.filtered <- spectacled.iso[names(spectacled.iso) %in% spectacled.ogs$Orthologs]
# spectacled.iso.filtered.og <- spectacled.iso.filtered
# names(spectacled.iso.filtered.og) <- paste('spectacled_',spectacled.ids[names(spectacled.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(spectacled.iso.filtered.og))
# 
# sloth.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/sloth.egapx.longestIso.cds.fna')  
# names(sloth.iso) <- str_split_fixed(names(sloth.iso),' ',2)[,1]
# names(sloth.iso) <- str_replace_all(names(sloth.iso),':','_')
# sloth.ids <- data.frame(row.names = sloth.ogs$Orthologs,val=sloth.ogs$Orthogroup)
# sloth.iso.filtered <- sloth.iso[names(sloth.iso) %in% sloth.ogs$Orthologs]
# sloth.iso.filtered.og <- sloth.iso.filtered
# names(sloth.iso.filtered.og) <- paste('sloth_',sloth.ids[names(sloth.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(sloth.iso.filtered.og))
# 
# asiatic.iso <- readAAStringSet('/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/cds/asiatic.egapx.longestIso.cds.fna')  
# names(asiatic.iso) <- str_split_fixed(names(asiatic.iso),' ',2)[,1]
# names(asiatic.iso) <- str_replace_all(names(asiatic.iso),':','_')
# asiatic.ids <- data.frame(row.names = asiatic.ogs$Orthologs,val=asiatic.ogs$Orthogroup)
# asiatic.iso.filtered <- asiatic.iso[names(asiatic.iso) %in% asiatic.ogs$Orthologs]
# asiatic.iso.filtered.og <- asiatic.iso.filtered
# names(asiatic.iso.filtered.og) <- paste('asiatic_',asiatic.ids[names(asiatic.iso.filtered.og),],sep='') # rename with species orthogroup IDs
# head(names(asiatic.iso.filtered.og))


# Export fasta with single copy orthogroup CDS's with simple headers
writeXStringSet(brown.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/brown.sc.ogids.apr23.fa')
# writeXStringSet(black.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/black.sc.ogids.apr23.fa')
# writeXStringSet(polar.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/polar.sc.ogids.apr23.fa')
# writeXStringSet(panda.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/panda.sc.ogids.apr23.fa')
# writeXStringSet(sun.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/sun.sc.ogids.apr23.fa')
# writeXStringSet(sloth.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/sloth.sc.ogids.apr23.fa')
# writeXStringSet(spectacled.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/spectacled.sc.ogids.apr23.fa')
# writeXStringSet(asiatic.iso.filtered.og,filepath = '/bearComparativeGenomics_Spring2025/hyphy/orthofinder_oldBrownAllNew_apr22/single_copy_fastas/asiatic.sc.ogids.apr23.fa')
# 










##
