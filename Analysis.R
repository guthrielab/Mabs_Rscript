-------------------------------------------------------------------------------
#### Mycobacterium abscessus genome analysis ####

-------------------------------------------------------------------------------

# This Rscript is a discription to reproduce some of the analysis for the M. abscessus
# project from the Montreal study
  
#### Mabs pipeline ####

# Clone pipeline from the hybrid branch of the repo: https://github.com/waglecn/mabs/tree/hybrid
# Follow the steps in the ReadMe to install
# Create a sample sheet file using the template of samples.test.yaml to point the config file to where samples are located
# FASTQs are stored in Superdome: /workspace/lab/guthrielab/1_projects/mabs/abscessus/raw_reads
# and also available in SRA under BioProject xxxx
# Make changes to the config.yaml file to reflect the output directory, GATK, and Kraken path
# Run the pipeline using snakemake following the ReadMe on the steps sequentially
#   Stage 1: QC
#       At the end of this step, you will have a QC_summary.csv file which allows you to exclude samples (using the config.yaml) that fail QC from other stages
#   Stage 2: Mapping and Variant Calling
#   Stage 3: Recombination and Phylogenetics

#### Data Dictionary for M. abscessus isolates from Montreal ####

# Patient ID.y = Unique patient identifier
# LSPQ = Sample identifier from LSPQ (Laboratoire de santé publique du Québec)
# Sex = Biological sex of the patient
# DOB = Date of birth of the patient
# CF status = Cystic Fibrosis status of the patient (0 = NO or not mentioned 1 = YES)
# Lab = Micorbiology Laboratory Hospital where samples were sent to (RVH 1, JGH 2, CHUM 3, HMR 4, StJu 5, other 6)
# DST = Drug susceptibility testing (0 = NO, 1 = YES)
# Specimen = Sample source ( 1 = respiratory, 2 = ear, 3 = other )
# Anatomic = Site of sample source
# Lung transplant (0 = NO, 1 = YES)

# Metadata file = new_metadata_merged.xlsx

library(dplyr)
library(readxl)
library(writexl)

epidata <- read_xlsx("/workspace/lab/guthrielab/1_projects/mabs/abscessus/1_epidata/Montreal_lab_isolates.xlsx")
epidata2 <- read_xlsx("/workspace/lab/guthrielab/1_projects/mabs/abscessus/1_epidata/Quebec_lab_isolates.xlsx")

# Merge both epidata by LSPQ
epi.patient <- epidata %>% 
  right_join(epidata2, by = c("LSPQ" = "samples received from LSPQ")) %>% 
  select(`Patient ID.y`, LSPQ, Sex, DOB, `CF status 0 = NO or not mentioned\r\n1 = YES`, `Lab (RVH 1, JGH 2, CHUM 3, HMR 4, StJu 5, other 6)`, `DST done 0 = NO\r\n1 = YES`, 
         `Specimen: 1 resp; 2 ear; 3 other`,'Anatomic', `Lung transplant 0 = NO \r\n1 = YES`)

write_xlsx(epi.patient, "new_metadata_merged.xlsx")



#### Figure 1 ####

library(lubridate)
library(ggplot2)
library(ggdist)
library(ggridges)
library(cowplot)

# Merge the QC_summary file with the 

mabs <- read.csv('mabs-count.csv', header = T)
mmas <- read.csv('mmas-count.csv', header = T)
mbol <- read.csv('mbol-count.csv', header = T)
mabs$sample_date <- as.Date(mabs$sample_date)
mmas$sample_date <- as.Date(mmas$sample_date)
mbol$sample_date <- as.Date(mbol$sample_date)
mabs$MRCA_ref <- as.factor(mabs$MRCA_ref)
mmas$MRCA_ref <- as.factor(mmas$MRCA_ref)
mbol$MRCA_ref <- as.factor(mbol$MRCA_ref)
mabs$date2 <- as.Date(cut(mabs$sample_date, breaks = '1 month', start.on.monday = F))
mmas$date2 <- as.Date(cut(mmas$sample_date, breaks = '1 month', start.on.monday = F))
mbol$date2 <- as.Date(cut(mbol$sample_date, breaks = '1 month', start.on.monday = F))
ggplot(mabs, aes(x = date2, y = MRCA_ref, height = stat(count))) +
  geom_density_ridges2(stat = 'binline', bins = 20, scale = 0.95, draw_baseline = T, aes(fill = MRCA_ref)) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "9 months") + scale_y_discrete() +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_classic() + theme(legend.position = 'none') #+
xlab("Sample collection date") + ylab("M. abscessus subsp.")


A <- ggplot(mabs, aes(x = date2, fill = MRCA_ref)) +
  geom_bar(stat = 'count', position = 'dodge', width = 60) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "9 months", limits = c(as.Date('2010-03-01'),as.Date('2018-01-01'))) + theme_classic() + theme(legend.position = 'none')+
  scale_fill_manual(values = '#00AFBB') + scale_y_continuous(limits = c(0,8)) +xlab("") +ylab("Frequency")

C <- ggplot(mmas, aes(x = date2, fill = MRCA_ref)) +
  geom_bar(stat = 'count', position = 'dodge', width = 60) +
  scale_x_date(date_labels = '%b\n%Y', date_breaks = '9 months', limits = c(as.Date('2010-03-01'),as.Date('2018-01-01'))) + theme_classic() + theme(legend.position = 'none') +
  scale_fill_manual(values = '#FC4E07')+ scale_y_continuous(limits = c(0,8)) +xlab("Sample collection date") +ylab("Frequency")

B <- ggplot(mbol, aes(x=date2, fill = MRCA_ref)) +
  geom_bar(stat = 'count', position = 'dodge', width = 60) +
  scale_x_date(date_labels = '%b\n%Y', date_breaks = '9 months', limits = c(as.Date('2010-03-01'),as.Date('2018-01-01'))) + theme_classic() + theme(legend.position = 'none') +
  scale_fill_manual(values= '#E7B800')+ scale_y_continuous(limits = c(0,8)) +xlab("") +ylab("Frequency")

plot_grid(A,B,C,nrow=3)

# DIversity plot
library(dplyr)
library(ggplot2)
library(cowplot)
library(matrixStats)

snpdist <- read.csv('mmassiliense.merge.gubbins.fasta.snpdists-new.csv', header=T, row.names = 1,)
# mmassilense patients
M07 <-  c('L00042522',	'MB081634',	'MB083272',	'MB085071',	'MB087124',	'MB088512',	'MB090324',	'MB091599',	'MB093159')
M19 <-  c('L00038032',	'L00040754')
M45 <- c('MB085750',	'MB088215',	'MB088425',	'MB089904',	'MB092343',	'MB092961')
M69 <- c('MB081436',	'MB084704',	'MB087316')
M71 <- c('L00003612',	'L00007189',	'L00016748')
M73 <- c('MB088685_1',	'MB088685_2')
M74 <- c('MB087103', 'MB087463')
M79 <- c('L00012354', 'L00042516', 'MB092254')
mmas.btw.pt <- c('MB084633','L00042522','L00038032','L00042538','MB087714','MB084195','MB092141',
                 'MB085750','L00037541','L00030602','L00036967','MB086252','MB081436','MB086213',
                 'L00003612','MB088685_1','MB087103','MB085218','MB085385','L00004793','L00012354',
                 'MB087203','L00003201','L00004856','L00036398','MB0153354','MB032793','MB045515','MB046809','MB092743','MB092800_1','MB092800_2')
snpdist.mbol <- read.csv('mbolletii.merge.gubbins.fasta.snpdists.csv', header=T, row.names = 1,)
# bolettii patients
mbol.within <- c('MB080688_ibis', 'MB080690_ibis','MB082744_ibis')

snpdist_m07 <- snpdist[rownames(snpdist) %in% M07, colnames(snpdist) %in% M07, ]
snpdist_m19 <- snpdist[rownames(snpdist) %in% M19, colnames(snpdist) %in% M19, ]
snpdist_m45 <- snpdist[rownames(snpdist) %in% M45, colnames(snpdist) %in% M45, ]
snpdist_m69 <- snpdist[rownames(snpdist) %in% M69, colnames(snpdist) %in% M69, ]
snpdist_m71 <- snpdist[rownames(snpdist) %in% M71, colnames(snpdist) %in% M71, ]
snpdist_m73 <- snpdist[rownames(snpdist) %in% M73, colnames(snpdist) %in% M73, ]
snpdist_m74 <- snpdist[rownames(snpdist) %in% M74, colnames(snpdist) %in% M74, ]
snpdist_m79 <- snpdist[rownames(snpdist) %in% M79, colnames(snpdist) %in% M79, ]
snpdist_btwn_pt <- snpdist[rownames(snpdist) %in% mmas.btw.pt, colnames(snpdist) %in% mmas.btw.pt, ]

mbol.snps.wihin <- snpdist.mbol[rownames(snpdist.mbol) %in% mbol.within, colnames(snpdist.mbol) %in% mbol.within, ]
mbol.patients.id <- c("L00006319_ibis","L00042507_ibis","MB082744_ibis")
mbol.snps.btw <- snpdist.mbol[rownames(snpdist.mbol)%in% mbol.patients.id, colnames(snpdist.mbol) %in% mbol.patients.id, ]

snpdist.mabs <- read.csv('mabscessus.merge.gubbins.fasta.snpdists-edit.csv', header=T, row.names = 1,)
# abscessus patients
M01 <- c('MB083519','MB088683','MB090456')
M04 <- c('MB084351','MB088330')
M09 <- c('MB087998','MB089698','MB091852')
M10 <- c('L00019464','MB080517','MB085874','MB093021')
M12 <- c('L00004052','MB081209','MB081274','MB085532')
M13 <- c('MB088931','MB090458','MB091794','MB092927_1')
M14 <- c('MB075002','MB080668','MB082680','MB087514')
M16 <- c('L00015420','MB089963')
M18 <- c('L00005920','L00011957','MB092409','MB093074')
M20 <- c('L00039460','MB082896','MB086151','MB089024','MB091099','MB093261')
M24 <- c('L00042518','MB083835','MB086490')
M30 <- c('MB078705','MB078716','MB080848')
M32 <- c('L00004048','MB084821')
M40 <- c('L00004074','L00012343','L00029190','L00040772','MB081208','MB085749','MB087072',
         'MB088060','MB088834','MB089866','MB091987','MB091987')
M41 <- c('MB086831','MB089146')
M43 <- c('L00017444','MB088826','MB090190','MB091939')
M50 <- c('MB087802','MB089257')
M51 <- c('L00007906','MB084982','MB088020','MB089718','MB091633')
M53 <- c('L00029597_1','L00029597_2')
M59 <- c('MB081484','MB083076','MB084806')
M62 <- c('MB081563','MB085338_1','MB085338_2','MB088150')
M63 <- c('MB081709','MB085526','MB086146')
M66 <- c('L00019475','MB085585','MB091027','MB093198')
M67 <- c('L00027380','L00028203','MB082123','MB083394','MB085135','MB087023',
         'MB088934','MB090595','MB092130')
M83 <- c('MB087710','MB089546')
mabs.btw.patients <- c('MB083519','MB084351','MB082781','L00014051','L00006438','MB087998',
                       'L00019464','MB087393','L00004052','MB088931','MB075002','MB085564','L00015420','L00011957',
                       'L00039460','L00009173','L00042518','L00005918','MB073717','L00038889','L00044208',
                       'MB078705','L00038890','L00004048','L00013157','MB080831','MB091990','MB081290',
                       'MB092665','MB088577','L00004074','MB086831','MB091679','L00017444','MB084946',
                       'L00018091','MB087802','L00007906','MB080919','L00029597_1','L00031343','L00036434',
                       'MB081484','MB086605','MB081563','MB081709','MB088646','MB091631','L00019475',
                       'L00027380','MB084356','MB092161','L00016581','MB089903','MB087512','MB087710',
                       'MB092520','MB084356','MB089963','MB087293', 'L00003204', 'L00018966', 'L00019656', 'L00027056',
                       'L00034794','MB019977','MB091565','MB091566','MB091830','MB091920','MB092084','MB092759','MB093049','MB093228','MB092927_2')

snpdist_m01 <- snpdist.mabs[rownames(snpdist.mabs) %in% M01, colnames(snpdist.mabs) %in% M01, ]
snpdist_m04 <- snpdist.mabs[rownames(snpdist.mabs) %in% M04, colnames(snpdist.mabs) %in% M04, ]
snpdist_m09 <- snpdist.mabs[rownames(snpdist.mabs) %in% M09, colnames(snpdist.mabs) %in% M09, ]
snpdist_m10 <- snpdist.mabs[rownames(snpdist.mabs) %in% M10, colnames(snpdist.mabs) %in% M10, ]
snpdist_m12 <- snpdist.mabs[rownames(snpdist.mabs) %in% M12, colnames(snpdist.mabs) %in% M12, ]
snpdist_m13 <- snpdist.mabs[rownames(snpdist.mabs) %in% M13, colnames(snpdist.mabs) %in% M13, ]
snpdist_m14 <- snpdist.mabs[rownames(snpdist.mabs) %in% M14, colnames(snpdist.mabs) %in% M14, ]
snpdist_m16 <- snpdist.mabs[rownames(snpdist.mabs) %in% M16, colnames(snpdist.mabs) %in% M16, ]
snpdist_m18 <- snpdist.mabs[rownames(snpdist.mabs) %in% M18, colnames(snpdist.mabs) %in% M18, ]
snpdist_m20 <- snpdist.mabs[rownames(snpdist.mabs) %in% M20, colnames(snpdist.mabs) %in% M20, ]
snpdist_m24 <- snpdist.mabs[rownames(snpdist.mabs) %in% M24, colnames(snpdist.mabs) %in% M24, ]
snpdist_m30 <- snpdist.mabs[rownames(snpdist.mabs) %in% M30, colnames(snpdist.mabs) %in% M30, ]
snpdist_m32 <- snpdist.mabs[rownames(snpdist.mabs) %in% M32, colnames(snpdist.mabs) %in% M32, ]
snpdist_m40 <- snpdist.mabs[rownames(snpdist.mabs) %in% M40, colnames(snpdist.mabs) %in% M40, ]
snpdist_m41 <- snpdist.mabs[rownames(snpdist.mabs) %in% M41, colnames(snpdist.mabs) %in% M41, ]
snpdist_m43 <- snpdist.mabs[rownames(snpdist.mabs) %in% M43, colnames(snpdist.mabs) %in% M43, ]
snpdist_m50 <- snpdist.mabs[rownames(snpdist.mabs) %in% M50, colnames(snpdist.mabs) %in% M50, ]
snpdist_m51 <- snpdist.mabs[rownames(snpdist.mabs) %in% M51, colnames(snpdist.mabs) %in% M51, ]
snpdist_m53 <- snpdist.mabs[rownames(snpdist.mabs) %in% M53, colnames(snpdist.mabs) %in% M53, ]
snpdist_m59 <- snpdist.mabs[rownames(snpdist.mabs) %in% M59, colnames(snpdist.mabs) %in% M59, ]
snpdist_m62 <- snpdist.mabs[rownames(snpdist.mabs) %in% M62, colnames(snpdist.mabs) %in% M62, ]
snpdist_m63 <- snpdist.mabs[rownames(snpdist.mabs) %in% M63, colnames(snpdist.mabs) %in% M63, ]
snpdist_m66 <- snpdist.mabs[rownames(snpdist.mabs) %in% M66, colnames(snpdist.mabs) %in% M66, ]
snpdist_m67 <- snpdist.mabs[rownames(snpdist.mabs) %in% M67, colnames(snpdist.mabs) %in% M67, ]
snpdist_m83 <- snpdist.mabs[rownames(snpdist.mabs) %in% M83, colnames(snpdist.mabs) %in% M83, ]
snpdist.mabs.btw <- snpdist.mabs[rownames(snpdist.mabs) %in% mabs.btw.patients, colnames(snpdist.mabs) %in% mabs.btw.patients, ]

within <- read.csv('within-snp.csv', header = T)
A <- ggplot(within, aes(x=Mabs_sp, y=SNPs, fill=Mabs_sp)) +
  stat_boxplot(geom = 'errorbar', width=0.5) + geom_boxplot(width=0.5, color='black') + scale_y_log10(limits=c(1,10000))+
  theme_classic() + theme(legend.position = "none") + xlab("Within patients") + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

between <- read.csv('between-snp.csv', header = T)
B <- ggplot(between, aes(x=Mabs_sp, y=SNPs, fill=Mabs_sp)) +
  stat_boxplot(geom = 'errorbar',width=0.5) + geom_boxplot(width=0.5, color='black') +scale_y_log10(limits=c(1,10000))+
  theme_classic() + theme(legend.position = 'none') + ylab('') + xlab("Between patients") + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))


plot_grid(A,B,nrow = 1,ncol=2)


#### Figure 2 ####
library(ggtree)
library(ape)
library(ggplot2)
library(cowplot)
library(tidytree)
library(treeio)
library(ggtreeExtra)
library(lubridate)
library(ggnewscale)
library(phangorn)
library(dplyr)
library(TreeDist)

# Treefile used in figure 2 and 3 are retrieved from recombination adjusted SNP sites from gubbins. 
# Pan-genome tree was retrieved from accessory genes using roary with default parameters
# Metadata file contains 8 columns: /sample/patient_ID/ST/CC/CF_status/date/hospital/important_ID/
# M. massiliense phylogenetic tree annotation with CF status, ST, CC, and DCC

mmas.tree <- read.newick('mmassiliense.merge.gubbins.fasta.treefile')
tree1.3 <- reorder(mmas.tree)
mmas.cc <- read.csv('mmas.mlst.csv', header=T)
mmas.cc$CC <- factor(mmas.cc$CC)
mmas.cc$ST <- factor(mmas.cc$ST)
mmas.pangenome <- read.newick('mmas-new_accessory_binary_genes.fa.newick')
mmas.pangenome <- reorder(mmas.pangenome)

tree.tips.mmas <- ggtree(midpoint(tree1.3)) %<+% mmas.cc +
  geom_tippoint(aes(colour=CF_status),size=2) + geom_tiplab(aes(label=factor(important_ID)), hjust=-.2) +
  scale_colour_manual(name='CF status', values = c("red","grey")) + theme(legend.position = 'none')

tree1.1 <- tree.tips.mmas + geom_fruit(geom=geom_tile, mapping=aes(fill=ST), width=0.01, offset = 0.15) + scale_fill_manual(values=c('#FC4E07','#a8e4a0','#7442c8','#ff48d0','#339624','#ff9baa','#b5674d','#1974d2'),na.value = '#D9D9D9') + 
  theme(legend.background = element_rect(colour = 'black',fill='white',linetype='solid')) + theme(legend.position = "top")
newmastree <- tree1.1 + geom_cladelabel(node=64, label = 'Clonal Complex 7 \n DCC7', align=T, angle=90, barsize=1, hjust='centre', offset = 0.01, offset.text = 0.003) +
  geom_cladelabel(node=101, label = 'Clonal Complex 3 \n DCC3', align=T, angle = 90, barsize = 1, hjust = 'centre', offset=0.01, offset.text = 0.003) +
  geom_cladelabel(node = 62, label = 'Clonal Complex 6 \n DCC6', align = F, barsize = 1, hjust = 'centre', offset.text = 0.02, offset=0.01) +
  geom_treescale(width = 0.05)

ggsave(filename = "mmas-revised-tree1.svg", plot = newmastree, width = 8, height = 20)
#tree1.1.1 <- tree1.1+ new_scale_fill() + geom_fruit(geom=geom_tile, mapping = aes(fill=CC), width=0.01, offset =0.06) +
#scale_fill_manual(name='CC', values = c('#AA0A3C','gold','green4')) + geom_treescale(width = 0.05)

# M. abscessus phylogenetic tree annotation with CF status, ST, CC, and DCC
mabs.tree <- read.newick('mabscessus.merge.gubbins.fasta.treefile')
tree2.3 <- reorder(mabs.tree)
mabs.cc <- read.csv('mabs.mlst.csv', header = T)
mabs.cc$CC <- factor(mabs.cc$CC)
mabs.cc$ST <- factor(mabs.cc$ST)
mabs.pangenome <- read.newick('mabs-new-accessory_binary_genes.fa.2.newick')
mabs.pangenome <- reorder(mabs.pangenome)

tree.tips.mabs <- ggtree(midpoint(tree2.3)) %<+% mabs.cc +
  geom_tippoint(aes(colour=CF_status),size=2) + geom_tiplab(aes(label=factor(important_ID)), hjust=-.2) +
  scale_colour_manual(name='CF status', values = c("red","grey")) + theme(legend.position = 'none')

tree2.2 <- tree.tips.mabs + geom_fruit(geom=geom_tile, mapping=aes(fill=ST), width=0.004, offset = 0.095) + scale_fill_manual(na.value = '#D9D9D9',values = c('#5318D9','#A571D9','#D90BD1','#E02676','#9E00A8',
                                                                                                                                                              '#00A82E','#75A800','#18A89F','#65EB0E','#5EEBB8',
                                                                                                                                                              '#EB3405','#EB9F1A','#EB7878','#8C0101','#E0DC55',
                                                                                                                                                              '#77ACF2','#2111F2','#1BA2F2','#00FFCD')) +
  theme(legend.background = element_rect(colour = 'black',fill='white',linetype='solid')) + theme(legend.position = "top")

newmabstree <- tree2.2 + geom_cladelab(node = 170, label='Clonal Complex 1 \n DCC1', align = F, angle=90, barsize=1, hjust= 'centre', offset=0.02,offset.text=0.002) +
  geom_cladelabel(node=193, label = 'Clonal Complex 2 \n DCC2', align=F, barsize=1, hjust='centre', offset=0.01,offset.text = 0.012) +
  geom_cladelabel(node=209, label='Clonal Complex 5 \n DCC5', align = F, angle=90, barsize=1, hjust='centre', offset = 0.003, offset.text = 0.001) +
  geom_cladelabel(node=242, label='Clonal Complex 4 \n DCC4', align=F, barsize = 1, hjust = 'centre', offset=0.01, offset.text = 0.012) +
  geom_treescale(width=0.05)

ggsave(filename = "mabs-revised-tree1.svg", plot = newmabstree, width = 8, height = 20)


saveme <- cowplot::plot_grid(newmabstree, newmastree, fakeplot,nrow = 1)
ggsave(filename = "figure2v2.svg", plot= saveme, width = 24, height=20)
#tree2.2 + new_scale_fill()+ geom_fruit(geom=geom_tile, mapping=aes(fill=CC),width=0.005, offset=0.05) +
#scale_fill_manual(name='CC', values = c('coral','brown','navy','#f5c9d9')) + geom_treescale(width = 0.05)


#### Figure 3 ####
# M. massiliense core and pan-genome tree comparison
tree1 <- ggtree(midpoint(tree1.3)) %<+% mmas.cc +
  geom_tiplab(aes(label=factor(patient_ID), colour=CF_status)) +
  scale_colour_manual(name='CF status', values = c('black','#AA0A3C')) + theme(legend.position = 'none')

mmas.pangenome1 <- ggtree(midpoint(mmas.pangenome)) %<+% mmas.cc +
  geom_tiplab(aes(label=factor(patient_ID), colour=CF_status), hjust=1) +
  scale_colour_manual(name='CF status', values = c('black','#AA0A3C')) + scale_x_reverse() + theme(legend.position = 'none')

d1 = tree1$data[tree1$data$isTip,]
d1$x[] = 1
d2 = mmas.pangenome1$data[mmas.pangenome1$data$isTip,]
d2$x[] = 2
TTcon <- rbind(d1,d2)

L1 = ggplot(TTcon, aes(x = x, y = y, colour = patient_ID, group = label)) + geom_line() +   
  theme_void() + theme(legend.position="none", plot.margin = unit(c(0,0,0,0),"cm"))
cowplot::plot_grid(tree1, L1, mmas.pangenome1, nrow = 1, align = 'hv')

# M. abscessus core and pan-genome tree comparision
tree2 <- ggtree(midpoint(tree2.3)) %<+% mabs.cc +
  geom_tiplab(aes(label=factor(patient_ID),colour=CF_status),size=3) +
  scale_colour_manual(name='CF status',values = c('black','#AA0A3C')) + theme(legend.position = 'none')

mabs.pangenome1 <- ggtree(midpoint(mabs.pangenome)) %<+% mabs.cc +
  geom_tiplab(aes(label=factor(patient_ID), colour=CF_status), hjust=1, size=3) +
  scale_colour_manual(name='CF status', values = c('black','#AA0A3C')) + scale_x_reverse() + theme(legend.position = 'none')

g1 = tree2$data[tree2$data$isTip,]
g1$x[] = 1
g2 = mabs.pangenome1$data[mabs.pangenome1$data$isTip,]
g2$x[] = 2
TTcon2 <- rbind(g1,g2)

L2 = ggplot(TTcon2, aes(x = x, y = y, colour = patient_ID, group = label)) + geom_line() +   
  theme_void() + theme(legend.position="none", plot.margin = unit(c(0,0,0,0),"cm"))
cowplot::plot_grid(tree2, L2, mabs.pangenome1, nrow = 1, align = 'hv')

#### Supplementary Data ####
# Supplementary figure 2 - Patients with multiple isolates
# Metadata file has 4 columns: /patient_ID/sample_date/Subspecies/sample_ID/
# Read Data

patients2 <- read.csv("~/Documents/abscesus/montreal/patient-timeline.csv", header=T)
patients2$sample_date <- as.Date(patients2$sample_date)
patients2$Subspecies <- as.factor(patients2$Subspecies)

# Plot timeline

ggplot(patients2, aes(x = sample_date, y = as.factor(patient_ID))) +
  geom_point(size = 3, aes(colour = Subspecies)) +
  geom_line(aes(group = patient_ID), linewidth = 0.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "6 month")+
  labs(title = "",
       x = "Sample collection date",
       y = "Patient ID") + 
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size = 14, face="bold")) +
  theme(axis.text.x = element_text(color="black", size = 12)) +
  theme(axis.title.y = element_text(color="black", size = 14, face="bold")) +
  theme(axis.text.y = element_text(color="black", size = 12)) +
  theme(legend.text = element_text(color='black', size = 12)) +
  theme(legend.position = "top")

