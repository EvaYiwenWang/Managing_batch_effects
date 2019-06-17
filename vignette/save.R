
# sponge

rawdata.sponge = read.csv("./datasets/Intrasp_variation_paper.csv",header = T,row.names = 1)
counts.sponge = rawdata.sponge[10:33]
dim(counts.sponge)

# batch effect
gel.batch = as.factor(rawdata.sponge$gel)
names(gel.batch) = rownames(rawdata.sponge)

# effect of interest 
Tissue.trt = as.factor(rawdata.sponge$ectcho)
names(Tissue.trt) = rownames(rawdata.sponge)
table(gel.batch,Tissue.trt)
sponge.tss <- counts.sponge/100
sponge.batch <- gel.batch
sponge.trt <- Tissue.trt

sponge.batch = sort(sponge.batch)
sponge.trt = sponge.trt[names(sponge.batch)]
sponge.count = sponge.count[names(sponge.batch),]


#### AD
rawdata.ad <- read.csv("./datasets/abundance.csv",header = T)
rownames(rawdata.ad) = rawdata.ad[,1]
counts.ad = rawdata.ad[,-c(1,2)]
counts.ad = t(counts.ad)
dim(counts.ad)

metadata.ad <- read.csv("./datasets/metadata_KA.csv", header = T)

# batch effect
time.batch = as.factor(metadata.ad$sequencing_run_date)
names(time.batch) = metadata.ad$analysis_name
time.batch = factor(time.batch, levels = c('09/04/2015', '14/04/2016', '01/07/2016', '14/11/2016', '21/09/2017'))

# effect of interest
metadata.ad$initial_phenol_concentration.regroup = NA
metadata.ad[c(which(metadata.ad$initial_phenol_concentration == '0'),which(metadata.ad$initial_phenol_concentration == '0.05'), which(metadata.ad$initial_phenol_concentration == '0.1'),which(metadata.ad$initial_phenol_concentration == '0.5')),]$initial_phenol_concentration.regroup = '0-0.5'
metadata.ad[c(which(metadata.ad$initial_phenol_concentration == '1'),which(metadata.ad$initial_phenol_concentration == '1.5'), which(metadata.ad$initial_phenol_concentration == '2')),]$initial_phenol_concentration.regroup = '1-2'

concentration.trt = as.factor(metadata.ad$initial_phenol_concentration.regroup)
names(concentration.trt) = metadata.ad$analysis_name

table(time.batch,concentration.trt)
ad.count <- counts.ad
ad.count <- ad.count[as.character(metadata.ad$analysis_name),]
ad.batch <- time.batch
ad.trt <- concentration.trt

#
ad.batch = sort(ad.batch)
ad.trt = ad.trt[names(ad.batch)]
ad.count = ad.count[names(ad.batch),]
ad.metadata <- metadata.ad
rownames(ad.metadata) <- as.character(ad.metadata$analysis_name)
ad.metadata <- ad.metadata[rownames(ad.count),]


###########
# hd data
##------------------
## Process 16S data 
##------------------
data.count = read.table('./datasets/Combined_data.tsv', row.names = 1, header = TRUE, sep = '\t')
#dim(data.count)
# trranspose data
data.count = t(data.count)

# Meta data
metadata = read.table('./datasets/Sample_metadata.tsv', row.names = 1, header = TRUE, sep = '\t')

# checking
table(metadata$Genotype, metadata$Cage)

# remove cage A,C,I
data.count.old = data.count
metadata.old = metadata

data.count = data.count.old[-c(which(metadata.old$Cage == 'A'),
                               which(metadata.old$Cage == 'C'),which(metadata.old$Cage == 'I')),]
metadata = metadata.old[-c(which(metadata.old$Cage == 'A'),
                           which(metadata.old$Cage == 'C'),which(metadata.old$Cage == 'I')),]

data.count.fem = data.count[metadata$Sex == 'F',]
metadata.fem = metadata[metadata$Sex == 'F',]

metadata.fem$Genotype = droplevels(metadata.fem$Genotype)
metadata.fem$Cage = droplevels(metadata.fem$Cage)

table(metadata.fem$Genotype, metadata.fem$Cage)

hd.count = data.count.fem
rownames(hd.count) = paste0('M',rownames(metadata.fem))
colnames(hd.count) = paste0('OTU',1:ncol(hd.count))
hd.batch = metadata.fem$Cage
hd.trt = metadata.fem$Genotype
names(hd.batch) = names(hd.trt) = paste0('M',rownames(metadata.fem))

hd.batch = sort(hd.batch)
hd.trt = hd.trt[names(hd.batch)]
hd.count = hd.count[names(hd.batch),]


save(sponge.tss,sponge.batch,sponge.trt,ad.count,ad.batch,ad.trt,ad.metadata, hd.count, hd.batch, hd.trt,
     file = './datasets/microbiome_datasets.RData')










