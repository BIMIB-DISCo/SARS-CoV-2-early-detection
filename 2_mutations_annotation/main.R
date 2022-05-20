# load data
load(file="final_data/annotations.RData")
load(file="final_data/human_cu.RData")
load(file="final_data/metadata.RData")
load(file="final_data/sars_cov2_fasta.RData")
load(file="final_data/variants.RData")

# process annotations
cds_start = as.numeric(annotations[which(annotations[,3]=="CDS"),4])
cds_end = as.numeric(annotations[which(annotations[,3]=="CDS"),5])
cds_name = c("ORF1ab","ORF1ab","ORF1ab","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")

# consider only variants spanning CDS regions
cds_positions = NULL
for(i in 1:length(cds_start)) {
    cds_positions = c(cds_positions,cds_start[i]:cds_end[i])
}
variants = variants[which(as.numeric(variants[,"POSITION"])%in%cds_positions),]

# remove variants in uncertain positions (different between REF-ANC and Wuhan1)
uncertain_positions = NULL
uncertain_variant = 8782
uncertain_variant_offset = uncertain_variant-cds_start[which(cds_start<=uncertain_variant&cds_end>=uncertain_variant)[1]]+1
if(uncertain_variant_offset%%3==1) { # first position of triplet
    uncertain_positions = c(uncertain_positions,c(uncertain_variant,(uncertain_variant+1),(uncertain_variant+2)))
}
if(uncertain_variant_offset%%3==2) { # second position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-1),uncertain_variant,c(uncertain_variant+1)))
}
if(uncertain_variant_offset%%3==0) { # third position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-2),(uncertain_variant-1),uncertain_variant))
}
uncertain_variant = 28144
uncertain_variant_offset = uncertain_variant-cds_start[which(cds_start<=uncertain_variant&cds_end>=uncertain_variant)[1]]+1
if(uncertain_variant_offset%%3==1) { # first position of triplet
    uncertain_positions = c(uncertain_positions,c(uncertain_variant,(uncertain_variant+1),(uncertain_variant+2)))
}
if(uncertain_variant_offset%%3==2) { # second position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-1),uncertain_variant,c(uncertain_variant+1)))
}
if(uncertain_variant_offset%%3==0) { # third position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-2),(uncertain_variant-1),uncertain_variant))
}
uncertain_variant = 29095
uncertain_variant_offset = uncertain_variant-cds_start[which(cds_start<=uncertain_variant&cds_end>=uncertain_variant)[1]]+1
if(uncertain_variant_offset%%3==1) { # first position of triplet
    uncertain_positions = c(uncertain_positions,c(uncertain_variant,(uncertain_variant+1),(uncertain_variant+2)))
}
if(uncertain_variant_offset%%3==2) { # second position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-1),uncertain_variant,c(uncertain_variant+1)))
}
if(uncertain_variant_offset%%3==0) { # third position of triplet
    uncertain_positions = c(uncertain_positions,c((uncertain_variant-2),(uncertain_variant-1),uncertain_variant))
}
variants = variants[which(!as.numeric(variants[,"POSITION"])%in%uncertain_positions),]

# perform analysis
dataset = array(NA,c(nrow(variants),15))
data_pos = 0
for(i in 1:nrow(variants)) {
    curr_variants = variants[i,,drop=TRUE]
    curr_position = which(cds_start<=as.numeric(curr_variants["POSITION"])&cds_end>=as.numeric(curr_variants["POSITION"]))
    curr_position = curr_position[1] # take first if I have multiple matches to CDS
    curr_cds = cds_name[curr_position]
    curr_offset = as.numeric(curr_variants["POSITION"])-cds_start[curr_position]+1
    curr_cds_codon = ceiling(curr_offset/3)
    if(curr_offset%%3==1) { # first position of triplet
        triplet_wt = sars_cov2_fasta[c(as.numeric(curr_variants["POSITION"]),(as.numeric(curr_variants["POSITION"])+1),(as.numeric(curr_variants["POSITION"])+2))]
        triplet_mut = c(as.character(curr_variants["VARIANT_ALLELE"]),triplet_wt[2],triplet_wt[3])
    }
    if(curr_offset%%3==2) { # second position of triplet
        triplet_wt = sars_cov2_fasta[c((as.numeric(curr_variants["POSITION"])-1),as.numeric(curr_variants["POSITION"]),(as.numeric(curr_variants["POSITION"])+1))]
        triplet_mut = c(triplet_wt[1],as.character(curr_variants["VARIANT_ALLELE"]),triplet_wt[3])
    }
    if(curr_offset%%3==0) { # third position of triplet
        triplet_wt = sars_cov2_fasta[c((as.numeric(curr_variants["POSITION"])-2),(as.numeric(curr_variants["POSITION"])-1),as.numeric(curr_variants["POSITION"]))]
        triplet_mut = c(triplet_wt[1],triplet_wt[2],as.character(curr_variants["VARIANT_ALLELE"]))
    }
    triplet_wt = paste0(triplet_wt,collapse="")
    triplet_wt = as.character(human_cu[which(human_cu[,1]==triplet_wt),])
    triplet_mut = paste0(triplet_mut,collapse="")
    triplet_mut = as.character(human_cu[which(human_cu[,1]==triplet_mut),])
    # if the variant is nonsynonymous then it is valid
    if(triplet_wt[2]!=triplet_mut[2]) {
        curr_res = as.character(c(curr_variants[c(3,4,5,6)],triplet_wt,triplet_mut,curr_variants[1]))
        curr_res = c(curr_res,metadata$CollectionDate[which(metadata$Run==curr_res[11])],metadata$Country[which(metadata$Run==curr_res[11])])
        data_pos = data_pos + 1
        dataset[data_pos,1:15] = c(curr_res,curr_cds,curr_cds_codon)
    }
    cat(i/nrow(variants),"\n")
}
dataset = dataset[1:data_pos,]
rownames(dataset) = 1:nrow(dataset)
colnames(dataset) = c("REF_POS","REF","ALT","VF","CWT","PWT","VWT","CMUT","PMUT","VMUT","PATIENT","DATE","COUNTRY","CDS_NAME","CDS_POS")
dataset = dataset[,c("CDS_NAME","CDS_POS","PATIENT","REF_POS","REF","ALT","VF","CWT","PWT","VWT","CMUT","PMUT","VMUT","DATE","COUNTRY")]

# post-processing
VARIANT_ID = as.character(paste0(dataset[,"CDS_NAME"],"_",dataset[,"CDS_POS"],"_",dataset[,"PMUT"]))
CDS_NAME = as.character(dataset[,"CDS_NAME"])
CDS_POS = as.numeric(dataset[,"CDS_POS"])
PATIENT_ID = as.character(dataset[,"PATIENT"])
REF_POS = as.numeric(dataset[,"REF_POS"])
REF_BASE = as.character(dataset[,"REF"])
ALT_BASE = as.character(dataset[,"ALT"])
VARIANT_FREQ = as.numeric(dataset[,"VF"])
CODON_REF = as.character(dataset[,"CWT"])
PROTEIN_REF = as.character(dataset[,"PWT"])
CODON_ALT = as.character(dataset[,"CMUT"])
PROTEIN_ALT = as.character(dataset[,"PMUT"])
COLLECTION_DATE = as.Date(as.character(dataset[,"DATE"]),format="%Y-%m-%d")
COUNTRY = as.character(dataset[,"COUNTRY"])
dataset = data.frame(VARIANT_ID=VARIANT_ID,CDS_NAME=CDS_NAME,CDS_POS=CDS_POS,PATIENT_ID=PATIENT_ID,REF_POS=REF_POS,REF_BASE=REF_BASE,ALT_BASE=ALT_BASE,VARIANT_FREQ=VARIANT_FREQ,CODON_REF=CODON_REF,PROTEIN_REF=PROTEIN_REF,CODON_ALT=CODON_ALT,PROTEIN_ALT=PROTEIN_ALT,COLLECTION_DATE=COLLECTION_DATE,COUNTRY=COUNTRY)

# save results
save(dataset,file="dataset.RData")
