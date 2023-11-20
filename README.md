### ID and SNP extraction

first extract a ID list, and screen out SNPs with imputation INFO score > 0.3
```shell
folders=("0410" "0411" "0412" "0413")

module load Plink/2.0

for fl in "${folders[@]}"; do

tmpdir=/PHShome/ys724/scratch/IDsubset/${fl}
mkdir ${tmpdir}

for chr in $(seq 1 22);do

cd /data/biobank_release/datasets/${fl}/vcf

plink2 \
  --vcf chr${chr}.dose.vcf.gz \
  dosage=DS \
  --make-pfile \
  --keep /PHShome/ys724/Documents/GMBI_endometriosis/all/Endo_ID \
  --out ${tmpdir}/chr${chr}_extractedID_extract

# imputation INFO score > 0.3
zcat chr${chr}.info.gz | awk '$7 > 0.3 {print $1}' >  ${tmpdir}/chr${chr}_info_0.3.snplist

cd ${tmpdir}

plink2 \
  --pfile chr${chr}_extractedID_extract \
  --extract chr${chr}_info_0.3.snplist \
  --out chr${chr}_extracted_extract

done
done
```

merge all chromosomes
```shell
mkdir /PHShome/ys724/scratch/IDsubset/merge/

folders=("0410" "0411" "0412" "0413")

for fl in "${folders[@]}"; do
cd /PHShome/ys724/scratch/IDsubset/${fl}
	ls chr*.pgen | sed -e 's/.pgen//' > merge-list_${fl}.txt
	/PHShome/ys724/software/plink2 --pmerge-list merge-list_${fl}.txt --make-bed --out /PHShome/ys724/scratch/IDsubset/merge/${fl}_merge    
done
```

merge different folders
```shell
#find common variants

cd /PHShome/ys724/scratch/IDsubset/merge/

folders=("0410" "0411" "0412" "0413")

for fl in "${folders[@]}"; do
awk '{print $3}' ${fl}_merge.pvar > ${fl}_snplist.txt
sort ${fl}_snplist.txt > ${fl}_snplist_sorted.txt
done

comm -12 0410_snplist_sorted.txt 0411_snplist_sorted.txt > common_0410_0411.txt
comm -12 common_0410_0411.txt 0412_snplist_sorted.txt > common_0410_0411_0412.txt
comm -12 common_0410_0411_0412.txt 0413_snplist_sorted.txt > common_all.txt

#extract common variants from different folders

module load Plink/2.0

for fl in "${folders[@]}"; do

plink2 \
  --pfile ${fl}_merge \
  --extract common_all.txt \
  --make-bed \
  --out ${fl}_common

done

#merge different folders

module load Plink/1.9

plink

```

QC
```
cd /PHShome/ys724/scratch/IDsubset/merge/

module load Plink/2.0

folders=("0410" "0411" "0412" "0413")

for fl in "${folders[@]}"; do

#Genotype QC: MAF > 1%, genotype missing rates < 5%, p-value of Hardy-Weinberg equilibrium > 1e-6
plink2 \
  --pfile ${fl}_merge \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-6 \
  --make-bed \
  --out /PHShome/ys724/scratch/IDsubset/merge/${fl}_geno_filtered

#Remove individuals with missing rates > 10%
plink2 \
  --bfile ${fl}_geno_filtered \
  --mind 0.1 \
  --make-bed
  --out /PHShome/ys724/scratch/IDsubset/merge/${fl}_geno_mind__filtered

#Remain only unrelated individuals with SNP-derived genetic relatedness < 0.05
plink2 \
  --bfile ${fl}_geno_mind__filtered \
  --genome \
  --min 0.05 \
  --out pihat_min0.05

awk '{if($8>0.9) print $0}' pihat_min0.2.genome  > zoom_pihat.genome

plink --bfile HapMap_3_r3_10 --filter-founders --make-bed --out HapMap_3_r3_11

plink --bfile reduced_missing --genome --min 0.05 --out ibd_result
XXXXXXXXXX



```
merge different folders
```
cd /PHShome/ys724/scratch/IDsubset/merge/

folders=("0410" "0411" "0412" "0413")

for fl in "${folders[@]}"; do
awk '{print $3}' ${fl}_merge.pvar > ${fl}_snplist.txt
sort ${fl}_snplist.txt > ${fl}_snplist_sorted.txt
done

comm -12 0410_snplist_sorted.txt 0411_snplist_sorted.txt > common_0410_0411.txt
comm -12 common_0410_0411.txt 0412_snplist_sorted.txt > common_0410_0411_0412.txt
comm -12 common_0410_0411_0412.txt 0413_snplist_sorted.txt > common_all.txt
```

### Project PC
```
cd /PHShome/ys724/scratch/pca/merge/
mkdir /PHShome/ys724/scratch/pca/projectpc

plink2 \
  --pfile 0410_merge \
  --score /PHShome/ys724/Documents/GMBI_endometriosis/pca/hgdp_tgp_pca_gbmi_snps_loadings.GRCh38.plink.tsv \
  variance-standardize \
  cols=-scoreavgs,+scoresums \
  list-variants \
  header-read \
  --score-col-nums 3-22 \
  --read-freq /PHShome/ys724/Documents/GMBI_endometriosis/pca/hgdp_tgp_pca_gbmi_snps_loadings.GRCh38.plink.afreq \
  --out /PHShome/ys724/scratch/pca/projectpc/0410_score
```

### Assign pop
```R
library(randomForest)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(data.table)

setwd("/Users/shiyue//Desktop/Ancestry/")
ref <- fread("gnomad_meta_hgdp_tgp_v1.txt") #4150 individuals
ref1 <- fread("hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz") # 3327 individuals
refs <- merge(x = ref, y = ref1, by = "s") #3105 individuals
colnames(refs)[colnames(refs) == "hgdp_tgp_meta.Genetic.region"] <- "superpop"
colnames(refs)[colnames(refs) == "hgdp_tgp_meta.Project"] <- "Project"
nvars <- nrow(fread("0410_score.sscore.vars", header = F))
proj <- fread("0410_score.sscore")
names(proj)[5:24] <- paste0("PC", 1:20)
proj[,superpop := "unknown"]
proj[,Project := "MGB"]
pcvecs <- paste0("PC", seq(10))
proj[, (pcvecs) := lapply(.SD, function(d) d/sqrt(nvars)), .SDcols = pcvecs]

pop_forest <- function(training_data, data, ntree=100, seed=2, pcs=1:npc) {
  set.seed(seed)
  form <- formula(paste('as.factor(superpop) ~', paste0('PC', pcs, collapse = '+' )))
  forest <- randomForest(form,
                         data = training_data,
                         importance = T,
                         ntree = ntree)
  print(forest)
  fit_data <- data.frame(predict(forest, data, type='prob'), sample = data$sample)
  fit_data %>%
    gather(predicted_pop, probability, -sample) %>%
    group_by(sample) %>%
    slice(which.max(probability))
}
trdat <- refs[,c("s", "superpop", paste0("PC", 1:20))]
names(trdat)[1] <- "sample"
tedat <- proj[,c("IID", paste0("PC", 1:20))]
names(tedat)[1] <- "sample"

npc=6
pop_pred <- as.data.table(pop_forest(training_data = trdat, data = tedat))
table(pop_pred$predicted_pop)
write.csv(pop_pred, "/Users/shiyue/Desktop/Ancestry/pop_pred.csv")
```

