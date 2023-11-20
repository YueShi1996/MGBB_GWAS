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

ls *_common.bim | sed -e 's/.bim//' > merge.txt

plink --merge-list merge.txt --allow-no-sex --make-bed --out MGBB_EUR
```

### Quality Control
```shell
cd /PHShome/ys724/scratch/IDsubset/merge/

module load Plink/2.0

folders=("0410" "0411" "0412" "0413")

for fl in "${folders[@]}"; do

#Genotype QC: MAF > 1%, genotype missing rates < 5%, p-value of Hardy-Weinberg equilibrium > 1e-6
plink2 \
  --pfile MGBB_EUR \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-6 \
  --make-bed \
  --out /PHShome/ys724/scratch/IDsubset/merge/MGBB_EUR_geno_filtered

#Remove individuals with missing rates > 10%
plink2 \
  --bfile MGBB_EUR_geno_filtered \
  --mind 0.1 \
  --make-bed \
  --out /PHShome/ys724/Documents/GMBI_endometriosis/MGBB/MGBB_EUR_final
```
### GWAS
```shell
module load Plink/2.0

cd /PHShome/ys724/Documents/GMBI_endometriosis/MGBB
mkdir /PHShome/ys724/Documents/GMBI_endometriosis/assoc/

folders=("W" "Wex" "SCNv1" "SCNv2" "PCNv1" "PCNv2")

for fl in "${folders[@]}"; do

plink2 \
  --bfile MGBB_EUR \
  --glm hide-covar\
  --allow-no-sex \
  --pheno /PHShome/ys724/Documents/GMBI_endometriosis/phenotype/${fl}_EUR_pheno_modified.txt \
  --covar /PHShome/ys724/Documents/GMBI_endometriosis/phenotype/${fl}_EUR_age_cleaned.txt \
  --out /PHShome/ys724/Documents/GMBI_endometriosis/assoc/MGBB_${fl}_EUR

done
```
plot manhattan and QQ plot, screen out variants with gwas statistically significance
```R
args <- commandArgs(TRUE)

assoc_unadj <- read.table(paste(args[1],".endo.glm.logistic",sep=""),header=TRUE)
assoc_unadj_final <- na.omit(assoc_unadj)
colnames(assoc_unadj_final) <- c( "CHR", "BP","SNP","REF","ALT","A1","TEST","OBS_CT","OR","LOG(OR)_SE","Z_STAT","P")
assoc_unadj_final_sig <- assoc_unadj_final[assoc_unadj_final[,12] < 0.000001,]

# Save assoc_unadj_final_sig to a text file
output_filename <- paste(args[1], "_sig.txt", sep="")
write.table(assoc_unadj_final_sig, file=output_filename, row.names=FALSE, sep="\t")

library(qqman)
png(paste(args[1],".manhattan.png",sep=""))
manhattan(assoc_unadj_final,main="Manhattan Plot")
dev.off()

png(paste(args[1],".qq.png",sep=""))
qq(assoc_unadj_final$P, main="Q-Q Plot")
dev.off()
```
```shell
folders=("MGBB_W_EUR" "MGBB_Wex_EUR" "MGBB_PCNv1_EUR" "MGBB_PCNv2_EUR" "MGBB_SCNv1_EUR" "MGBB_SCNv2_EUR")

module load R/3.5.1-foss-2018b

cd /PHShome/ys724/Documents/GMBI_endometriosis/assoc

for fl in "${folders[@]}"; do

Rscript glm_result.R ${fl} 2> error_${fl}_log.txt

done
```
