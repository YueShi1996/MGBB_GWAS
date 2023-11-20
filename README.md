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

QC
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
  --out /PHShome/ys724/scratch/IDsubset/merge/MGBB_EUR_geno_mind__filtered
```
