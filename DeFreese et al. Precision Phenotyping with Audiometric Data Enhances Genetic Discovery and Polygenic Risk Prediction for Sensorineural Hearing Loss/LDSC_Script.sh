###############################################
# AudBase Project - Validation analyses of the SNHL GWAS signals
# Project Leaders: AD and TRM
# Author of the Script: Tanguy Rubat du Mérac
# Last edited on: April 3rd, 2026
###############################################

# Make sure paths are anonymized for public sharing
path2GWAS_AudBase_GWASes=/your_path/AudBase/
path2GWAS_HL_Trpchevska_metaGWAS=/your_path/AudBase/Published_HL_Sumstats/
path2GWAS_hearing_aid_use_Liu=/your_path/AudBase/Published_HL_Sumstats/hearing_aid_use_liu_cleaned.sumstats
path2GWAS_SPiN_Liu=/your_path/AudBase/Published_HL_Sumstats/speech_in_noise_liu_cleaned.sumstats
path2GWAS_PTA_Nagtegaal=/your_path/AudBase/PTA_Nagtegaal_2019/

wdir=/your_path/AudBase/LDSC/AudBase_Manuscript
cd $wdir 

# Set up LDSC env
ml StdEnvACCRE/2023 miniconda3/23.9.0-0 scipy-stack

# Change directory to LDSC folder:
path2LDSC=/your_path/AudBase/ldsc_updt

# Activate virtual environment
source activate ldsc

##########################################
# Munge Summary Statistics Files
##########################################

# AudBase ICD-based EUR postlingual GWAS
python $path2LDSC/munge_sumstats.py \
    --out $wdir/SNHL_EUR_ICD_postlingual_REGENIE_with_rsID_reformat.txt \
    --N 61468 \
    --a1 A1 \
    --a2 A2 \
    --snp rsID \
    --p P \
    --chunksize 500000 \
    --merge-alleles /data/p_gordon_lab/users/rubatdtb/LDSC_analyses_CROISER/eur_w_ld_chr/w_hm3.snplist \
    --sumstats $path2GWAS_AudBase_GWASes/SNHL_EUR_ICD_postlingual_REGENIE_with_rsID_reformat.txt

# AudBase PTA-based EUR postlingual GWAS
python $path2LDSC/munge_sumstats.py \
    --out $wdir/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt \
    --N 16057 \
    --a1 A1 \
    --a2 A2 \
    --snp rsID \
    --p P \
    --chunksize 500000 \
    --merge-alleles /data/p_gordon_lab/users/rubatdtb/LDSC_analyses_CROISER/eur_w_ld_chr/w_hm3.snplist \
    --sumstats $path2GWAS_AudBase_GWASes/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt

# Hearing Loss meta-GWAS in EUR
python $path2LDSC/munge_sumstats.py \
    --out $wdir/Trpchevska_et_al_2022_Hearing_Loss_MetaGWAS_munged.txt \
    --N 723266 \
    --a1 Allele1 \
    --a2 Allele2 \
    --snp SNP \
    --p P.value \
    --chunksize 500000 \
    --merge-alleles $path2LDSC/eur_w_ld_chr/w_hm3.snplist \
    --sumstats $path2GWAS_HL_Trpchevska_metaGWAS

# PTA GWAS in EUR from Nagtegaal et al. 2019 (Low Frequency only)
python $path2LDSC/munge_sumstats.py \
    --out $wdir/Nagtegaal_et_al_2019_PTA_EUR_LowFreq_munged.txt \
    --N 9675 \
    --a1 EFFECT_ALLELE \
    --a2 OTHER_ALLELE \
    --snp SNPID \
    --p PVAL \
    --chunksize 500000 \
    --ignore P \
    --merge-alleles $path2LDSC/eur_w_ld_chr/w_hm3.snplist \
    --sumstats $path2GWAS_PTA_Nagtegaal/CLEANED_META_LOW_20160107.with_rsid.txt.gz

##########################################
# Calculate h2 on liability scale for AudBase GWASes
##########################################

# Compute the heritability of the AudBase GWAS PTA in EUR for postlingual 
python $path2LDSC/ldsc.py \
    --ref-ld-chr $path2LDSC/eur_w_ld_chr/ \
    --out $wdir/AudBase_GWAS_PTA_postlingual_EUR_h2 \
    --h2  $wdir/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt.sumstats.gz \
    --w-ld-chr $path2LDSC/eur_w_ld_chr/

# Compute the heritability of the AudBase GWAS ICD in EUR for postlingual
python $path2LDSC/ldsc.py \
    --ref-ld-chr $path2LDSC/eur_w_ld_chr/ \
    --out $wdir/AudBase_GWAS_ICD_postlingual_EUR_h2 \
    --h2  $wdir/SNHL_EUR_ICD_postlingual_REGENIE_with_rsID_reformat.txt.sumstats.gz \
    --w-ld-chr $path2LDSC/eur_w_ld_chr/

##########################################
# Run rg between AudBase PTA Postlingual EUR and HL Meta-GWAS
##########################################

# Run rg between AudBase PTA Postlingual EUR and HL Meta-GWAS
python $path2LDSC/ldsc.py \
    --ref-ld-chr $path2LDSC/eur_w_ld_chr/ \
    --out $wdir/AudBase_PTA_Postlingual_EUR_vs_HL_MetaGWAS_rg \
    --rg $wdir/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt.sumstats.gz,$path2GWAS_HL_Trpchevska_metaGWAS \
    --w-ld-chr $path2LDSC/eur_w_ld_chr/

# Run rg between AudBase PTA Postlingual EUR and Liu Hearing Aid Use GWAS
python $path2LDSC/ldsc.py \
    --ref-ld-chr $path2LDSC/eur_w_ld_chr/ \
    --out $wdir/AudBase_PTA_Postlingual_EUR_vs_Liu_Hearing_Aid_Use_rg \
    --rg $wdir/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt.sumstats.gz,$path2GWAS_hearing_aid_use_Liu \
    --w-ld-chr $path2LDSC/eur_w_ld_chr/

# Run rg between AudBase PTA Postlingual EUR and Liu Speech in Noise GWAS
python $path2LDSC/ldsc.py \
    --ref-ld-chr $path2LDSC/eur_w_ld_chr/ \
    --out $wdir/AudBase_PTA_Postlingual_EUR_vs_Liu_Speech_in_Noise_rg \
    --rg $wdir/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt.sumstats.gz,$path2GWAS_SPiN_Liu \
    --w-ld-chr $path2LDSC/eur_w_ld_chr/

# Run rg between AudBase PTA Postlingual EUR and Nagtegaal PTA GWAS (Low Frequency only)
python $path2LDSC/ldsc.py \
    --ref-ld-chr $path2LDSC/eur_w_ld_chr/ \
    --out $wdir/AudBase_PTA_Postlingual_EUR_vs_Nagtegaal_PTA_EUR_LowFreq_rg \
    --rg $wdir/SNHL_PTA_EUR_postlingual_REGENIE_with_rsID_reformat.txt.sumstats.gz,$wdir/Nagtegaal_et_al_2019_PTA_EUR_LowFreq_munged.txt.sumstats.gz \
    --w-ld-chr $path2LDSC/eur_w_ld_chr/
