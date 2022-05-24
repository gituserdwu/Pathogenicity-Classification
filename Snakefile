#!/usr/bin/python

import glob
import sys
import os
import shutil
import codecs
import pysam
from snakemake.utils import R



proj = 'pediatric_sarcoma' #'CPG'

# pheno_file = '/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/Eric/pheno/OS.exome.pheno.csv'

#inVcf = '/DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/Tongwu_germline/variants_InterVar_annotated.vcf.gz'
#inVcf = '/DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/TCGA_Breast_cancer/variants_InterVar_annotated.vcf.gz'
#inVcf = '/DCEG/Projects/CCSS/steve/shareXfer/karlinsCV/variants_gene_lists_plusMissing_REVEL_InterVar_ClinVar.extCV.byDate.vcf.gz'
inVcf = '/DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/variants_InterVar_annotated.vcf.gz'
#inBed = '/DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/Bed_file/Exome_v3_UTR_merged.bed'
#inBed = '/DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/TCGA_Breast_cancer/new_annotation_revel_clinvar/new_variant_classification/no_chr_coordinated_gene_list_20181012_merged_sorted.bed'
#inBed = '/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/Eric/gene_files/253cancerPredisGenes.sorted.merged.canonical.bed'

inBed = '/DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/seqcap_EZ_Exome_v3_v3utr_HG19_vcrome2.1_intersect_correct.bed'

#vcf_fields_file = '/DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/TCGA_Breast_cancer/new_annotation_revel_clinvar/new_variant_classification/vcf_fields_to_keep.txt'
#vcf_fields_file = '/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/Eric/vcf_fields_to_keep_v4.txt'
vcf_fields_file = '/DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/vcf_fields_to_keep_v4_updated.txt'

#gene_file = '/DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/TCGA_Breast_cancer/annotation_revel_clinvar/variant_classification/gene_list_of_185genes.txt'
#gene_file = '/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/Eric/gene_files/245CPG_HGNC_InherMode.txt'
gene_file = '/DCEG/Projects/Exome/SequencingData/BED_FILES/NimbleGen_V3_download/sorted_unique_exome_gene.txt'

rule all:
    input:
        'subset_vcf/pediatric_sarcoma.genesOfInterest.rare.nonsyn.vcf',
        'subset_vcf/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt',
        'impact_filtered/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt',
        'classify_variant/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt',
        'delivery/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt'
        #'delivery/HKBC.genesOfInterest.rare.nonsyn.txt'
        #'subset_vcf/HKBC.genesOfInterest.rare.nonsyn.txt'


rule subset_vcf:
    input:
        vcf = inVcf,
        bed = inBed
    output:
        'subset_vcf/pediatric_sarcoma.genesOfInterest.rare.nonsyn.vcf'
    params:
        script = '/DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/scripts/subsetVcfMultiSplit_v2.py'
   
    shell:
        'module load python/2.7.5;python {params.script} {input.vcf} {input.bed} {output}'


rule vcf_to_txt:
    input:
        vcf = 'subset_vcf/pediatric_sarcoma.genesOfInterest.rare.nonsyn.vcf',
        fields_to_keep = vcf_fields_file
    output:
        'subset_vcf/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt'
    params:
        gatk = '/DCEG/Projects/Exome/SequencingData/GATK_binaries/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar',
        #gatk = '/DCEG/Projects/Exome/SequencingData/GATK_binaries/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar',
        #ref = '/CGF/Resources/Data/genome/reference_wgs_b37/human_g1k_v37.fasta'
        ref = '/CGF/Resources/Data/genome/hg19_canonical_correct_chr_order.fa'
    run:
        fields = ''
        with open(input.fields_to_keep) as f:
            for line in f:
                if line.strip():
                    fields += ' ' + line.rstrip()
        shell('module load jdk/1.8.0_111;java -Xmx10g -jar {params.gatk} -R {params.ref} -T VariantsToTable -raw -V {input.vcf}' + fields + ' -o {output}')
        #shell('module load jdk/0.7.17;java -Xmx10g -jar {params.gatk} -R {params.ref} -T VariantsToTable -raw -AMD -V {input.vcf}' + fields + ' -o {output}')

'''
rule update_txt:
    input:
        txt = 'subset_vcf/OsteoExome.genesOfInterest.rare.nonsyn.txt',
        pheno = pheno_file
    output:
        'add_new_columns/OsteoExome.genesOfInterest.rare.nonsyn.txt'
    params:
        script = '/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/Eric/makeOsteoTxtSplitMulti.py'
    shell:
        'module load python/2.7.5;python {params.script} {input.txt} {input.pheno} {output}'
'''

rule impact_filter:
    input:
        txt = 'subset_vcf/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt',
        #txt = 'add_new_columns/OsteoExome.genesOfInterest.rare.nonsyn.txt',
        gene = gene_file
    output:
        'impact_filtered/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt'
    params:
        script = '/DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/scripts/impactFilterAddCols_v3.py' 
        #script = '/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/Eric/impactFilterAddCols_v3.py'
    shell:
        'module load python/2.7.5;python {params.script} {input.txt} {output} {input.gene}'

rule classify_variants:
    input:
        txt = 'impact_filtered/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt',
        gene = gene_file
    output:
        'classify_variant/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt'
    params:
        script = '/DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/scripts/OSTEOSARCOMA_pathogenicity_v2.py'
        #script = '/DCEG/CGF/Bioinformatics/Production/Eric/scripts/OS_pathogenicity_v2.py'
    shell:
        'module load python/2.7.5;python {params.script} {input.txt} {output} {input.gene}'


rule remove_dups:
    input:
        'classify_variant/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt'
    output:
        'delivery/pediatric_sarcoma.genesOfInterest.rare.nonsyn.txt'
    run:
        varToSetDict = {}
        with open(input[0]) as f:
            head = f.readline()
            head_list = head.rstrip().split('\t')
            setCol = None
            for i in range(len(head_list)):
                if head_list[i] == 'set':
                    setCol = i
            if setCol == None:
                print('no column called "set" in input file')
                sys.exit(1)
            line = f.readline()
            while line != '':
                line_list = line.rstrip().split('\t')
                (chrom, pos, snp, ref, alt) = line_list[:5]
                varId = '_'.join([chrom, pos, ref, alt])
                varSet = line_list[setCol]
                if not varToSetDict.get(varId):
                    varToSetDict[varId] = [varSet]
                else:
                    varToSetDict[varId].append(varSet)
                line = f.readline()
        dupVarDict = {}
        for key in varToSetDict.keys():
            setList = varToSetDict[key]
            if len(setList) > 1:
                dupVarDict[key] = setList
        with open(input[0]) as f, open(output[0], 'w') as out:
            head = f.readline()
            out.write(head.rstrip() + '\tisMulti\n')
            head_list = head.rstrip().split('\t')
            setCol = None
            filtCol = None
            indCol = None
            annCol = None
            for i in range(len(head_list)):
                if head_list[i] == 'set':
                    setCol = i
                elif head_list[i] == 'FILTER':
                    filtCol = i
                elif head_list[i] == 'IND':
                    indCol = i
                elif head_list[i] == 'ANN':
                    annCol = i
            if setCol == None or filtCol == None or indCol == None or annCol == None:
                print('Need columns "set" and "FILTER" in input file')
                sys.exit(1)
            line = f.readline()
            while line != '':
                line_list = line.rstrip().split('\t')
                oldId = line_list[5]
                if len(oldId.split('_')) > 4:
                    multi = 'Y'
                else:
                    multi = 'N'
                (chrom, pos, snp, ref, alt) = line_list[:5]
                ind = line_list[indCol]
                if len(ind) > 30000:
                    line_list[indCol] = 'TooManyToShow'
                ann = line_list[annCol]
                if len(ann) > 30000:
                    line_list[annCol] = ann[:30000]
                line_list.append(multi)
                varId = '_'.join([chrom, pos, ref, alt])
                varSet = line_list[setCol]
                if dupVarDict.get(varId) and varSet == 'fb':
                    print('Found dup for ' + line[:100])
                elif dupVarDict.get(varId) and varSet != 'fb':
                    line_list[filtCol] = 'PASS'
                    out.write('\t'.join(line_list) + '\n')
                else:
                    out.write('\t'.join(line_list) + '\n')
                line = f.readline()
