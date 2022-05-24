#!/bin/bash
threads=4
#cd /DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/Tongwu_germline/variant_classification
cd /DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification
#/DCEG/Projects/Exome/builds/build_Osteosarcoma_2017_18794/Results/intervar/ExAC_nonTCGA/CPG #/DCEG/Projects/Exome/builds/build_Guatemala_2017_19009/Results/Eric/CPG_v1
module load sge
module load python3/3.5.1
module load R/3.4.0
mkdir -p logs
#snakemake --cluster "qsub -q research.q -pe by_node ${threads} -o /DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/Tongwu_germline/variant_classification/logs/ -e /DCEG/Projects/Exome/builds/build_Hong-Kong-Breast-Cancer_Exome_2017_18839/Tongwu_germline/variant_classification/logs/" --jobs 6000 --latency-wait 300

snakemake --cluster "qsub -q seq-calling2.q -pe by_node ${threads} -o /DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/logs/ -e /DCEG/Projects/Exome/builds/build_pediatric_sarcoma_exome_2021_29739/variant_classification/logs/" --jobs 6000 --latency-wait 100
