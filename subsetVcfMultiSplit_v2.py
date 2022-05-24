import sys
import pysam
import os

popFreqCols = '''KG_AMR_AF
KG_EUR_AF
KG_EAS_AF
KG_SAS_AF
KG_AF
KG_AFR_AF'''.split()


popFreqColDict = {}
for pop in popFreqCols:
    popFreqColDict[pop] = 1

def keepLine(line, maxPopFreq = .01):
    line_list = line.split()
    alt = line_list[4]
    info = line_list[7]
    info_list = info.split(';')
    PopMax = 0.0
    impacts = []
    for i in info_list:
        if i.startswith('ANN='):
            snpeff = i.split('=')[1]
            if 'HIGH' not in snpeff and 'MODERATE' not in snpeff and 'LOW' not in snpeff:
                return False
        elif popFreqColDict.get(i.split('=')[0]):
            freqStr = i.split('=')[1]
            if ',' in freqStr:
                print('there is a "," in ' + i.split('=')[0] + ' in line:  ' + line[:100])
                sys.exit(1)
            if freqStr != '.':
                if float(freqStr) > 0.5:
                    maf = 1.0 - float(freqStr)
                else:
                    maf = float(freqStr)
                if maf > PopMax:
                    PopMax = maf
    if PopMax < maxPopFreq:
        return True
    return False

def UpdateOutVcf(outVCF):
    #outVCF = "/home/lixin/lxwg/Data/ad-hoc/VCFClusification/OSTEOSARCOMA.genesOfInterest.rare.nonsyn.vcf"
    outVCFTmp = outVCF + ".tmp"
    outVCFOrg = outVCF + ".org"
    if os.path.exists(outVCF):
        if os.path.exists(outVCFTmp):
            CMD = "rm " + outVCFTmp
            os.system(CMD)
        if os.path.exists(outVCFOrg):
            CMD = "rm " + outVCFOrg
            os.system(CMD)
        CMD = "sed '/#CHROM/q' " + outVCF + " >> " + outVCFTmp
        os.system(CMD)
        CMD = "sed -n '/#CHROM/,$p' " + outVCF + " | tail -n +2 " + "| sort -k1,1V -k2,2n | uniq >> " + outVCFTmp
        os.system(CMD)
        CMD = "mv " + outVCF + " " + outVCFOrg
        os.system(CMD)
        CMD = "mv " + outVCFTmp + " " + outVCF
        os.system(CMD)

def outputSubset(inVcf, inBed, outVcf):
    vcf = pysam.TabixFile(inVcf)
    with open(inBed) as f, open(outVcf, 'w') as out:
        for line in vcf.header:
            out.write(line + '\n')
        for bed_line in f:
            bed_line_list = bed_line.split()
            chrom = bed_line_list[0]
            start = int(bed_line_list[1])
            end = int(bed_line_list[2])
            if chrom in vcf.contigs:
                for line in vcf.fetch(chrom, start, end):
                    if keepLine(line) == True:
                        out.write(line + '\n')
    vcf.close()
    # Update outVcf: sort based on position and remove the duplicate lines.
    UpdateOutVcf(outVcf)


def main():
    args = sys.argv[1:]
    if len(args) != 3:
        print ("error: usage: python subsetVcf.py /path/to/in.vcf.gz /path/to/in.bed /path/to/out.vcf")
        sys.exit(1)
    else:
        outputSubset(args[0], args[1], args[2])


if __name__ == "__main__":
    main()
