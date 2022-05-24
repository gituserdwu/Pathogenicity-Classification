#!/usr/bin/env python

import sys


popFreqColsNames = '''ExAC_nonTCGA_AF_AMR
ExAC_nonTCGA_AF_AFR
ExAC_nonTCGA_AF_OTH
ExAC_nonTCGA_AF_Adj
ExAC_nonTCGA_AF_SAS
ExAC_nonTCGA_AF_FIN
ExAC_nonTCGA_AF_EAS
ExAC_nonTCGA_AF_NFE
ESP_EA_MAF
ESP_EA_AF
ESP_ALL_MAF
ESP_AA_MAF
ESP_AA_AF
KG_AMR_AF
KG_EUR_AF
KG_EAS_AF
KG_SAS_AF
KG_AF
KG_AFR_AF'''.split()

popColNameDict = {}
for pop in popFreqColsNames:
    popColNameDict[pop] = 1


otherColNames = '''highest_impact
effect
SnpEffGene
badged_NH_TBbD_call
HGMD_CLASS
InterVarSnpEffGene
REVEL_SCORE
dbNSFP_CADD_phred
dbNSFP_MetaSVM_pred'''.split()
# AC_NFE'''.split()
# num_cases
# num_controls'''.split()

otherNameDict = {}
for name in otherColNames:
    otherNameDict[name] = 1

def makeGeneArDict(geneTxt):
    geneDict = {}
    with open(geneTxt) as f:
        for line in f:
            (gene, ar) = line.split()
            geneDict[gene] = ar
    return geneDict



def getColDict(colNameDict, allColNameList):
    '''
    (dict) -> dict
    return a dict where the keys are the column numbers for the names in the list in that order
    '''
    colDict = {}
    for i in range(len(allColNameList)):
        colName = allColNameList[i]
        if colNameDict.get(colName):
            colDict[colName] = i
    for colName in colNameDict.keys():
        if not colDict.get(colName):
            print(colName + ' not found in header')
            sys.exit(1)
    return colDict


def getMaxPopFreq(popFreqColDict, popColNameDict, line_list):
    PopMax = 0.0
    for colName in popColNameDict.keys():
        col = popFreqColDict[colName]       
        freqStr = line_list[col]
        if ',' in freqStr:
            print('there is a "," in ' + colName  + ' in line:  ' + '\t'.join(line_list[:6]))
            sys.exit(1)
        if freqStr != '.' and freqStr != 'NA':
            if float(freqStr) > 0.5:
                maf = 1.0 - float(freqStr)
            else:
                maf = float(freqStr)
            if maf > PopMax:
                PopMax = maf
    return PopMax


def classifyVariant(PopMax, PopMaxThresh, InterVar, revelScore, revelThresh, ClinVar, HGMD, highest_impact, SnpEffAnn, numCases, numControls, totCountThresh, caddScore, caddThresh, metaSVMpred):
    '''
    (float, float, str, float, float, str, str, str, str, int, int, int) -> (str, str)                                                                             
    Return (variantClasiffication, Reason) using the input annotations.                                                                                            
    '''
    hgmd_list = HGMD.split(',')
    cadd_list = []
    if caddScore == 'NA':
        cadd_list.append(None)
    else:
        for caddStr in caddScore.split(','):
            cadd_list.append(float(caddStr))
    metaSVM_list = metaSVMpred.split(',')
    prediction = 'VUS_Unc'
    reason = 'Does not satisfy other rules'
    if ClinVar != 'NA':
        prediction = ClinVar
        reason = 'CinVar ' + prediction
    elif InterVar == 'Pathogenic':
        prediction = 'P'
        reason = 'ClinVar NA and InterVar P'
    elif InterVar == 'Likely_pathogenic':
        prediction = 'LP'
        reason = 'ClinVar NA and InterVar LP'
    elif InterVar == 'Likely_benign':
        prediction = 'LB'
        reason = 'ClinVar NA and InterVar LB'
    elif InterVar == 'Benign':
        prediction = 'B'
        reason = 'ClinVar NA and InterVar B'
    elif InterVar== 'Uncertain_significance':
        prediction = 'VUS'
        reason = 'ClinVar NA and InterVar U'
    else:
        prediction = 'VUS'
        reason = 'ClinVar NA and InterVar ' + InterVar
    if (prediction == 'P' or prediction == 'LP'):
        if PopMax > PopMaxThresh:
            prediction = 'VUS_Unc'
            reason += ' and PopMax > threshold'
        elif (numCases + numControls) > totCountThresh:
            prediction = 'VUS_Unc'
            reason += ' and Total Count > threshold'
        else:
            reason += ' and PopMax <= threshold and Total Count <= threshold'
    elif prediction == 'VUS':
        if PopMax <= PopMaxThresh and (numCases + numControls) <= totCountThresh:
            if 'DM' in hgmd_list and highest_impact == 'HIGH' and ('frameshift_variant' in SnpEffAnn or 'stop_gained' in SnpEffAnn or 'stop_lost' in SnpEffAnn or 'start_lost' in SnpEffAnn or 'splice_acceptor_variant' in SnpEffAnn or 'splice_donor_variant' in SnpEffAnn):
                prediction = 'LP'
                reason += ' and Freq <= threshold and Classic high impact and HGMD DM'
            elif highest_impact == 'HIGH' and ('frameshift_variant' in SnpEffAnn or 'stop_gained' in SnpEffAnn or 'stop_lost' in SnpEffAnn or 'start_lost' in SnpEffAnn or 'splice_acceptor_variant' in SnpEffAnn or 'splice_donor_variant' in SnpEffAnn):
                prediction = 'VUS_Pop'
                reason += ' and Freq <= threshold and Classic high impact'
            elif 'DM' in hgmd_list:
                prediction = 'VUS_Pop'
                reason += ' and Freq <= threshold and HGMD DM'
            elif revelScore >= revelThresh and max(cadd_list) >= caddThresh and 'D' in metaSVM_list:
                prediction = 'VUS_Pop'
                reason += ' and Freq <= threshold and in silico check > threshold'
            else:
                prediction = 'VUS_Unc'
                reason += ' and Freq <= threshold and fails other checks'
        else:
            prediction = 'VUS_Unc'
            reason += ' and Freq > threshold'
    return (prediction, reason)


def outputTxt(inputTxt, outputTxt, geneFile, revelThresh = 0.5, ArFreqThresh = 0.005, otherGeneFreqThresh = 0.001, totCountThresh = 10, caddThresh = 20):
    ArDict = makeGeneArDict(geneFile)
    with open(inputTxt) as f, open(outputTxt, 'w') as output:
        head = f.readline()
        output.write(head.rstrip('\n') + '\tnum_total\tPopMax\tAutoClassify\tAutoReason\n')
        head_list = head.rstrip('\n').split('\t')
        popFreqColDict = getColDict(popColNameDict, head_list)
        otherColDict = getColDict(otherNameDict, head_list)
        geneCol = otherColDict['SnpEffGene']
        impactCol = otherColDict['highest_impact']
        effCol = otherColDict['effect']
        clnSigCol = otherColDict['badged_NH_TBbD_call']
        hgmdCol = otherColDict['HGMD_CLASS']
        interCol = otherColDict['InterVarSnpEffGene']
        revelCol = otherColDict['REVEL_SCORE']
        caddCol = otherColDict['dbNSFP_CADD_phred']
        metaSvmCol = otherColDict['dbNSFP_MetaSVM_pred']
        numCaseCol = -1; #No column number, bypass this filter. #For ExAC, no cases. #otherColDict['num_cases']
        numContCol = -1; #No column number, bypass this filter. #otherColDict['AC_NFE'] #Consider NFE here.  otherColDict['num_controls']
        line = f.readline()
        while line != '':
            line_list = line.rstrip('\n').split('\t')
            gene = line_list[geneCol]
            numCases = 0 #always smaller than 10, bypass this filter. #For ExAC int(line_list[numCaseCol])
            numControls = 0 #always smaller than 10, bypyass this filter. #int(line_list[numContCol])
            numTotal = numCases + numControls
            if gene == 'NA':
                PopMax = 'NA'
                variantClass = 'NA'
            else:
                PopMax = getMaxPopFreq(popFreqColDict, popColNameDict, line_list)
                if ArDict[gene] == 'AR':
                    PopMaxThresh = ArFreqThresh
                else:
                    PopMaxThresh = otherGeneFreqThresh
                InterVar = line_list[interCol]
                revelScoreStr = line_list[revelCol]
                if revelScoreStr == 'NA' or revelScoreStr == '.':
                    revelScore = None
                else:
                    revelScore = float(revelScoreStr)
                caddScore = line_list[caddCol]
                metaSVMpred = line_list[metaSvmCol]
                ClinVar = line_list[clnSigCol]
                HGMD = line_list[hgmdCol]
                highest_impact = line_list[impactCol]
                SnpEffEffect = line_list[effCol]
                (variantClass, reason) = classifyVariant(PopMax, PopMaxThresh, InterVar, revelScore, revelThresh, ClinVar, HGMD, highest_impact, SnpEffEffect, numCases, numControls, totCountThresh, caddScore, caddThresh, metaSVMpred)
            output.write(line.rstrip('\n'))
            output.write('\t' + str(numTotal) + '\t' + str(PopMax))
            output.write('\t' + variantClass + '\t' + reason + '\n')
            line = f.readline()



def main():
    args = sys.argv[1:]
    if len(args) != 3:
        print ("error: usage: python OS_pathogenicity_v1.py /path/to/in.txt /path/to/out.txt /path/to/AR_gene_file.txt")
        sys.exit(1)
    else:
        outputTxt(args[0], args[1], args[2])
    




if __name__ == "__main__":
    main()
