import sys

##update v2 by also adding a column IsClassicHighImpact


def makeRefGeneDict(refGeneFile = '/DCEG/CGF/Bioinformatics/Production/Dongjing/refFiles/refGeneToKnownGene.txt'):
    refGeneToKnownListDict = {}
    with open(refGeneFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            (known, ref) = line.rstrip('\n').split('\t')
            ref = ref.strip(",")
            if not refGeneToKnownListDict.get(ref):
                refGeneToKnownListDict[ref] = [known]
            else:
                refGeneToKnownListDict[ref].append(known)
            line = f.readline()
    return refGeneToKnownListDict



def getNewHead(head):
    head_list = head.rstrip().split('\t')
    newHead_list = []
    snpEffCol = None
    InterVarCol = None
    hgmdCol = None
    ClinVarCol = None
    for i in range(len(head_list)):
        field = head_list[i]
        newHead_list.append(field)
        if field == 'ANN':
            snpEffCol = i
            newHead_list.append('highest_impact')
            newHead_list.append('num_impact_transcripts')
            newHead_list.append('total_transcripts')
            newHead_list.append('HGVS.c')
            newHead_list.append('HGVS.p')
            newHead_list.append('effect')
            newHead_list.append('IsClassicHighImpact')
            newHead_list.append('SnpEffGene')
        elif field == 'InterVar':
            InterVarCol = i
            newHead_list.append('InterVarSnpEffGene')
        elif field == 'HGMD_CLASS':
            hgmdCol = i
        elif field == 'CLNSIG':
            ClinVarCol = i
    if snpEffCol == None or InterVarCol == None or hgmdCol == None:
        print('SnpEff or InterVar or HGMD_CLASS not found in header')
        sys.exit(1)
    newHead = '\t'.join(newHead_list) + '\n'
    return (newHead, snpEffCol, InterVarCol, hgmdCol, ClinVarCol)


def keepHGMD(hgmd_class):
    hgmd_class_list = hgmd_class.split(',')
    if 'DM' in hgmd_class_list or 'DM?' in hgmd_class_list:
        return True
    return False


def keepClinVar(clinsig):
    keep_list = '''Conflicting_interpretations_of_pathogenicity
Likely_pathogenic
Pathogenic
Pathogenic/Likely_pathogenic
Uncertain_significance'''.split()
    clinsig_list = clinsig.split(',')
    for prediction in keep_list:
        if prediction in clinsig_list:
            return True
    return False



def checkImpact(snpEff, geneDict):
    if snpEff == 'NA':
        return ('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')
    trans_list = snpEff.split(',')
    tot_trans = len(trans_list)
    impactDict = {}
    for trans in trans_list:
        ann_list = trans.split('|')
        effect = ann_list[1]
        impact = ann_list[2]
        gene = ann_list[3]
        hgvs_c = ann_list[9]
        hgvs_p = ann_list[10]
        if geneDict.get(gene):
            if not impactDict.get(impact):
                impactDict[impact] = [1, [hgvs_c], [hgvs_p], [effect], [gene]]
            else:
                impactDict[impact][0] += 1
                impactDict[impact][1].append(hgvs_c)
                impactDict[impact][2].append(hgvs_p)
                impactDict[impact][3].append(effect)
                impactDict[impact][4].append(gene)
    if 'HIGH' in impactDict.keys():
        highest_impact = 'HIGH'
    elif 'MODERATE' in impactDict.keys():
        highest_impact = 'MODERATE'
    elif 'LOW' in impactDict.keys():
        highest_impact = 'LOW'
    elif 'MODIFIER' in impactDict.keys():
        highest_impact = 'MODIFIER'
    else:
        highest_impact = 'GeneNotFound'
    if highest_impact == 'GeneNotFound':
        (num_impact, tot_trans, hgvs_c_all, hgvs_p_all, effect_all, foundGene) = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
    else:
        (num_impact, hgvs_c_list, hgvs_p_list, effect_list, foundGenes) = impactDict[highest_impact]
        hgvs_c_all = ';'.join(hgvs_c_list)
        hgvs_p_all = ';'.join(hgvs_p_list)
        foundGene = ';'.join(set(foundGenes))
        effect_all = ';'.join(set(effect_list))
    return (highest_impact, num_impact, tot_trans, hgvs_c_all, hgvs_p_all, effect_all, foundGene)


def getInterPredict(InterVar, SnpEffGene, refGeneToKnownDict):
    if InterVar == 'NA':
        return 'NA'
    predictions = []
    inter_var_list = InterVar.split(',')
    for inter in inter_var_list:
        gene_list = inter.split('|')
        for x in gene_list:
            (gene, predict) = x.split(':')
            if gene == SnpEffGene:
                return predict
            elif refGeneToKnownDict.get(gene) and SnpEffGene in refGeneToKnownDict[gene]:
                return predict
    return 'GeneNotFound'


def makeGeneDict(geneTxt):
    geneDict = {}
    with open(geneTxt) as f:
        for line in f:
            gene = line.split()[0].strip()
            geneDict[gene] = 1
    return geneDict

def outputTxt(inTxt, outTxt, geneTxt):
    geneDict = makeGeneDict(geneTxt)
    refGeneToKnownDict = makeRefGeneDict()
    with open(inTxt) as f, open(outTxt, 'w') as output:
        head = f.readline()
        (newHead, snpEffCol, InterVarCol, hgmdCol, ClinVarCol) = getNewHead(head)
        output.write(newHead)
        line = f.readline()
        while line != '':
            line_list = line.rstrip().split('\t')
            snpEff = line_list[snpEffCol]
            (highest_impact, num_impact, tot_trans, hgvs_c, hgvs_p, effect_all, foundGene) = checkImpact(snpEff, geneDict)
            if highest_impact == 'NA' or highest_impact == 'HIGH' or highest_impact == 'MODERATE' or ((highest_impact == 'LOW' and keepHGMD(line_list[hgmdCol])) or (highest_impact == 'LOW' and keepClinVar(line_list[ClinVarCol]))):
                if highest_impact == 'HIGH' and ('frameshift_variant' in effect_all or 'stop_gained' in effect_all or 'stop_lost' in effect_all or 'start_lost' in effect_all or 'splice_acceptor_variant' in effect_all or 'splice_donor_variant' in effect_all):
                    classicHigh = 'Y'
                else:
                    classicHigh = 'N'
                output.write('\t'.join(line_list[:snpEffCol + 1]) + '\t')
                output.write('\t'.join([highest_impact, str(num_impact), str(tot_trans), hgvs_c, hgvs_p, effect_all, classicHigh, foundGene]) + '\t')
                output.write('\t'.join(line_list[snpEffCol + 1:InterVarCol + 1]) + '\t')
                output.write(getInterPredict(line_list[InterVarCol], foundGene, refGeneToKnownDict) + '\t')
                output.write('\t'.join(line_list[InterVarCol + 1:]) + '\n')
            line = f.readline()


def main():
    args = sys.argv[1:]
    if len(args) != 3:
        print ("error: usage: python impactFilterAddCols.py /path/to/in.txt /path/to/out.txt /path/to/genes.txt")
        sys.exit(1)
    else:
        outputTxt(args[0], args[1], args[2])


if __name__ == "__main__":
    main()

