import csv

builtData = csv.DictReader(open("built_with_priors.tsv", "r"), delimiter="\t")


"""Counting things in sean's file to compare to built.tsv
Only problem is you can't expect them to really add up because there are different variants in BRCAEx vs. HCI
Thus, below counts are useless"""
noCount = 0
yesCount = 0

variantData = csv.DictReader(open("seanFile.txt", "r"), delimiter="\t") #named variantData to resemble calcVarPriors.py
for item in variantData:
    if item['protein_prior'] == '0.03' and item['removesKeyDomain'] == 'Y':
        yesCount += 1
    if item['protein_prior'] == '0.03'  and item['removesKeyDomain'] == 'N':
        noCount += 1
#print 'Y: ' + str(yesCount)
#print 'N: ' + str(noCount)

"""Below begins the counts to check built.tsv (which has priors data) for ?????? """




def checkVCsForProtPriorAmbig():
    """Checks Var consequences for substitution_variant or missense_variant to see how successful the ambiguity resolution using VCs would be"""
    countVCFails = 0
    VCFails = []
    total = 0
    countIntron = 0
    builtData = csv.DictReader(open("built_with_priors.tsv", "r"), delimiter="\t")
    for item in builtData:
        if item['proteinPrior'] ==  '0.02':
            total += 1
            if 'synonymous_variant' in item['varConsequences']:
                continue
                #print item['varConsequences']
            elif 'missense_variant' in item['varConsequences']:
                continue
                #print item['varConsequences']
            #elif '-' in item['pyhgvs_cDNA']:
                #print item['pyhgvs_cDNA']
                #countIntron += 1
            #elif '+' in item['pyhgvs_cDNA']:
                #print item['pyhgvs_cDNA']
                #countIntron += 1
            else:
                 countVCFails += 1
                 VCFails.append(item['pyhgvs_cDNA'])
    #print countIntron
    print countVCFails
    print total
    percentError = float(countVCFails) / float(total)
    print percentError
    for v in VCFails:
        print v


def checkVCsForDNAmbig():
    """Checks Var consequences for substitution_variant or missense_variant to see how successful the ambiguity resolution using VCs would be"""
    builtData = csv.DictReader(open("built_with_priors.tsv", "r"), delimiter="\t")
    total = 0
    for item in builtData:
        if item['deNovoDonorPrior'] ==  '0.02':
            total += 1
            print item['varConsequences']




def checkCI():
    countCI = 0
    countNotCIButKey = 0
    for item in builtData:
        if item['varLoc'] == 'CI_domain_variant' or item['varLoc'] == 'CI_splice_acceptor_variant' or item['varLoc'] == 'CI_splice_donor_variant':
            if item['proteinPrior'] == '0.03': #random choice to look into specifics...
                countCI += 1
        elif item['proteinPrior'] == '0.03': #counts the number of times a 0.03 (presumably key domain) protein prior is attached to a non CI domain variant
            countNotCIButKey += 1
    print str(countCI) + ' CI variants with proteinPrior = 0.03 in built.tsv'
    print str(countNotCIButKey) + ' nonCI variants with proteinPrior = 0.03 in built.tsv'
    print str(countCI + countNotCIButKey) + ' total proteinPrior = 0.03'


def writeFrameshifts():
    fs1 = []
    fs0 = []
    fsNA = []
    fsNull = []
    count1 = 0
    count0 = 0
    countNA = 0
    countNull = 0
    countIndel = 0
    fsArray = []
    for item in builtData:
        if item['varType'] == 'insertion' or item['varType'] == 'deletion' or item['varType'] == 'delins':
            countIndel += 1
            # counts number of insertion, deletion, or delins variants
        for cDNA in pyHGVSTesters:
            newcDNA = item['HGVS_cDNA']
            if cDNA == newcDNA and newcDNA is not '-':
                testers.append(item)
                # print len(cDNA)
                # if cDNA == item['HGVS_cDNA']:
                #     print item['HGVS_cDNA']
                #     print cDNA #NOT WORKING
        if item['frameshiftFlag'] == '1':
            if item['varType'] == 'substitution':
                count1 += 1
                fs1.append(item)
                # appends all substitution variants with fs = 1 to eventually be written to a file (below ifs do same for 0, NA, -)
        if item['frameshiftFlag'] == '0':
            if item['varType'] == 'substitution':
                count0 += 1
                fs0.append(item)
        if item['frameshiftFlag'] == 'N/A':
            if item['varType'] == 'substitution':
                countNA += 1
                fsNA.append(item)
        if item['frameshiftFlag'] == '-':
            if item['varType'] == 'substitution':
                countNull += 1
                fsNull.append(item)
    print str(count0) + ' substitutions with framehsiftFlag = 0'
    print str(count1) + ' substitutions with framehsiftFlag = 1'
    print str(countNA) + ' substitutions with framehsiftFlag = N/A'
    print str(countNull) + ' substitutions with framehsiftFlag = -'
    print str(countIndel) + ' indel variants total'

    with open('fs1variants.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=builtData.fieldnames)
        writer.writeheader()
        for variant in fs1:
            writer.writerow(variant)

    with open('fs0variants.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=builtData.fieldnames)
        writer.writeheader()
        for variant in fs0:
            writer.writerow(variant)

    with open('fsNAvariants.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=builtData.fieldnames)
        writer.writeheader()
        for variant in fsNA:
            writer.writerow(variant)

    with open('fsNullvariants.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=builtData.fieldnames)
        writer.writeheader()
        for variant in fsNull:
            writer.writerow(variant)


def writeTesters():
    testData = csv.DictReader(open("TestFile1.csv", "r"))
    print testData.fieldnames
    pyHGVSTesters = []
    for item in testData:
        pyHGVSTesters.append(item['HGVS_cDNA'])
    # print len(pyHGVSTesters)
    count = 1
    testers = []
    print testers
    with open('testers.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=builtData.fieldnames)
        writer.writeheader()
        for variant in testers:
            writer.writerow(variant)

def checkNAvarCs():
    fsNAData = csv.DictReader(open('fsNAvariants.csv', "r"))
    countBlankConsequences = 0
    for item in fsNAData:
        if item['varConsequences'] == 'stop_gained':
            countBlankConsequences += 1
    print countBlankConsequences

def pickAnyVar(builtDataDict, cdot):
    "allows you to print prior data given any portion of a variant's cDNA nomenclature"
    builtData = builtDataDict
    for item in builtData:
        if cdot in item['pyhgvs_cDNA']:
            print item['pyhgvs_cDNA'], ', AP: ' + item['applicablePrior'], ', PP: ' + item['proteinPrior'], ',RDP: ', item['refDonorPrior'], ',DNP: ', item['deNovoDonorPrior'],  ',RAP: ', item['refAccPrior'], ',DNAP: ', item['deNovoAccPrior']

def checkVarInRange(builtDataDict):
    """Checks if a variant is in a hardocded range (varRange)"""
    builtData = builtDataDict
    varRange = range(9331, 9539)
    varStrRange = [] #becomes a 'numberline' of varRange
    for v in varRange:
        varStrRange.append(str(v)) #makes varRange strings because numbers in pyhgvs_cDNA are strings
    for item in builtData:
        for v in varStrRange:
            if v in item['pyhgvs_cDNA']:
                if item['proteinPrior'] == '0.99':
                    print item['pyhgvs_cDNA']

def checkStrangePriors(builtDataDict):
    """Checks for applicable priors that seem to not come from any preliminary prior. These are priors which have 
    applicable priors assigned as values, but N/As as all other preliminary prior values"""
    builtData = builtDataDict
    print builtData.fieldnames
    priorValueList =  ['0.02', '0.29', '0.66', '0.81', '0.99', '0.04', '0.34', '0.64''0.97']
    strangePriorCount = 0
    strangePriors = []
    fieldnames = ['pyhgvs_cDNA', 'applicablePrior', 'proteinPrior', 'refDonorPrior', 'deNovoDonorPrior', 'refAccPrior',
                  'deNovoAccPrior']
    strangePriors.append(fieldnames)
    for item in builtData:
        #print type(item)
        if  item['applicablePrior'] in priorValueList and item['proteinPrior'] == 'N/A' and item['deNovoDonorPrior'] == 'N/A' and item['refDonorPrior'] == 'N/A'  and item['refAccPrior'] == 'N/A':
            strangePriorCount += 1
            strangePriors.append(str(item['pyhgvs_cDNA']))



    f = open('strangePriors.csv', 'a')
    fnames = str(fieldnames)
    fnames = fnames[1:-1]
    #f.append(fnames)
    #f.write(str(fieldnames)
    with open('applicableTrueRestNA.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile)
        #writer.writeheader()
         #for variant in strangePriors:
        writer.writerows(strangePriors)
    # print strangePriorCount

def concordanceCaseCreate(sp, op):
    """Takes seans prior, our prior as argument to then append the name (pyhgvs_cDNA  from built) of the variant to a 
    list (caseList), returns list of variant names depending on arguments provided"""
        # 17 0.02       0.81
        # 6 0.96       0.03
        # 6 0.81       0.02
        # 3 0.97       0.34
        # 3 0.66       0.02
        # 2 0.96       0.29
        #
        #
        # 12 0.02       0.29
        # 11 0.02       0.66
        # 10 0.29       0.02
        # 5 0.3        0.03
        # 5 0.02       0.3
        #
        # 5 0.96       0.81
        # 4 0.5        0.99
        # 4 0.5        0.97
    countCase = 0
    concordanceData = csv.DictReader(open("priors_concordance.tsv", "r"), delimiter="\t")
    caseList = []
    for row in concordanceData:
        sPrior = row['applicable_prior']
        ourPrior =  row['applicablePrior']
        if sPrior == sp and ourPrior ==  op:
            caseList.append(row['pyhgvs_cDNA'])
    return caseList #list of pyhgvs_cDNAs

def getCaseDataFromBuilt(file, caseList):
    """obtains more data from built_with_priors.tsv given a case"""
    f = file
    builtData = csv.DictReader(open(f, "r"), delimiter="\t")
    #you have to read the file new again every time
    caseDataCollection = []
    #caseDict.keys = ['pyhgvs_cDNA', 'Genomic_Coordinate_hg38', 'applicablePrior', 'proteinPrior', 'refDonorPrior', 'deNovoDonorPrior', 'refAccPrior', 'deNovoAccPrior']
    #caseList = [] #NEEDS TO BE A CSV DICT...?
    #print builtData.fieldnames
    #print 'cDNA, hg38, applicablePrior, proteinPrior, refDonorPrior, deNovoDonorPrior, refAccPrior, deNovoAccPrior'
    for item in builtData:
        caseDict = {'pyhgvs_cDNA': '', 'Genomic_Coordinate_hg38': '', 'applicablePrior': '', 'proteinPrior': '',
                    'refDonorPrior': '', 'deNovoDonorPrior': '', 'refAccPrior': '', 'deNovoAccPrior': ''}
        if item['pyhgvs_cDNA'] in caseList:
            caseDict['pyhgvs_cDNA'] = item['pyhgvs_cDNA']
            caseDict['Genomic_Coordinate_hg38'] = item['Genomic_Coordinate_hg38']
            caseDict['applicablePrior'] = item['applicablePrior']
            caseDict['proteinPrior']  =  item['proteinPrior']
            caseDict['refDonorPrior'] = item['refDonorPrior']
            caseDict['deNovoDonorPrior'] = item['deNovoDonorPrior']
            caseDict['refAccPrior']  = item['refAccPrior']
            caseDict['deNovoAccPrior'] = item['deNovoAccPrior']
            caseDataCollection.append(caseDict)
    return caseDataCollection #list of dictionaries with priors headers as keys

def getCaseDataFromMupitBuilt(file, caseList):
    f = file
    caseList = caseList
    builtData = csv.DictReader(open(f, "r"), delimiter="\t")
    fn = builtData.fieldnames
    caseDataCollection =  []
    for var in builtData:
        caseDict = {"Combined_prior_probablility_exLOVD": "", "Genomic_Coordinate_hg38":"","pyhgvs_cDNA" :""}
        if var["pyhgvs_cDNA"] in caseList:
            caseDict["Combined_prior_probablility_exLOVD"] = var["Combined_prior_probablility_exLOVD"]
            caseDict["pyhgvs_cDNA"] = var["pyhgvs_cDNA"]
            caseDict["Genomic_Coordinate_hg38"] = var["Genomic_Coordinate_hg38"]
            caseDataCollection.append(caseDict)
    return caseDataCollection
    #for item in builtData

def checkProtPriorsFile(specVar):
    variantData = csv.DictReader(open("seanFile.txt", "r"), delimiter="\t")  # named variantData to resemble calcVarPriors.py
    for varEntry in variantData:
        caseDict = {'nthgvs': '', 'gene':'', 'protein_prior': '', 'de_novo_prior':  '', 'applicable_prior': ''}
        if specVar == varEntry['nthgvs']:
            caseDict['nthgvs'] = varEntry['nthgvs']
            caseDict['gene'] = varEntry['gene']
            caseDict['protein_prior'] = varEntry['protein_prior']
            caseDict['de_novo_prior'] = varEntry['de_novo_prior']
            caseDict['applicable_prior'] = varEntry['applicable_prior']
            break
    return caseDict


def getCaseDataFromSean(caseList):
    caseDataCollection = []
    for var in caseList:
        splitVar = var.split(':')
        caseDict = checkProtPriorsFile(splitVar[1])
        caseDataCollection.append(caseDict)
    #print 'hello'
    return caseDataCollection

