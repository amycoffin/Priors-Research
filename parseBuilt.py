import builtParsing as bp
import csv

"""This module prints cases of priors discordances and outputs (discordanceCases.py) a caseDirectory as well as multiple 
caseData objects (arrays of dictionaries) that correspond to multiple discordance data criteria. caseDict provides
 the case number as keys and the a tuple of priors criteria arguments (seans, brcaexchange) as values. caseData objects 
 come in the form of 'seanPriorsData' or 'builtPriorsData', each numbered according to the case that they correspond 
 with.
 
 caseLists originate from priors_concordance.tsv (built_with_mupit.tsv was used to make priors_concordance.tsv)
 caseBuiltCollect (builtPriorsData) objects orginate from built_with_priors.tsv
 caseSeanCollect (seanPriorsData) objects originate from seanFile.txt ('mod_res_dn_brca20160525.txt' in calcVarPriors)"""

caseDirectory = {} #will tell you which cases correspond with which applicable priors arguments

"""Likely to be resolved cases"""

"""Case 1 is not in a key domain and should be 0.02"""
caseList1 = bp.concordanceCaseCreate('0.02', '0.81') #returns a list that meets the seans prior, our prior args criteria
caseDirectory['case1'] = ('0.02', '0.81')
caseBuiltCollect1 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList1)
caseSeanCollect1 = bp.getCaseDataFromSean(caseList1)

"""Case 2 is a start codon mutation, and should be 0.96"""
caseList2 = bp.concordanceCaseCreate('0.96', '0.03') #originates from priors_concordance.tsv
#args are sean's applicable prior, our applicable prior
caseDirectory['case2'] = ('0.96', '0.03')
caseBuiltCollect2 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList2) #origniates from built.tsv
caseSeanCollect2 = bp.getCaseDataFromSean(caseList2) #originates from seanFile.txt (weird file name in brcaex repo)

"""Variants in caseList3 are within exon 25 of BRCA 2, within the DNA binding domain. 
Protein prior should be 0.81"""
caseList3 = bp.concordanceCaseCreate('0.81', '0.02')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case3'] = ('0.81', '0.02')
caseBuiltCollect3 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList3)
caseSeanCollect3 = bp.getCaseDataFromSean(caseList3)

"""There are 3 case 4 variants. 2/3 of them do not have data because they are intronic; the 3rd shows to have \
protein_prior = 0.02 in source file used by calcVarPriors.py, and an applicable prior of 0.97 with no explanation 
 (seanFile.txt AKA mod_res_dn_brca20160525.txt)"""
caseList4 = bp.concordanceCaseCreate('0.97','0.34')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case4'] = ('0.97','0.34')
caseBuiltCollect4 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList4)
caseSeanCollect4 = bp.getCaseDataFromSean(caseList4)

""""""
caseList5 = bp.concordanceCaseCreate('0.66', '0.02')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case5'] = ('0.66', '0.02')
caseBuiltCollect5 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList5)
caseSeanCollect5 = bp.getCaseDataFromSean(caseList5)

caseList6 = bp.concordanceCaseCreate('0.96', '0.29')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case6'] = ('0.96', '0.29')
caseBuiltCollect6 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList6)
caseSeanCollect6 = bp.getCaseDataFromSean(caseList6)

"""Concerning Cases"""

caseList7 = bp.concordanceCaseCreate('0.02', '0.29')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case7'] = ('0.02', '0.29')
caseBuiltCollect7 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList7)
caseSeanCollect7 = bp.getCaseDataFromSean(caseList7)

caseList8 = bp.concordanceCaseCreate('0.02', '0.66')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case8'] = ('0.02', '0.66')
caseBuiltCollect8 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList8)
caseSeanCollect8 = bp.getCaseDataFromSean(caseList8)

caseList9 = bp.concordanceCaseCreate('0.29', '0.02')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case9'] = ('0.29', '0.02')
caseBuiltCollect9 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList9)
caseSeanCollect9 = bp.getCaseDataFromSean(caseList9)

caseList10 = bp.concordanceCaseCreate('0.3', '0.03')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case10'] = ('0.3', '0.03')
caseBuiltCollect10 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList10)
caseSeanCollect10 = bp.getCaseDataFromSean(caseList10)

caseList11 = bp.concordanceCaseCreate('0.02', '0.3')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case11'] = ('0.02', '0.3')
caseBuiltCollect11 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList11)
caseSeanCollect11 = bp.getCaseDataFromSean(caseList11)

"""Likely-to-be-resolved-cases"""
caseList12 = bp.concordanceCaseCreate('0.96', '0.81')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case12'] = ('0.96', '0.81')
caseBuiltCollect12 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList12)
caseSeanCollect12 = bp.getCaseDataFromSean(caseList12)

caseList13 = bp.concordanceCaseCreate('0.5', '0.99')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case13'] = ('0.5', '0.99')
caseBuiltCollect13 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList13)
caseSeanCollect13 = bp.getCaseDataFromSean(caseList13)

caseList14 = bp.concordanceCaseCreate('0.5', '0.97')
#args are sean's applicable prior first, our applicable prior second
caseDirectory['case14'] = ('0.5', '0.97')
caseBuiltCollect14 = bp.getCaseDataFromBuilt("built_with_priors.tsv", caseList14)
caseSeanCollect14 = bp.getCaseDataFromSean(caseList14)

f = open("discordanceCases.py", 'w')
f.write('caseDirectory = ' + str(caseDirectory) + '\n' +
        'caseList1 = ' + str(caseList1) + '\n' + 'seanPriorsData1 = ' + str(caseSeanCollect1) + '\n' + 'builtPriorsData1 = ' + str(caseBuiltCollect1) + '\n' +
        'caseList2 = ' + str(caseList2) + '\n' + 'seanPriorsData2 = ' + str(caseSeanCollect2) + '\n' + 'builtPriorsData2 = ' + str(caseBuiltCollect2) + '\n' +
        'caseList3 = ' + str(caseList3) + '\n' + 'seanPriorsData3 = ' + str(caseSeanCollect3) + '\n' + 'builtPriorsData3 = ' + str(caseBuiltCollect3) + '\n' +
        'caseList4 = ' + str(caseList4) + '\n' + 'seanPriorsData4 = ' + str(caseSeanCollect4) + '\n' + 'builtPriorsData4 = ' + str(caseBuiltCollect4) + '\n' +
        'caseList5 = ' + str(caseList5) + '\n' + 'seanPriorsData5 = ' + str(caseSeanCollect5) + '\n' + 'builtPriorsData5 = ' + str(caseBuiltCollect5) + '\n' +
        'caseList6 = ' + str(caseList6) + '\n' + 'seanPriorsData6 = ' + str(caseSeanCollect6) + '\n' + 'builtPriorsData6 = ' + str(caseBuiltCollect6) + '\n' +
        'caseList7 = ' + str(caseList7) + '\n' + 'seanPriorsData7 = ' + str(caseSeanCollect7) + '\n' + 'builtPriorsData7 = ' + str(caseBuiltCollect7) + '\n' +
        'caseList8 = ' + str(caseList8) + '\n' + 'seanPriorsData8 = ' + str(caseSeanCollect8) + '\n' + 'builtPriorsData8 = ' + str(caseBuiltCollect8) + '\n'+
        'caseList9 = ' + str(caseList9) + '\n' + 'seanPriorsData9 = ' + str(caseSeanCollect9) + '\n' + 'builtPriorsData9 = ' + str(caseBuiltCollect9) + '\n'+
        'caseList10 = ' + str(caseList10) + '\n' + 'seanPriorsData10 = ' + str(caseSeanCollect10) + '\n' + 'builtPriorsData10 = ' + str(caseBuiltCollect10) + '\n'+
        'caseList11 = ' + str(caseList11) + '\n' + 'seanPriorsData11 = ' + str(caseSeanCollect11) + '\n' + 'builtPriorsData11 = ' + str(caseBuiltCollect11) + '\n'+
        'caseList12 = ' + str(caseList12) + '\n' + 'seanPriorsData12 = ' + str(caseSeanCollect12) + '\n' + 'builtPriorsData12 = ' + str(caseBuiltCollect12) + '\n'+
        'caseList13 = ' + str(caseList13) + '\n' + 'seanPriorsData13 = ' + str(caseSeanCollect13) + '\n' + 'builtPriorsData13 = ' + str(caseBuiltCollect13) + '\n'+
        'caseList14 = ' + str(caseList14) + '\n' + 'seanPriorsData14 = ' + str(caseSeanCollect14) + '\n' + 'builtPriorsData14 = ' + str(caseBuiltCollect14) + '\n')




# dnb = range(32356433, 32396954)
# startError  = 32396954 - 32394763
# endError = 32396954 - 32396935
