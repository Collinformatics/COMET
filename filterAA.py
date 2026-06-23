# PURPOSE: This code will load in your extracted substrates for processing

# IMPORTANT: Process all of your data with extractSubstrates before using this script


from functions import getFileNames, NGS
from itertools import product
import os
import sys



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inEnzymeName = 'Mpro2'
inPathFolder = os.path.join('Enzymes', inEnzymeName)
inSaveFigures = True
inSetFigureTimer = False

# Input 2: Computational Parameters
inFixResidues = True
inFixedResidue = 'Q'
inFixedPosition = 4
inExcludeResidues = False
inExcludedResidue = [['C','I','V'],['C','I','V'],['C','I','V'],['C','I','V']] # ['Y','Y','Y','Y','Y','Y','Y']
inExcludedPosition = [1,2,3,4,6,7,8,9]
inMinimumSubstrateCount = 1
inMinDeltaS = 0.6
inPrintFixedSubs = True
inShowSampleSize = True
inUseEnrichmentFactor = False
inCodonSequence = 'NNS' # Baseline probs of degenerate codons (can be N, S, or K)
inUseCodonProb = False # Use AA prob from inCodonSequence to calculate enrichment
inAvgInitialProb = False
inDropResidue = ['R9'] # To drop: inDropResidue = ['R9'], For nothing: inDropResidue = []

# Input 3: Making Figures
inBlockFigures = False
inPlotEntropy = True
inPlotEnrichmentMap = True
inPlotEnrichmentMapScaled = False
inPlotLogo = True
inPlotWeblogo = True
inPlotMotifEnrichment = True
inPlotWordCloud = True
inPlotAADistribution = False
inPlotBarGraphs = False
inPlotPCA = False
inPlotCounts = False
inPlotPositionalProbDist = False # For understanding shannon entropy
if inBlockFigures:
    inPlotEntropy = False
    inPlotEnrichmentMap = False
    inPlotEnrichmentMapScaled = False
    inPlotLogo = False
    inPlotWeblogo = False
    inPlotMotifEnrichment = False
    inPlotWordCloud = False # Word cloud
    inPlotMotifEnrichment = False
    inPlotAADistribution = False
    inPlotBarGraphs = False
    inPlotPCA = False
    inPlotCounts = False

# Input 4: Inspecting The data
inPrintNumber = 10
inFindSequences = False
inFindSeq = ['AVLQS', 'VILQS','VILQT','VILQS','VILHS','VIMQS','VPLQS','NILQS']
inFindSeq = ['CC', 'CI', 'CV', 'VV', 'II', 'IC', 'IV', 'VV', 'VC', 'VI']

# Input 5: CSV
inSaveCSV = False # Save substrates in a csv file
inMinSubsCSV = 10000 # Minimum counts for saved substrates
inSubLengthCSV = False # If: False, use full seq, If: 6, use 6 AA seq
inUseBgSubs = True
inExcludeSeq = 'LQ'
inMaxBgSubstrateCount = 20
inModulo = 2500 # Increase to select fewer Bg substrates
inScaleModulo = False # Continually increase modulus value to keep fewer low count bg subs

# Input 6: Plot Heatmap
inShowEnrichmentScores = True
inShowEnrichmentAsSquares = False
inYLabelEnrichmentMap = 2 # 0 for full Residue name, 1 for 3-letter code, 2 for 1 letter

# Input 7: Plot Sequence Motif
inBigLettersOnTop = False
inLimitYAxis = False

# Input 8: Word Cloud
inLimitWords = True
inTotalWords = 50

# Input 9: Bar Graphs
inNSequences = 50

# Input 10: PCA
inPCAMotif = False
inNumberOfPCs = 2
inTotalSubsPCA = 10000
inIncludeSubCountsESM = False
inExtractPopulations = False
inPlotEntropyPCAPopulations = False
inAdjustZeroCounts = False # Prevent counts of 0 in PCA EM & Motif

# Input 11: Predict Activity
inPredictActivity = False
inPredictionTag = 'pp1a/b Substrates'
inPredictSubstrates = ['AVLQSGFR', 'VTFQSAVK', 'ATVQSKMS', 'ATLQAIAS',
                       'VKLQNNEL', 'VRLQAGNA', 'PMLQSADA', 'TVLQAVGA',
                       'ATLQAENV', 'TRLQSLEN', 'PKLQSSQA']
# inPredictionTag = 'Substrates'
# inPredictSubstrates = ['AVLQSGFR', 'VILQAGFR', 'VILQAPFR', 'LVLQSNDL',
#                        'ATLQGLMI', 'TVLQAAML', 'VSLQSTYK', 'VSLQGAEL']
# inPredictionTag = 'FP14-18'
# inPredictSubstrates = ['AVLQSGFR', 'TVLQAAMH', 'VLLQGCVH',
#                        'WVLQAKLL', 'AILQCMLG', 'VLLQGVVH']
# inPredictionTag = 'FP19-23'
# inPredictSubstrates = ['AVLQSGFR', 'CILQAVFH', 'VVLQAVMH',
#                        'SILQCVLM', 'VMLQAVFH', 'PLLQAILM']
inPredictionTag = 'Heatmap Substrates'
inPredictSubstrates = ['AVLQSGFR', 'VILQSGFR', 'VILQSPFR', 'VILHSGFR', 'VIMQSGFR',
                       'VPLQSGFR', 'NILQSGFR', 'VILQTGFR', 'PILQSGFR', 'PIMQSGFR']
inRankScores = False
inScalePredMatrix = False # Scale EM by ΔS

# Input 12: Optimal Substrates
inEvaluateOS = False
inPrintOSNumber = 10
inMaxResidueCount = 4

# Input 13: Evaluate Substrate Enrichment
inEvaluateSubstrateEnrichment = False # ============= Fix: Load Initial Subs =============
inSaveEnrichedSubstrates = False
inNumberOfSavedSubstrates = 10**6

# Input 14: Evaluate Positional Preferences
inPlotPosProb = False # Plot RF distributions of a given AA
inCompairAA = 'L' # Select AA of interest (different A than inFixedResidue)



# =================================== Setup Parameters ===================================
# Colors:
white = '\033[38;2;255;255;255m'
greyDark = '\033[38;2;144;144;144m'
purple = '\033[38;2;189;22;255m'
magenta = '\033[38;2;255;0;128m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'

# Load: Dataset labels
enzymeName, filesInitial, filesFinal, labelAAPos = getFileNames(enzyme=inEnzymeName)

inPlotMotifEnrichmentNBars = True



# =================================== Initialize Class ===================================
ngs = NGS(
    enzyme=inEnzymeName, enzymeName=enzymeName, substrateLength=len(labelAAPos),
    filterSubs=inFixResidues, fixedAA=inFixedResidue, fixedPosition=inFixedPosition,
    excludeAAs=inExcludeResidues, excludeAA=inExcludedResidue,
    excludePosition=inExcludedPosition, minCounts=inMinimumSubstrateCount,
    minEntropy=inMinDeltaS, figEMSquares=inShowEnrichmentAsSquares,
    xAxisLabels=labelAAPos, printNumber=inPrintNumber, showNValues=inShowSampleSize,
    bigAAonTop=inBigLettersOnTop, findMotif=False, folderPath=inPathFolder,
    filesInit=filesInitial, filesFinal=filesFinal, useEF=inUseEnrichmentFactor,
    plotPosS=inPlotEntropy, plotFigEM=inPlotEnrichmentMap,
    plotFigEMScaled=inPlotEnrichmentMapScaled, plotFigLogo=inPlotLogo,
    plotFigWebLogo=inPlotWeblogo, plotFigMotifEnrich=inPlotMotifEnrichment,
    plotFigWords=inPlotWordCloud, wordLimit=inLimitWords, wordsTotal=inTotalWords,
    plotFigBars=inPlotBarGraphs, NSubBars=inNSequences, plotFigPCA=inPlotPCA,
    numPCs=inNumberOfPCs, NSubsPCA=inTotalSubsPCA, plotSuffixTree=False,
    saveFigures=inSaveFigures, setFigureTimer=inSetFigureTimer
)



# ====================================== Load data =======================================
# Get dataset tag
ngs.getDatasetTag(useCodonProb=inUseCodonProb, codon=inCodonSequence)

# Load: Counts
countsInitial, countsInitialTotal = ngs.loadCounts(filter=False, fileType='Initial Sort',
                                                   dropColumn=inDropResidue)

# Load: Substrates
if inFindSequences or inEvaluateSubstrateEnrichment or inUseBgSubs:
    substratesInitial, totalSubsInitial = ngs.loadUnfilteredSubs(loadInitial=True)

# Calculate: Initial sort probabilities
rfInitial = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
                              fileType='Initial Sort', calcAvg=inAvgInitialProb)

substratesInitial = None
loadUnfilteredSubs = False
filePathFixedCountsFinal, filePathFixedSubsFinal = None, None
substratesFinal, countsFinal, countsFinalTotal = None, None, None
if inFixResidues or inExcludeResidues:
    filePathFixedSubsFinal, filePathFixedCountsFinal = (
        ngs.getFilePath(datasetTag=ngs.datasetTagSave, sortType='FinalSort'))

    # Verify that the file exists
    if (os.path.exists(filePathFixedSubsFinal) and
            os.path.exists(filePathFixedCountsFinal)):

        # Load: Counts
        countsFinal, countsFinalTotal = ngs.loadCounts(filter=True,
                                                       fileType='Final Sort',
                                                       datasetTag=ngs.datasetTagSave,
                                                       dropColumn=inDropResidue)

        # Load: Substrates
        substratesFinal, totalSubsFinal = ngs.loadSubstratesFiltered()

        if countsFinalTotal != totalSubsFinal:
            print(f'{orange}ERROR: '
                  f'The total number of Loaded Counts ({cyan}{countsFinalTotal:,}'
                  f'{orange}) =/= number of Total Substrates '
                  f'({cyan}{totalSubsFinal:,}{orange})\n')
            sys.exit()
    else:
        loadUnfilteredSubs = True

        # Load: Substrates
        substratesFinal, totalSubsFinal = ngs.loadUnfilteredSubs(loadFinal=True)
else:
    substratesFinal, totalSubsFinal = ngs.loadUnfilteredSubs(loadFinal=True)

    # Load: Counts
    countsFinal, countsFinalTotal = ngs.loadCounts(fileType='Final Sort', filter=False,
                                                   dropColumn=inDropResidue)



# ================================== Evaluate The Data ===================================
if inFixResidues:
    if loadUnfilteredSubs:
        # Fix AA
        substratesFinal, countsFinalTotal = ngs.fixResidue(
            substrates=substratesFinal, fixedString=ngs.datasetTag,
            printRankedSubs=inPrintFixedSubs, sortType='Final Sort')

        # Count fixed substrates
        countsFinal, countsFinalTotal = ngs.countResidues(substrates=substratesFinal,
                                                          datasetType='Final Sort')

        # Save the data
        ngs.saveData(substrates=substratesFinal, counts=countsFinal)

        # Filter counts matrix
        if inDropResidue:
            countsFinal = ngs.dropColumnsFromMatrix(countMatrix=countsFinal,
                                                    datasetType='Final Sort',
                                                    dropColumn=inDropResidue)

    if inEvaluateSubstrateEnrichment:
        fixedSubsInitial, countsInitialFixedTotal = ngs.fixResidue(
            substrates=substratesInitial, fixedString=ngs.datasetTag,
            printRankedSubs=inPrintFixedSubs, sortType='Initial Sort')

if inExcludeResidues and not inFixResidues:
    if loadUnfilteredSubs:
        # Exclude AAs
        # Fix AA
        substratesFinal, countsFinalTotal = ngs.exclResidue(
            substrates=substratesFinal, fixedString=ngs.datasetTag,
            printRankedSubs=inPrintFixedSubs, sortType='Final Sort')

        # Count fixed substrates
        countsFinal, countsFinalTotal = ngs.countResidues(substrates=substratesFinal,
                                                          datasetType='Final Sort')

        # Save the data
        ngs.saveData(substrates=substratesFinal, counts=countsFinal)


# Display current sample size
ngs.recordSampleSize(NInitial=countsInitialTotal, NFinal=countsFinalTotal,
                     NFinalUnique=len(substratesFinal.keys()))

# Calculate: RF
rfFinal = ngs.calculateRF(counts=countsFinal, N=countsFinalTotal,
                            fileType='Final Sort')

if inPlotAADistribution:
    # Plot: AA probabilities in initial and final sorts
    ngs.plotLibraryAADist(rfInitial=rfInitial, rfFinal=rfFinal,
                            codonType=inCodonSequence, datasetTag=ngs.datasetTag)

    # Evaluate: Degenerate codon probabilities
    probCodon = ngs.calculateProbCodon(codonSeq=inCodonSequence)

    # Plot: Codon probabilities
    ngs.plotLibraryAADist(rfInitial=rfFinal, rfFinal=probCodon,
                            codonType=inCodonSequence, datasetTag=inCodonSequence,
                            skipInitial=True)

# Calculate: Positional entropy
entropy = ngs.calculateEntropy(rf=rfFinal)

# Calculate: Enrichment scores
if inUseCodonProb:
    # Evaluate: Degenerate codon probabilities
    probCodon = ngs.calculateProbCodon(codonSeq=inCodonSequence)
    enrichmentScores = ngs.calculateEnrichment(rfInitial=probCodon,
                                               rfFinal=rfFinal)
else:
    enrichmentScores = ngs.calculateEnrichment(rfInitial=rfInitial,
                                               rfFinal=rfFinal)

# Create csv
if inSaveCSV:
    if inUseBgSubs:
        ngs.saveSubstrateCSV(
            seqs=substratesFinal, initialRF=rfInitial, finalRF=rfFinal,
            minCounts=inMinSubsCSV, seqsBg=substratesInitial, excludeAA=inExcludeSeq,
            maxCountsBg=inMaxBgSubstrateCount, mod=inModulo, modScale=inScaleModulo,
            chopSeq=inSubLengthCSV
        )
    else:
        ngs.saveSubstrateCSV(
            seqs=substratesFinal, initialRF=rfInitial,
            finalRF=rfFinal, minCounts=inMinSubsCSV, chopSeq=inSubLengthCSV
        )

# Evaluate: Sequences
if inUseEnrichmentFactor:
    motifs = ngs.processSubstrates(
        subsInit=substratesInitial, subsFinal=substratesFinal, motifs=substratesFinal,
        subLabel=labelAAPos)

if inPlotPCA:
    finalSubsMotif = ngs.getMotif(substrates=substratesFinal)
    if inUseEnrichmentFactor:
        if substratesInitial is None:
            print(f'{red}Write code to load in the initial substrates{resetColor}\n')
            sys.exit()

    # # Plot: Work cloud
    # if inPlotWordCloud:
    #     if inFixResidues:
    #         ngs.plotWordCloud(substrates=finalSubsMotif, plotEF=inUseEnrichmentFactor)
    #     else:
    #         ngs.plotWordCloud(substrates=substratesFinal, plotEF=inUseEnrichmentFactor)

    # Plot: PCA
    if inPlotPCA:
        if inPCAMotif:
            datasetTag = ngs.datasetTagMotif
            saveTag = ngs.saveTagMotif
            labelPos = ngs.xAxisLabelsMotif
            subsPCA = finalSubsMotif
        else:
            datasetTag =  ngs.datasetTag
            saveTag = datasetTag
            labelPos = labelAAPos
            subsPCA = substratesFinal

        # Convert substrate data to numerical
        tokensESM, subsESM, subCountsESM = ngs.ESM(
            substrates=subsPCA, collectionNumber=int(inTotalSubsPCA),
            useSubCounts=inIncludeSubCountsESM, subPositions=labelPos,
            datasetTag=saveTag)

        # Cluster substrates
        subPopulations = ngs.plotPCA(substrates=finalSubsMotif, data=tokensESM,
                                     indices=subsESM, numberOfPCs=inNumberOfPCs,
                                     N=subCountsESM, fixedSubs=inFixResidues,
                                     datasetTag=datasetTag, saveTag=saveTag)

        # Plot: Substrate clusters
        if subPopulations is not None:
            clusterCount = len(subPopulations)
            for index, subCluster in enumerate(subPopulations):
                # Plot data
                ngs.plotSubstratePopulations(
                    substrates=subCluster, clusterIndex=index, numClusters=clusterCount,
                    datasetTag=datasetTag, saveTag=saveTag)

            print(f'Debug PCA')
            sys.exit()

# Predict substrate activity
if inPredictActivity:
    ngs.predictActivityHeatmap(
        predSubstrates=inPredictSubstrates, predModel=ngs.datasetTag,
        predLabel=inPredictionTag, rankScores=inRankScores, scaleEMap=inScalePredMatrix)

# Plot: Word cloud
if inPlotWordCloud and not inUseEnrichmentFactor:
    ngs.plotWordCloud(substrates=substratesFinal)

# Plot: Bar graphs
if inPlotBarGraphs and not inUseEnrichmentFactor:
    ngs.plotBarGraph(substrates=substratesFinal, dataType='Counts')
    ngs.plotBarGraph(substrates=substratesFinal, dataType='Relative Frequency')

# Find sequences
if inFindSequences:
    # ngs.findSequence(substrates=substratesInitial,
    #                  sequence=inFindSeq,
    #                  sortType='Initial Sort')
    #sys.exit()
    ngs.findSequence(substrates=substratesFinal,
                     sequence=inFindSeq,
                     sortType='Final Sort')


# ========================================================================================
if inEvaluateOS:
    print('============================== Evaluate Optimal Substrates '
          '==============================')
    if inFixResidues:
        # Determine the OS
        combinations = 1
        optimalAA = []
        substratesOS = {}
        for indexColumn, column in enumerate(ngs.eMap.columns):
            # Find the best residues at this position
            optimalAAPos = ngs.eMap[column].nlargest(inMaxResidueCount)

            # Filter the data
            for rank, (AA, ES) in (
                    enumerate(zip(optimalAAPos.index, optimalAAPos.values), start=1)):
                if ES <= 0:
                    optimalAAPos = optimalAAPos.drop(index=AA)
            optimalAA.append(optimalAAPos)
        print(f'Optimal Residues: {purple}{inEnzymeName} - {ngs.datasetTag}{resetColor}')
        for index, data in enumerate(optimalAA, start=1):
            # Determine the number of variable residues at this position
            numberAA = len(data)
            combinations *= numberAA

            # Define substrate position
            positionSub = labelAAPos[index-1]
            print(f'Preferred Residues:\n'
                  f'{greenLight}{data}{resetColor}\n')
        print(f'Possible Substrate Combinations: {red}{combinations:,}'
              f'{resetColor}\n\n')

        # Make all possible substrate sequences
        residueChoices = [series.index.tolist() for series in optimalAA]
        subCombos = list(product(*residueChoices))
        substrates = [''.join(combo) for combo in subCombos]

        # Predict activity
        substratesOS = {}
        for substrate in substrates:
            score = 0
            for index in range(len(substrate)):
                AA = substrate[index]
                pos = ngs.eMap.columns[index]
                score += ngs.eMap.loc[AA, pos]
            substratesOS[substrate] = score
        substratesOS = dict(sorted(substratesOS.items(),
                                   key=lambda x: x[1], reverse=True))

        print(f'Predicted Optimal Substrates:')
        for index, (substrate, ES) in enumerate(substratesOS.items()):
            print(f'     {pink} {substrate}{resetColor}, '
                  f'ES:{red} {ES:.3f}{resetColor}')
            if index >= inPrintOSNumber:
                sys.exit()

        print(f'\nNumber of substrates:{greyDark} {len(substratesOS):,}{resetColor}\n\n')



if inEvaluateSubstrateEnrichment:
    ngs.substrateEnrichment(initialSubs=substratesInitial, finalSubs=substratesFinal,
                            NSubs=inNumberOfSavedSubstrates,
                            saveData=inSaveEnrichedSubstrates)

    if inFixResidues:
        ngs.substrateEnrichment(initialSubs=substratesInitial, finalSubs=substratesFinal,
                                NSubs=inNumberOfSavedSubstrates,
                                saveData=inSaveEnrichedSubstrates)

if inPlotCounts:
    # Plot the data
    ngs.plotCounts(countedData=countsFinal, totalCounts=countsFinalTotal,
                   fileName=ngs.datasetTag)

if inPlotPositionalProbDist:
    ngs.plotPositionalProbDist(probability=rfFinal, entropyScores=entropy,
                               sortType='Final Sort', datasetTag=ngs.datasetTag)
