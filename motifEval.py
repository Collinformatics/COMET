from functions import getFileNames, NGS
import os
import pandas as pd
import sys



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inEnzymeName = 'Den'
inPathFolder = os.path.join('Enzymes', inEnzymeName)
inSaveFigures = True
inSetFigureTimer = False

# Input 2: Experimental Parameters
inMotifPositions = ['P3','P2','P1','P1\'','P2\'','P3\'']
# inMotifPositions = ['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
# inMotifPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
inIndexNTerminus = 1 # Define the index if the first AA in the motif

# Input 3: Computational Parameters
inFixedResidue = ['R', ['A', 'G', 'S']]
inFixedPosition = [[4,5],[5,6]]
inExcludeResidues = True
inExcludedResidue = ['A']
inExcludedPosition = [9,10]
inMinimumSubstrateCount = 1
inCodonSequence = 'NNS' # Baseline probs of degenerate codons (can be N, S, or K)
inUseCodonProb = False # Use AA prob from inCodonSequence to calculate enrichment
inAvgInitialProb = True
inDropResidue = ['R9'] # To drop: inDropResidue = ['R9'], For nothing: inDropResidue = []

# Input 4: Figures
# inPlotPCA = False # PCA plot of an individual fixed frame
inBlockFigures = False
inPlotEntropy = True
inPlotEnrichmentMap = True
inPlotEnrichmentMapScaled = False
inPlotLogo = True
inPlotWeblogo = True
inPlotMotifEnrichment = True
inPlotWordCloud = True
inPlotStats = True
inPlotBarGraphs = True
inPlotPCA = False # PCA plot of the combined set of motifs
inPlotSuffixTree = False
inPredictSubstrateActivity = False
inPredictSubstrateActivityPCA = False
inPlotBinnedSubstrateES = False
inPlotBinnedSubstratePrediction = False
inPlotCounts = False
inPlotFilteredSubs = True
inShowSampleSize = True # Include the sample size in your figures
if inBlockFigures:
    inPlotEntropy = False
    inPlotEnrichmentMap = False
    inPlotEnrichmentMapScaled = False
    inPlotLogo = False
    inPlotWeblogo = False
    inPlotMotifEnrichment = False
    inPlotWordCloud = False # Word cloud
    inPlotStats = False
    inPlotBarGraphs = False
    inPlotPCA = False
    inPlotSuffixTree = False
    inPredictSubstrateActivity = False
    inPredictSubstrateActivityPCA = False
    inPlotBinnedSubstrateES = False
    inPlotBinnedSubstratePrediction = False
    inPlotCounts = False
    inPlotFilteredSubs = False

# Input 5: CSV
inSaveCSV = False # Save substrates in a csv file
inMinSubsCSV = 1000 # Minimum counts for saved substrates
inSubLengthCSV = 6 # If: False, use full seq. If: 6, use 6 AA seq
inUseBgSubs = True
inExcludeSeq = 'LQ'
inMaxBgSubstrateCount = 30
inModulo = 40000 # Increase to select fewer Bg substrates
inScaleModulo = True # Continually increase modulus value to keep fewer low count bg subs

# Input 6: Printing The Data
inPrintLoadedSubs = True
inPrintSampleSize = True
inPrintCounts = True
inPrintRF = True
inPrintES = True
inPrintEntropy = True
inPrintMotifData = True
inPrintNumber = 10

# Input 7: Find Protein Sequences
inFindSequences = False
inFindSeq = ['AVLQSG', 'VILQSG','VILQTG','VILQSP','VILHSG','VIMQSG','VPLQSG','NILQSG']
inFindAAInSequence = False
inFindAA = ['A', 'F', 'W']
inAAPos = 4

# Input 8: Plot Heatmap
inShowEnrichmentScores = True
inShowEnrichmentAsSquares = False

# Input 9: Plot Sequence Motif
inNormLetters = False  # Normalize fixed letter heights
inPlotWeblogoMotif = False
inShowWeblogoYTicks = True
inAddHorizontalLines = False
inPlotNegativeWeblogoMotif = False
inBigLettersOnTop = False

# Input 10: Motif Enrichment
inPlotNBars = 50

# Input 11: Word Cloud
inLimitWords = True
inTotalWords = inPlotNBars

# Input 12: PCA
inNumberOfPCs = 2
inTotalSubsPCA = int(5*10**4)
inIncludeSubCountsESM = True
inPlotEntropyPCAPopulations = False

# Input 13: Predict Activity
inPredictActivity = True
inPredictSubstrates = []
inUseNaturalSubs = False
if inUseNaturalSubs:
    inPredictionLabel = 'pp1a/b Substrates'
    inPredictSubstrates = ['AVLQSGFR', 'VTFQSAVK', 'ATVQSKMS', 'ATLQAIAS',
                           'VKLQNNEL', 'VRLQAGNA', 'PMLQSADA', 'TVLQAVGA',
                           'ATLQAENV', 'TRLQSLEN', 'PKLQSSQA']
    inSubstrateActivity = {}
    for substrate in inPredictSubstrates:
        inSubstrateActivity[substrate] = 50.0
    inErrorBars = []
else:
    inPredictionLabel = '30 Min'
    inSubstrateActivity = {
        'AVLQSG': 55,  # 60,
        'VILQSG': 66,  # 70,
        'VILQTG': 34,  # 6,
        'VILQSP': 0,  # 0,
        'VILHSG': 26,  # 15,
        'VIMQSG': 61,  # 50,
        'VPLQSG': 0,  # 0,
        'NILQSG': 22,  # 6,
    }
    inSubstrateActivity = {
        'PLALWR': 0,
        'PMALVV': 0,
        'PLALVV': 0,
        'PMELVV': 0,
        'PMVLVV': 0,
        'TMALVV': 0,
        'FMALVV': 0,
        'PMTLVV': 0,
        'PMALPV': 0,
        'PMALVP': 0,
        'PAALVV': 0,
        'PTALVV': 0,
        'PMAMVV': 0,
        'PMATVV': 0,
        'VMALVV': 0,
        'PMAIVV': 0,
        'PMLLVV': 0,
        'PMMLVV': 0,
        'MMALVV': 1
    }
    # inSubstrateActivity = {
    #     'PMCMELVV': 4.06 * 10 ** -8,
    #     'PMCMALVV': 3.91 * 10 ** -8,
    #     'PMVMELVV': 3.95 * 10 ** -8,
    #     'PMVMALVV': 3.80 * 10 ** -8,
    #     'PVLALMLM': 5.39 * 10 ** -9,
    #     'ASQGLLDR': 2.83 * 10 ** -11,
    #     'RDDTTWPP': 1.21 * 10 ** -16
    # }
    inPredictionLabel = ''
    inSubstrateActivity = {
        'CMELVV': 1,
        'CMALVV': 0,
        'VMELVV': 0,
        'VMALVV': 0,
        'VLALML': 0,
        'QGLLDR': 0,
        'DTTWPP': 0
    }
    inErrorBars = [0.01, 0.058, 0.025, 0.0, 0.027, 0.044, 0.0, 0.033], # Avg stdev
inEMapStartIndex = 0  # Sub: ACDEFGHI, if idx = 0 start at A
inRankScores = False
inScalePredMatrix = False  # Scale EM by ΔS

# Input 14: Codon Enrichment
inPredictCodonsEnrichment = False

# Input 15: Evaluate Known Substrates
inNormalizePredictions = True
inYMaxPred = 1.05
inYMinPred, inYMinPredScaled, inYMinPredAI = 0, 0, -0.25
inYTickMinPred, inYTickMinScaled, inYTickMinAI = inYMinPred, inYMinPredScaled, -0.4
inSubsPredict = ['VVLQSGFR', 'VVLQSPFR', 'VYLQSGFR', 'VVLQAGFR', 'VVMQSGFR',
                 'IVLQSGFR', 'VVLHSGFR', 'VGLQSGFR', 'VVLMSGFR', 'VVVQSGFR',
                 'VVLQIGFR', 'VVGQSGFR', 'KVLQSGFR', 'VVLQNGFR', 'VVLYSGFR']
inSubsPredictStartIndex = 0
inKnownTarget = ['nsp4/5', 'nsp5/6', 'nsp6/7', 'nsp7/8', 'nsp8/9', 'nsp9/10',
                 'nsp10/12', 'nsp12/13', 'nsp13/14', 'nsp14/15', 'nsp15/16']
inBarWidth = 0.75
inBarColor = '#BF5700' # Burnt orange
inEdgeColor = 'black'
inEdgeColorOrange = '#F8971F' # Orange
inDatapointColor = []
for _ in inSubsPredict:
    inDatapointColor.append(inBarColor)

# Input 16: Evaluate Binned Substrates
inPlotEnrichedSubstrateFrame = False
inPrintLoadedFrames = True
inPlotBinnedSubNumber = 30
inPlotBinnedSubProb = True
inPlotBinnedSubYMax = 0.07

# Input 17: Predict Binned Substrate Enrichment
inEvaluatePredictions = False
inPrintPredictions = False
inBottomParam = 0.16
inPredictionDatapointColor = '#BF5700'
inMiniumSubstrateScoreLimit = False
inMiniumSubstrateScore = -55
inNormalizeValues = False
inPlotSubsetOfSubstrates = False
inPrintPredictionAccuracy = False
inInspectExperimentalES = True
inExperimentalESUpperLimit = 3.6
inExperimentalESLowerLimit = 3.0
inInspectPredictedES = False
inPredictedESUpperLimit = 3.5
inPredictedESLowerLimit = 2.5
inSetAxisLimits = False
inPlotSubstrateText = False
inTestBinnedSubES = True
inSaveBinnedSubES = False



# ==================================== Set Parameters ====================================
# Colors:
white = '\033[38;2;255;255;255m'
greyDark = '\033[38;2;144;144;144m'
purple = '\033[38;2;189;22;255m'
magenta = '\033[38;2;255;0;128m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
blue = '\033[38;5;51m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'

# Print options
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.float_format', '{:,.3f}'.format)

# Load: Dataset labels
enzymeName, filesInitial, filesFinal, labelAAPos = getFileNames(enzyme=inEnzymeName)
print(f'Label: {labelAAPos}')
# inMotifPositions = labelAAPos
motifLen = len(inMotifPositions)
motifFramePos = [inIndexNTerminus, inIndexNTerminus + motifLen]


# =================================== Initialize Class ===================================
ngs = NGS(
    enzyme=inEnzymeName, enzymeName=enzymeName, substrateLength=len(labelAAPos),
    filterSubs=True, fixedAA=inFixedResidue, fixedPosition=inFixedPosition,
    excludeAAs=inExcludeResidues, excludeAA=inExcludedResidue,
    excludePosition=inExcludedPosition, minCounts=inMinimumSubstrateCount,
    minEntropy=None, figEMSquares=inShowEnrichmentAsSquares, xAxisLabels=labelAAPos,
    xAxisLabelsMotif=inMotifPositions, printNumber=inPrintNumber,
    showNValues=inShowSampleSize, bigAAonTop=inBigLettersOnTop,
    findMotif=False, folderPath=inPathFolder, filesInit=filesInitial,
    filesFinal=filesFinal, releasedCounts=True, plotPosS=inPlotEntropy,
    plotFigEM=inPlotEnrichmentMap, plotFigEMScaled=inPlotEnrichmentMapScaled,
    plotFigLogo=inPlotLogo, plotFigWebLogo=inPlotWeblogo,
    plotFigMotifEnrich=inPlotMotifEnrichment, plotFigWords=inPlotWordCloud,
    wordLimit=inLimitWords, wordsTotal=inTotalWords, plotFigBars=inPlotBarGraphs,
    NSubBars=inPlotNBars, plotFigPCA=inPlotPCA, numPCs=inNumberOfPCs,
    NSubsPCA=inTotalSubsPCA, plotSuffixTree=inPlotSuffixTree,
    saveFigures=inSaveFigures, setFigureTimer=inSetFigureTimer
)



# ====================================== Load Data =======================================
# Load: Counts
countsInitial, countsInitialTotal = ngs.loadCounts(filter=False, fileType='Initial Sort',
                                                   dropColumn=inDropResidue)

# Calculate: Initial sort probabilities
if inUseCodonProb:
    # Evaluate: Degenerate codon probabilities
    rfInitial = ngs.calculateProbCodon(codonSeq=inCodonSequence)
else:
    rfInitial = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
                                fileType='Initial Sort', calcAvg=inAvgInitialProb)
    # if len(labelAAPos) == len(inMotifPositions):
    #     rfInitial = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
    #                                    fileType='Initial Sort')
    # else:
    #     rfInitial = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
    #                                    fileType='Initial Sort', calcAvg=True)



# ===================================== Run The Code =====================================
# Get dataset tag
ngs.getDatasetTag(combinedMotifs=True, useCodonProb=inUseCodonProb, codon=inCodonSequence)

# Load: Substrates
if inSaveCSV and inUseBgSubs or inFindSequences:
    substratesInitial, totalSubsInitial = ngs.loadUnfilteredSubs(loadInitial=True)

# Load: Substrate motifs
motifs, motifsCountsTotal, substratesFiltered = ngs.loadMotifSeqs(
    motifLabel=inMotifPositions, motifIndex=motifFramePos
)

# Display current sample size
ngs.recordSampleSize(
    NInitial=countsInitialTotal, NFinal=motifsCountsTotal, NFinalUnique=len(motifs.keys())
)

# Evaluate dataset
combinedMotifs = False
if len(ngs.motifIndexExtracted) > 1:
    combinedMotifs = True

# Load: Motif counts
countsMotifs, countsRelCombined, countsRelCombinedTotal = ngs.loadMotifCounts(
        motifLabel=inMotifPositions, motifIndex=motifFramePos, returnList=True
)

# Calculate: RF
rfCombinedReleasedMotif = ngs.calculateRFCombinedMotif(
    countsCombinedMotifs=countsRelCombined
)

# Calculate: Positional entropy
ngs.calculateEntropy(rf=rfCombinedReleasedMotif, combinedMotifs=combinedMotifs)

# Evaluate statistics
if inPlotStats and len(countsMotifs) > 1:
    ngs.fixedMotifStats(
        countsList=countsMotifs, initialRF=rfInitial,
        motifFrame=inMotifPositions, datasetTag=ngs.datasetTag
    )

# Calculate enrichment scores
ngs.calculateEnrichment(
    rfInitial=rfInitial, rfFinal=rfCombinedReleasedMotif, combinedMotifs=combinedMotifs
)

# Create csv
if inSaveCSV:
    if motifLen != len(labelAAPos):
        substrates = motifs
    else:
        substrates = substratesFiltered
    if inUseBgSubs:
        ngs.saveSubstrateCSV(
            seqs=substrates, initialRF=rfInitial, finalRF=rfCombinedReleasedMotif,
            minCounts=inMinSubsCSV, seqsBg=substratesInitial, excludeAA=inExcludeSeq,
            maxCountsBg=inMaxBgSubstrateCount, mod=inModulo, modScale=inScaleModulo,
            chopSeq=inSubLengthCSV
        )
    else:
        ngs.saveSubstrateCSV(
            seqs=substrates, initialRF=rfInitial, finalRF=rfCombinedReleasedMotif,
            minCounts=inMinSubsCSV, chopSeq=inSubLengthCSV
        )


# Find sequences
if inFindSequences:
    ngs.findSequence(substrates=substratesInitial, sequence=inFindSeq,
                     sortType='Initial Sort')
    ngs.findSequence(substrates=substratesFiltered, sequence=inFindSeq,
                     sortType='Final Sort')
    if combinedMotifs:
        ngs.findSequence(substrates=motifs, sequence=inFindSeq,
                         sortType='Motifs', combinedMotifs=combinedMotifs)
if inFindAAInSequence:
    # ngs.findAAInSequence(substrates=substratesInitial, AA=inFindAA, idxPos=inAAPos,
    #                  sortType='Initial Sort')
    # ngs.findAAInSequence(substrates=substratesFiltered, AA=inFindAA, idxPos=inAAPos,
    #                  sortType='Final Sort')
    if combinedMotifs:
        ngs.findAAInSequence(substrates=motifs, AA=inFindAA, idxPos=inAAPos,
                         sortType='Motifs', combinedMotifs=combinedMotifs)


# Predict substrate activity
if inPredictActivity:
    ngs.predictActivity(
        activityExp=inSubstrateActivity, errorBars=inErrorBars,
        finalRF=rfCombinedReleasedMotif, initialRF=rfInitial, predModel=ngs.datasetTag,
        predLabel=inPredictionLabel, combinedMotifs=combinedMotifs)

    # # Square Root predictions
    # ngs.predictActivity(
    #     activityExp=inSubstrateActivity, errBars=inErrorBars,
    #     finalRF=rfCombinedReleasedMotif, initialRF=rfInitial, predModel=ngs.datasetTag,
    #     predLabel=f'{inPredictionLabel} - Square Root', combinedMotifs=combinedMotifs)


    if not inPredictActivity:
        ngs.predictActivityHeatmap(predSubstrates=inPredictSubstrates,
                                   predModel=ngs.datasetTag, predLabel=inPredictionLabel,
                                   RF=rfCombinedReleasedMotif, rankScores=inRankScores,
                                   scaleEMap=inScalePredMatrix)

if inPlotFilteredSubs or inPlotWordCloud or inPlotBarGraphs:
    # Plot count related figures
    ngs.processSubstrates(subsInit=substratesInitial, subsFinal=substratesFiltered,
                          motifs=motifs, subLabel=inMotifPositions,
                          combinedMotifs=combinedMotifs)


    # # Evaluate: Motif Sequences
    # Count fixed substrates
    motifCountsFinal, motifsCountsTotal = ngs.countResidues(substrates=motifs,
                                                            datasetType='Final Sort')

    # Calculate: RF
    rfMotif = ngs.calculateRF(counts=motifCountsFinal, N=motifsCountsTotal,
                                fileType='Final Sort')

    # Calculate: Positional entropy
    ngs.calculateEntropy(rf=rfMotif, combinedMotifs=combinedMotifs)

    # Calculate: AA Enrichment
    ngs.calculateEnrichment(rfInitial=rfInitial, rfFinal=rfMotif,
                            combinedMotifs=combinedMotifs)


if inPredictCodonsEnrichment:
# Evaluate codon
    rfInitial = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
                                     fileType='Initial Sort', calcAvg=True)
    probCodon = ngs.calculateProbCodon(codonSeq=inCodonSequence)
    ngs.codonPredictions(codon=inCodonSequence, codonProb=probCodon, substrates=motifs)