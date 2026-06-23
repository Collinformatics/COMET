from functions import getFileNames, NGS
import os
import pandas as pd
import sys



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inEnzymeName = 'Mpro2'
inPathFolder = os.path.join('Enzymes', inEnzymeName)
inSaveFigures = True
inSetFigureTimer = False

# Input 2: Experimental Parameters
inMotifPositions = ['P4','P3','P2','P1','P1\'','P2\'']
# inMotifPositions = ['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
# inMotifPositions = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
inIndexNTerminus = 0 # Define the index if the first AA in the motif

# Input 3: Computational Parameters
inFixedResidue = 'Q' # ['R', ['A', 'G']]
inFixedPosition = [4,5,6]
inExcludeResidues = False
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
inPlotBinnedSubstrateES = False
inPlotBinnedSubstratePrediction = False
inPlotCounts = False
inPlotFilteredSubs = False
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
    inPredictionLabel = ''
    inActivityMpro = [0.54, 0.66, 0.34, 0.0, 0.26, 0.61, 0.0, 0.22]
    inSubstrateActivity = {
        'AVLQSG': inActivityMpro[0],
        'VILQSG': inActivityMpro[1],
        'VILQTG': inActivityMpro[2],
        'VILQSP': inActivityMpro[3],
        'VILHSG': inActivityMpro[4],
        'VIMQSG': inActivityMpro[5],
        'VPLQSG': inActivityMpro[6],
        'NILQSG': inActivityMpro[7],
    }
    inActivityMpro2 = [0.76, 0.82, 0.38, 0.0, 0.32, 0.66, 0.0, 0.23]
    inSubstrateActivity = {
        'AVLQSG': inActivityMpro2[0],
        'VILQSG': inActivityMpro2[1],
        'VILQTG': inActivityMpro2[2],
        'VILQSP': inActivityMpro2[3],
        'VILHSG': inActivityMpro2[4],
        'VIMQSG': inActivityMpro2[5],
        'VPLQSG': inActivityMpro2[6],
        'NILQSG': inActivityMpro2[7],
    }
    inErrorBars = [0.01, 0.058, 0.025, 0.0, 0.027, 0.044, 0.0, 0.033], # Avg stdev
    inErrorBars = []
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
if inSaveCSV and inUseBgSubs or inFindSequences or inPlotFilteredSubs:
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

if inPlotFilteredSubs:
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