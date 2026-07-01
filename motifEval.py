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
inFixedResidue = [['L','M'],'L']
inFixedPosition = [[4,5], [6,7]]
# inFixedResidue = ['R',['A','G','S']] # ['R',['A','G']]
# inFixedPosition = [[3,4,5,6],[4,5,6,7]]
inFixedResidue = 'Q'
inFixedPosition = [4,5,6]
inExcludeResidues = False
inExcludedResidue = ['A']
inExcludedPosition = [9,10]
inMinimumSubstrateCount = 1
inCodonSequence = 'NNS' # Baseline probs of degenerate codons (can be N, S, or K)
inUseCodonProb = False # Use AA prob from inCodonSequence to calculate enrichment
inAvgInitialProb = True
inDropResidue = [''] # To drop: inDropResidue = ['R9'], For nothing: inDropResidue = []

# Input 4: Figures
# inPlotPCA = False # PCA plot of an individual fixed frame
inBlockFigures = True
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
inPlotCounts = False
inShowSampleSize = True # Include the sample size in your figures
if inBlockFigures:
    inPlotEntropy = False
    inPlotEnrichmentMap = False
    inPlotEnrichmentMapScaled = False
    inPlotLogo = False
    inPlotWeblogo = False
    inPlotMotifEnrichment = False
    # inPlotWordCloud = False
    inPlotStats = False
    # inPlotBarGraphs = False
    inPlotPCA = False
    inPlotSuffixTree = False
    inPlotCounts = False

# Input 5: CSV
inSaveCSV = False # Save substrates in a csv file
inMinSubsCSV = 1000 # Minimum counts for saved substrates
inSubLengthCSV = 6 # If: False, use full seq. If: 6, use 6 AA seq
inUseBgSubs = True
inExcludeSeq = 'LQ'
inMaxBgSubstrateCount = 30
inModulo = 40000 # Increase to select fewer Bg substrates
inScaleModulo = True # Continually increase modulus value to keep fewer low count bg subs

# Input 6: Find Protein Sequences
inFindSequences = False
inFindSeq = ['AVLQSG', 'VILQSG','VILQTG','VILQSP','VILHSG','VIMQSG','VPLQSG','NILQSG']
inFindAAInSequence = False
inFindAA = ['A', 'F', 'W']
inAAPos = 4

# Input 7: Predict Activity
inPredictActivity = True
inRotateXLabel = 0
inPredictSubstrates = []
inUseNaturalSubs = False
if inUseNaturalSubs:
    inPredictionLabel = 'pp1a/b Substrates'
    inPredictSubstrates = ['AVLQSGFR', 'VTFQSAVK', 'ATVQSKMS', 'ATLQAIAS',
                           'VKLQNNEL', 'VRLQAGNA', 'PMLQSADA', 'TVLQAVGA',
                           'ATLQAENV', 'TRLQSLEN', 'PKLQSSQA']
    inSubstrateActivity = {}
    for substrate in inPredictSubstrates:
        sub = substrate[inIndexNTerminus:inIndexNTerminus+len(inMotifPositions)]
        inSubstrateActivity[sub] = 50.0
    inErrorBars = [] # Avg st. dev.
elif 'mpro1' in inEnzymeName.lower() or 'mpro' == inEnzymeName.lower():
    inPredictionLabel = ''
    inActivityMpro = [32.1, 39.1, 14.9, 0.0, 16.0, 36.5, 0.0, 15.6]
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
    inErrorBars = [] # Avg st. dev.
elif 'mpro2' in inEnzymeName.lower():
    inPredictionLabel = ''
    inActivityMpro2 = [46.1, 49.5, 14.5, 0.0, 13.1, 37.0, 0.0, 16.1]
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
    inErrorBars = [0.01, 0.058, 0.025, 0.0, 0.027, 0.044, 0.0, 0.033] # Avg st. dev.
    inErrorBars = []
elif 'mmp7' in inEnzymeName.lower():
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
else:
    inPredictActivity = False
    inPredictionLabel = ''
    inSubstrateActivity = {}
inErrorBars = []
inEMapStartIndex = 0  # Sub: ACDEFGHI, if idx = 1 start at C
inRankScores = False
inScalePredMatrix = False  # Scale EM by ΔS

# Input 8: Heatmap
inShowEnrichmentScores = True # Both cannot be True
inShowEnrichmentAsSquares = False # Both cannot be True

# Input 9: Enrichment Logo
inBigLettersOnTop = False

# Input 10: Bar Graphs
inPlotNBars = 50

# Input 11: Word Cloud
inLimitWords = True
inTotalWords = inPlotNBars

# Input 12: PCA
inNumberOfPCs = 2
inTotalSubsPCA = int(5*10**4)

# Input 13: Printing The Data
inPrintLoadedSubs = True
inPrintSampleSize = True
inPrintCounts = True
inPrintRF = True
inPrintES = True
inPrintEntropy = True
inPrintMotifData = True
inPrintNumber = 10



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
    NSubsPCA=inTotalSubsPCA, plotSuffixTree=False,
    saveFigures=inSaveFigures, setFigureTimer=inSetFigureTimer
)


# ===================================== Run The Code =====================================
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

# Get dataset tag
ngs.getDatasetTag(combinedMotifs=True, useCodonProb=inUseCodonProb, codon=inCodonSequence)

# Load: Substrates
substratesInitial, totalSubsInitial = None, None
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

# Plot count related figures
ngs.processSubstrates(motifs=motifs, subLabel=inMotifPositions,
                      combinedMotifs=combinedMotifs)

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
        finalRF=rfCombinedReleasedMotif, initialRF=rfInitial,
        predModel=ngs.datasetTag, predLabel=inPredictionLabel,
        combinedMotifs=combinedMotifs, rotateLabel=inRotateXLabel
    )