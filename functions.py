import base64
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
from concurrent.futures import ProcessPoolExecutor
from counter import counter
import gzip
import io
from itertools import batched
import logomaker
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os.path
import pandas as pd
import pickle as pk
import queue
import seaborn as sns
import secrets
import sys
import threading
import time
from wordcloud import WordCloud




defaultResidues = (
    ('Alanine', 'Ala', 'A'), ('Arginine', 'Arg', 'R'), ('Asparagine', 'Asn', 'N'),
    ('Aspartic Acid', 'Asp', 'D'), ('Cysteine', 'Cys', 'C'), ('Glutamic Acid', 'Glu', 'E'),
    ('Glutamine', 'Gln', 'Q'), ('Glycine', 'Gly', 'G'), ('Histidine', 'His ', 'H'),
    ('Isoleucine', 'Ile', 'I'), ('Leucine', 'Leu', 'L'), ('Lysine', 'Lys', 'K'),
    ('Methionine', 'Met', 'M'), ('Phenylalanine', 'Phe', 'F'), ('Proline', 'Pro', 'P'),
    ('Serine', 'Ser', 'S'), ('Threonine', 'Thr', 'T'), ('Tryptophan', 'Typ', 'W'),
    ('Tyrosine', 'Tyr', 'Y'), ('Valine', 'Val', 'V')
)


# Generate figures entirely in memory without opening a window
matplotlib.use('Agg')  # Use a non-interactive backend for servers


class WebApp:
    def __init__(self):
        # Params: Job
        self.jobDone = False
        self.jobParams = {}
        self.datasetTag = ''
        self.datasetTagMotif = ''

        # Params: Dataset
        self.enzymeName = ''
        self.seqLength = False
        self.minCounts = 1
        self.printN = 10
        self.roundVal = 3
        self.xAxisLabel = []
        self.entropy = pd.DataFrame(0.0, index=[], columns=['∆S'])
        self.entropyMax = None
        self.subsExp = {}
        self.subsExpAll = {}
        self.countsExp = pd.DataFrame()
        self.countExpTotal = 0
        self.countExpUnique = 0
        self.rfExp = None
        self.rfExpScaled = None
        self.saveTagExp = {}
        self.subsBg = {}
        self.countsBg = pd.DataFrame()
        self.countBgTotal = 0
        self.countBgUnique = 0
        self.rfBg = None
        self.eMap = None
        self.eMapScaled = None
        self.eMapReleased = None
        self.eMapReleasedScaled = None
        self.saveTagBg = {}
        self.saveTagFig = {}
        self.datasetTypes = {'Exp': 'Experimental', 'Bg': 'Background'}

        # Params: COMET
        self.iteration = 0
        self.minS = 0.65
        self.minES = 0
        self.minESRel = -1
        self.motifFilter = False
        self.motifLen = 0
        self.motifPos = pd.DataFrame()
        self.substrateProfile = pd.DataFrame()

        # Params: Files
        self.queueLog = queue.Queue()
        self.fileExp = []
        self.fileExpRev = []
        self.fileBg = []
        self.fileBgRev = []
        self.pathDir = ''
        self.pathData = ''
        self.pathFigs = ''
        self.pathLog = ''

        # Params: Process dna
        self.seq5Prime = ''
        self.seq3Prime = ''
        self.minPhred = False

        # Params: Filter Dataset
        self.filterPos = False
        self.fixAA = {}
        self.exclAA = {}

        # Params: Figures
        self.numSamples = 50
        self.figureResolution = 600
        self.titleCombined = ''
        self.titleReleased = ''
        self.titleWeblogo = ''
        self.titleWeblogoCombined = ''
        self.titleWeblogoReleased = ''
        self.titleWords = ''
        self.titleWordsCombined = ''
        self.figEMSquares = False
        self.figSize = (9.5, 8) # (width, height)
        self.figSizeSq = (5, 8)
        self.figSizeMini = (self.figSize[0], 6)
        self.residueLabelType = 2  # 0 = full AA name, 1 = 3-letter code, 2 = 1 letter
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.colorsAA = self.residueColors()
        self.residues = defaultResidues
        self.AA = [residue[2] for residue in self.residues]
        self.bigAAonTop = False
        self.figures = {}

        # Colors
        self.orange = '#FA8128'
        self.orangeBurnt = '#BF5700'

        pd.set_option('display.max_columns', None)
        pd.options.display.float_format = '{:,.3f}'.format


    @staticmethod
    def createCustomColorMap(colorType):
        colorType = colorType.lower()
        if colorType == 'counts':
            useGreen = True
            if useGreen:
                # Green
                colors = ['#FFFFFF','#ABFF9B','#39FF14','#2E9418','#2E9418','#005000']
            else:
                # Orange
                colors = ['white','white','#FF76FA','#FF50F9','#FF00F2','#CA00DF','#BD16FF']
        elif colorType == 'stdev':
            colors = ['white','white','#FF76FA','#FF50F9','#FF00F2','#CA00DF','#BD16FF']
        elif colorType == 'word cloud':
            # ,'#F2A900','#2E8B57','black'
            colors = ['#CC5500','#CC5500','#F79620','#FAA338',
                      '#00C01E','#1D680D','#003000','black']
        elif colorType == 'em':
            colors = ['navy','royalblue','dodgerblue','lightskyblue','white',
                      'white','lightcoral','red','firebrick','darkred']
        else:
            print(f'ERROR: Cannot create colormap. '
                  f'Unrecognized colorType parameter: {colorType}\n')
            sys.exit(1)

        # Create colormap
        if len(colors) == 1:
            colorList = [(0, colors[0]), (1, colors[0])]
        else:
            colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        return LinearSegmentedColormap.from_list('custom_colormap', colorList)


    @staticmethod
    def residueColors():
        color = ['darkgreen', 'firebrick', 'deepskyblue', 'pink', 'navy', 'black', 'gold']
        # Aliphatic, Acidic, Basic, Hydroxyl, Amide, Aromatic, Sulfur

        return {
            'A': color[0],
            'R': color[2],
            'N': color[4],
            'D': color[1],
            'C': color[6],
            'E': color[1],
            'Q': color[4],
            'G': color[0],
            'H': color[2],
            'I': color[0],
            'L': color[0],
            'K': color[2],
            'M': color[6],
            'F': color[5],
            'P': color[0],
            'S': color[3],
            'T': color[3],
            'W': color[5],
            'Y': color[5],
            'V': color[0]
        }


    @staticmethod
    def getKey(app):  # required for CSRF
        app.config['SECRET_KEY'] = secrets.token_hex(nbytes=32)
        # print(f'Key: {app.config['SECRET_KEY']}\n'
        #       f'Len: {len(app.config['SECRET_KEY'])}\n')


    @staticmethod
    def pressButton(message):
        print(f'Received data: {message}')
        return {'key': 'Returned data'}


    def encodeFig(self, fig):
        # Save to a memory buffer instead of disk
        buffer = io.BytesIO()
        plt.savefig(
            buffer, format='png', bbox_inches='tight', dpi=self.figureResolution
        )
        buffer.seek(0)

        # Encode as base64 for embedding in HTML
        figBase64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        plt.close(fig)  # or fig.clear() to release memory
        buffer.close()

        return figBase64


    def getDatasetTag(self, subProfile=False, log=False):
        self.fixAA = dict(sorted(self.fixAA.items()))
        self.exclAA = dict(sorted(self.exclAA.items()))
        tagFix = 'Fix '
        tagExcl = 'Excl '

        # Evaluate filters
        if self.exclAA:
            for index, (pos, AA) in enumerate(self.exclAA.items()):
                if len(AA) > 1:
                    tagExcl += f'[{','.join(AA)}]@{pos} '
                else:
                    tagExcl += f'{AA}@{pos} '
            tagExcl = tagExcl[:-1]
        if self.fixAA:
            for index, (pos, AA) in enumerate(self.fixAA.items()):
                if len(AA) > 1:
                    tagFix += f'[{','.join(AA)}]@{pos} '
                else:
                    tagFix += f'{AA}@{pos} '
            tagFix = tagFix[:-1]
        if tagExcl != 'Excl ' and tagFix != 'Fix ':
            self.datasetTag = f'{tagExcl} {tagFix}'
        elif tagFix != 'Fix ':
            self.datasetTag = tagFix
        elif tagExcl != 'Excl ':
            self.datasetTag = tagExcl
        else:
            self.datasetTag = 'Unfiltered'
        if self.motifFilter and  not self.datasetTagMotif:
            self.datasetTagMotif = f'Motif - {self.datasetTag}'
            self.jobParams['Dataset Tag'] = self.datasetTagMotif
        else:
            self.jobParams['Dataset Tag'] = self.datasetTag
        if log:
            self.log(f'Filter: {self.datasetTag}')

        # Initialize: Save tags
        self.saveTagExp = {
            'subs': f'{self.enzymeName}-Subs_Exp-{self.datasetTag}-'
                    f'MinCounts_{self.minCounts}-{self.seqLength}AA.pkl',
            'counts': f'{self.enzymeName}-AA_Counts_Exp-{self.datasetTag}-'
                      f'MinCounts_{self.minCounts}-{self.seqLength}AA.csv'
        }
        self.saveTagBg = {
            'subs': f'{self.enzymeName}-Subs_Bg-{self.datasetTag}-'
                    f'MinCounts_{self.minCounts}-{self.seqLength}AA.pkl',
            'counts': f'{self.enzymeName}-AA_Counts_Bg-{self.datasetTag}-'
                      f'MinCounts_{self.minCounts}-{self.seqLength}AA.csv',
        }
        if subProfile:
            datasetTag = self.datasetTagMotif.replace('Motif - ', '')
            for tag, path in self.saveTagExp.items():
                path = path.replace(
                    f'{self.enzymeName}', f'{self.enzymeName}-SubProfile'
                ).replace(
                    f'{self.datasetTag}', f'{datasetTag}'
                )
                self.saveTagExp[tag] = path
            print(f'Tag:\n{self.saveTagExp}')

        self.saveTagFig = (f'{self.enzymeName}-Fig-{self.datasetTag}-'
                           f'Min_Counts_{self.minCounts}-{self.seqLength}AA')


    def getSaveTag(self, saveTag=''):
        # Evaluate filters
        if self.motifFilter:
            tag = 'COMET'
        else:
            tag = self.datasetTag
            tag = tag.replace(' Fix ', '-Fix_').replace('Fix ', 'Fix_')
            tag = tag.replace('Excl ', 'Excl_').replace(' ', '_')
        if saveTag:
            return saveTag.replace(self.datasetTag, tag)
        else:
            return tag


    def getFilter(self, data):
        self.fixAA = {}
        self.exclAA = {}
        if 'filterPos' in data.keys():
            self.filterPos = data['filterPos']
            # Get filter params
            for key, value in data.items():
                if 'fix' in key:
                    self.fixAA[key.replace('fix', '')] = value
                if 'excl' in key:
                    self.exclAA[key.replace('excl', '')] = value
        self.getDatasetTag(log=True)


    def initDataStructures(self):
        # Initialize data structures
        self.fileExp = []
        self.fileExpRev = []
        self.fileBg = []
        self.fileBgRev = []
        self.subsExp = {}
        self.subsBg = {}
        self.xAxisLabel = [f'R{index}' for index in range(1, self.seqLength + 1)]
        self.countsExp = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.countsBg = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)


    def jobInit(self, form, job):
        self.jobParams['Job ID'] = form['jobID']
        self.jobParams['Job'] = job

        # Initialize directories
        # self.pathDir = os.path.join('dset', f"{form['enzymeName']}-{form['jobID']}")
        self.pathDir = os.path.join('dset', form['enzymeName'])
        self.pathData = os.path.join(self.pathDir, 'data')
        self.pathFigs = os.path.join(self.pathDir, 'figures')
        self.pathLog = os.path.join(self.pathDir, 'log.txt')
        if self.pathDir is not None:
            if not os.path.exists(self.pathDir):
                os.makedirs(self.pathDir, exist_ok=True)
        if self.pathData is not None:
            if not os.path.exists(self.pathData):
                os.makedirs(self.pathData, exist_ok=True)
        if self.pathFigs is not None:
            if not os.path.exists(self.pathFigs):
                os.makedirs(self.pathFigs, exist_ok=True)
            # else:
            #     import shutil
            #     # Remove everything inside the directory
            #     for filename in os.listdir(self.pathFigs):
            #         path = os.path.join(self.pathFigs, filename)
            #         if os.path.isfile(path) or os.path.islink(path):
            #             os.unlink(path)  # delete file or link
            #         elif os.path.isdir(path):
            #             shutil.rmtree(path)  # delete subdirectory
            #     # time.sleep(5)
            #     os.makedirs(self.pathFigs, exist_ok=True)

        self.log() # Clear the log
        self.log('================================ Job Summary '
                 '=================================')

        # Record job params
        self.log(f'Job ID: {self.jobParams['Job ID']}')
        self.log(f'Job: {job}')
        self.enzymeName = form['enzymeName']
        self.jobParams['Enzyme Name'] = self.enzymeName
        self.log(f'Enzyme: {self.enzymeName}')
        self.seqLength = int(form['seqLength'])
        self.jobParams['Substrate Length'] = self.seqLength
        self.log(f'Substrate Length: {self.seqLength}')

        # Initialize params
        self.initDataStructures()
        self.figures = {}
        self.motifFilter = False


        def addFile(data, value):
            if isinstance(value, list):
                for f in value:
                    data.append(f)
            else:
                data.append(value)

        # Get the files
        for key, value in form.items():
            if 'fileExpRev' in key:
                addFile(self.fileExpRev, value)
            elif 'fileExp' in key:
                addFile(self.fileExp, value)
            elif 'fileBgRev' in key:
                addFile(self.fileBgRev, value)
            elif 'fileBg' in key:
                addFile(self.fileBg, value)

        # Job dependant parameters
        if job == 'Process DNA':
            self.seq5Prime = form['seq5Prime']
            self.seq3Prime = form['seq3Prime']
            self.minCounts = int((round(float(form['minCounts']))))
            self.jobParams['Minimum Count'] = self.minCounts
            self.minPhred = int(round(float(form['minPhred']))) if form['minPhred'] != '' else 0
            self.log(f'5\' Sequence: {self.seq5Prime}\n'
                     f'3\' Sequence: {self.seq3Prime}\n'
                     f'Min Phred Score: {self.minPhred}')
        elif job == 'Filter AA':
            print(f'Job: {job}')
        elif job == 'Filter Motif':
            print(f'Job: {job}')
            self.iteration = 0
            self.motifFilter = True
        elif job == 'Combine Motifs':
            print(f'Job: {job}')
            self.motifLen = form['motifLength']
        else:
            print('ERROR: What Script Is Running')
            sys.exit()

        # Add: Min counts
        if self.fileExp:
            print(f'File Exp: {type(self.fileExp)}\n'
                  f'* {self.fileExp}')
        if self.fileExpRev:
            print(f'File Exp Rev: {type(self.fileExpRev)}\n'
                  f'* {self.fileExpRev}')
        if self.fileBg:
            print(f'File Bg: {type(self.fileBg)}\n'
                  f'* {self.fileBg}')
        if self.fileBgRev:
            print(f'File Bg Rev: {type(self.fileBgRev)}\n'
                  f'* {self.fileBgRev}')

        # Get the filter and initialize the data structures
        self.getFilter(form)
        # self.fixAA['R4'] = ['C', 'G', 'H', 'K', 'T'] ## Delete me
        # self.getDatasetTag() ## Delete me
        # print(f'\nFix AA Tag: fixR4: {self.fixAA}')


    def log(self, txt=None):
        if txt is None:
            with open(self.pathLog, 'w'):
                pass
        else:
            with open(self.pathLog, 'a') as log:
                log.write(f'{txt}\n')


    def logInQueue(self, logQueue):
        with open(self.pathLog, 'a') as log:
            while not logQueue.empty():
                log.write(logQueue.get() + '\n')


    def logErrorFn(self, function, msg, getStr=False):
        if getStr:
            return (f'\n========================================='
                    f'========================================\n'
                     f'========================================='
                     f'========================================\n'
                     f'***** ERROR: {function} *****\n\n'
                     f'{msg}\n'
                     f'========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================')
        else:
            self.log(f'\n========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n'
                     f'***** ERROR: {function} *****\n\n'
                     f'{msg}\n'
                     f'========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================')
            sys.exit(1)


    def processSubs(self, substrates, datasetType, filteredAA):
        self.log('\n\n================================= Substrates '
                 '=================================')
        self.log(f'Dataset: {datasetType}')

        # Inspect sequences
        if not filteredAA:
            filteredSubs = {}
            for substrate, count in substrates.items():
                for AA in substrate:
                    if AA not in self.AA:
                        filteredSubs[substrate] = count
            if filteredSubs:
                self.log(f'\nFiltering Substrates:\n'
                         f'     If a substrate contains an '
                         f'unaccented AA it will be removed.\n'
                         f'     Accepted: {self.AA}\n\n'
                         f'     Removed Substrates:')
                for substrate, count in filteredSubs.items():
                    substrates.pop(substrate, count)
                    self.log(f'          {substrate}: {count}')

        # Sort data
        substrates = dict(sorted(substrates.items(),
                                 key=lambda item: item[1], reverse=True))


        # Count AAs
        countMatrix = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.log('\nSubstrate Totals:')
        if datasetType == self.datasetTypes['Exp']:
            self.subsExp = substrates
            self.countExpTotal = sum(substrates.values())
            self.countExpUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countExpTotal:,}\n'
                     f'    Unique Substrates: {self.countExpUnique:,}\n')

            # Record job params
            self.jobParams['Total Experimental Substrates'] = f'{self.countExpTotal:,}'
            self.jobParams['Unique Experimental Substrates'] = f'{self.countExpUnique:,}'
        elif datasetType == self.datasetTypes['Bg']:
            self.subsBg = substrates
            self.countBgTotal = sum(substrates.values())
            self.countBgUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countBgTotal:,}\n'
                     f'    Unique Substrates: {self.countBgUnique:,}\n')

            # Record job params
            self.jobParams['Total Background Substrates'] = f'{self.countBgTotal:,}'
            self.jobParams['Unique Background Substrates'] = f'{self.countBgUnique:,}'
        else:
            self.logErrorFn(function='processSubs()',
                            msg=f'Unknown dataset type: {datasetType}')

        self.log(f'Top {self.printN:,} {datasetType} Sequences')
        for index, (sub, count) in enumerate(substrates.items()):
            if index >= self.printN:
                break
            self.log(f'     {sub}: {count}')

        # Save data
        self.saveSubstrates(substrates=substrates, datasetType=datasetType)

        # Count AAs
        countMatrix = self.countAA(substrates=substrates, countMatrix=countMatrix,
                                   datasetType=datasetType)
        if datasetType == self.datasetTypes['Bg']:
            self.countsBg = countMatrix
        else:
            self.countsExp = countMatrix

    def loadDNA(self, path, datasetType, queueData, queueLog, reverseRead):
        translate = True
        filename = path.filename if hasattr(path, 'filename') else path.name
        try:
            # Open the file
            if filename.endswith('.gz'):
                file_handle = gzip.open(path, 'rt')
            else:
                path.seek(0)  # ensure at start
                file_handle = io.StringIO(path.read().decode('utf-8'))
            data = None

            if path.filename.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                data = SeqIO.parse(file_handle, 'fastq')
            elif path.filename.endswith(('.fasta', '.fa', '.fasta.gz', '.fa.gz')):
                data = SeqIO.parse(file_handle, 'fasta')
            else:
                queueLog.put(self.logErrorFn(
                    function='loadDNA()',
                    msg=f'Unrecognized file\n     {path}',
                    getStr=True
                    )
                )
                translate = False

            # Translate the dna
            if translate:
                substrates = self.translate(
                    data, filename, datasetType, queueLog, reverseRead
                )
                queueData.put(substrates) # Put the substrates in the queue
        except Exception as e:
            queueLog.put(self.logErrorFn(
                function='loadDNA()',
                msg=f'Failed to load file:\n     {path}\n     {e}',
                getStr=True))


    def translate(self, data, fileName, datasetType, queueLog, revRead):
        queueLog.put('\n\n================================ Translate DNA '
                     '===============================')
        data = list(data)
        self.printN = 10

        queueLog.put(f'File Name: {fileName}')
        if revRead:
            queueLog.put('Read Type: Reverse Read')
        else:
            queueLog.put('Read Type: Forward Read')
        queueLog.put(f'  Dataset: {datasetType}')

        # Inspect the file
        useQS = False
        for datapoint in data:
            if 'phred_quality' in datapoint.letter_annotations:
                useQS = True
            break
        queueLog.put(f'  Eval QS: {useQS}')
        if useQS:
            queueLog.put(f'Min Phred: {self.minPhred}')
        queueLog.put('\n')

        # Inspect the datasetType parameter
        if (datasetType != self.datasetTypes['Exp']
                and datasetType != self.datasetTypes['Bg']):
            self.logErrorFn(function='translate()',
                            msg=f'Unknown dataset type: {datasetType}')


        def reverseComplement(seq):
            """Returns the reverse complement of a DNA sequence."""
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
            return ''.join(complement[base] for base in reversed(seq))


        def extractionEfficiency(totalSeqsDNA, totalSubsExtracted, fullSet=False):
            if totalSeqsDNA == 0:
                totalSeqsDNA = 1
            perExtracted = (totalSubsExtracted / totalSeqsDNA) * 100
            if fullSet:
                queueLog.put(f'\nExtraction Efficiency: - All Sequences')
            else:
                queueLog.put(f'\nExtraction Efficiency: {self.printN} Substrates')
            queueLog.put(f'     Evaluated DNA Sequences: {totalSeqsDNA:,}\n'
                         f'        Extracted Substrates: {totalSubsExtracted:,}\n'
                         f'       Extraction Efficiency: {round(perExtracted, 3)} %')


        def logDNA(datapoint, totalSeqs, totalSubsExtracted, useQS):
            """
                Scan dataset and log datapoints
            """
            totalSeqs += 1
            dna = str(datapoint.seq)
            # print(f'DNA: {dna}')
            if revRead:
                dna = reverseComplement(dna)
            queueLog.put(f'DNA Seq: {dna}')

            # Inspect full dna seq
            queueLog.put(f'Tags: {self.seq5Prime}, {self.seq5Prime in dna} - '
                  f'{self.seq3Prime}, {self.seq3Prime in dna}')
            if (self.seq5Prime in dna and self.seq3Prime in dna or
                    not self.seq5Prime and self.seq3Prime in dna or
                    self.seq5Prime in dna and not self.seq3Prime):

                # Find: Substrate indices
                start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                end = dna.find(self.seq3Prime)

                # Extract substrate dna seq
                if end:
                    substrateDNA = dna[start:end].strip()
                else:
                    substrateDNA = dna[start:].strip()
                queueLog.put(f'* Sub DNA: {substrateDNA}')
                if len(substrateDNA) == self.seqLength * 3:
                    substrate = str(Seq.translate(substrateDNA))
                    queueLog.put(f'Sub Seq: {substrate}')

                    # Inspect substrate seq: log only
                    if 'X' not in substrate and '*' not in substrate:
                        if useQS:
                            qs = datapoint.letter_annotations['phred_quality']
                            if end:
                                qs = qs[start:end]
                            else:
                                qs = qs[start:]
                            queueLog.put(f'     QS: {qs}')
                            if all(score < self.minPhred for score in qs):
                                print(f'Low QS: {qs}')
                                return totalSeqs, totalSubsExtracted
                        queueLog.put(f'* Keep Substrate')
                        totalSubsExtracted += 1
            queueLog.put('')

            return totalSeqs, totalSubsExtracted

        # Translate and log DNA
        totalSeqs = 0
        totalSubsExtracted = 0
        for idx, datapoint in enumerate(data):
            if totalSubsExtracted >= self.printN:
                break
            totalSeqs, totalSubsExtracted = logDNA(
                datapoint, totalSeqs, totalSubsExtracted, useQS
            )
        extractionEfficiency(totalSeqs, totalSubsExtracted)


        def processDNA(datapoint, totalSeqs, totalSubsExtracted, substrates, useQS):
            # Process datapoint
            totalSeqs += 1
            dna = str(datapoint.seq)
            if revRead:
                # print(f'Reverse DNA: {dna}')
                dna = reverseComplement(dna)

            # Inspect full dna seq
            if (self.seq5Prime in dna and self.seq3Prime in dna or
                    not self.seq5Prime and self.seq3Prime in dna or
                    self.seq5Prime in dna and not self.seq3Prime):

                # Find: Substrate indices
                start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                end = dna.find(self.seq3Prime)

                # Extract substrate dna seq
                if end:
                    substrateDNA = dna[start:end].strip()
                else:
                    substrateDNA = dna[start:].strip()
                if len(substrateDNA) == self.seqLength * 3:
                    substrate = str(Seq.translate(substrateDNA))

                    # Inspect substrate seq
                    if 'X' not in substrate and '*' not in substrate:
                        if useQS:
                            qs = datapoint.letter_annotations['phred_quality']
                            if end:
                                qs = qs[start:end]
                            else:
                                qs = qs[start:]
                            if all(score < self.minPhred for score in qs):
                                return substrates, totalSeqs, totalSubsExtracted
                        if substrate in substrates.keys():
                            substrates[substrate] += 1
                        else:
                            substrates[substrate] = 1
                        totalSubsExtracted += 1

            return substrates, totalSeqs, totalSubsExtracted

        # Translate DNA and record substrates
        substrates = {}
        totalSeqs = 0
        totalSubsExtracted = 0
        for datapoint in data:
            substrates, totalSeqs, totalSubsExtracted = processDNA(
                datapoint, totalSeqs, totalSubsExtracted, substrates, useQS
            )
            # if len(substrates.keys()) > 100: ## Delete me
            #     break
        extractionEfficiency(totalSeqs, totalSubsExtracted, fullSet=True)
        return substrates


    def evalDNA(self, form):
        self.jobInit(form, job='Process DNA')

        # Load the data
        threads = []
        queuesExp = []
        queuesExpLog = []
        queuesExpRev = []
        queuesExpRevLog = []
        queuesBg = []
        queuesBgLog = []
        queuesBgRev = []
        queuesBgRevLog = []
        if self.fileExp:
            for file in self.fileExp:
                queueExp = queue.Queue()
                queueLog = queue.Queue()
                queuesExp.append(queueExp)
                queuesExpLog.append(queueLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Exp'],
                          queueExp, queueLog, False,)
                )
                thread.start()
                threads.append(thread)
        if self.fileExpRev:
            for file in self.fileExpRev:
                queueExpRev = queue.Queue()
                queueExpRevLog = queue.Queue()
                queuesExpRev.append(queueExpRev)
                queuesExpRevLog.append(queueExpRevLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Exp'],
                          queueExpRev, queueExpRevLog, True,)
                )
                thread.start()
                threads.append(thread)
        if self.fileBg:
            for file in self.fileBg:
                queueBg = queue.Queue()
                queueLog = queue.Queue()
                queuesBg.append(queueBg)
                queuesBgLog.append(queueLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Bg'],
                          queueBg, queueLog, False,)
                )
                thread.start()
                threads.append(thread)
        if self.fileBgRev:
            for file in self.fileBgRev:
                queueBgRev = queue.Queue()
                queueBgRevLog = queue.Queue()
                queuesBgRev.append(queueBgRev)
                queuesBgRevLog.append(queueBgRevLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Bg'],
                          queueBgRev, queueBgRevLog, True,)
                )
                thread.start()
                threads.append(thread)
        # Wait for all threads to finish
        for thread in threads:
            thread.join()

        # Log the output
        if queuesExpLog:
            for log in queuesExpLog:
                self.logInQueue(log)
        if queuesExpRevLog:
            for log in queuesExpRevLog:
                self.logInQueue(log)
        if queuesBgLog:
            for log in queuesBgLog:
                self.logInQueue(log)
        if queuesBgRevLog:
            for log in queuesBgRevLog:
                self.logInQueue(log)

        def addSubs(data, subs):
            for sub, count in subs.items():
                if count >= self.minCounts:
                    if sub in data.keys():
                        data[sub] += count
                    else:
                        data[sub] = count

        # Get results from queue
        if self.fileExp:
            subs = {}
            for queueData in queuesExp:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in subs.keys():
                        subs[substrate] += count
                    else:
                        subs[substrate] = count
            subs = dict(sorted(subs.items(), key=lambda item: item[1], reverse=True))
            addSubs(self.subsExp, subs)
        if self.fileExpRev:
            subs = {}
            for queueData in queuesExpRev:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in subs.keys():
                        subs[substrate] += count
                    else:
                        subs[substrate] = count
            subs = dict(sorted(subs.items(), key=lambda item: item[1], reverse=True))
            addSubs(self.subsExp, subs)
        if self.fileBg:
            subs = {}
            for queueData in queuesBg:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in subs.keys():
                        subs[substrate] += count
                    else:
                        subs[substrate] = count
            subs = dict(sorted(subs.items(), key=lambda item: item[1], reverse=True))
            addSubs(self.subsBg, subs)
        if self.fileBgRev:
            subs = {}
            for queueData in queuesBgRev:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in subs.keys():
                        subs[substrate] += count
                    else:
                        subs[substrate] = count
            subs = dict(sorted(subs.items(), key=lambda item: item[1], reverse=True))
            addSubs(self.subsBg, subs)

        # Make figures
        if self.subsExp:
            # Sort substrates and count AA
            self.processSubs(substrates=self.subsExp,
                             datasetType=self.datasetTypes['Exp'],
                             filteredAA=False)

            # Plot counts
            self.figures['exp_counts'] = (
                self.plotCounts(
                    countedData=self.countsExp, totalCounts=self.countExpTotal,
                    datasetType=self.datasetTypes['Exp']
                )
            )
        if self.subsBg:
            # Sort substrates and count AA
            self.processSubs(substrates=self.subsBg,
                             datasetType=self.datasetTypes['Bg'],
                             filteredAA=False)

            # Plot counts
            self.figures['bg_counts'] = (
                self.plotCounts(
                    countedData=self.countsBg, totalCounts=self.countBgTotal,
                    datasetType=self.datasetTypes['Bg']
                )
            )

        if self.subsExp and self.subsBg:
            self.calculateRF()
            self.calculateEntropy()
            self.evalEnrichment()
        self.jobDone = True


    def loadSubstrates(self, path, queueData, queueLog):
        # with open(path, 'rb') as openedFile:  # Open file
        #     data = pk.load(openedFile)  # Access the data
        # queueData.put(data)
        # queueLog.put(f'     {path}')
        try:
            path.seek(0)  # ensure at start
            data = pk.load(path)  # BytesIO is already file-like
            queueData.put(data)
            queueLog.put(f'     {path}')
        except Exception as e:
            queueLog.put(self.logErrorFn(
                function='loadSubstrates()',
                msg=f'Failed to load file:\n     {path}\n     {e}',
                getStr=True))
            sys.exit(0)


    def loadCounts(self, path, queueData, queueLog):
        try:
            # Load file
            data = pd.read_csv(path, index_col=0)
            data = data.astype(int)
            queueData.put(data)
            queueLog.put(f'     {path}\n')
        except Exception as e:
            queueLog.put(self.logErrorFn(
                function='loadCounts()',
                msg=f'Failed to load file:\n     {path}\n     {e}',
                getStr=True))


    def evalSubs(self, form, filterMotifs=False):
        if filterMotifs:
            print('\nFilter Motif')
            self.jobInit(form, job='Filter Motif')
        else:
            self.jobInit(form, job='Filter AA')
        self.log('\n\n================================== Load Data '
                 '=================================')

        # Load the data
        threads = []
        queuesExp = []
        queuesExpLog = []
        queuesBg = []
        queuesBgLog = []
        for file in self.fileExp:
            queueExp = queue.Queue()
            queueLog = queue.Queue()
            queuesExp.append(queueExp)
            queuesExpLog.append(queueLog)
            thread = threading.Thread(
                target=self.loadSubstrates, args=(file, queueExp, queueLog)
            )
            thread.start()
            threads.append(thread)
        for file in self.fileBg:
            queueBg = queue.Queue()
            queueLog = queue.Queue()
            queuesBg.append(queueBg)
            queuesBgLog.append(queueLog)
            thread = threading.Thread(
                target=self.loadCounts, args=(file, queueBg, queueLog)
            )
            thread.start()
            threads.append(thread)

        # Wait for all threads to finish
        for thread in threads:
            thread.join()

        # Process input
        if queuesExpLog:
            self.log('Loading Substrates: Experimental')
            # for log in queuesExpLog:
            #     self.logInQueue(log)
        if queuesExp:
            for q in queuesExp:
                substrates = q.get()
                for substrate, count in substrates.items():
                    if substrate in self.subsExp.keys():
                        self.subsExp[substrate] += count
                    else:
                        self.subsExp[substrate] = count
            self.log('\nSubstrates:')
            if self.subsExp:
                i = 0
                for substrate, count in self.subsExp.items():
                    i += 1
                    self.log(f'    {substrate}, {count:,}')
                    if i >= self.printN:
                        break
                self.log('\n')
            else:
                self.logErrorFn(function='loadSubstrates()',
                                msg='No experimental substrates were loaded')
        if queuesBgLog:
            self.log('Loading Substrate Counts: Background')
        if queuesBg:
            for idx, q in enumerate(queuesBg):
                counts = q.get()
                self.countsBg += counts
                if len(queuesBg) > 1:
                    self.log(f'\nLoaded Count Set: {idx}\n{counts}\n')
            self.countBgTotal = sum(self.countsBg.iloc[:, 0])
            if self.countBgTotal == 0:
                print(f'DataFrame consists entirely of zeros.\n{self.countsBg}')
                self.logErrorFn(function='loadCounts()',
                                msg='No background substrates were loaded')
            self.log(f'\nBackground Counts:\n{self.countsBg}')
        self.motifLen = len(next(iter(self.subsExp)))
        print(f'Motif Length: {self.motifLen}')

        # Filter AAs
        self.filterSubs(plotEntropy=True)

        # Plot figures
        if filterMotifs:
            self.subsExpAll = self.subsExp
            # Count AAs #
            self.countAA(substrates=self.subsExp, countMatrix=self.countsExp,
                         datasetType=self.datasetTypes['Exp'])
            self.selectMotifPos()
        else:
            self.evalEnrichment()
        self.jobDone = True

    def filterSubs(self, plotEntropy=False, saveData=True):
        self.log('\n\n============================== Filter Substrates '
                 '=============================')
        self.log(f'Dataset tag: {self.datasetTag}')
        if self.fixAA or self.exclAA:
            totalSubs = 0
            for count in self.subsExp.values():
                totalSubs += count
            totalSubsUnique = len(self.subsExp.keys())
            self.log(f'\nUnfiltered Substrates:\n'
                     f'     Total: {totalSubs:,}\n'
                     f'    Unique: {totalSubsUnique:,}')

            if self.motifFilter:
                substrates = self.subsExpAll
            else:
                substrates = self.subsExp

            subs = {}
            print(f'Filter: {self.fixAA}')
            # print(f'\nFilters:\n* Fix: {self.fixAA}\n* Excl: {self.exclAA}\n')
            for substrate, count in substrates.items():
                # print(f'Substrate: {substrate}')
                keepSub = True
                for posExcl, exclAA in self.exclAA.items():
                    if not isinstance(exclAA, list):
                        exclAA = list(exclAA)
                    idx = int(posExcl.replace('R', '')) - 1
                    if substrate[idx] in exclAA:
                        keepSub = False
                        break
                if not keepSub:
                    continue

                for posFix, fixAA in self.fixAA.items():
                    if not isinstance(fixAA, list):
                        fixAA = list(fixAA)
                    idx = int(posFix.replace('R', '')) - 1
                    # print(f'* Pos: {posFix} {substrate[idx]} {fixAA}')
                    if substrate[idx] not in fixAA:
                        keepSub = False
                        # print(f'  Drop - {posFix} {substrate[idx]}')
                        break
                if keepSub:
                    # print('  keep')
                    subs[substrate] = count
            self.subsExp = subs
        self.countExpTotal = sum(self.subsExp.values())
        self.countExpUnique = len(self.subsExp.keys())
        self.subsExp = dict(sorted(
            self.subsExp.items(), key=lambda item: item[1], reverse=True)
        )

        # Log substrates
        self.log(f'\nFiltered Substrates:\n'
                 f'     Total: {self.countExpTotal:,}\n'
                 f'    Unique: {self.countExpUnique:,}\n')
        self.log('Substrates:')
        for i, (substrate, count) in enumerate(self.subsExp.items()):
            self.log(f'    {substrate}, {count:,}')
            i += 1
            if i >= self.printN:
                break

        # Save data
        if saveData:
            self.saveSubstrates(substrates=self.subsExp,
                                datasetType=self.datasetTypes['Exp'])

        # Count AAs
        self.countsExp = self.countAA(
            substrates=self.subsExp, countMatrix=self.countsExp,
            datasetType=self.datasetTypes['Exp'], saveData=saveData
        )
        self.calculateRF()
        if saveData:
            self.calculateEntropy(plotFig=plotEntropy)


    def filterMotifs(self, form):
        self.jobDone = False
        self.evalEnrichment()
        self.log('\n\n================================ Filter Motif '
                 '================================')
        self.minS = float(form['minS'])
        self.log(f'Minimum ∆S: {self.minS}')
        self.minES = float(form['minES'])
        self.log(f'Minimum ES: {self.minES}')
        self.minESRel = float(form['minESRel'])
        self.log(f'Minimum ES Release: {self.minESRel}')
        self.selectMotifPos()  ##
        self.log('\nRecognition Sites:')
        self.log(pd.DataFrame.from_dict(self.motifPos, orient='index', columns=['∆S']))

        def evalAAs(position, minES):
            self.fixAA[position] = []
            for aa in self.eMap.index:
                if self.eMap.loc[aa, position] >= minES:
                    self.fixAA[position].append(aa)
            # print(f'Filter (minES={minES}): {position} - {self.fixAA[position]}')


        # Apply Filter
        for pos in self.motifPos.keys():
            if pos not in self.fixAA.keys():
                evalAAs(pos, self.minES)
                self.filterSubs()
                self.evalEnrichment()


        # Refine Filter
        for posRel in self.motifPos.keys():
            self.fixAA = {}
            filter = []
            for pos in self.motifPos.keys():
                if pos != posRel:
                    filter.append(pos)

            # Release
            for posFix in filter:
                evalAAs(posFix, self.minESRel)
            self.getDatasetTag()
            self.filterSubs()
            self.evalEnrichment()

            # Refilter
            evalAAs(posRel, self.minESRel)
            self.getDatasetTag()
            self.filterSubs()
            self.evalEnrichment()
        # print(f'* Fixing AA: {self.datasetTag}')
        # for k, v in self.fixAA.items():
        #     print(f'    {k}: {v}')


        # Release filter
        self.fixAA = {}
        eMap = self.eMap.copy()
        self.substrateProfile = pd.DataFrame(0.0, index=self.countsBg.index,
                                             columns=self.countsBg.columns)
        motifFilter = self.fixAA
        print(f'Filter:\n* {motifFilter}')
        print('Substrate Profile:')
        for posRel in self.motifPos.keys():
            print(f'* Release: {posRel}')
            filter = []
            for pos in self.motifPos.keys():
                if pos != posRel:
                    filter.append(pos)
                for posFix in filter:
                    evalAAs(posFix, self.minESRel)

            self.filterSubs(saveData=False)
            self.substrateProfile.loc[:, posRel] = self.eMap.loc[:, posRel]
        for pos in eMap.columns:
            if pos not in self.motifPos.keys():
                self.substrateProfile.loc[:, pos] = eMap.loc[:, pos]
        self.figures['subProfile'] = self.plotEnrichmentScores(dataType='Enrichment',
                                                               subProfile=True)
        self.jobDone = True


    def combineProfiles(self, form):
        self.jobInit(form, job='Combine Motifs')


    def selectMotifPos(self):
        self.motifPos = {}
        entropy = self.entropy.sort_values(by=self.entropy.columns[0], ascending=False)
        for pos in entropy.index:
            S = round(float(self.entropy.loc[pos, self.entropy.columns[0]]), 2)
            if S >= self.minS:
                self.motifPos[pos] = S
        # print(f'Motif Pos:\n{self.motifPos}')


    def saveSubstrates(self, substrates, datasetType):
        saveTag = None
        if self.datasetTag is None:
            print(f'Dont save, dataset tag: {self.datasetTag}\n')
            sys.exit()
        if datasetType == self.datasetTypes['Exp']:
            saveTag = self.getSaveTag(self.saveTagExp['subs'])
        elif datasetType == self.datasetTypes['Bg']:
            saveTag = self.getSaveTag(self.saveTagBg['subs'])
        else:
            self.logErrorFn(function='saveSubstrates()',
                            msg=f'Unknown dataset type: {datasetType}')

        # Save the substrates
        path = os.path.join(self.pathData, saveTag)
        self.log(f'\nSaving Substrates:\n     {path}')
        with open(path, 'wb') as file:
            pk.dump(substrates, file)


    def countAA(self, substrates, countMatrix, datasetType, saveData=True):
        self.log('\n\n================================== Count AA '
                 '==================================')
        self.log(f'Dataset: {datasetType}\n')
        countMatrix.loc[:, :] = 0
        totalCounts = pd.DataFrame(0, index=self.xAxisLabel, columns=['Sum'])

        def splitData(data, matrix):
            cores = os.cpu_count()
            size = max(1, int(np.ceil(len(data) / cores)))
            batches = list(batched(data.items(), size))
            print(f'Cores: {cores}, Batches: {len(batches)}, Seq/Batch: {size:,}')
            args = [(dict(batch), matrix.columns.tolist(), matrix.index.tolist())
                    for batch in batches]

            with ProcessPoolExecutor(max_workers=cores) as executor:
                results = list(executor.map(counter, args))
            return sum(results)

        start = time.time()
        countMatrix = splitData(substrates, countMatrix)
        end = time.time()
        runtime = (end - start) / 60
        runtime = round(runtime, 3)
        self.log(f'Counted {self.countExpUnique:,} unique substrates in: {runtime} min')
        print(f'Counted {self.countExpUnique:,} unique substrates in: {runtime} min')
        self.log(f'Counts:\n{countMatrix}')

        for pos in countMatrix.columns:
            counts = sum(countMatrix.loc[:, pos])
            totalCounts.loc[pos, 'Sum'] = counts
        self.log(f'\nCount Totals:\n{totalCounts}')

        # File path
        saveTag = None
        if self.datasetTag is None:
            print(f'Dont save, dataset tag: {self.datasetTag}\n')
            sys.exit()
        if datasetType == self.datasetTypes['Exp']:
            saveTag = self.getSaveTag(self.saveTagExp['counts'])
        elif datasetType == self.datasetTypes['Bg']:
            saveTag = self.getSaveTag(self.saveTagBg['counts'])
        else:
            self.logErrorFn(function='countAA()',
                            msg=f'Unknown dataset type: {datasetType}')

        # Save the counts
        print(f'Saving Counts: {saveData}\n     {saveTag}')
        if saveData:
            path = os.path.join(self.pathData, saveTag)
            # self.log(f'\nSaving Counts:\n     {path}')
            countMatrix.to_csv(path)

        return countMatrix


    def plotCounts(self, countedData, totalCounts, datasetType):
        countedData.index = self.AA

        # Create color map
        cMapCustom = self.createCustomColorMap(colorType='Counts')

        # Set figure title
        title = f'{self.enzymeName}\n{datasetType} Substrates\nN={totalCounts:,}'

        # Define the yLabel
        if self.residueLabelType == 0:
            countedData.index = [residue[0] for residue in self.residues]
        elif self.residueLabelType == 1:
            countedData.index = [residue[1] for residue in self.residues]
        elif self.residueLabelType == 2:
            countedData.index = [residue[2] for residue in self.residues]

        # Define color bar limit
        maxCount = np.max(countedData.values)
        mag = math.floor(math.log10(maxCount)) - 1
        maxNorm = maxCount / (10 ** mag)
        maxRound = (math.ceil(maxNorm)) * (10**mag)
        # print(f'Max: {maxCount}\n'
        #       f'Magnitude: {mag}\n'
        #       f'Norm: {maxNorm}\n'
        #       f'Round: {maxRound}\n')


        # Plot the heatmap with numbers centered inside the squares
        fig, ax = plt.subplots(figsize=self.figSize)
        sns.heatmap(countedData, annot=True, fmt=',d', cmap=cMapCustom,
                    cbar=False, linewidths=self.lineThickness-1, square=False,
                    linecolor='black', center=None, vmax=maxRound,
                    annot_kws={'fontweight': 'bold'}, cbar_kws={'pad': 0.02})
        ax.set_xlabel('Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        fig.tight_layout()
        fig.set_size_inches(self.figSize)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        ax.tick_params(axis='y', labelrotation=0)

        # Set x-ticks
        xTicks = np.arange(len(countedData.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(countedData.columns)

        # Set y-ticks
        yTicks = np.arange(len(countedData.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(countedData.index)

        # Colorbar
        normalize = Normalize(vmin=0, vmax=maxRound)  # Normalize the entropy values
        sm = plt.cm.ScalarMappable(norm=normalize, cmap=cMapCustom)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # File path
        figName = f'counts-{self.enzymeName}-{datasetType}.png'
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving Figure: {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def calculateRF(self):
        self.log('\n\n=============================== Calculate: RF '
                 '================================')
        self.log(f'Filter: {self.datasetTag}\n')
        self.rfExp = pd.DataFrame(
            0.0, index=self.countsExp.index, columns=self.countsExp.columns
        )
        for pos in self.countsExp.columns:
            self.rfExp.loc[:, pos] = self.countsExp[pos] / sum(self.countsExp[pos])
        self.log(f'RF Experimental:\n{self.rfExp}')

        if self.rfBg is None:
            self.rfBg = pd.DataFrame(
                0.0, index=self.countsBg.index, columns=self.countsBg.columns
            )
            for pos in self.countsBg.columns:
                totalCounts = sum(self.countsBg[pos])
                for AA in self.countsBg.index:
                    count = self.countsBg.loc[AA, pos]
                    if count == 0:
                        count = 1
                    self.rfBg.loc[AA, pos] = count / totalCounts
            self.log(f'\nRF Background:\n{self.rfBg}')


    def calculateEntropy(self, plotFig=True):
        self.log('\n\n============================= Calculate: Entropy '
                 '=============================')
        self.log(f'Filter: {self.datasetTag}\n')
        self.entropyMax = np.log2(len(self.rfExp.index))
        for indexColumn in self.rfExp.columns:
            S = 0
            for indexRow, probRatio in self.rfExp.iterrows():
                prob = probRatio[indexColumn]
                if prob == 0:
                    continue
                else:
                    S += -prob * np.log2(prob)
            self.entropy.loc[indexColumn, self.entropy.columns[0]] = self.entropyMax - S
        self.log(f'{self.entropy}\n\nMax Entropy: {self.entropyMax.round(6)}')

        if plotFig:
            self.figures['entropy'] = self.plotEntropy()


    def calculateWebLogo(self):
        # Evaluate weblogo
        self.rfExpScaled = pd.DataFrame(
            0.0, index=self.countsExp.index, columns=self.countsExp.columns
        )
        for pos in self.entropy.index:
            self.rfExpScaled.loc[:, pos] = (
                    self.rfExp.loc[:, pos] *
                    self.entropy.loc[pos, self.entropy.columns[0]]
            )
        self.figures['wLogo'] = self.plotWebLogo()


    def evalEnrichment(self, releasedCounts=False):
        self.log('\n\n======================== Calculate: Enrichment Score '
                 '=========================')
        self.log(f'Enrichment Scores:\n'
                 f'     log₂(RF Experimental / RF Background)\n')

        def evalMatrix(data):
            stacks = pd.DataFrame(0.0, index=data.columns,
                                  columns=['+Stack', '-Stack'])
            for pos in data.columns:
                totalPos = 0
                totalNeg = 0
                for value in data.loc[:, pos]:
                    if value > 0:
                        totalPos += value
                    if value < 0:
                        totalNeg += value
                stacks.loc[pos, '+Stack'] = totalPos
                stacks.loc[pos, '-Stack'] = totalNeg
            self.log(f'Stack Heights:\n{stacks}')


        # Calculate: Enrichment scores
        matrix = pd.DataFrame(0.0, index=self.rfExp.index,
                              columns=self.rfExp.columns)
        if len(self.rfBg.columns) == 1:
            # Eval: ES
            for pos in self.rfExp.columns:
                for AA in self.rfExp.index:
                    rf = self.rfExp.loc[AA, pos]
                    if rf == 0:
                        matrix.loc[AA, pos] = -np.inf
                    else:
                        rfBg = self.rfBg.loc[AA, self.rfBg.columns[0]]
                        matrix.loc[AA, pos] = np.log2(rf / rfBg)
        else:
            if len(self.rfBg.columns) != len(self.rfExp.columns):
                self.log(f'ERROR: The number of columns in the Initial Sort '
                      f'({len(self.rfBg.columns)}) needs to equal to the '
                      f'number of columns in the Final Sort '
                      f'({len(self.rfExp.columns)})\n'
                      f'     Initial: {self.rfBg.columns}\n'
                      f'       Final: {self.rfExp.columns}\n\n')
                sys.exit(1)

            # Eval: ES
            for pos in self.rfExp.columns:
                for AA in self.rfExp.index:
                    rf = self.rfExp.loc[AA, pos]
                    if rf == 0:
                        matrix.loc[AA, pos] = -np.inf
                    else:
                        matrix.loc[AA, pos] = np.log2(rf / self.rfBg.loc[AA, pos])
        self.eMap = matrix
        if releasedCounts:
            self.log(f'Enrichment Score: Released Counts\n'
                  f'{matrix.round(self.roundVal)}\n')
            self.log(f'RF Experimental:\n{self.rfExp}\n\n'
                     f'RF Background:\n{self.rfBg}\n')
        else:
            self.log(f'Enrichment Score: {self.datasetTag}\n'
                  f'{matrix.round(self.roundVal)}\n')

        # Evaluate stack heights
        evalMatrix(matrix.replace([np.inf, -np.inf], 0))


        self.log('\n\n===================== Calculate: Scaled Enrichment Score '
                 '=====================')
        if releasedCounts:
            self.log(f'Scale Enrichment Scores: Released Counts\n'
                  f'     Enrichment Scores * ΔS\n')
        else:
            self.log(f'Scale Enrichment Scores:\n'
                  f'     Enrichment Scores * ΔS\n')

        # Calculate: Letter heights
        heights = pd.DataFrame(0, index=matrix.index,
                               columns=matrix.columns, dtype=float)
        for indexColumn in heights.columns:
            heights.loc[:, indexColumn] = (matrix.loc[:, indexColumn] *
                                           self.entropy.loc[indexColumn,
                                           self.entropy.columns[0]])

        # Calculate: Max positive
        columnTotals = []
        for pos in heights.columns:
            totalPos = 0
            for value in heights.loc[:, pos]:
                if value > 0:
                    totalPos += value
            columnTotals.append(totalPos)
        yMax = max(columnTotals)

        # Adjust values
        for column in heights.columns:
            if heights.loc[:, column].isna().any():
                nValues = heights[column].notna().sum()
                if nValues > 0:
                    self.log(f'{len(heights[column]) - nValues} NaN values at: {column}')
                heights.loc[heights[column].notna(), column] = yMax / nValues
                heights.loc[:, column] = heights.loc[:, column].fillna(0)
        self.eMapScaled = heights.copy()
        self.log(f'Residue Heights: {self.datasetTag}\n'
                 f'{heights}\n')

        # Evaluate stack heights
        evalMatrix(heights.replace([np.inf, -np.inf], 0))

        x = {
            'eMap': False, 'eLogo': False, 'eLogoMin': False, 'eMapSc': False,
            'wLogo': False, 'words': False, 'barCounts': False, 'barRF': False,
            'exp_counts': False, 'bg_counts': False
        }

        # Plot: Enrichment Map
        self.figures['eMap'] = (
            self.plotEnrichmentScores(dataType='Enrichment')
        )
        self.figures['eMapSc'] = (
            self.plotEnrichmentScores(dataType='Scaled Enrichment')
        )
        
        # Plot: Enrichment Logo
        self.plotEnrichmentLogo()

        # Plot: Weblogo
        self.calculateWebLogo()

        # Plot: Wordcloud
        self.figures['words'] = self.plotWordCloud(self.subsExp)

        if self.motifFilter:
            self.iteration += 1


    def plotEntropy(self):
        # Set figure title
        title = self.enzymeName

        # Figure parameters
        yMax = self.entropyMax + 0.2
        xMax = len(self.entropy.iloc[:, 0])

        # Map entropy values to colors using the colormap
        colors = [(0, 'navy'),
                  (0.3/self.entropyMax, 'navy'),
                  (0.7/self.entropyMax, 'dodgerblue'),
                  (0.97/self.entropyMax, 'white'),
                  (0.98/self.entropyMax, 'white'),
                  (1.0/self.entropyMax, 'white'),
                  (1.65/self.entropyMax, 'red'),
                  (3/self.entropyMax, 'firebrick'),
                  (1, 'darkred')]
        colorBar = LinearSegmentedColormap.from_list('custom_colormap', colors)

        # Map entropy values to colors using the colormap
        normalize = Normalize(vmin=0, vmax=yMax) # Normalize the entropy values
        cMap = [colorBar(normalize(value)) for value in self.entropy[self.entropy.columns[0]].astype(float)]

        # Plotting the entropy values as a bar graph
        fig, ax = plt.subplots(figsize=self.figSize)
        if self.motifFilter:
            plt.hlines(y=[self.minS], xmin=-0.5, xmax=xMax, colors=[self.orange], zorder=0)
        plt.bar(self.entropy.index, self.entropy[self.entropy.columns[0]], color=cMap,
                edgecolor='black', linewidth=self.lineThickness, width=0.8)
        plt.xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        plt.ylabel(self.entropy.columns[0], fontsize=self.labelSizeAxis, rotation=0, labelpad=15)
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.tight_layout()

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks)

        # Set x-ticks
        xTicks = np.arange(0, xMax, 1)
        ax.set_xticks(xTicks)
        ax.set_xticklabels(self.entropy.index, rotation=0, ha='center')
        for tick in ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Set y-ticks
        yTicks = range(0, 5)
        yTickLabels = [f'{tick:.0f}' if tick != 0 else f'{int(tick)}' for tick in yTicks]
        ax.set_yticks(yTicks)
        ax.set_yticklabels(yTickLabels)
        for tick in ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness) # Set tick width

        # Set the edge thickness
        for spine in ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # Set axis limits
        ax.set_xlim(-0.5, xMax-0.5)
        ax.set_ylim(0, yMax)

        # Colorbar
        sm = plt.cm.ScalarMappable(norm=normalize, cmap=colorBar)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02)
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # File path
        figName = f'entropy-{self.enzymeName}-{self.getSaveTag()}-{self.motifLen}AA.png'
        if self.motifFilter:
            figName = figName.replace('entropy', 'entropyMin')
        path = os.path.join(self.pathFigs, figName)

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotEnrichmentScores(self, dataType, subProfile=False):
        # Select: Dataset
        scaleData = False
        if subProfile:
            scores = self.substrateProfile
        else:
            if 'scaled' in dataType.lower():
                scaleData = True
                scores = self.eMapScaled
            else:
                scores = self.eMap

        # Define: Figure title
        title = f'{self.enzymeName}'

        # Create heatmap
        cMapCustom = self.createCustomColorMap(colorType='EM')

        # Define the yLabel
        if self.residueLabelType == 0:
            scores.index = [residue[0] for residue in self.residues]
        elif self.residueLabelType == 1:
            scores.index = [residue[1] for residue in self.residues]
        elif self.residueLabelType == 2:
            scores.index = [residue[2] for residue in self.residues]

        # Define color bar limits
        if np.max(scores) >= np.min(scores):
            cBarMax = np.ceil(np.max(scores) * 10) / 10
            cBarMin = -1 * cBarMax
        else:
            cBarMin = np.floor(np.min(scores) * 10) / 10
            cBarMax = -1 * cBarMin

        # Plot the heatmap with numbers centered inside the squares
        if self.figEMSquares:
            fig, ax = plt.subplots(figsize=self.figSizeSq)
            heatmap = sns.heatmap(
                scores, annot=False, cmap=cMapCustom, cbar=True,
                linewidths=self.lineThickness - 1, linecolor='black',
                square=True, center=None, vmax=cBarMax, vmin=cBarMin,
                cbar_kws={'pad': 0.02}
            )
        else:
            fig, ax = plt.subplots(figsize=self.figSize)
            heatmap = sns.heatmap(
                scores, annot=True, fmt='.3f', cmap=cMapCustom, cbar=True,
                linewidths=self.lineThickness - 1, linecolor='black',
                square=False, center=None, vmax=cBarMax, vmin=cBarMin,
                annot_kws={'fontweight': 'bold'}, cbar_kws={'pad': 0.02}
            )
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        ax.set_xlabel('Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        fig.tight_layout()
        fig.set_size_inches(self.figSize)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', rotation=0, length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set x-ticks
        xTicks = np.arange(len(scores.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(scores.columns)

        # Set y-ticks
        yTicks = np.arange(len(scores.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(scores.index)

        # Set invalid values to grey
        cmap = plt.cm.get_cmap(cMapCustom)
        cmap.set_bad(color='lightgrey')

        # Colorbar
        normalize = Normalize(vmax=cBarMax, vmin=cBarMin)
        sm = plt.cm.ScalarMappable(norm=normalize, cmap=cMapCustom)
        sm.set_array([])
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        for tick in cbar.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # File path
        figName = f'eMap-{self.enzymeName}-{self.getSaveTag()}-{self.motifLen}AA.png'
        if self.motifFilter and not subProfile:
            figName = figName.replace('eMap', f'eMap_{self.iteration}')
        if scaleData:
            figName = figName.replace('eMap', 'eMap_Scaled')
        if subProfile:
            figName = figName.replace('eMap', 'eMap_SubProfile')
        path = os.path.join(self.pathFigs, figName)
        self.log(f'\nSaving EM at:\n   {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName


    def plotEnrichmentLogo(self):
        # Define: Figure title
        title = f'{self.enzymeName}'

        # Set parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'

        # Rename columns for logomaker script
        data = self.eMapScaled.copy().replace([np.inf, -np.inf], 0)
        xticks = data.columns
        data.columns = range(len(data.columns))

        # Calculate: Max and min
        columnTotals = [[], []]
        for indexColumn in data.columns:
            totalPos = 0
            totalNeg = 0
            for value in data.loc[:, indexColumn]:
                if value > 0:
                    totalPos += value
                elif value < 0:
                    totalNeg += value
            columnTotals[0].append(totalPos)
            columnTotals[1].append(totalNeg)
        yMax = max(columnTotals[0])
        yMin = min(columnTotals[1])


        def plotLogo(matrix, limitYAxis=False):
            # Plot the sequence motif
            fig, ax = plt.subplots(figsize=self.figSize)
            motif = logomaker.Logo(matrix.transpose(), ax=ax, color_scheme=self.colorsAA,
                                   width=0.95, stack_order=stackOrder)
            motif.ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
            fig.tight_layout()
            fig.set_size_inches(self.figSize)

            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks)

            # Set borders
            motif.style_spines(visible=False)
            motif.style_spines(spines=['left', 'bottom'], visible=True)
            for spine in motif.ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            # Set x-ticks
            motif.ax.set_xticks([pos for pos in range(len(xticks))])
            motif.ax.set_xticklabels(xticks, fontsize=self.labelSizeTicks,
                                     rotation=0, ha='center')

            # Set y-ticks
            yTicks = [yMin, 0, yMax]
            yTickLabels = [f'{tick:.2f}' if tick != 0 else f'{int(tick)}'
                           for tick in yTicks]
            motif.ax.set_yticks(yTicks)
            motif.ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
            motif.ax.set_ylim(yMin, yMax)

            # Set tick width
            for tick in motif.ax.xaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness)
            for tick in motif.ax.yaxis.get_major_ticks():
                tick.tick1line.set_markeredgewidth(self.lineThickness)

            # Label the axes
            motif.ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
            motif.ax.set_ylabel('Scaled Enrichment', fontsize=self.labelSizeAxis)

            # Set horizontal line
            motif.ax.axhline(y=0, color='black', linestyle='-',
                             linewidth=self.lineThickness)

            # File path
            if self.datasetTag is None:
                print(f'Dont save, dataset tag: {self.datasetTag}\n')
                sys.exit()
            figName = f'eLogo-{self.enzymeName}-{self.getSaveTag()}-{self.motifLen}AA.png'
            if self.motifFilter:
                figName = figName.replace('eLogo', f'eLogo_{self.iteration}')
            if limitYAxis:
                figName = figName.replace('eLogo-', 'eLogo_yMin-')
            path = os.path.join(self.pathFigs, figName)
            # self.log(f'\nSaving Enrichment Logo:\n     {path}')

            # Encode the figure
            figBase64 = self.encodeFig(fig)
            with open(path, "wb") as file:
                file.write(base64.b64decode(figBase64))

            # Close the figure to free memory
            plt.close(fig)

            return figName

        # Plot figure
        self.figures['eLogo'] = plotLogo(data) # Full y-axis

        # Adjust yMin to fit the largest negative AA
        yMin = 0
        for col in data.columns:
            for row in data.index:
                if data.loc[row, col] < yMin:
                    yMin = data.loc[row, col]
        self.figures['eLogoMin'] = plotLogo(data, limitYAxis=True) # Limited y-axis


    def plotWordCloud(self, substrates):
        # Limit the number of words
        subs = {}
        iteration = 0
        for substrate, count in substrates.items():
            subs[substrate] = count
            iteration += 1
            if iteration >= self.numSamples:
                break
        substrates = subs
        totalWords = len(substrates)

        # Define: Figure title
        title = self.enzymeName

        # Create word cloud
        cmap = self.createCustomColorMap(colorType='Word Cloud')
        wordcloud = (WordCloud(
            width=950,
            height=800,
            background_color='white',
            min_font_size=10, # Minimum font size
            max_font_size=120, # Maximum font size
            scale=5,  # Increase scale for larger words
            colormap=cmap
        ).generate_from_frequencies(substrates))

        # Display the word cloud
        fig = plt.figure(figsize=self.figSize, facecolor='white')
        plt.imshow(wordcloud, interpolation='bilinear')
        plt.title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axis('off')
        fig.tight_layout()
        fig.set_size_inches(self.figSize)

        # File path
        if self.datasetTag is None:
            print(f'Dont save, dataset tag: {self.datasetTag}\n')
            sys.exit()
        figName = (f'wordcloud-{self.enzymeName}-{self.getSaveTag()}-{self.motifLen}AA'
                   f'{totalWords}_Words.png')
        if self.motifFilter:
            figName = figName.replace('wordcloud',
                                      f'wordcloud_{self.iteration}')
        path = os.path.join(self.pathFigs, figName)
        # self.log(f'\nSaving Wordcloud:\n     {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName

    def plotWebLogo(self):
        # Define: Figure title
        title = f'{self.enzymeName}'

        # Set parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'

        # Rename columns for logomaker script
        data = self.rfExpScaled.copy().replace([np.inf, -np.inf], 0)
        xticks = data.columns
        data.columns = range(len(data.columns))

        # Plot the sequence motif
        fig, ax = plt.subplots(figsize=self.figSize)
        motif = logomaker.Logo(data.transpose(), ax=ax, color_scheme=self.colorsAA,
                               width=0.95, stack_order=stackOrder)
        motif.ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        fig.tight_layout()
        fig.set_size_inches(self.figSize)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks)

        # Set borders
        motif.style_spines(visible=False)
        motif.style_spines(spines=['left', 'bottom'], visible=True)
        for spine in motif.ax.spines.values():
            spine.set_linewidth(self.lineThickness)

        # Set x-ticks
        motif.ax.set_xticks([pos for pos in range(len(xticks))])
        motif.ax.set_xticklabels(xticks, fontsize=self.labelSizeTicks,
                                 rotation=0, ha='center')

        # Set y-ticks
        yMax = self.entropyMax
        yTicks = range(0, 5)
        yTickLabels = [f'{tick:.0f}' if tick != yMax else f'{yMax:.2f}' for tick in
                       yTicks]
        # yTicks.append(4.32)
        # yTickLabels.append('')
        motif.ax.set_yticks(yTicks)
        motif.ax.set_yticklabels(yTickLabels, fontsize=self.labelSizeTicks)
        motif.ax.set_ylim(0, yMax)

        # Set tick width
        for tick in motif.ax.xaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)
        for tick in motif.ax.yaxis.get_major_ticks():
            tick.tick1line.set_markeredgewidth(self.lineThickness)

        # Label the axes
        motif.ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        motif.ax.set_ylabel('Bits', fontsize=self.labelSizeAxis)

        # Set horizontal line
        motif.ax.axhline(y=0, color='black', linestyle='-',
                         linewidth=self.lineThickness)

        # File path
        if self.datasetTag is None:
            print(f'Dont save, dataset tag: {self.datasetTag}\n')
            sys.exit()
        figName = f'webLogo-{self.enzymeName}-{self.getSaveTag()}-{self.motifLen}AA.png'
        if self.motifFilter:
            figName = figName.replace('webLogo', f'webLogo_{self.iteration}')
        path = os.path.join(self.pathFigs, figName)
        # self.log(f'\nSaving WebLogo:\n     {path}')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName
