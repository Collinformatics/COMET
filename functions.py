import base64
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import gzip
import hashlib
import io
# import logomaker
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
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
import warnings
# from wordcloud import WordCloud


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
        # Params: Dataset
        self.datasetTag = 'Unfiltered'
        self.datasetTagMotif = None
        self.motifFilter = None
        self.jobParams = {}
        self.jobID = ''
        self.jobHash = ''
        self.enzymeName = ''
        self.seqLength = False
        self.minCounts = 1
        self.printN = 10
        self.roundVal = 3
        self.xAxisLabel = False
        self.entropy = None
        self.subsExp = {}
        self.countsExp = 'Initialize me'
        self.countExpTotal = 0
        self.countExpUnique = 0
        self.rfExp = None
        self.saveTagExp = {}
        self.subsBg = {}
        self.countsBg = 'Initialize me'
        self.countBgTotal = 0
        self.countBgUnique = 0
        self.rfBg = None
        self.eMap = None
        self.eMapScaled = None
        self.eMapReleased = None
        self.eMapReleasedScaled = None
        self.saveTagBg = {}
        self.saveTagFig = {}
        self.fixMotif = False
        self.datasetTypes = {'Exp': 'Experimental', 'Bg': 'Background'}

        # Params: Files
        self.figures = {}
        self.queueLog = queue.Queue()
        self.fileExp = []
        self.seqExp = None
        self.fileBg = []
        self.seqBg = None
        self.pathData = ''
        self.pathSeqs = ''
        self.pathFigs = ''
        self.pathLog = ''

        # Params: Process dna
        self.seq5Prime = False
        self.seq3Prime = False
        self.minPhred = False

        # Params: Filter Dataset
        self.filterPos = False
        self.fixAA = {}
        self.exclAA = {}

        # Params: Figures
        self.figureResolution = 600
        self.saveFigureIteration = None
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

        # # Params:
        # self. = False
        # self. = False
        # self. = False

        pd.options.display.float_format = '{:,.3f}'.format



    @staticmethod
    def createCustomColorMap(colorType):
        colorType = colorType.lower()
        if colorType == 'counts':
            useGreen = True
            if useGreen:
                # Green
                colors = ['#FFFFFF', '#ABFF9B', '#39FF14', '#2E9418', '#2E9418',
                          '#005000']
            else:
                # Orange
                colors = ['white', 'white', '#FF76FA', '#FF50F9', '#FF00F2',
                          '#CA00DF', '#BD16FF']
        elif colorType == 'stdev':
            colors = ['white', 'white', '#FF76FA', '#FF50F9', '#FF00F2', '#CA00DF',
                      '#BD16FF']
        elif colorType == 'word cloud':
            # ,'#F2A900','#2E8B57','black'
            colors = ['#CC5500', '#CC5500', '#F79620', '#FAA338',
                      '#00C01E', '#1D680D', '#003000', 'black']
            # colors = ['#008631','#39E75F','#CC5500','#F79620','black']
        elif colorType == 'em':
            colors = ['navy', 'royalblue', 'dodgerblue', 'lightskyblue', 'white', 'white',
                      'lightcoral', 'red', 'firebrick', 'darkred']
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
    def hashStr(string):
        return hashlib.sha256(string.encode('utf-8')).hexdigest()


    @staticmethod
    def pressButton(self, message):
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



    def filterAA(self):
        self.log()  # Clear the log

        print(f'Code: filterAA()')
        sys.exit()

        return {'AA': 'VEHTVALKQNR'}



    def filterMotif(self):
        self.log()  # Clear the log

        print(f'Code: filterAA()')
        sys.exit()

        return {'Motif': 'TVALK'}



    def getDatasetTag(self):
        tagFix = 'Fix '
        tagExcl = 'Excl '
        fixNPos = len(self.fixAA)

        if self.exclAA:
            print('Exclude AA:')
            for index, (pos, AA) in enumerate(self.exclAA.items()):
                if len(AA) > 1:
                    tag = f'[{','.join(AA)}]@{pos.replace('fix', '')}'
                else:
                    tag = f'{AA}@{pos.replace('fix', '')}'

                if fixNPos > 1 and index != fixNPos - 1:
                    tagExcl += f'{tag}_'
                else:
                    tagExcl += tag
                print(f'Tag: {tagExcl}')
            self.datasetTag = tagExcl

        if self.fixAA:
            print('Fix AA')
            for index, (pos, AA) in enumerate(self.fixAA.items()):
                print(index, pos, AA)
                if len(AA) > 1:
                    tag = f'[{','.join(AA)}]@{pos.replace('fix', '')}'
                else:
                    tag = f'{AA}@{pos.replace('fix', '')}'

                if fixNPos > 1 and index != fixNPos - 1:
                    tagFix += f'{tag}_'
                else:
                    tagFix += tag
                print(f'Tag: {tagFix}')
            self.datasetTag = tagFix
        self.jobParams['Dataset Tag'] = self.datasetTag
        self.log(f'Dataset Filter: {self.datasetTag}')

        # Initialize: Save tags
        self.saveTagExp = {
            'subsRaw': f'{self.enzymeName} - Subs Exp - {self.datasetTag} - '
                    f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'countsRaw': f'{self.enzymeName} - AA Counts Exp - {self.datasetTag} - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'subs': f'{self.enzymeName} - Subs Exp - {self.datasetTag} - '
                    f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'counts': f'{self.enzymeName} - AA Counts Exp - {self.datasetTag} - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA'
        }
        self.saveTagBg = {
            'subsRaw': f'{self.enzymeName} - Subs Bg - {self.datasetTag} - '
                   f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'countsRaw': f'{self.enzymeName} - AA Counts Bg - {self.datasetTag} - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'subs': f'{self.enzymeName} - Subs Bg - {self.datasetTag} - '
                    f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'counts': f'{self.enzymeName} - AA Counts Bg - {self.datasetTag} - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA',
        }
        if self.fixMotif:
            for tag, path in self.saveTagExp:
                self.saveTagExp[tag] = (path.replace(f'{self.enzymeName}',
                                                     f'{self.enzymeName} - Motif'))

        self.saveTagFig = (f'{self.enzymeName} - Fig - {self.datasetTag} - '
                           f'Min Counts {self.minCounts} - {self.seqLength} AA')



    def initDataStructures(self):
        # Initialize data structures
        self.subsExp = {}
        self.subsBg = {}
        self.xAxisLabel = [f'R{index}' for index in range(1, self.seqLength + 1)]
        self.countsExp = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.countsBg = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)



    def jobInit(self, form, evalDNA=False, fixAA=False, fixMotif=False):
        # Directories
        self.pathData = os.path.join('data', form['enzymeName'])
        self.pathSeqs = os.path.join(self.pathData, 'sequences')
        self.pathFigs = os.path.join(self.pathData, 'figures')
        self.pathLog = os.path.join(self.pathData, 'log.txt')
        if self.pathData is not None:
            if not os.path.exists(self.pathData):
                os.makedirs(self.pathData, exist_ok=True)
        if self.pathSeqs is not None:
            if not os.path.exists(self.pathSeqs):
                os.makedirs(self.pathSeqs, exist_ok=True)
        if self.pathFigs is not None:
            if not os.path.exists(self.pathFigs):
                os.makedirs(self.pathFigs, exist_ok=True)
            else:
                import shutil
                # Remove everything inside the directory
                for filename in os.listdir(self.pathFigs):
                    file_path = os.path.join(self.pathFigs, filename)
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)  # delete file or link
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)  # delete subdirectory
                time.sleep(5)
                os.makedirs(self.pathFigs, exist_ok=True)

        print('Job Parameters:')
        for key, value in form.items():
            print(f'     {key}: {value}')
        print()

        # Get the files
        for key, value in form.items():
            if 'fileExp' in key:
                self.fileExp.append(value)
            elif 'fileBg' in key:
                self.fileBg.append(value)

        self.log()  # Clear the log
        self.log('================================ Job Summary '
                 '=================================')

        # Record job params
        self.enzymeName = form['enzymeName']
        self.jobParams['Enzyme Name'] = self.enzymeName
        self.log(f'Enzyme: {self.enzymeName}')
        self.seqLength = int(form['seqLength'])
        self.jobParams['Substrate Length'] = self.seqLength
        self.log(f'Substrate Length: {self.seqLength}')
        self.jobID = f'{self.enzymeName}_{self.seqLength}'

        ## Placeholder for files
        self.fileExp = ['data/validation/variantsExp.fastq'] # , 'data/validation/variantsExp2.fastq'
        self.fileBg = ['data/validation/variantsBg.fasta'] # , 'data/validation/variantsBg2.fasta'
        print(f'File Exp: {type(self.fileExp)}\n'
              f'{self.fileExp}\n')
        print(f'File Bg: {type(self.fileBg)}\n'
              f'{self.fileBg}\n')

        # Job dependant parameters
        if evalDNA:
            self.log('Job: Evaluate DNA')
            self.seq5Prime = form['seq5Prime']
            self.seq3Prime = form['seq3Prime']
            self.minPhred = int(form['minPhred']) if form['minPhred'] != '' else 0
            self.log(f'5\' Sequence: {self.seq5Prime}\n'
                     f'3\' Sequence: {self.seq3Prime}\n'
                     f'Min Phred Score: {self.minPhred}')
            self.jobID += f'{self.seq5Prime}_{self.seq3Prime}_Phred-{self.minPhred}'
        elif fixAA:
            print('Fix AA')
        elif fixMotif:
            print('Fix Motif')
        else:
            print('ERROR: What Script Is Running')
            sys.exit()

        # Complete initialization
        self.getFilter(form)
        self.initDataStructures()
        self.jobID += (f'_{self.datasetTag.replace(' ','-')}'
                          f'_Exp-{"-".join(self.fileExp)}'
                          f'_Bg-{"-".join(self.fileBg)}'
                          f'_{time.ctime().replace(' ', '-')}')
        self.jobHash = self.hashStr(self.jobID)
        # print(f'Job:\n'
        #       f'* Label: {self.jobID}\n'
        #       f'* Hash: {self.jobHash}\n')



    def getFilter(self, data):
        # print('Filter:')
        # for key, value in data.items():
        #     print(f'  {key}: {value}')
        # print()

        if 'filterPos' in data.keys():
            self.filterPos = data['filterPos']
            self.fixAA = {}
            self.exclAA = {}

            # Get filter params
            print('Filter:')
            for key, value in data.items():
                if 'fix' in key:
                    self.fixAA[key] = value
                    print(f'  {key}: {value}')
                if 'excl' in key:
                    self.exclAA[key] = value
                    print(f'  {key}: {value}')

            # Sort filter params
            self.fixAA = dict(sorted(self.fixAA.items()))
            self.exclAA = dict(sorted(self.exclAA.items()))
        self.getDatasetTag()



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
                     f'========================================\n\n'
                     f'ERROR: {function}\n'
                     f'{msg}\n\n'
                     f'========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n')
        else:
            self.log(f'\n========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n\n'
                     f'ERROR: {function}\n'
                     f'{msg}\n\n'
                     f'========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n')
            sys.exit(1)



    def processSubs(self, substrates, datasetType, filteredAA):
        self.log('\n\n================================= Substrates '
                 '=================================')
        self.log(f'Dataset: {datasetType}\n')

        # Inspect sequences
        if not filteredAA:
            filteredSubs = {}
            for substrate, count in substrates.items():
                for AA in substrate:
                    if AA not in self.AA:
                        filteredSubs[substrate] = count
            if filteredSubs:
                self.log(f'Filtering Substrates:\n'
                         f'     If a substrate contains an '
                         f'unaccented AA it will be removed.\n'
                         f'     Accepted: {self.AA}\n\n'
                         f'     Removed Substrates:')
                for substrate, count in filteredSubs.items():
                    substrates.pop(substrate, count)
                    self.log(f'          {substrate}: {count}')
                self.log('')

        # Sort data
        substrates = dict(sorted(substrates.items(),
                                 key=lambda item: item[1], reverse=True))

        # Count AAs
        countMatrix = None
        self.log('Substrate Totals:')
        if datasetType == self.datasetTypes['Exp']:
            self.subsExp = substrates
            countMatrix = self.countsExp
            self.countExpTotal = sum(substrates.values())
            self.countExpUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countExpTotal:,}\n'
                     f'    Unique Substrates: {self.countExpUnique:,}\n')

            # Record job params
            self.jobParams['Total Experimental Substrates'] = f'{self.countExpTotal:,}'
            self.jobParams['Unique Experimental Substrates'] = f'{self.countExpUnique:,}'
        elif datasetType == self.datasetTypes['Bg']:
            self.subsBg = substrates
            countMatrix = self.countsBg
            self.countBgTotal = sum(substrates.values())
            self.countBgUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countBgTotal:,}\n'
                     f'    Unique Substrates: {self.countBgUnique:,}\n')

            # Record job params
            self.jobParams['Total Background Substrates'] = f'{self.countBgTotal:,}'
            self.jobParams['Unique Background Substrates'] = f'{self.countBgUnique:,}'
        else:
            self.logErrorFn(function='sampleSize()',
                            msg=f'Unknown dataset type: {datasetType}')

        self.log(f'Top {self.printN:,} {datasetType} Sequences')
        for index, (sub, count) in enumerate(substrates.items()):
            if index >= self.printN:
                break
            self.log(f'     {sub}: {count}')
        self.log('')

        # Save data
        self.saveSubstrates(substrates=substrates,
                            datasetType=datasetType,
                            filteredAA=filteredAA)

        # Count AAs
        self.countAA(substrates=substrates, countMatrix=countMatrix,
                     datasetType=datasetType, filteredAA=filteredAA)



    def evalDNA(self, form):
        self.jobInit(form, evalDNA=True)

        # Load the data
        threads = []
        queuesExp = []
        queuesExpLog = []
        queuesBg = []
        queuesBgLog = []
        if self.fileExp:
            for file in self.fileExp:
                queueExp = queue.Queue()
                queueLog = queue.Queue()
                queuesExp.append(queueExp)
                queuesExpLog.append(queueLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Exp'], queueExp, queueLog, True,))
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
                    args=(file, self.datasetTypes['Bg'], queueBg, queueLog, True,))
                thread.start()
                threads.append(thread)

        # Wait for all threads to finish
        for thread in threads:
            thread.join()

        # Log the output
        if queuesExpLog:
            for log in queuesExpLog:
                self.logInQueue(log)
        if queuesBgLog:
            for log in queuesBgLog:
                self.logInQueue(log)

        # Get results from queue
        if self.fileExp:
            for queueData in queuesExp:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in self.subsExp.keys():
                        self.subsExp[substrate] += count
                    else:
                        self.subsExp[substrate] = count
        if self.fileBg:
            for queueData in queuesBg:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in self.subsBg.keys():
                        self.subsBg[substrate] += count
                    else:
                        self.subsBg[substrate] = count

        # Make figures
        self.figures = {
            'eMap': False, 'eMapSc': False, 'eLogo': False, 'wLogo': False,
            'words': False, 'barCounts': False, 'barRF': False,
            'exp_counts': False, 'bg_counts': False
        }
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
            self.calculateEnrichment()



    def loadDNA(self, path, datasetType, queueData, queueLog, forwardRead):
        translate = True

        # Open the file
        openFn = gzip.open if path.endswith('.gz') else open  # Define open function
        with openFn(path, 'rt') as file:  # 'rt' = read text mode
            if '.fastq' in path or '.fq' in path:
                data = SeqIO.parse(file, 'fastq')
                warnings.simplefilter('ignore', BiopythonWarning)
            elif '.fasta' in path or '.fa' in path:
                data = SeqIO.parse(file, 'fasta')
                warnings.simplefilter('ignore', BiopythonWarning)
            else:
                queueLog.put(self.logErrorFn(
                    function='loadDNA()',
                    msg=f'Unrecognized file\n     {path}',
                    getStr=True))
                translate = False

            # Translate the dna
            if translate:
                substrates = self.translate(data, path, datasetType,
                                            queueLog, forwardRead)
                queueData.put(substrates) # Put the substrates in the queue



    def translate(self, data, fileName, datasetType, queueLog, forwardRead):
        queueLog.put('\n\n================================ Translate DNA '
                     '===============================')
        data = list(data)
        substrates = {}
        totalSubsExtracted = 0
        totalSeqsDNA = 0
        self.printN = 10

        queueLog.put(f'File Name: {fileName}')
        if forwardRead:
            queueLog.put('Read Type: Forward Read')
        else:
            queueLog.put('Read Type: Reverse Read')
        queueLog.put(f'  Dataset: {datasetType}')

        # Inspect the file
        useQS = False
        for datapoint in data:
            if 'phred_quality' in datapoint.letter_annotations:
                useQS = True
            break
        queueLog.put(f'  Eval QS: {useQS}\n\n')

        # Inspect the datasetType parameter
        if (datasetType != self.datasetTypes['Exp']
                and datasetType != self.datasetTypes['Bg']):
            self.logErrorFn(function='translate()',
                            msg=f'Unknown dataset type: {datasetType}')


        def extractionEfficiency(fullSet=False):
            perExtracted = (totalSubsExtracted / totalSeqsDNA) * 100
            if fullSet:
                queueLog.put('- All Sequences')
            else:
                queueLog.put(f'\nExtraction Efficiency: {datasetType}\n'
                         f'- First {self.printN} Sequences')
            queueLog.put(f'     Evaluated DNA Sequences: {totalSeqsDNA:,}\n'
                         f'        Extracted Substrates: {totalSubsExtracted:,}\n'
                         f'       Extraction Efficiency: {round(perExtracted, 3)} %')


        # Translate DNA - Sample Set
        if useQS:
            for index, datapoint in enumerate(data):
                if totalSubsExtracted >= self.printN: # Exit loop
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                queueLog.put(f'DNA Seq: {dna}')


                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:
                    qs = datapoint.letter_annotations['phred_quality']

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    queueLog.put(f'Sub DNA: {substrateDNA}')
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        queueLog.put(f'Sub Seq: {substrate}')

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            qs = qs[start:end]
                            queueLog.put(f'     QS: {qs}')
                            if all(score >= self.minPhred for score in qs):
                                queueLog.put(f'Keep Substrate\n')
                                if substrate in substrates.keys():
                                    substrates[substrate] += 1
                                else:
                                    substrates[substrate] = 1
                                totalSubsExtracted += 1
                            else:
                                queueLog.put('')
        else:
            for index, datapoint in enumerate(data):

                if totalSubsExtracted >= self.printN:
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                queueLog.put(f'DNA Seq: {dna}')

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    queueLog.put(f'Sub DNA: {substrateDNA}')

                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        queueLog.put(f'Sub Seq: {substrate}')

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            queueLog.put(f'Keep Substrate\n')
                            if substrate in substrates.keys():
                                substrates[substrate] += 1
                            else:
                                substrates[substrate] = 1
                            totalSubsExtracted += 1
                        else:
                            queueLog.put('')
        extractionEfficiency() # Evaluate data quality


        # Translate DNA - Full Set
        substrates = {}
        totalSeqsDNA = 0
        totalSubsExtracted = 0
        if useQS:
            for index, datapoint in enumerate(data):
                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:
                    qs = datapoint.letter_annotations['phred_quality']

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            qs = qs[start:end]
                            if all(score >= self.minPhred for score in qs):
                                if substrate in substrates.keys():
                                    substrates[substrate] += 1
                                else:
                                    substrates[substrate] = 1
                                totalSubsExtracted += 1
        else:
            for index, datapoint in enumerate(data):
                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            if substrate in substrates.keys():
                                substrates[substrate] += 1
                            else:
                                substrates[substrate] = 1
                            totalSubsExtracted += 1
        extractionEfficiency(fullSet=True)  # Evaluate data quality

        return substrates



    def saveSubstrates(self, substrates, datasetType, filteredAA):
        saveTag = None
        if filteredAA:
            if self.datasetTag is None:
                print(f'Dont save, dataset tag: {self.datasetTag}\n')
                sys.exit()
            if datasetType == self.datasetTypes['Exp']:
                saveTag = self.saveTagExp['subs']
            elif datasetType == self.datasetTypes['Bg']:
                saveTag = self.saveTagBg['subs']
            else:
                self.logErrorFn(function='saveSubstrates()',
                                msg=f'Unknown dataset type: {datasetType}')
        else:
            if datasetType == self.datasetTypes['Exp']:
                saveTag = self.saveTagExp['subsRaw']
            elif datasetType == self.datasetTypes['Bg']:
                saveTag = self.saveTagBg['subsRaw']
            else:
                self.logErrorFn(function='saveSubstrates()',
                                msg=f'Unknown dataset type: {datasetType}')

        # Save the substrates
        path = os.path.join(self.pathSeqs, saveTag)
        self.log(f'Saving Substrates:\n     {path}')
        with open(path, 'wb') as file:
            pk.dump(substrates, file)



    def countAA(self, substrates, countMatrix, datasetType, filteredAA):
        self.log('\n\n================================== Count AA '
                 '==================================')
        self.log(f'Dataset: {datasetType}\n')
        totalCounts = pd.DataFrame(0, index=self.xAxisLabel, columns=['Sum'])
        for substrate, count in substrates.items():
            for index, AA in enumerate(substrate):
                countMatrix.loc[AA, self.xAxisLabel[index]] += count
        for pos in countMatrix.columns:
            counts = sum(countMatrix.loc[:, pos])
            totalCounts.loc[pos, 'Sum'] = counts
        self.log(f'Counts:\n{countMatrix}\n\n{totalCounts}\n')

        # File path ##
        saveTag = None
        if filteredAA:
            if self.datasetTag is None:
                print(f'Dont save, dataset tag: {self.datasetTag}\n')
                sys.exit()
            if datasetType == self.datasetTypes['Exp']:
                saveTag = self.saveTagExp['counts']
            elif datasetType == self.datasetTypes['Bg']:
                saveTag = self.saveTagBg['counts']
            else:
                self.logErrorFn(function='countAA()',
                                msg=f'Unknown dataset type: {datasetType}')
        else:
            if datasetType == self.datasetTypes['Exp']:
                saveTag = self.saveTagExp['countsRaw']
            elif datasetType == self.datasetTypes['Bg']:
                saveTag = self.saveTagBg['countsRaw']
            else:
                self.logErrorFn(function='countAA()',
                                msg=f'Unknown dataset type: {datasetType}')

        # Save the counts
        path = os.path.join(self.pathSeqs, saveTag)
        self.log(f'Saving Counts:\n     {path}')
        countMatrix.to_csv(path)



    def plotCounts(self, countedData, totalCounts, datasetType):
        # Remove commas from string values and convert to float
        # countedData = countedData.applymap(lambda x:
        #                                    float(x.replace(',', ''))
        #                                    if isinstance(x, str) else x)
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
        heatmap = sns.heatmap(countedData, annot=True, fmt=',d', cmap=cMapCustom,
                              cbar=True, linewidths=self.lineThickness-1,
                              linecolor='black', square=False, center=None,
                              annot_kws={'fontweight': 'bold'}, vmax=maxRound)
        ax.set_xlabel('Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        fig.tight_layout()

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

        # Modify the colorbar
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        # File path
        figName = f'counts - {self.enzymeName} - {datasetType}.png'
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving Fig: {datasetType}\n     {path}\n')

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
        self.log(f'Dataset: {self.datasetTag}\n')
        self.rfExp = pd.DataFrame(
            0.0, index=self.countsExp.index, columns=self.countsExp.columns
        )
        for pos in self.countsExp.columns:
            self.rfExp.loc[:, pos] = self.countsExp[pos] / sum(self.countsExp[pos])
        self.log(f'RF Experimental:\n{self.rfExp}\n')

        self.rfBg = pd.DataFrame(
            0.0, index=self.countsBg.index, columns=self.countsBg.columns
        )
        for pos in self.countsBg.columns:
            self.rfBg.loc[:, pos] = self.countsBg[pos] / sum(self.countsBg[pos])
        self.log(f'RF Background:\n{self.rfBg}')



    def calculateEntropy(self):
        self.log('\n\n============================= Calculate: Entropy '
                 '=============================')
        self.log(f'Dataset: {self.datasetTag}\n')

        self.entropy = pd.DataFrame(0.0, index=self.rfExp.columns, columns=['ΔS'])
        self.entropyMax = np.log2(len(self.rfExp.index))
        for indexColumn in self.rfExp.columns:
            S = 0
            for indexRow, probRatio in self.rfExp.iterrows():
                prob = probRatio[indexColumn]
                if prob == 0:
                    continue
                else:
                    S += -prob * np.log2(prob)
            self.entropy.loc[indexColumn, 'ΔS'] = self.entropyMax - S
        self.log(f'{self.entropy}\n\nMax Entropy: {self.entropyMax.round(6)}')



    def calculateEnrichment(self, releasedCounts=False, combinedMotifs=False,
                            posFilter=False, relFilter=False, releasedIteration=False):
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
        if len(self.rfBg.columns) == 1: ##
            # Eval: ES
            for pos in self.rfExp.columns:
                for AA in self.rfExp.index:
                    rf = self.rfBg.loc[AA, self.rfBg.columns[0]]
                    if rf == 0:
                        rf = 1
                    matrix.loc[AA, pos] = np.log2(self.rfExp.loc[AA, pos] / rf)
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
                    rf = self.rfBg.loc[AA, pos]
                    if rf == 0:
                        rf = 1
                    matrix.loc[AA, pos] = np.log2(self.rfExp.loc[AA, pos] / rf)
        if releasedCounts:
            self.log(f'Enrichment Score: Released Counts\n'
                  f'{matrix.round(self.roundVal)}\n')
            self.log(f'RF Experimental:\n{self.rfExp}\n\n'
                     f'RF Background:\n{self.rfBg}\n')
        else:
            self.log(f'Enrichment Score: {self.datasetTag}\n'
                  f'{matrix.round(self.roundVal)}\n')

            # Evaluate stack heights
            evalMatrix(matrix)

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
                                           self.entropy.loc[indexColumn, 'ΔS'])

        # Record values
        if releasedCounts:
            self.eMapReleased = matrix.copy()
            self.eMapReleasedScaled = heights.copy()
        else:
            self.eMap = matrix.copy()
            self.eMapScaled = heights.copy()

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
        heights = heights.replace([np.inf, -np.inf], 0)
        self.heights = heights
        self.log(f'\nResidue Heights: {self.datasetTag}\n'
                 f'{self.heights}\n')

        # Evaluate stack heights
        evalMatrix(self.heights)


        # Plot: Enrichment Map
        x = {
            'eMap': False, 'elogo': False, 'eMapSc': False, 'wLogo': False,
            'words': False, 'barCounts': False, 'barRF': False,
            'exp_counts': False, 'bg_counts': False
        }
        self.figures['eMap'] = (
            self.plotEnrichmentScores(
                dataType='Enrichment', releasedCounts=releasedCounts,
                posFilter=posFilter, relFilter=relFilter
            )
        )
        self.figures['eMapSc'] = (
            self.plotEnrichmentScores(
                dataType='Scaled Enrichment', releasedCounts=releasedCounts,
                posFilter=posFilter, relFilter=relFilter
            )
        )
        # Plot: Enrichment Logo
        # self.figures['eLogo'] = (
        #     self.plotEnrichmentLogo(
        #         releasedCounts=releasedCounts,
        #         posFilter=posFilter, relFilter=relFilter
        #     )
        # )

        # # Calculate & Plot: Weblogo
        # self.calculateWeblogo(probability=self.rfExp, releasedCounts=releasedCounts,
        #                       combinedMotifs=combinedMotifs)



    def plotEnrichmentScores(self, dataType, releasedCounts=False,
                             posFilter=False, relFilter=False):
        # Select: Dataset
        scaleData = False
        if 'scaled' in dataType.lower():
            scaleData = True
            if releasedCounts:
                scores = self.eMapReleasedScaled
            else:
                scores = self.eMapScaled
        else:
            if releasedCounts:
                scores = self.eMapReleased
            else:
                scores = self.eMap

        # Define: Figure title
        datasetType = self.datasetTag
        if releasedCounts:
            title = f'{self.enzymeName}'
            datasetType = self.datasetTagMotif # <----- Do we need this -----
        else:
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
            cBarMax = np.max(scores)
            cBarMin = -1 * cBarMax
        else:
            cBarMin = np.min(scores)
            cBarMax = -1 * cBarMin

        # Plot the heatmap with numbers centered inside the squares
        if self.figEMSquares:
            fig, ax = plt.subplots(figsize=self.figSizeSq)
            heatmap = sns.heatmap(
                scores, annot=False, cmap=cMapCustom, cbar=True,
                linewidths=self.lineThickness - 1, linecolor='black',
                square=True, center=None, vmax=cBarMax, vmin=cBarMin
            )
        else:
            fig, ax = plt.subplots(figsize=self.figSize)
            heatmap = sns.heatmap(
                scores, annot=True, fmt='.3f', cmap=cMapCustom, cbar=True,
                linewidths=self.lineThickness - 1, linecolor='black',
                square=False, center=None, vmax=cBarMax, vmin=cBarMin,
                annot_kws={'fontweight': 'bold'}
            )
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        ax.set_xlabel('Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        fig.tight_layout()

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

        # Modify the colorbar
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        # File path
        if scaleData:
            figName = f'eMap - {self.enzymeName} - {datasetType}.png'
        else:
            figName = f'eMap Scaled - {self.enzymeName} - {datasetType}.png'
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving Fig: {datasetType}\n     {path}\n')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig)

        return figName



    def plotEnrichmentLogo(self, combinedMotifs=False, releasedCounts=False,
                           posFilter=False, relFilter=False, relIteration=False):
        print('============================= Plot: Enrichment Logo '
              '=============================')
        # Define: Figure title
        datasetType = self.datasetTag
        if releasedCounts or combinedMotifs or len(self.motifIndexExtracted) > 1:
            title = self.titleReleased
            datasetType = self.datasetTagMotif
        else:
            title = self.title
            datasetType = self.datasetTagMotif
        if len(self.datasetTag.replace('[', '').replace(']', '').replace('-', '')) > 40:
            title = title.replace('Register ', 'Register\n')

        # Print: data
        print(f'Dataset: {self.datasetTag}\n'
              f'Unique Substrates: {self.countExpUnique:,}')
        if self.motifFilter:
            print(f'Figure Number: '
                  f'{self.saveFigureIteration}')
        if posFilter:
            if relFilter:
                print(f'Releasing Filter: {posFilter}')
            else:
                print(f'Applying Filter: {posFilter}')
        print(f'\nResidue heights:\n'
              f'{self.heights}\n')

        # Set local parameters
        if self.bigAAonTop:
            stackOrder = 'big_on_top'
        else:
            stackOrder = 'small_on_top'

        # Calculate: Max and min
        columnTotals = [[], []]
        for indexColumn in self.heights.columns:
            totalPos = 0
            totalNeg = 0
            for value in self.heights.loc[:, indexColumn]:
                if value > 0:
                    totalPos += value
                elif value < 0:
                    totalNeg += value
            columnTotals[0].append(totalPos)
            columnTotals[1].append(totalNeg)
        yMax = max(columnTotals[0])
        yMin = min(columnTotals[1])
        print(f'y Max: {np.round(yMax, 4)}\n'
              f'y Min: {np.round(yMin, 4)}\n')

        # Rename columns for logomaker script
        data = self.heights.copy()
        data.columns = range(len(data.columns))


        def plotLogo(limitYAxis=False):
            # Plot the sequence motif
            fig, ax = plt.subplots(figsize=self.figSize)
            motif = logomaker.Logo(data.transpose(), ax=ax, color_scheme=self.colorsAA,
                                   width=0.95, stack_order=stackOrder)
            motif.ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
            fig.tight_layout()

            # Set tick parameters
            ax.tick_params(axis='both', which='major', length=self.tickLength,
                           labelsize=self.labelSizeTicks)

            # Set borders
            motif.style_spines(visible=False)
            motif.style_spines(spines=['left', 'bottom'], visible=True)
            for spine in motif.ax.spines.values():
                spine.set_linewidth(self.lineThickness)

            # Set x-ticks
            motif.ax.set_xticks([pos for pos in range(len(self.heights.columns))])
            motif.ax.set_xticklabels(self.heights.columns, fontsize=self.labelSizeTicks,
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

            # Evaluate dataset for fixed residues
            spacer = np.diff(motif.ax.get_xticks())  # Find the space between each tick
            spacer = spacer[0] / 2

            # Use the spacer to set a gray background to fixed residues
            for index, position in enumerate(self.xAxisLabels):
                if position in self.fixedPos:
                    # Plot gray boxes on each side of the xtick
                    motif.ax.axvspan(index - spacer, index + spacer,
                                     facecolor='darkgrey', alpha=0.2)

            # File path
            if limitYAxis:
                figName = f'eLogo yMin - {self.enzymeName} - {datasetType}.png'
            else:
                figName = f'eLogo - {self.enzymeName} - {datasetType}.png'
            path = os.path.join(self.pathFigs, figName)
            print(f'Saving Fig: {datasetType}\n     {path}\n')

            # Encode the figure
            figBase64 = self.encodeFig(fig)
            with open(path, "wb") as file:
                file.write(base64.b64decode(figBase64))

            # Close the figure to free memory
            plt.close(fig)

            return figName

        # Plot figure
        plotLogo() # Full y-axis

        # # Adjust yMin to fit the largest negative AA
        # yMin = 0
        # for col in self.heights.columns:
        #     for row in self.heights.index:
        #         if self.heights.loc[row, col] < yMin:
        #             yMin = self.heights.loc[row, col]
        # print(f'Adjusting Y Min:\n'
        #       f'y Max: {np.round(yMax, 4)}\n'
        #       f'y Min: {np.round(yMin, 4)}\n\n')
        # plotLogo(limitYAxis=True) # Limited y-axis
