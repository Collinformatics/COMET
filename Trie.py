class TrieNode:
    # Create nodes
    def __init__(self):
        self.children = {}
        self.endOfWord = False


class Trie:
    def __init__(self):
        # Create the root as an empty node
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        # Add words to the trie
        current =  self.root

        for AA in word:
            if AA not in current.children:
                current.children[AA] = TrieNode()
            current = current.children[AA]

    def search(self, word):
        # Search for a word in the trie
        current = self.root

        for AA in word:
            if AA not in current.children:
                return False
            current = current.children[AA]
        return True

    def startsWith(self, prefix):
        # Find the prefix in the trie
        current = self.root

        for AA in prefix:
            if AA not in current.children:
                return False
        return True



def evaluateSubtrees(self, trie, motifTrie):
    print('============================= Evaluate Suffix Tree '
          '==============================')
    print(f'Datapoints: {len(motifTrie.keys())}')

    def subtreeTable(subtreeFreq):
        # Sort motifs by length
        sortedMotifs = sorted(subtreeFreq.keys(), key=len)

        # Organize motifs by their length and sort by frequency (highest first)
        motifGroups = {}
        for motif in sortedMotifs:
            length = len(motif)
            if length not in motifGroups:
                motifGroups[length] = []
            motifGroups[length].append((motif, subtreeFreq[motif]))

        # Sort motifs in each length group by frequency (descending)
        for length in motifGroups:
            motifGroups[length].sort(key=lambda x: x[1], reverse=True)

        # Convert motifs back to formatted strings
        for length in motifGroups:
            motifGroups[length] = [f"{motif}: {round(freq, 5)}"
                                   for motif, freq in motifGroups[length]]

        # Find the max number of motifs in any length group
        maxRows = max(len(motifs) for motifs in motifGroups.values())

        # Construct the table row by row
        tableData = []
        for i in range(maxRows):
            row = []
            for length in sorted(motifGroups.keys()):
                motifs = motifGroups[length]
                row.append(motifs[i] if i < len(
                    motifs) else "") # Fill missing values with empty strings
            tableData.append(row)

        # Convert to DataFrame
        motifTable = pd.DataFrame(tableData,
                                  index=range(1, len(motifTrie.keys()) + 1),
                                  columns=[str(length)
                                           for length in sorted(motifGroups.keys())])
        print(f'{motifTable}\n\n')

        return motifTable



    def printTrie(node, level=0, path=""):
        # Recursively print the Trie structure
        if node is None:
            return

        # Print: Current node's path and level
        print("  " * level + f"Level {level}: {path}")

        # Recursively print all children of the current node
        for char, nodeChild in node.children.items():
            printTrie(nodeChild, level + 1, path + char)

    motifsTotal = 0
    for motif, count in motifTrie.items():
        motifsTotal += count
    print(f'Total Motifs: {motifsTotal:,}\n')

    # Evaluate: Partial sequence counts
    subtreeCount = {}
    motifLength = len(next(iter(motifTrie)))
    for index in range(motifLength):
        for motif, count in motifTrie.items():
            subSeq = motif[0:index +1]
            if subSeq in subtreeCount.keys():
                subtreeCount[subSeq] += count
            else:
                subtreeCount[subSeq] = count

    # Evaluate: Partial sequence frequency
    subtreeFreq = {}
    for subSeq, count in subtreeCount.items():
        subtreeFreq[subSeq] = count / motifsTotal
    prevSeqLen = 1
    for subSeq, count in subtreeFreq.items():
        if len(subSeq) != prevSeqLen:
            prevSeqLen = len(subSeq)
    motifTable = subtreeTable(subtreeFreq)


    # Plot the trie
    printTrie(trie.root)
    print('\n')

    return motifTable



def suffixTree(self, substrates, N, entropySubFrame, indexSubFrame, entropyMin,
               datasetTag, dataType):
    print('================================== Suffix Tree '
          '==================================')
    if datasetTag is None:
        print(f'Dataset: {purple}{self.enzymeName} - Unfixed{resetColor}\n')
    else:
        print(f'Dataset: {purple}{self.enzymeName} {dataType} {datasetTag}'
              f'{resetColor}\n')


    trie = Trie() # Initialize Trie
    motifs = {}
    indexStart = min(indexSubFrame)
    indexEnd = max(indexSubFrame)

    # Print: Substrates
    iteration = 0
    substrates = dict(sorted(substrates.items(), key=lambda item: item[1],
                             reverse=True))
    for substrate, count in substrates.items():
        iteration += 1
        print(f'     {pink}{substrate}{resetColor}, '
              f'Counts: {red}{count:,}{resetColor}')
        if iteration >= self.printNumber:
            break
    print('\n')

    # Find motif positions based on the entropy threshold
    indexPos = []
    for index in entropySubFrame.index:
        posEntropy = entropySubFrame.loc[index, 'ΔS']
        if posEntropy >= entropyMin:
            indexPos.append(int(index.replace('R', '')) - 1)
    print(f'Index Pos: {indexPos}')

    motifTrie = {}
    countsMotif = 0
    def addMotif(motif, count):
        # Extract important AAs from the motif
        motif = ''.join(motif[index] for index in indexPos)

        # Add motif to the trie
        if motif in motifTrie.keys():
            motifTrie[motif] += count
        else:
            motifTrie[motif] = count
            trie.insert(motif)


    # Extract the motifs
    motifCount = 0
    for substrate, count in substrates.items():
        motif = substrate[indexStart:indexEnd + 1]
        if motif in motifs:
            motifs[motif] += count
        else:
            motifs[motif] = count
            motifCount += 1

        # Add the motif to the tree
        addMotif(motif, count)
        countsMotif = len(motifTrie.keys())
        if countsMotif >= N:
            break
    motifs = dict(sorted(motifs.items(), key=lambda item: item[1], reverse=True))
    motifTrie = dict(sorted(motifTrie.items(), key=lambda item: item[1],
                            reverse=True))

    # Print: Motifs
    print(f'Extracted Motifs:')
    for index, (motif, count) in enumerate(motifs.items()):
        print(f'{index+1}:{blue} {motif}{resetColor} '
              f'Count:{red} {count:,}{resetColor}')
    print('\n')

    # Print: Trie
    print(f'Extracted Trie:')
    for index, (seq, count) in enumerate(motifTrie.items()):
        print(f'{index + 1}:{pink} {seq}{resetColor} '
              f'Count:{red} {count:,}{resetColor}')
    print('\n')

    # Calculate: RF
    motifTable = self.evaluateSubtrees(trie=trie, motifTrie=motifTrie)

    # Plot the Trie
    self.plotTrie(trie=trie, motifTable=motifTable, countsMotif=countsMotif,
                  datasetTag=datasetTag)



def plotTrie(self, trie, motifTable, countsMotif, datasetTag):
    print('=============================== Plot: Suffix Tree '
          '===============================')
    import networkx as nx

    inOffset = 2000
    inNodeSizeMax = 800
    inNodeSizeMin = 100
    inFontSize = 10
    inScaleX = 2
    inScaleY = 1

    # Calculate: Node size
    nodeSizes = pd.DataFrame('',
                             index=motifTable.index,
                             columns=motifTable.columns)
    for col in motifTable.columns:
        for index, entry in enumerate(motifTable[col].dropna()):
            if ": " in entry:
                motif, rf = entry.split(": ")
                nodeSize = inNodeSizeMax - (inNodeSizeMax * (1 - float(rf)))
                if nodeSize < 100:
                    nodeSize = inNodeSizeMin
                if len(motif) > 2:
                    motif = motif[-2:]
                nodeSizes.loc[index+1, col] = f'{motif}: {nodeSize:.2f}'
    print(f'Node Size:\n{nodeSizes}\n')


    def addNodesToGraph(node, graph, scaleX, scaleY, offset=inOffset,
                        nodeSizesDF=None):
        pos = {}
        nodeSizes = {}
        nodeCountLevel = {}

        # Track node index separately
        queue = [(node, None, '', '', 0,
                  1)]  # (currentNode, parentID, char, fullMotif, level, index)

        while queue:
            nodeCurrent, parent, char, motifSoFar, level, index = queue.pop(0)
            nodeID = f"{char}-{level}-{id(nodeCurrent)}"

            if level not in nodeCountLevel:
                nodeCountLevel[level] = []
            nodeCountLevel[level].append(
                (nodeCurrent, parent, char, nodeID, motifSoFar, index))

            fullMotif = motifSoFar + char  # Build full motif sequence

            # Assign node size from nodeSizesDF if available
            nodeSize = inNodeSizeMin  # Default size
            if nodeSizesDF is not None and level in nodeSizesDF.columns:
                entry = nodeSizesDF.iloc[index - 1, level]  # Use the tracked index
                if isinstance(entry, str) and ": " in entry:
                    _, size = entry.split(": ")
                    nodeSize = float(size)

            graph.add_node(nodeID, label=char, size=nodeSize)
            nodeSizes[nodeID] = nodeSize

            if parent is not None:
                graph.add_edge(parent, nodeID, arrowstyle='->')

            # Track child nodes with incremented index
            childIndex = 1
            for child_char, nodeChild in nodeCurrent.children.items():
                queue.append(
                    (nodeChild, nodeID, child_char, fullMotif, level + 1, childIndex))
                childIndex += 1  # Ensure a unique index for each child

        return pos, nodeSizes


    # Build the graph
    graph = nx.DiGraph()
    pos, nodeSizes = addNodesToGraph(trie.root, graph, scaleX=inScaleX,
                                     scaleY=inScaleY, offset=inOffset,
                                     nodeSizesDF=nodeSizes)
    finalNodeSizes = [graph.nodes[node]["size"] for node in graph.nodes]

    # Get node labels
    labels = {node: data['label'] for node, data in graph.nodes(data=True)}

    # Print: Dataset tag
    if datasetTag is None:
        figLabel = f'Suffix Tree-{self.enzymeName}-{countsMotif}-Unfixed'

    else:
        figLabel = f'Suffix Tree-{self.enzymeName}-{countsMotif}-{datasetTag}'


    # Plot the data
    fig, ax = plt.subplots(figsize=self.figSize)
    fig.canvas.mpl_connect('key_press_event', pressKey)

    # Draw graph
    nx.draw(graph, pos, with_labels=True, labels=labels, node_size=finalNodeSizes,
            node_color="#F18837", font_size=inFontSize, font_weight="bold",
            edge_color="#101010", ax=ax, arrows=False)
    plt.title(f'{self.enzymeName}: {datasetTag}\nTop {countsMotif:,} Motifs',
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.show()

    # Save the figure
    if self.saveFigures:
        # Define: Save location
        figLabel += '.png'
        saveLocation = os.path.join(self.pathSaveFigs, figLabel)

        # Save figure
        if os.path.exists(saveLocation):
            print(f'{yellow}The figure was not saved\n\n'
                  f'File was already found at path:\n'
                  f'     {saveLocation}{resetColor}\n\n')
        else:
            print(f'Saving figure at path:\n'
                  f'     {greenDark}{saveLocation}{resetColor}\n\n')
            fig.savefig(saveLocation, dpi=self.figureResolution)