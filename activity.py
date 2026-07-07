from functions import pressKey
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import sys


# Set options
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
pd.set_option('display.float_format', '{:,.3f}'.format)


"""
    If a key in inData is "% St Dev", table will exclude it.
    The "% Product" values will be converted from decimal to percentage.
    The order of the keys in inData will have the same order in the table.
"""

# Input: Data
inPlotBoth = False # Add the secondary enzyme to the figures
inEnzyme = f'M{"ᵖʳᵒ"}2'
inEnzyme2 = f'M{"ᵖʳᵒ"}' # Secondary enzyme
inSigFigs = 0
inRoundVal = 3
inNatLog = False
inSubstrates = ['AVLQSGFR', 'VILQSGFR', 'VILQTGFR', 'VILQSPFR',
                'VILHSGFR', 'VIMQSGFR', 'VPLQSGFR', 'NILQSGFR']
inPredActivity = [0.595, 1.0, 0.008, 0.004, 0.055, 0.417, 0.003, 0.049]
inPredActivity2 = [0.748, 1.0, 0.007, 0.009, 0.03, 0.453, 0.005, 0.027]
inEmptyList = [0 for i in range(1, len(inSubstrates) + 1)]
inData = {
    'Substrates': inSubstrates,

    f'% Product {inEnzyme}': [46.1, 49.5, 14.5, 0.0, 13.1, 37.0, 0.0, 16.1],
    f'Activity Z {inEnzyme}': inEmptyList,
    f'Activity Rank {inEnzyme}': inEmptyList,
    f'% St Dev {inEnzyme}': [0.1, 0.09, 0.02, 0, 0.06, 0.09, 0, 0.05],
    f'Predicted {inEnzyme}': inPredActivity,
    f'Predicted Z {inEnzyme}': inEmptyList,
    f'Predicted Rank {inEnzyme}': inEmptyList,

    f'% Product {inEnzyme2}': [32.1, 39.1, 14.9, 0.0, 16.0, 36.5, 0.0, 15.6],
    f'Activity Z {inEnzyme2}': inEmptyList,
    f'Activity Rank {inEnzyme2}': inEmptyList,
    f'% St Dev {inEnzyme2}': [0.01, 0.058, 0.025, 0.0, 0.027, 0.044, 0.0, 0.033],
    f'Predicted {inEnzyme2}': inPredActivity2,
    f'Predicted Z {inEnzyme2}': inEmptyList,
    f'Predicted Rank {inEnzyme2}': inEmptyList,
}

# inActivityMpro2 = [61, 29, 5, 28, 13, 37, 22, 22, 22, 28, 33]
# inPredictSubstrates = ['AVLQSGFR', 'VTFQSAVK', 'ATVQSKMS', 'ATLQAIAS',
#                        'VKLQNNEL', 'VRLQAGNA', 'PMLQSADA', 'TVLQAVGA',
#                        'ATLQAENV', 'TRLQSLEN', 'PKLQSSQA']
# inSubstrates = ['AVLQSGFR', 'VILQSGFR', 'VILQTGFR', 'VILQSPFR',
#                 'VILHSGFR', 'VIMQSGFR', 'VPLQSGFR', 'NILQSGFR',
#
#                 'VTFQSAVK', 'ATVQSKMS', 'ATLQAIAS', 'VKLQNNEL',
#                 'VRLQAGNA', 'PMLQSADA', 'TVLQAVGA', 'ATLQAENV',
#                 'TRLQSLEN', 'PKLQSSQA'
#                 ]
# inEmptyList = [0 for i in range(1, len(inSubstrates) + 1)]
# inExpActivity = [1.0, 1.07, 0.315]
# inData = {
#     'Substrates': inSubstrates,
#
#     f'% Product {inEnzyme}': [46.1, 49.5, 14.5, 0.0, 13.1, 37.0, 0.0, 16.1,
#                               29, 5, 28, 13, 37, 22, 22, 22, 28, 33],
#     f'Activity Z {inEnzyme}': inEmptyList,
#     f'Activity Rank {inEnzyme}': inEmptyList,
#     f'% St Dev {inEnzyme}': [0.1, 0.09, 0.02, 0, 0.06, 0.09, 0, 0.05,
#                              0,0,0,0,0,
#                              0,0,0,0],
#     f'Predicted {inEnzyme}': [0.595, 1.0, 0.008, 0.004, 0.055, 0.417, 0.003, 0.049,
#                               0.258, 0.005, 0.151, 0.053, 0.342, 0.149, 0.073,
#                               0.708, 0.124, 0.051],
#     f'Predicted Z {inEnzyme}': inEmptyList,
#     f'Predicted Rank {inEnzyme}': inEmptyList,
#
#     # f'% Product {inEnzyme2}': inEmptyList,
#     # f'Activity Z {inEnzyme2}': inEmptyList,
#     # f'Activity Rank {inEnzyme2}': inEmptyList,
#     # f'% St Dev {inEnzyme2}': inEmptyList,
#     # f'Predicted {inEnzyme2}': inEmptyList,
#     # f'Predicted Z {inEnzyme2}': inEmptyList,
#     # f'Predicted Rank {inEnzyme2}': inEmptyList,
# }

# Input: Figures
inPlotBarGraph = False
inPlotTable = False
inFigSaveTag = '' # Add label to saved figures
inSavePath = '/Users/ca34522/Documents/Papers/COMET/Figures/'
inFigTitle = f'Enzyme Activity'
inColor1 = '#BF5700'
inColor2 = '#F8971F'
inFigResolution = 600

# Input: Figure Params
inFigSize = (9.5, 8)
inTickLength = 4
inLinewidth = 1.5
inTitleSize = 20
inLabelSize = 17
inLabelTickSize = 15

# Input: Table
if inPlotBoth:
    inTableCols = [
        f'% Product {inEnzyme2}', f'Activity Rank {inEnzyme2}',
        f'Predicted Rank {inEnzyme2}',
        f'% Product {inEnzyme}', f'Activity Rank {inEnzyme}',
        f'Predicted Rank {inEnzyme}'
    ]

    # Input: Z-Scores
    inCalcZScores = [ # Format: (x, y), use "x" values to calc the Z-Score saved in "y"
        (f'% Product {inEnzyme2}', f'Activity Z {inEnzyme2}'),
        (f'Predicted {inEnzyme2}', f'Predicted Z {inEnzyme2}'),
        (f'% Product {inEnzyme}', f'Activity Z {inEnzyme}'),
        (f'Predicted {inEnzyme}', f'Predicted Z {inEnzyme}')
    ]
else:
    inTableCols = [
        f'% Product {inEnzyme}', f'Activity Rank {inEnzyme}', f'Predicted Rank {inEnzyme}'
    ]

    # Input: Z-Scores
    inCalcZScores = [  # Format: (x, y), use "x" values to calc the Z-Score saved in "y"
        (f'% Product {inEnzyme}', f'Activity Z {inEnzyme}'),
        (f'Predicted {inEnzyme}', f'Predicted Z {inEnzyme}')
    ]

def normalizeData():
    for key in inData.keys():
        if inEnzyme2 in key:
            continue
        if '% Product' in key and 'Rank' not in key and 'ln' not in key:
            maxValue = max(inData[key])
            norm = []
            for val in inData[key]:
                norm.append(val / maxValue)
            inData[key] = norm


def convertNum(data, key):
    keyStDev = '% St Dev'
    if key != '% St Dev':
        keyStDev = f'% St Dev{key.replace("% Product", "")}'

    activity = [] # Convert to %
    for idx, substrate in enumerate(data['Substrates']):
        act = data[key][idx] * 100
        stdev = int(data[keyStDev][idx] * 100)
        if inSigFigs == 0 or not inSigFigs:
            act = int(act)
        else:
            act = round(act, inSigFigs)
        if stdev == 0.0:
            activity.append(act)
        else:
            activity.append(f'{act} ± {stdev}')
    return activity


def plotTable():
    data = {}
    for col in inData.keys():
        # print(f'\nCOL: {col}')
        if '% st dev' in col.lower():
            continue
        elif 'substrates' in col.lower():
            data[col] = inData[col]
        elif '% product' in col.lower() and 'rank' not in col.lower():
            data[col] = convertNum(inData, col)
        elif ('predicted' in col.lower() and
              'rank' not in col.lower() and
              'ln' not in col.lower()):
            # print(f'Col: {col}')
            key = col.replace('Rank', 'Z')
            vals = []
            for v in inData[col]:

                vals.append(round(v, inRoundVal))
                # if inSigFigs == 0 or not inSigFigs:
                #     vals.append(int(v * 100))
                # else:
                #     vals.append(str(round((v * 100), 0)))
            # data[col] = [str(v) for v in vals]
            data[col] = vals
            # print(f'* Pred: {data[col]}\n')
        elif 'rank' in col.lower():
            key = col.replace('Rank', 'Z')
            print(f'* Rank: {col}, {key}')
            # print(inData[key])
            rank = list(pd.Series(inData[key]).rank(ascending=False, method='min'))
            l = []
            for val in rank:
                l.append(f'{str(int(val))}')
            data[col] = l
        else:
            x = []
            for i in range(len(inData[col])):
                x.append(round(inData[col][i], inRoundVal))
            data[col] = x
            # print(f'* Data: {data[col]}')
    # print()
    # for k, v in data.items():
    #     print(f'* {k}\n{v}\n')
    # print()

    # Create table
    if 'Substrates' not in inTableCols:
        inTableCols.insert(0, 'Substrates')
    table = pd.DataFrame('', index=[], columns=inTableCols)
    for col in table.columns:
        table.loc[:, col] = data[col]

    # table.drop(f'Predicted {inEnzyme2}', axis=1, inplace=True)
    # table.drop(f'Predicted {inEnzyme}', axis=1, inplace=True)
    print(f'Table:\n{table}\n')

    print()
    for col in table.columns:
        x = table.loc[:, col]
        print(col)
        for i in x:
            print(i)
        print()
    # sys.exit()

    # Make figure
    h = (len(table.index) / 2) - 1
    w = len(table.columns) * 2
    fig, ax = plt.subplots(figsize=(1.2*w, 1.2*h))
    ax.axis('off')  # hide axes
    tbl = plt.table(cellText=table.values,
                    colLabels=table.columns,
                    loc='center',
                    cellLoc='center',
                    bbox=(0, 0, 1, 1)
                    )
    tbl.scale(1, 2)

    # Bold headers
    for i in range(len(table.columns)):
        tbl[0, i].set_text_props(weight='bold')

    tbl.auto_set_font_size(False)
    tbl.set_fontsize(16)

    for (row, col), cell in tbl.get_celld().items():
        text = cell.get_text()
        text.set_fontname('Times New Roman')
        text.set_verticalalignment('top')

    fig.tight_layout(pad=0.5)
    fig.canvas.mpl_connect('key_press_event', pressKey)
    plt.show()

    if inSavePath:
        figName = 'enzActivity_table.png'
        if inFigSaveTag:
            figName = figName.replace('.png', f'_{inFigSaveTag}.png')
        path = os.path.join(inSavePath, figName)
        fig.savefig(path, dpi=inFigResolution)
        print(f'Saving figure at path:\n'
              f'     {path}\n\n')
    else:
        print(f'The figure was not saved\n\n')


def zScore(tag):
    # print(f'Calculate Z Score: {tag}')
    values = inData[tag]
    avg = np.mean(values)
    stdev = np.std(values)
    z = []
    for x in values:
        z.append((x - avg) / stdev)
    for i in range(len(values)):
        x, y = round(float(values[i]), inRoundVal), round(float(z[i]), inRoundVal)
        # print(f'* Value: {x}, Z Score: {y}')
    # print()
    return z


def ln(vals):
    x = []
    for val in vals:
        l = np.log(val)
        x.append(l)
        print(f'*: ln({round(val, inRoundVal)}) -> {round(l, inRoundVal)}')
    print()
    return x


def fnExp(x, a, b, c):
    return a * np.exp(b * x) + c


def fitData(x, y):
    # Fit the curve
    popt, pcov = curve_fit(fnExp, x, y, p0=[1, 1, 0], maxfev=10000)
    a, b, c = popt

    # Generate smooth curve for plotting
    xFit = np.linspace(x.min(), x.max(), 300)
    yFit = fnExp(xFit, *popt)

    # R² for the exponential fit
    yPred = fnExp(x, *popt)
    ss_res = np.sum((y - yPred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    return xFit, yFit, r2


def plotBars(data, e1, e2, barWidth=0.35):
    substrates = data['Substrates']
    xTicks = np.arange(len(substrates))
    y1 = data[f'% Product {e1}']
    y2 = data[f'% Product {e2}']
    d = pd.DataFrame(0.0, index=substrates, columns=[e1, e2])
    d.loc[substrates, e1] = y1
    d.loc[substrates, e2] = y2
    print(f'Bar Graph: Normalized Activity\n'
          f'{d}\n\n')

    # Labels
    l1 = e1
    l2 = e2
    if f'M{"ᵖʳᵒ"}' in e1 or f'M{"ᵖʳᵒ"}' in e2:
        l1 = f'SARS-CoV {e1}'
        l2 = f'SARS-CoV-2 {e1}'


    # Plot bar graph
    fig, ax = plt.subplots(figsize=inFigSize)
    ax.bar(xTicks - barWidth / 2, y1, barWidth, label=l1,
           color=inColor2, edgecolor='black', linewidth=inLinewidth)
    ax.bar(xTicks + barWidth / 2, y2, barWidth, label=l2,
           color=inColor1, edgecolor='black', linewidth=inLinewidth)
    plt.title(inFigTitle, fontsize=inTitleSize, fontweight='bold')
    ax.set_ylabel('Normalized Activity', fontsize=inLabelSize)

    # Set the thickness of the figure border
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(inLinewidth)

    # Legend
    legend_props = {
        'size': inLabelTickSize-2,
        'weight': 'bold'
    }
    ax.legend(edgecolor='black', prop=legend_props, loc='best')

    # Set xticks
    ax.set_xticks(xTicks)
    ax.set_xticklabels(substrates, rotation=45)

    # Set yticks
    ax.set_ylim([0, 1.1])

    # Set tick parameters
    ax.tick_params(axis='both', which='major', length=inTickLength,
                   labelsize=inLabelTickSize, width=inLinewidth)

    fig.canvas.mpl_connect('key_press_event', pressKey)
    plt.tight_layout()
    plt.show()

    if inSavePath:
        figName = 'enzActivity_bars.png'
        if inFigSaveTag:
            figName = figName.replace('.png', f'_{inFigSaveTag}.png')
        path = os.path.join(inSavePath, figName)
        fig.savefig(path, dpi=inFigResolution)
        print(f'Saving figure at path:\n'
              f'     {path}\n\n')
    else:
        print(f'The figure was not saved\n\n')


# ========================================================================================
normalizeData()

if inPlotBarGraph:
    plotBars(data=inData, e1=inEnzyme2, e2=inEnzyme)

# Calculate Z-Scores
for tags in inCalcZScores:
    inData[tags[1]] = zScore(tags[0])
    # inData[tags[1]] = ln(inData[tags[0]])

if inPlotTable:
    plotTable()
# sys.exit()

# Evaluate the natural log
if inNatLog:
    inData[f'Predicted {inEnzyme2}'] = ln(inData[f'Predicted {inEnzyme2}'])
    inData[f'Predicted {inEnzyme}'] = ln(inData[f'Predicted {inEnzyme}'])


# ========================================================================================
# Evaluate data
data = pd.DataFrame(0.0, index=inData['Substrates'], columns=[])
data[f'Activity {inEnzyme}'] = inData[f'% Product {inEnzyme}']
data[f'Activity Z {inEnzyme}'] = inData[f'Activity Z {inEnzyme}']
# print(f'L1: {len(inData['Substrates'])}\n'
#       f'L2: {len(inData[f'Predicted {inEnzyme}'])}')
data[f'Pred {inEnzyme}'] = inData[f'Predicted {inEnzyme}']
data[f'Predicted Z {inEnzyme}'] = inData[f'Predicted Z {inEnzyme}']
x = ['AVLQSG', 'VILQSG', 'VILQTG', 'VILQSP',
     'VILHSG', 'VIMQSG', 'VPLQSG', 'NILQSG']
actR2 = round(r2_score(data[f'Activity {inEnzyme}'],
                       data[f'Pred {inEnzyme}']), inRoundVal)
actZR2 = round(r2_score(data[f'Activity Z {inEnzyme}'],
                        data[f'Predicted Z {inEnzyme}']), inRoundVal)
print(f'{data}\n\n')

# Figure labels
l1, l2 = f'{inEnzyme}', f'{inEnzyme2}'
l1, l2 = f'SARS-CoV-2 {inEnzyme.replace('2', '')}', f'SARS-CoV {inEnzyme2}'

# Plot data
mkr1, mkr2, edgeWidth = 'D', 'o', 1
fig, ax = plt.subplots(figsize=inFigSize)
plt.title(inFigTitle, fontsize=inTitleSize, fontweight='bold')
x, y = f'Activity Z {inEnzyme}', f'Predicted Z {inEnzyme}'
x_fit, y_fit, fitCurve = fitData(x=data[x].values, y=data[y].values)
l1 += f' R² = {fitCurve:.3f}'
data.plot(x=x, y=y, ax=ax, color=inColor1, marker=mkr1, linestyle='none',
          markeredgecolor='black', markeredgewidth=edgeWidth, legend=l1)
ax.plot(x_fit, y_fit, color=inColor1, linestyle='-', linewidth=inLinewidth)


# Evaluate the secondary set
if inPlotBoth:
    data2 = pd.DataFrame(0.0, index=inData['Substrates'], columns=[])
    data2[f'Activity {inEnzyme2}'] = inData[f'% Product {inEnzyme2}']
    data2[f'Activity Z {inEnzyme2}'] = inData[f'Activity Z {inEnzyme2}']
    data2[f'Pred {inEnzyme2}'] = inData[f'Predicted {inEnzyme2}']
    data2[f'Predicted Z {inEnzyme2}'] = inData[f'Predicted Z {inEnzyme2}']
    print(f'{data2}\n\n')
    actEnzR2 = round(
        r2_score(data2[f'Activity {inEnzyme2}'],
                 data2[f'Pred {inEnzyme2}']), inRoundVal
    )
    actEnzZR2 = round(
        r2_score(data2[f'Activity Z {inEnzyme2}'],
                 data2[f'Predicted Z {inEnzyme2}']), inRoundVal
    )

    x, y = f'Activity Z {inEnzyme2}', f'Predicted Z {inEnzyme2}'
    x_fit2, y_fit2, fitCurve2 = fitData(x=data2[x].values, y=data2[y].values)
    l2 += f' R² = {fitCurve2:.3f}'
    data2.plot(x=x, y=y, ax=ax, color=inColor2, marker=mkr2, linestyle='none',
               markeredgecolor='black', markeredgewidth=edgeWidth, legend=l2)
    ax.plot(x_fit2, y_fit2, color=inColor2, linestyle='-', linewidth=1.5)
    ax.legend(
        prop=FontProperties(size=inLabelTickSize-2, weight='bold'), loc='upper left',
        handles=[
            Line2D([], [], color=inColor2, marker=mkr2,
                   markeredgecolor='black', markeredgewidth=edgeWidth,
                   linewidth=inLinewidth, label=l2),
            Line2D([], [], color=inColor1, marker=mkr1,
                   markeredgecolor='black', markeredgewidth=edgeWidth,
                   linewidth=inLinewidth, label=l1),
        ]
    )
else:
    ax.legend(prop=FontProperties(size=inLabelTickSize-2, weight='bold'), loc='upper left',
              handles=[Line2D([], [], linestyle='None', marker='None',
                              markeredgecolor='black', markeredgewidth=edgeWidth,
                              label=l1)], handletextpad=0, handlelength=0
              )

# Set the thickness of the figure border
for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(inLinewidth)

# X Axis
ax.set_xlabel('Experimental Activity Z Scores', fontsize=inLabelSize)
ax.set_xticks(ax.get_xticks())
ax.set_xticklabels(ax.get_xticklabels(), ha='center', fontsize=inLabelSize-2)

# Y Axis
ax.set_ylabel('Predicted Activity Z Scores', fontsize=inLabelSize)
ax.set_yticks(ax.get_yticks())
# ax.set_ylim([-1.0, 2.5]) # Set yticks

# Set tick parameters
ax.tick_params(axis='both', which='major', length=inTickLength,
                   labelsize=inLabelTickSize, width=inLinewidth)

plt.tight_layout()
fig.canvas.mpl_connect('key_press_event', pressKey)
plt.show()

if inSavePath:
    figName = 'enzActivity.png'
    if inFigSaveTag:
        figName = figName.replace('.png', f'_{inFigSaveTag}.png')
    path = os.path.join(inSavePath, figName)
    fig.savefig(path, dpi=inFigResolution)
    print(f'Saving figure at path: 1\n'
          f'     {path}\n\n')
else:
    print(f'The figure was not saved\n\n')