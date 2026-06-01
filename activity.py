import os.path

from functions import pressKey
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
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
inEnzyme = f'M{"ᵖʳᵒ"}2'
inEnzyme2 = f'M{"ᵖʳᵒ"}' # Secondary enzyme
inSigFigs = 0
inRoundVal = 3
inNatLog = False
inSubstrates = ['AVLQSGFR', 'VILQSGFR', 'VILQTGFR', 'VILQSPFR',
                'VILHSGFR', 'VIMQSGFR', 'VPLQSGFR', 'NILQSGFR']
inEmptyList = [0 for i in range(1, len(inSubstrates) + 1)]
inData = {
    'Substrates': inSubstrates,

    f'% Product {inEnzyme}': [0.76, 0.82, 0.38, 0.0, 0.32, 0.66, 0.0, 0.23],
    f'Activity Z {inEnzyme}': inEmptyList,
    f'Activity Rank {inEnzyme}': inEmptyList,
    f'% St Dev {inEnzyme}': [0.1, 0.09, 0.02, 0, 0.06, 0.09, 0, 0.05],
    f'Predicted {inEnzyme}': [0.976, 1.0, 0.035, 0.005, 0.053, 0.362, 0.006, 0.01],
    f'Predicted Z {inEnzyme}': inEmptyList,
    f'Predicted Rank {inEnzyme}': inEmptyList,

    f'% Product {inEnzyme2}': [0.54, 0.66, 0.34, 0.0, 0.26, 0.61, 0.0, 0.22],
    f'Activity Z {inEnzyme2}': inEmptyList,
    f'Activity Rank {inEnzyme2}': inEmptyList,
    f'% St Dev {inEnzyme2}': [0.01, 0.058, 0.025, 0.0, 0.027, 0.044, 0.0, 0.033],
    f'Predicted {inEnzyme2}': [0.572, 1.0, 0.015, 0.01, 0.038, 0.51, 0.034, 0.066],
    f'Predicted Z {inEnzyme2}': inEmptyList,
    f'Predicted Rank {inEnzyme2}': inEmptyList,
}

# Input: Table
inPlotTable = False
inTableEnz = inEnzyme
inTableCols = [f'% Product {inTableEnz}',
               f'Activity Z {inTableEnz}', f'Activity Rank {inTableEnz}',
               f'Predicted Z {inTableEnz}', f'Predicted Rank {inTableEnz}']

# Input: Z-Scores
inCalcZScores = [ # Format: (x, y), use "x" values to calc the Z-Score saved in "y"
    (f'% Product {inEnzyme2}', f'Activity Z {inEnzyme2}'),
    (f'Predicted {inEnzyme2}', f'Predicted Z {inEnzyme2}'),
    (f'% Product {inEnzyme}', f'Activity Z {inEnzyme}'),
    (f'Predicted {inEnzyme}', f'Predicted Z {inEnzyme}')
]

# Input: Scatter Plot
inPlotBoth = True # Add the secondary enzyme to the scatter plot

# Input: Save Figure
inSavePath = ''
inFigResolution = 600

# Input: Figure Params
inFigSize = (9.5, 8)
inTitle = 'Enzyme Activity'
inTickLength = 4
inLinewidth = 1.5
inTitleSize = 18
inLabelSize = 16
inLabelTickSize = 13


def normalizeData():
    for key in inData.keys():
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
        path = os.path.join(inSavePath, 'activity_table.png')
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

    # Plot bar graph
    fig, ax = plt.subplots(figsize=inFigSize)
    ax.bar(xTicks - barWidth / 2, y1, barWidth, label=e1,
           color='#F8971F', edgecolor='black', linewidth=inLinewidth)
    ax.bar(xTicks + barWidth / 2, y2, barWidth, label=e2,
           color='#BF5700', edgecolor='black', linewidth=inLinewidth)
    plt.title(inTitle, fontsize=inTitleSize, fontweight='bold')
    ax.set_ylabel('Normalized Activity', fontsize=inLabelSize)
    ax.legend(edgecolor='black', fontsize=inLabelTickSize)

    # Set xticks
    ax.set_xticks(xTicks)
    ax.set_xticklabels(substrates, rotation=0)

    # Set yticks
    ax.set_ylim([0, 1.1])

    # Set tick parameters
    ax.tick_params(axis='both', which='major', length=inTickLength,
                   labelsize=inLabelTickSize)
    ax.tick_params(axis='x', labelsize=inLabelTickSize)

    fig.canvas.mpl_connect('key_press_event', pressKey)
    plt.tight_layout()
    plt.show()

    if inSavePath:
        path = os.path.join(inSavePath, 'enzActivity_bars.png')
        fig.savefig(path, dpi=inFigResolution)
        print(f'Saving figure at path:\n'
              f'     {path}\n\n')
    else:
        print(f'The figure was not saved\n\n')


# ========================================================================================
normalizeData()
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
data[f'Pred {inEnzyme}'] = inData[f'Predicted {inEnzyme}']
data[f'Predicted Z {inEnzyme}'] = inData[f'Predicted Z {inEnzyme}']
actR2 = round(r2_score(data[f'Activity {inEnzyme}'],
                       data[f'Pred {inEnzyme}']), inRoundVal)
actZR2 = round(r2_score(data[f'Activity Z {inEnzyme}'],
                        data[f'Predicted Z {inEnzyme}']), inRoundVal)
print(f'{data}\n\n')


# Evaluate the secondary set
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


# Plot data
c1, c2, = '#101010', '#BF5700'
fig, ax = plt.subplots(figsize=inFigSize)
x, y = f'Activity Z {inEnzyme}', f'Predicted Z {inEnzyme}'
x_fit, y_fit, fitCurve = fitData(x=data[x].values, y=data[y].values)
data.plot(
    x=x, y=y, ax=ax, marker='D', linestyle='none', color=c2, legend=f'R² = {fitCurve:.3f}'
)
ax.plot(x_fit, y_fit, color=c2, linestyle='-', linewidth=inLinewidth)
plt.title(f'Enzyme Activity', fontsize=inTitleSize,
             fontweight='bold')

# X Axis
ax.set_xlabel('Experimental Activity Z Scores', fontsize=inLabelSize)
ax.set_xticks(ax.get_xticks()) # Fix the tick positions
ax.set_xticklabels(ax.get_xticklabels(), ha='center', fontsize=inLabelSize-2)

# Y Axis
ax.set_ylabel('Predicted Z Scores', fontsize=inLabelSize)
ax.tick_params(axis='y', labelsize=inLabelSize-2)

# Legend
if inPlotBoth:
    x, y = f'Activity Z {inEnzyme2}', f'Predicted Z {inEnzyme2}'
    x_fit4, y_fit4, fitCurve4 = fitData(x=data2[x].values, y=data2[y].values)
    data2.plot(
        x=x, y=y, ax=ax, marker='o', linestyle='none', color=c1,
        legend=f'R² = {fitCurve4:.3f}'
    )
    ax.plot(x_fit4, y_fit4, color=c1, linestyle='-', linewidth=1.5)

    ax.legend(
        prop=FontProperties(size=10, weight='bold'),
        handles=[
            Line2D([], [], color=c1, marker='o', linewidth=1.5,
                   label=f'{inEnzyme2} R² = {fitCurve4:.3f}'),
            Line2D([], [], color=c2, marker='D', linewidth=inLinewidth,
                   label=f'{inEnzyme} R² = {fitCurve:.3f}'),
        ]
    )
else:
    ax.legend(prop=FontProperties(size=10, weight='bold'), handles=[Line2D(
        [], [], linestyle='None', marker='None',
        label=f'R² {inEnzyme} = {fitCurve:.3f}')], handletextpad=0, handlelength=0
               )

plt.tight_layout()
fig.canvas.mpl_connect('key_press_event', pressKey)
plt.show()

if inSavePath:
    path = os.path.join(inSavePath, 'enzActivity.png')
    fig.savefig(path, dpi=inFigResolution)
    print(f'Saving figure at path:\n'
          f'     {path}\n\n')
else:
    print(f'The figure was not saved\n\n')


# ========================================================================================


# Plot data
c1, c2, = '#A800FF', '#BF5700'
inTitleSize = 16
inLabelSize = 14
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
x, y = f'Activity {inEnzyme}', f'Pred {inEnzyme}'
x_fit, y_fit, fitCurve = fitData(x=data[x].values, y=data[y].values)
data.plot(
    x=x, y=y, ax=ax1, marker='^',linestyle='none', color=c2, label=f'R² = {fitCurve:.3f}'
)
ax1.plot(x_fit, y_fit, color=c2, linestyle='-', linewidth=inLinewidth)

x, y = f'Activity Z {inEnzyme}', f'Predicted Z {inEnzyme}'
x_fit2, y_fit2, fitCurve2 = fitData(x=data[x].values, y=data[y].values)
data.plot(
    x=x, y=y, ax=ax2, marker='D', linestyle='none', color=c2, legend=f'R² = {fitCurve2:.3f}'
)
ax2.plot(x_fit2, y_fit2, color=c2, linestyle='-', linewidth=inLinewidth)
plt.suptitle(inTitle, fontsize=inTitleSize, fontweight='bold')

# X Axis
ax1.set_xlabel('Experimental Activity', fontsize=inLabelSize)
ax2.set_xlabel('Experimental Activity Z Scores', fontsize=inLabelSize)
ax1.set_xticks(ax1.get_xticks()) # Fix the tick positions
ax1.set_xticklabels(ax1.get_xticklabels(), ha='center', fontsize=inLabelSize-2)
ax2.set_xticks(ax2.get_xticks()) # Fix the tick positions
ax2.set_xticklabels(ax2.get_xticklabels(), ha='center', fontsize=inLabelSize-2)

# Y Axis
ax1.set_ylabel('Predicted Activity', fontsize=inLabelSize)
ax2.set_ylabel('Predicted Z Scores', fontsize=inLabelSize)
ax1.tick_params(axis='y', labelsize=inLabelSize-2)
ax2.tick_params(axis='y', labelsize=inLabelSize-2)

# Legend
if inPlotBoth:
    x, y = f'Activity {inEnzyme2}', f'Pred {inEnzyme2}'
    x_fit3, y_fit3, fitCurve3 = fitData(x=data2[x].values, y=data2[y].values)
    data2.plot(
        x=x, y=y, ax=ax1, marker='*', linestyle='none', color=c1,
        legend=f'R² = {fitCurve3:.3f}'
    )
    ax1.plot(x_fit3, y_fit3, color=c1, linestyle='--', linewidth=inLinewidth)

    x, y = f'Activity Z {inEnzyme2}', f'Predicted Z {inEnzyme2}'
    x_fit4, y_fit4, fitCurve4 = fitData(x=data2[x].values, y=data2[y].values)
    data2.plot(
        x=x, y=y, ax=ax2, marker='o', linestyle='none', color=c1,
        legend=f'R² = {fitCurve4:.3f}'
    )
    ax2.plot(x_fit4, y_fit4, color=c1, linestyle='--', linewidth=inLinewidth)

    ax1.legend(prop=FontProperties(size=10, weight='bold'), handles=[Line2D(
        [], [], linestyle='None', marker='None',
        label=f'R² {inEnzyme2}  = {fitCurve3:.3f}\nR² {inEnzyme} = {fitCurve:.3f}')],
               handletextpad=0, handlelength=0
               )
    ax2.legend(prop=FontProperties(size=10, weight='bold'), handles=[Line2D(
        [], [], linestyle='None', marker='None',
        label=f'R² {inEnzyme2}  = {fitCurve4:.3f}\nR² {inEnzyme} = {fitCurve2:.3f}')],
               handletextpad=0, handlelength=0
               )
else:
    ax1.legend(prop=FontProperties(size=10, weight='bold'), handles=[Line2D(
        [], [], linestyle='None', marker='None',
        label=f'R² {inEnzyme} = {fitCurve:.3f}')], handletextpad=0, handlelength=0
               )
    ax2.legend(prop=FontProperties(size=10, weight='bold'), handles=[Line2D(
        [], [], linestyle='None', marker='None',
        label=f'R² {inEnzyme} = {fitCurve2:.3f}')], handletextpad=0, handlelength=0
               )

plt.tight_layout()
fig.canvas.mpl_connect('key_press_event', pressKey)
plt.show()

if inSavePath:
    path = os.path.join(inSavePath, 'enzActivity_norm+zscores.png')
    fig.savefig(path, dpi=inFigResolution)
    print(f'Saving figure at path:\n'
          f'     {path}\n\n')
else:
    print(f'The figure was not saved\n\n')