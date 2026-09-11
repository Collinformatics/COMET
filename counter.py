import pandas as pd

def counter(args):
    batch, columns, index = args
    matrixLocal = pd.DataFrame(0, index=index, columns=columns)
    for k, v in batch.items():  # batch is a dict
        for i, aa in enumerate(k):
            col = matrixLocal.columns[i]
            if aa in matrixLocal.index:
                matrixLocal.loc[aa, col] += v
    return matrixLocal