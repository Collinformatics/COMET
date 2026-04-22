import pandas as pd

def counter(args):
    batch, columns, index = args
    matrixLocal = pd.DataFrame(0, index=index, columns=columns)
    for k, v in batch.items():  # batch is a dict
        for i, a in enumerate(k, start=1):
            col = f'R{i}'
            if col in matrixLocal.columns and a in matrixLocal.index:  # was matrix.index
                matrixLocal.loc[a, col] += v
    return matrixLocal