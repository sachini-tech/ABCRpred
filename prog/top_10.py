import pandas as pd
import numpy as np
import sys
df1 = pd.read_csv(sys.argv[1])
cc = list(df1.feature)
df2 = pd.read_csv(sys.argv[2])
df3 = df2[cc]
df3.to_csv(sys.argv[3],index=None)
