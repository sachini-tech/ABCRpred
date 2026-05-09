import sys
import pandas as pd
import numpy as np
df1 = pd.read_csv(sys.argv[1], header=None)
thr = float(sys.argv[2])
aa = []
for i in range(0,len(df1)):
    if df1[0][i] >= thr:
        aa.append("Resistant")
    else:
        aa.append("Susceptible")
df2 = pd.DataFrame(aa)
df2.to_csv(sys.argv[3],index=None, header=False)
