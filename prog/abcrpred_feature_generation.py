#!/usr/bin/python
from __future__ import print_function
import getopt
import sys
import os
import numpy as np
import pandas as pd
import math
import itertools
from itertools import repeat
import argparse
import csv
from collections import Counter
import re
import glob
import time
from time import sleep
from argparse import RawTextHelpFormatter
import uuid
import warnings
warnings.filterwarnings("ignore") 
std = list("ACDEFGHIKLMNPQRSTVWY")
PCP= pd.read_csv('Data/PhysicoChemical.csv', header=None)
AAindices = 'Data/aaind.txt'
AAIndex = pd.read_csv('Data/aaindex.csv',index_col='INDEX');
AAIndexNames = pd.read_csv('Data/AAIndexNames.csv',header=None);
dir_1 = sys.argv[3]
###########################BTC###############
def bond(file) :
    tota = []
    hy = []
    Si = []
    Du = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []
    bb = pd.DataFrame()
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    bonds=pd.read_csv("Data/bonds.csv", sep = ",")
    for i in range(0,len(df)) :
        tot = 0
        h = 0
        S = 0
        D = 0
        tota.append([i])
        hy.append([i])
        Si.append([i])
        Du.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    tot = tot + bonds.iloc[:,1][k]
                    h = h + bonds.iloc[:,2][k]
                    S = S + bonds.iloc[:,3][k]
                    D = D + bonds.iloc[:,4][k]
        tota[i].append(tot)
        hy[i].append(h)
        Si[i].append(S)
        Du[i].append(D)
    for m in range(0,len(df)) :
        b1.append(tota[m][1])
        b2.append(hy[m][1])
        b3.append(Si[m][1])
        b4.append(Du[m][1])

    bb["BTC_T"] = b1
    bb["BTC_H"] = b2
    bb["BTC_S"] = b3
    bb["BTC_D"] = b4

    bb.to_csv(dir_1+'/sam_allcomp.btc', index=None, encoding="utf-8")
#################################DDOR####################################
def DDOR(file) :
    df = pd.read_csv(file, header = None)
    df1 = pd.DataFrame(df[0].str.upper())
    f = open(dir_1+'/sam_allcomp.ddr','w')
    sys.stdout = f
    for i in std:
        print('DDR_'+i, end=",")
    print("")
    for i in range(0,len(df1)):
        s = df1[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            print("%.2f"%zz2,end=",")
        print("")
    f.truncate()
######################################AAC###################################
def ctd(file):
    attr=pd.read_csv("Data/aa_attr_group.csv", sep="\t")
    filename, file_extension = os.path.splitext(file)
    df1 = pd.read_csv(file, header = None)
    df = pd.DataFrame(df1[0].str.upper())
    n = 0
    stt1 = []
    m = 1
    for i in range(0,len(attr)) :
        st =[]
        stt1.append([])
        for j in range(0,len(df)) :
            stt1[i].append([])
            for k in range(0,len(df[0][j])) :
                while m < 4 :
                    while n < len(attr.iloc[i,m]) :
                        if df[0][j][k] == attr.iloc[i,m][n] :
                            st.append(m)
                            stt1[i][j].append(m)
                        n += 2
                    n = 0
                    m += 1
                m = 1
#####################Composition######################
    f = open(dir_1+"/compout_1", 'w')
    sys.stdout = f
    std = [1,2,3]
    print("1,2,3,")
    for p in range (0,len(df)) :
        for ii in range(0,len(stt1)) :
            #for jj in stt1[ii][p]:
            for pp in std :
                count = 0
                for kk in stt1[ii][p] :
                    temp1 = kk
                    if temp1 == pp :
                        count += 1
                    composition = (count/len(stt1[ii][p]))*100
                print("%.2f"%composition, end = ",")
            print("")
    f.truncate()

#################################Transition#############
    tt = []
    tr=[]
    kk =0
    for ii in range(0,len(stt1)) :
        tt = []
        tr.append([])
        for p in range (0,len(df)) :
            tr[ii].append([])
            while kk < len(stt1[ii][p]) :
                if kk+1 <len(stt1[ii][p]):
                #if  stt1[ii][p][kk] < stt1[ii][p][kk+1] or stt1[ii][p][kk] > stt1[ii][p][kk+1]: # condition for adjacent values
                    tt.append(stt1[ii][p][kk])
                    tt.append(stt1[ii][p][kk+1])
                    tr[ii][p].append(stt1[ii][p][kk])
                    tr[ii][p].append(stt1[ii][p][kk+1])

                kk += 1
            kk = 0

    pp = 0
    xx = []
    xxx = []
    for mm in range(0,len(tr)) :
        xx = []
        xxx.append([])
        for nn in range(0,len(tr[mm])):
            xxx[mm].append([])
            while pp < len(tr[mm][nn]) :
                xx .append(tr[mm][nn][pp:pp+2])
                xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                pp+=2
            pp = 0

    f1 = open(dir_1+"/compout_2", 'w')
    sys.stdout = f1
    std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
    print("1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->3,")
    for rr in range(0,len(df)) :
        for qq in range(0,len(xxx)):
            for tt in std1 :
                count = 0
                for ss in xxx[qq][rr] :
                    temp2 = ss
                    if temp2 == tt :
                        count += 1
                print(count, end = ",")
            print("")
    f1.truncate()

    #################################Distribution#############
    c_11 = []
    c_22 = []
    c_33 = []
    zz = []
    #print("0% 25% 50% 75% 100%")
    for x in range(0,len(stt1)) :
        #c_11.append([])
        c_22.append([])
        #c_33.append([])
        yy_c_1 = []
        yy_c_2 = []
        yy_c_3 = []
        ccc = []

        k = 0
        j = 0
        for y in range(0,len(stt1[x])):
            #c_11[x].append([])
            c_22[x].append([])
            for i in range(1,4) :
                cc = []
                c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                c_22[x][y].append(c1)
    cc = []
    for ss in range(0,len(df)):
        for uu in range(0,len(c_22)):
            for mm in range(0,3):
                for ee in range(0,101,25):
                    k = (ee*(len(c_22[uu][ss][mm])))/100
                    cc.append(math.floor(k))
    f2 = open(dir_1+'/compout_3', 'w')
    sys.stdout = f2
    print("0% 25% 50% 75% 100%")
    for i in range (0,len(cc),5):
        print(*cc[i:i+5])
    f2.truncate()
    head = []
    header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
    for i in header1:
        for j in range(1,4):
            head.append(i+str(j))
    df11 = pd.read_csv(dir_1+"/compout_1")
    df_1 = df11.iloc[:,:-1]
    zz = pd.DataFrame()
    for i in range(0,len(df_1),7):
        zz = zz.append(pd.DataFrame(pd.concat([df_1.loc[i],df_1.loc[i+1],df_1.loc[i+2],df_1.loc[i+3],df_1.loc[i+4],df_1.loc[i+5],df_1.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    zz.columns = head
    #zz.to_csv(filename+".ctd_comp", index=None, encoding='utf-8')
    head2 = []
    header2 = ['CeTD_11','CeTD_12','CeTD_13','CeTD_21','CeTD_22','CeTD_23','CeTD_31','CeTD_32','CeTD_33']
    for i in header2:
        for j in ('HB','VW','PO','PZ','CH','SS','SA'):
            head2.append(i+'_'+str(j))
    df12 = pd.read_csv(dir_1+"/compout_2")
    df_2 = df12.iloc[:,:-1]
    ss = pd.DataFrame()
    for i in range(0,len(df_2),7):
        ss = ss.append(pd.DataFrame(pd.concat([df_2.loc[i],df_2.loc[i+1],df_2.loc[i+2],df_2.loc[i+3],df_2.loc[i+4],df_2.loc[i+5],df_2.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
    ss.columns = head2
    #ss.to_csv(filename+".ctd_trans", index=None, encoding='utf-8')
    head3 = []
    header3 = ['CeTD_0_p','CeTD_25_p','CeTD_50_p','CeTD_75_p','CeTD_100_p']
    header4 = ['HB','VW','PO','PZ','CH','SS','SA']
    for j in range(1,4):
        for k in header4:
            for i in header3:
                head3.append(i+'_'+k+str(j))
    df_3 = pd.read_csv(dir_1+"/compout_3", sep=" ")
    rr = pd.DataFrame()
    for i in range(0,len(df_3),21):
        rr = rr.append(pd.DataFrame(pd.concat([df_3.loc[i],df_3.loc[i+1],df_3.loc[i+2],df_3.loc[i+3],df_3.loc[i+4],df_3.loc[i+5],df_3.loc[i+6],df_3.loc[i+7],df_3.loc[i+8],df_3.loc[i+9],df_3.loc[i+10],df_3.loc[i+11],df_3.loc[i+12],df_3.loc[i+13],df_3.loc[i+14],df_3.loc[i+15],df_3.loc[i+16],df_3.loc[i+17],df_3.loc[i+18],df_3.loc[i+19],df_3.loc[i+20]],axis=0)).transpose()).reset_index(drop=True)
    rr.columns = head3
    cotrdi= pd.concat([zz,ss,rr],axis=1)
    cotrdi.to_csv(dir_1+'/sam_allcomp.ctd', index=None, encoding='utf-8')
##################################QSO#############################################################
def qos(file):
    w = 0.1
    gap = 3
    ff = []
    filename, file_extension = os.path.splitext(file)
    df = pd.read_csv(file, header = None)
    df2 = pd.DataFrame(df[0].str.upper())
    for i in range(0,len(df2)):
        ff.append(len(df2[0][i]))
    if min(ff) < gap:
        print("Error: All sequences' length should be higher than :", gap)
    else:
        mat1 = pd.read_csv("Data/Schneider-Wrede.csv", index_col = 'Name')
        mat2 = pd.read_csv("Data/Grantham.csv", index_col = 'Name')
        s1 = []
        s2 = []
        for i in range(0,len(df2)):
            for n in range(1, gap+1):
                sum1 = 0
                sum2 = 0
                for j in range(0,(len(df2[0][i])-n)):
                    sum1 = sum1 + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                s1.append(sum1)
                s2.append(sum2)
        zz = pd.DataFrame(np.array(s1).reshape(len(df2),gap))
        zz["sum"] = zz.sum(axis=1)
        zz2 = pd.DataFrame(np.array(s2).reshape(len(df2),gap))
        zz2["sum"] = zz2.sum(axis=1)
        c1 = []
        c2 = []
        c3 = []
        c4 = []
        h1 = []
        h2 = []
        h3 = []
        h4 = []
        for aa in std:
            h1.append('QSO'+str(gap)+'_SC_' + aa)
        for aa in std:
            h2.append('QSO'+str(gap)+'_G_' + aa)
        for n in range(1, gap+1):
            h3.append('SC' + str(n))
        h3 = ['QSO'+str(gap)+'_'+sam for sam in h3]
        for n in range(1, gap+1):
            h4.append('G' + str(n))
        h4 = ['QSO'+str(gap)+'_'+sam for sam in h4]
        for i in range(0,len(df2)):
            AA = {}
            for j in std:
                AA[j] = df2[0][i].count(j)
                c1.append(AA[j] / (1 + w * zz['sum'][i]))
                c2.append(AA[j] / (1 + w * zz2['sum'][i]))
            for k in range(0,gap):
                c3.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
                c4.append((w * zz[k][i]) / (1 + w * zz['sum'][i]))
        pp1 = np.array(c1).reshape(len(df2),len(std))
        pp2 = np.array(c2).reshape(len(df2),len(std))
        pp3 = np.array(c3).reshape(len(df2),gap)
        pp4 = np.array(c4).reshape(len(df2),gap)
        zz5 = round(pd.concat([pd.DataFrame(pp1, columns = h1),pd.DataFrame(pp2,columns = h2),pd.DataFrame(pp3, columns = h3),pd.DataFrame(pp4, columns = h4)], axis = 1),4)
        zz5.to_csv(dir_1+'/sam_allcomp.qso', index = None, encoding = 'utf-8')

bond(sys.argv[1])
DDOR(sys.argv[1])
ctd(sys.argv[1])
qos(sys.argv[1])
df5 = pd.read_csv(dir_1+"/sam_allcomp.btc")
df9 = pd.read_csv(dir_1+"/sam_allcomp.ddr")
df14 = pd.read_csv(dir_1+"/sam_allcomp.ctd")
df15 = pd.read_csv(dir_1+"/sam_allcomp.qso")
df19 = pd.concat([df5,df9.iloc[:,:-1],df14,df15],axis=1)
df19.to_csv(sys.argv[2], index=None)
filelist=glob.glob(dir_1+"/sam_allcomp*")
for file_2 in filelist:
    os.remove(file_2)
