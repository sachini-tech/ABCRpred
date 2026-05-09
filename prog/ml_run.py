import re
import sys
import csv
import pandas as pd
import getopt
import sys
import numpy
import pathlib
import os
import numpy as np
import pickle
from sklearn.metrics import auc
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.cluster import KMeans
from sklearn.model_selection import GridSearchCV
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils import shuffle
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
a=[]
file_name = sys.argv[1]
file_name1 = sys.argv[2]
file_name2 = sys.argv[3]
clf = pickle.load(open(file_name2,'rb'))
#X_test, y_test = load_svmlight_file(file_name)
data_test = pd.read_csv(file_name)
#X_test = np.loadtxt(file_name, delimiter=',')
X_test = data_test
#y_test = data_test[:, -1].astype(np.int)
y_p_score1=clf.predict_proba(X_test)
y_p_s1=y_p_score1.tolist()
a.extend(y_p_s1)
df = pd.DataFrame(a)
df1 = df.iloc[:,-1].round(2)
df1.to_csv(file_name1, index=None, header=False)
