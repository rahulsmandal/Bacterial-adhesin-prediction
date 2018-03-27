import numpy as np
import sys
from scipy import interp
import pylab as pl
import pandas as pd
from sklearn.cross_validation import train_test_split
from pydpi.pypro import PyPro # from Prediction_protocol.py
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.metrics import *
from sklearn import preprocessing
import matplotlib.image as mpimg
from sklearn.ensemble import RandomForestClassifier

Selected_Descriptor =['APAAC9','APAAC5','_PolarityD2001','413','T','652','APAAC16','_SecondaryStrD2001','_HydrophobicityD3050','_ChargeC3','_SolventAccessibilityC2','GearyAuto_Hydrophobicity2','PAAC23','_PolarityC1','_SolventAccessibilityD1025','R','S','_SecondaryStrC1']
class DC_CLASS(object):
    
    def full_prot(self, pro_in): 
        self.sequences = []
        self.hdr = []
        self.seq = None
        self.pro_in = pro_in

        pf = open(self.pro_in)
        lines = pf.readlines()

        for line in lines:
            line = line.rstrip()

            if line.startswith(">"):
            	self.hdr.append(line)
                if self.seq: 
                    self.sequences.append(self.seq)
                    
                self.seq = ''
            else:
                self.seq += line

        if self.seq:
            self.sequences.append(self.seq)
        return self.sequences, self.hdr
                
    def Decriptor_generator(self, ps):

        protein = PyPro()
        protein.ReadProteinSequence(ps)
        DS_1 = protein.GetAAComp()
        DS_2 = protein.GetDPComp()
        DS_4 = protein.GetTriad()
        DS_5 = protein.GetPAAC(lamda=5,weight=0.5)
        DS_6 = protein.GetAPAAC(lamda=5,weight=0.5)
        DS_7 = protein.GetCTD()
        DS_8 = protein.GetGearyAuto()
        DS_9 = protein.GetMoranAuto()
        DS_10 = protein.GetMoreauBrotoAuto()
        DS_11 = protein.GetQSO()
        DS_12 = protein.GetSOCN()

        DS_ALL = {}
        
        for DS in (DS_1,DS_2,DS_4,DS_5,DS_6,DS_7,DS_8,DS_9,DS_10,DS_11,DS_12):
            DS_ALL.update(DS)
        return DS_ALL

    def Return_DF(self,f_file):
        
        self.f_file = f_file
        values = []
        seql, _ = DC_CLASS().full_prot(self.f_file)
        for i, seq in enumerate(seql):
            try:
                value = DC_CLASS().Decriptor_generator(seq)
            except:
				print "error in sequence: "
				pass 
            values.append(value)
        
        df1 = pd.DataFrame(values) 
        return df1, _
	
        
    def main_p(self,f):
        self.f = f 
        df1, _ =  DC_CLASS().Return_DF(f)
        new_df = df1[Selected_Descriptor]
        return new_df, _	

class PRED_CLASS(object):

    def data_gen(self,csv_path):
        df = pd.read_csv(csv_path)
        clm_list = []
        clm_list = list(df.columns)
        X_data = df[clm_list[0:len(clm_list)-1]].values
        y_data = df[clm_list[len(clm_list)-1]].values
        return X_data, y_data, clm_list   
    
    def Prediction(self,xdata,ydata,xtest,alg):
        alg.fit(xdata,ydata)  
        l = pd.read_csv('reduced_Data_for_Query_norm.csv')
        full_df = l.append(xtest)
        xtest_p=preprocessing.scale(full_df) #Normalization of input test descriptors
        xtest_p=pd.DataFrame(xtest_p)
        xtest = xtest_p[300:300+xtest.shape[0]]
        b = alg.predict(xtest) 
        pb = alg.predict_proba(xtest) 
        return b, pb
    
    def main_process(self,query_Seq):
        RF = RandomForestClassifier(n_estimators=10)
        LR = LogisticRegression(C = 1e10, solver = 'lbfgs')
        GNB = GaussianNB()
        KNB = KNeighborsClassifier()
        DT = DecisionTreeClassifier()
        SV = SVC(probability=True)
        train_x, train_y,train_l = PRED_CLASS().data_gen("train.csv") 
        test_x, _ = DC_CLASS().main_p(query_Seq)
        preO, proO = PRED_CLASS().Prediction(train_x, train_y, test_x,SV)
        return _ , preO,proO

if __name__=="__main__":

    _, predicted, probability = PRED_CLASS().main_process("sequence.fasta") #input sequence file name

    for x, i in enumerate(range(0,len(predicted))):
        print predicted[i],"---->", round(probability[i][1],2), _[x]
