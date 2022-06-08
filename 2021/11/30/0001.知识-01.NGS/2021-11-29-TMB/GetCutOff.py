import pandas as pd
import numpy as np

Data = pd.read_excel("tTMB-2019.09.15-2022.01.27.xlsx")

CancerType = Data["癌种"].drop_duplicates().tolist()

data = {'Name':['cancer', 'Sample_Num', 'Percentile_90', 'Percentile_80']}

Result = ""# pd.DataFrame(data)
for cancer in CancerType:
    cancer_data = Data[Data["癌种"]==cancer]
    TMB =list(map(float,cancer_data['TMB']))
    Percentile_90 = np.percentile(TMB,90)
    Percentile_80 = np.percentile(TMB,80)
    Result = Result + cancer + "\t" + str(len(cancer_data))+"\t"+str(Percentile_90)+"\t"+str(Percentile_80) + "\n"
    print (cancer + "\t" + str(len(cancer_data))+"\t"+str(Percentile_90)+"\t"+str(Percentile_80))
    #Result.loc[len(Result)+1] = [cancer,len(cancer_data),Percentile_90,Percentile_80]
#np.percentile()

print(Result)