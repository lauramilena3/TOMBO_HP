import pandas as pd
import numpy as np
import csv  
import sys
inFile = sys.argv[1]
name=inFile.split(".")[0]
with open(inFile, 'r') as read_obj:
    csv_reader = csv.reader(read_obj,delimiter =" ")
    list_of_rows = list(csv_reader)
    data=(list_of_rows[2:])
df=pd.DataFrame(list(data))
df.columns=["a", "b"]
df['a']=df['a'].astype(int)
df_pos=pd.DataFrame(list(range(1,(df["a"].max()))))
df_pos.columns=["a"]
df['a']=df['a'].astype(int)
df_pos['a']=df_pos['a'].astype(int)
ones=[-1] * len(list(range(1,df["a"].max())))
df_pos['b']=ones
df3 = pd.concat([df,df_pos])
df3.drop_duplicates(subset=['a'], inplace=True, keep='first')
df3=df3.sort_values('a')
file_object = open(name+"_corrected.wig", 'w')
file_object.write(" ".join(list_of_rows[0]))
file_object.write(" ".join(list_of_rows[1]))
file_object.write("\n")
file_object.close()
df3.to_csv(name+"_corrected.wig",sep=' ', index=False, header=False, mode='a')
