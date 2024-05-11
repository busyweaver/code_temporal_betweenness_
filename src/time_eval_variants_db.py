#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import operator
import math
cost_type = [["shortest","passive"],["shortestfastest","passive"]]
names = [ ["bison_dire", True, "crimson","bis"] , ["cattle_dire", True, "green","cat"] , ["primate",False, "blue","pri"], ["racoon", False, "brown","rac"], ["sheep_dire",True,"olive","she"], ["weaver", False, "pink","wea"]]
#names = [["employees", False, "red","emp"], ["voles", False, "black", "vol"]]


# In[ ]:


import subprocess
import time
from pathlib import Path
import os
fold ="time_computation/" 
Path(fold).mkdir(parents=True, exist_ok=True)
which=["dij", "bell", "bfs"]
path = "../build.build/"
path_fold = "databases/"
names.reverse()
d = { e[0] : { tuple(c):{ w:0 for w in which  }  for c in cost_type  }  for e in names   }
for e in names:
    print(e)
    for c in cost_type:    
        for w in which:
            current = e[0]+"_"+c[0]+"_"+c[1]+"_"+w
            if not os.path.isfile(fold+current):
                st = time.time()
                if w == "dij":
                    l = [path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-o" ]
                    if e[1] == True:
                        l.append("-d")
                    subprocess.run(l) 
                elif w == "bell":
                    l = [path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-o", "-z" ]
                    if e[1] == True:
                        l.append("-d")
                    subprocess.run(l)
                else:
                    l = [path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-b" ]
                    if e[1] == True:
                        l.append("-d")
                    subprocess.run(l) 
                ed = time.time()
                tot = (ed - st)
                d[e[0]][tuple(c)][w] = tot
                #with open(fold+current,"w") as f:
                #    f.write(str(tot))
            else:
                print("already here")


# In[ ]:


import pickle

def save_dic(d,s):
    with open(s+'.pickle', 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)


# In[ ]:


save_dic(d,"time_comput")


# In[ ]:




