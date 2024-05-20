#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import operator
import math
cost_type = [["shortest","passive"],["shortestfastest","passive"]]
names = [ ["ht09", False, "black", "ht"],  ["dblp", False, "crimson","dbl"]  , ["primate",False, "blue","pri"], ["racoon", False, "brown","rac"], ["weaver", False, "pink","wea"], ["employees", False, "red","emp"], ["voles", False, "black", "vol"], ["workplace_2013", False, "green", "wp" ] ]


# In[ ]:


import subprocess
import time
from pathlib import Path
import os
fold ="time_computation/" 
Path(fold).mkdir(parents=True, exist_ok=True)
which=["dij", "bfs", "bell"]
path = "../build.build/"
path_fold = "databases/"
names.reverse()
d = { e[0] : { tuple(c):{ w:0 for w in which  }  for c in cost_type  }  for e in names   }
for e in names:
    print(e)
    for w in which:
        for c in cost_type:
            current = e[0]+"_"+c[0]+"_"+c[1]+"_"+w
            if not os.path.isfile(current):
                st = time.time()
                if w == "dij":
                    l = ["timeout", "2h",  path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-o" ]
                    if e[1] == True:
                        l.append("-d")
                    subprocess.run(l) 
                elif w == "bell":
                    l = ["timeout", "2h",  path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-o", "-z" ]
                    if e[1] == True:
                        l.append("-d")
                    subprocess.run(l)
                else:
                    l = ["timeout", "2h",  path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-b" ]
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




