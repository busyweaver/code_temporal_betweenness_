#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import operator
import math
cost_type = [["shortest","passive"],["shortestfastest","passive"], ["shortest","active"],["shortestfastest","active"] ]
names = [ ["grenoble", False, "crimson","gre"] , ["rennes", False, "green","ren"] , ["belfast", False, "red","bel"], ["kuopio",False, "blue","kuo"], ["primaryschool", False, "brown","prim"], ["highschool_2011",False,"olive","hs11"], ["highschool_2012", False, "pink","hs12"], ["hospital_ward", False, "black","hp"], ["ht09", False, "orange","ht"], ["workplace_2013",False,"purple","wp"] ]


# In[ ]:


import subprocess
import time
from pathlib import Path
import os
fold ="time_computation" 
Path(fold).mkdir(parents=True, exist_ok=True)
which=["dij", "bell", "bfs"]
path = "../build.build/"
path_fold = "examples/"
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
                    subprocess.run([path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-o" ]) 
                elif w == "bell":
                    subprocess.run([path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-o", "-z" ])
                else:
                    subprocess.run([path+"btwBenchmark", "-f", path_fold+e[0]+".csv","-c",c[0],"-y",c[1], "-b" ]) 
                ed = time.time()
                tot = (ed - st)
                d[e[0]][tuple(c)][w] = tot
                with open(fold+current,"w") as f:
                    f.write(str(tot))
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




