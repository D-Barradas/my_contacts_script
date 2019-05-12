
# coding: utf-8

# In[1]:


import os, sys
from Bio.PDB import PDBParser
from optparse import OptionParser


# In[2]:
##Parse the options
usage = "USAGE: python my_contacts_script.py --f1 FirstPDB --c1 firts chains  --c2 second chains \n"
parser = OptionParser(usage=usage)

##options
parser.add_option("--f1",help="First molecule pdb", dest="f1")
parser.add_option("--c1",help="First chains", dest="c1")
parser.add_option("--c2",help="Second chains", dest="c2")

(options, args) = parser.parse_args()
if (options.f1 and options.c1 and options.c2):
    print "%s"%(options.f1.split(".")[0]),

else:

    print "Not enough input arguments supplied"
    print usage

#str_1 = PDBParser().get_structure('first_one', 'F_1N2C_9998_I.pdb') # load your molecule
str_1 = PDBParser(QUIET=True).get_structure('first_one', options.f1) # load your molecule


# In[3]:


CHAIN_A , CHAIN_B = [],[]
for model in str_1:
    for chain in model:
       # print chain.id
        for residue in chain:
            for atom in residue:
                #if chain.id == "A":
                if chain.id in options.c1 :
                    x,y,z =atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]
                    #print x,y,z
                    CHAIN_A.append((x,y,z))
                #if chain.id == "B":
                if chain.id in options.c2 :

                    x,y,z =atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]
                    #print x,y,z
                    CHAIN_B.append((x,y,z))
                    
                    #vector1 = atom1.get_vector()


# In[4]:


import numpy as np


# In[5]:


def dist(i,j):
    dx = i[0] -j[0]
    dy = i[1] - j[1]
    dz = i[2] - j[2]
    d = np.sqrt(dx**2 + dy**2 + dz**2)
    return d


# In[6]:


counter = 0
for w in CHAIN_A:
    for f in CHAIN_B:
        my_d = dist(w,f)
        if my_d <=5 and my_d >= 3:
           # print my_d 
            counter +=1 


# In[7]:


print counter

