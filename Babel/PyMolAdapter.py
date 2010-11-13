# -*- coding: utf-8 -*-

from pymol import cmd, selector, stored
import pprint
import math
import copy
import operator

from csm import computeCsm

import pickle

def csmFunc( sel, csm_type): 
  """
  Compute and print the csm of the selected protein
  """

  print "Running CSM"
  
  # Read the molecular coordinates
  mol = cmd.identify(sel, 1)[0][0]
    
  stored.elems = []
  stored.mol = []
  
  cmd.iterate_state(1, mol, "stored.mol.append((x,y,z))")  
  cmd.iterate(mol, "stored.elems.append(elem)")

  options = ["-findperm"]
#  f = open("ttt.ppp","w")
#  params = {}
#  params["mol"] = stored.mol
#  params["elems"] = stored.elems
#  params["csm_type"] = csm_type
#  params["options"] = options
#  pickle.dump(params,f)
#  f.close()
  
  
  # call the c function to compute csm and return data
  # first version includes no connectivity data
  [csmVal] = computeCsm(stored.mol, stored.elems, csm_type, options)

#  print "Fake Run"
  print csmVal

cmd.extend("csm", csmFunc)


