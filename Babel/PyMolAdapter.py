# -*- coding: utf-8 -*-

from pymol import cmd, selector, stored
import pprint
import math
import copy
import operator

def csm( sel, csm_type, options): 
  """
  Compute and print the csm of the selected protein
  """
  
  # Read the molecular coordinates
  mol = cmd.identify(sel1, 1)[0][0]
    
  stored.elem = []
  stored.mol = []
  
  cmd.iterate_state(1, mol, "stored.mol.append([x,y,z])")  
  cmd.iterate(1, mol, "stored.elem.append((elem))")  
  
  # call the c function to compute csm and return data
  # first version includes no connectivity data
  [csmVal] = computeCsm(stored.mol, stored.elem, csm_type, options)
  
  print csm

cmd.extend("csm", csm)
  

