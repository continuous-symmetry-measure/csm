# -*- coding: utf-8 -*-

from pymol import cmd, selector, stored
from pymol.cgo import *    # get constants
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

  #store stuff
  stored.csm_sel = sel
  stored.csm_type = csm_type
  stored.elems = []
  stored.mol = []
  stored.chain = []
  stored.bonds = []
  
  cmd.iterate_state(1, mol, "stored.mol.append((x,y,z))")  
  cmd.iterate(mol, "stored.elems.append(elem)")
  cmd.iterate(mol, "stored.chain.append(ID)")

  # prepare bond structure
  for i in xrange(len(stored.elems)):
    stored.bonds.append([])

  atomModel = cmd.get_model(sel)
  for b in atomModel.bond:
    stored.bonds[b.index[0]].append(b.index[1])
    stored.bonds[b.index[1]].append(b.index[0])

  options = []
  if (len(stored.elems) > 20):
    options.append("-findperm")

  #  f = open("ttt.ppp","w")
#  params = {}
#  params["mol"] = stored.mol
#  params["elems"] = stored.elems
#  params["csm_type"] = csm_type
#  params["options"] = options
#  params["bonds"] = stored.bonds
#  pickle.dump(params,f)
#  f.close()
  # call the c function to compute csm and return data
  # first version includes no connectivity data
  [csmVal, atomicCSM,symElement, outAtomPos] = computeCsm(stored.mol, stored.elems, stored.bonds, csm_type, options)
  
  stored.csm_center = [0.0,0.0,0.0]
  stored.csm_center[0] = float(sum([x[0] for x in stored.mol])) / len(outAtomPos)
  stored.csm_center[1] = float(sum([x[1] for x in stored.mol])) / len(outAtomPos)
  stored.csm_center[2] = float(sum([x[2] for x in stored.mol])) / len(outAtomPos)
  stored.csm_norm = float(sum([((x[0]-stored.csm_center[0])**2 +
                          (x[1]-stored.csm_center[1])**2 +
                          (x[2]-stored.csm_center[2])**2) for x in stored.mol])) / len(outAtomPos)
  stored.csm_norm = math.sqrt(stored.csm_norm)
  for i in xrange(len(outAtomPos)):
    newTup = (outAtomPos[i][0] + stored.csm_center[0],
              outAtomPos[i][1] + stored.csm_center[1],
              outAtomPos[i][2] + stored.csm_center[2])
    outAtomPos[i] = newTup

  source_obj = sel
  new_object = source_obj+"_"+csm_type+"_sym"
  cmd.copy(new_object, source_obj)
  
  stored.atomicCSM = atomicCSM
  stored.outAtomPos = outAtomPos
  stored.symElement = symElement

  cmd.alter_state(1, new_object, "(x,y,z)=stored.outAtomPos.pop(0)")  
  stored.outAtomPos = outAtomPos
  
  
  print csmVal

def csm_map():
  normalizedCsm = [(x/max(stored.atomicCSM)) for x in stored.atomicCSM]
  for i in xrange(len(normalizedCsm)):
    color_name = "color_" + str(i)
    cmd.set_color(color_name, [0.89+(normalizedCsm[i]/10),1-normalizedCsm[i],1-normalizedCsm[i]])
#    cmd.set_color(color_name, [normalizedCsm[i],0.0,1-normalizedCsm[i]])
#    print [normalizedCsm[i],0.0,1-normalizedCsm[i]]
    cmd.color(color_name, "("+ stored.csm_sel + " and ID " + str(i + 1) + ")")    

def csm_element(r=1.0, g=0.01, b=0, width=5.0):
  i = stored.symElement[0]
  j = stored.symElement[1]
  k = stored.symElement[2]
  x = stored.csm_center[0]
  y = stored.csm_center[1]
  z = stored.csm_center[2]  
  length = stored.csm_norm 
  x2,y2,z2 = x+i*length,y+j*length,z+k*length  
  obj = [
    LINEWIDTH, width,
    BEGIN, LINES,
    COLOR,   r,  g,  b,
    VERTEX, x, y, z,
    VERTEX, x2, y2, z2,
    END
    ]
  cmd.load_cgo(obj,stored.csm_sel + "_" + stored.csm_type + '_axis')
  
  
cmd.extend("csm", csmFunc)
cmd.extend("csm_map", csm_map)
cmd.extend("csm_element", csm_element)



