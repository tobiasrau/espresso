#
# Copyright (C) 2013,2014 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
include "myconfig.pxi"
from highlander import *
import numpy as np
# Non-bonded interactions



@highlander
class ElectrostaticInteraction(object):

  def __init__(self, *args, **kwargs):
    """
    Represents an instance of a Electrostatic interaction, such as P3M.
    Required aguments: lb
    """
    self._bjerrum_length=0
    self._params={}
    
    # Bjerrum length and/or tuneing as argument
    second_arg = ""
    if len(args)<=2:
      if len(args)==1:
        self._bejrrum_length=args[0]
      elif len(args)==2:
        if args[0]>0:
          self._bjerrum_length=args[0]
          second_arg=args[1]
        elif args[1]>0:
          self._bjerrum_length=args[1]
          second_arg=args[0]
        else:
          raise ValueError("Invalid value for Bjerrum length.")
          
      if coulomb_set_bjerrum(self._bjerrum_length):
        raise ValueError("Bjerrum_length should be a positive double")

      self._params=self.defaultParams()
 

      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
        self._params[k]=kwargs[k]


      for k in kwargs:
        if k in self.validKeys():
          self._params[k] = kwargs[k]
        else:
          raise KeyError("%s is not a vaild key" %k)


      # Validation of parameters
      self.validateParams()
        
      if second_arg == "tune" or second_arg == "tune_alpha":
        self._tune()
      elif second_arg == "":
        pass
      else:
        raise ValueError("Invalid argument %s" %second_arg)


    else:
      raise Exception("At least the Bjerrum length has to be given as argument.")

   
  def isValid(self):
    """Check, if the data stored in the instance still matches what is in Espresso"""

    # check, if the bond parameters saved in the class still match those saved in Espresso
    tempParams =self._getParamsFromEsCore()
    if self._params != tempParams:
      return False
    
    # If we're still here, the instance is valid
    return True
 
 
  def getParams(self):
    """Get interaction parameters"""
    # If this instance refers to an actual interaction defined in the es core, load
    # current parameters from there
    update=self._getParamsFromEsCore()
    self._params.update(update)    
    return self._bjerrum_length,self._params

  def setParams(self,**p):
    """Update parameters. Only given """
    # Check for Bjerrum length
    if "bjerrum_length" in p.keys():
      self._bjerrum_length = p["bjerrum_length"]
      p.pop("bjerrum_length",None)

    # Check, if any key was passed, which is not known
    for k in p.keys():
      if k not in self.validKeys():
        raise ValueError("Only the following keys are supported: "+self.validKeys().__str__())

    # When an interaction is newly activated, all required keys must be given
    if not self.isActive():
      for k in self.requiredKeys():
        if k not in p:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())

    
    
    self._params.update(p)
    # vaidate updated parameters
    self.validateParams()
    # Put in values given by the user
    self._setParamsInEsCore()
    if coulomb_set_bjerrum(self._bjerrum_length):
      raise ValueError("Bjerrum_length should be a positive double")
   


  def validateParams(self):
    return True

  def _getParamsFromEsCore(self):
    raise Exception("Subclasses of ElectrostaticInteraction must define the _getParamsFromEsCore() method.")
  
  def _setParamsInEsCore(self):
    raise Exception("Subclasses of ElectrostaticInteraction must define the _setParamsFromEsCore() method.")

  def _tune(self):
    raise Exception("Subclasses of ElectrostaticInteraction must define the _tune() method or chosen method does not support tuning.")

  def defaultParams(self):
    raise Exception("Subclasses of ElectrostaticInteraction must define the defaultParams() method.")
  
  def isActive(self):
    # If this instance refers to an actual interaction defined in the es core, load
    raise Exception("Subclasses of ElectrostaticInteraction must define the isActive() method.")

  
  def typeName(self): 
    raise Exception("Subclasses of ElectrostaticInteraction must define the typeName() method.")
  
  def validKeys(self): 
    raise Exception("Subclasses of ElectrostaticInteraction must define the validKeys() method.")
  
  def requiredKeys(self): 
    raise Exception("Subclasses of ElectrostaticInteraction must define the requiredKeys() method.")

  def tuneKeys(self):
    raise Exception("Subclasses of ElectrostaticInteraction must define the tuneKeys() method.")


class P3M(ElectrostaticInteraction):
  def typeNumber(self):
    return 0

  def typeName(self): 
    return "P3M"

  def validateParams(self):
    default_params=self.defaultParams()
    if not (self._params["r_cut"]>=0 or self._params["r_cut"]==default_params["r_cut"]):
      raise ValueError("P3M r_cut has to be >=0")
    if not (isinstance(self._params["mesh"],int) or len(self._params["mesh"])):
      raise ValueError("P3M mesh has to be an integer or integer list of length 3")
    if (isinstance(self._params["mesh"],basestring) and len(self._params["mesh"]) == 3):
      if (self._params["mesh"][0]%2 != 0 and self._params["mesh"][0] != -1) or (self._params["mesh"][1]%2 != 0 and self._params["mesh"][1] != -1) or (self._params["mesh"][2]%2 != 0 and self._params["mesh"][2] != -1):
        raise ValueError("P3M requires an even number of mesh points in all directions")
    if not (isinstance(self._params["cao"],int) and self._params["cao"] >= -1 and self._params["cao"] <= 7):
      raise ValueError("P3M cao has to be an integer between -1 and 7")
    if not (self._params["accuracy"] > 0):
      raise ValueError("P3M accuracy has to be positive")
    if self._params["epsilon"] == "metallic":
      self._params = 0.0
    if p3m_set_eps(self._params["epsilon"]):
      raise ValueError("epsilon should be a double or 'metallic'")
    if self._params["n_interpol"] != default_params["n_interpol"]:
      if p3m_set_ninterpol(self._params["n_interpol"]):
        raise ValueError("n_interpol should be a positive integer")
    if self._params["mesh_off"] != default_params["mesh_off"]:
      if python_p3m_set_mesh_offset(self._params["mesh_off"]):
        raise ValueError("mesh_off should be a list of length 3 and values between 0.0 and 1.0")

    return True

  def validKeys(self):
    return "alpha_L","r_cut_iL","mesh","mesh_off","cao","inter","accuracy","epsilon","cao_cut","a","ai","alpha","r_cut","inter2","cao3","additional_mesh","n_interpol"

  def requiredKeys(self): 
    return ["accuracy"]

  def defaultParams(self):
    return {"cao":-1,\
            "n_interpol":-1,\
            "r_cut":-1,\
            "accuracy":-1,\
            "mesh":[-1,-1,-1],\
            "epsilon":0.0,\
            "mesh_off":[-1,-1,-1]}

  def _getParamsFromEsCore(self):
    cdef p3m_data_struct p3m
    p3m = p3m_get_params()
    return p3m.params

  def _setParamsInEsCore(self):
    python_p3m_set_params(self._params["r_cut"],self._params["mesh"],self._params["cao"],self._params["alpha"],self._params["accuracy"],self._params["n_interpol"])
    

  def _tune(self):
    python_p3m_set_tune_params(self._params["r_cut"], self._params["mesh"], self._params["cao"], -1.0, self._params["accuracy"], self._params["n_interpol"])
    resp,log=python_p3m_adaptive_tune()
    if resp:
      raise Exception("failed to tune P3M parameters to required accuracy")
    print log
    p3m_params = self._getParamsFromEsCore()
    self._params.update(p3m_params)

  def tuneKeys(self):
    return "r_cut","mesh","cao","alpha","accuracy","n_interpol"

  def isActive(self):
    return True

# class P3M_GPU(ElectrostaticInteraction):
#   def __init__(self):
#     print "P3M_GPU"


# @highlander
# class HydrodynamicInteraction(object):
#   def __init__(self):
#     print "Hydrodynamic"


# class LB(HydrodynamicInteraction):
#   def __init__(self):
#     print "LB"


# class LB_GPU(HydrodynamicInteraction):
#   def __init__(self):
#     print "LB_GPU"
