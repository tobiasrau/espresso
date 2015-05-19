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



cdef class NonBondedInteraction(object):
  
  cdef public object _partTypes
  cdef object _params

  def __init__(self, *args, **kwargs):
    """Represents an instance of a non-bonded interaction, such as lennard jones
    Either called with two particle type id, in which case, the interaction 
    will represent the bonded interaction as it is defined in Espresso core
    Or called with keyword arguments describing a new interaction."""
    
    
    # Interaction id as argument
    if len(args)==2 and isinstance(args[0],int) and isinstance(args[1],int):
      self._partTypes=args
      
      # Load the parameters currently set in the Espresso core
      self._params=self._getParamsFromEsCore()
    
    # Or have we been called with keyword args describing the interaction
    elif len(args)==0:
      # Initialize default values
      self._params=self.defaultParams()
      self._partTypes=[-1,-1]
      
      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
      
      self._params = kwargs
      
      # Validation of parameters
      self.validateParams()
      
    else: 
      raise Exception("The constructor has to be called either with two particle type ids (as interger), or with a set of keyword arguments describing a new interaction")

   
      
   


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
    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._params=self._getParamsFromEsCore()
    
    return self._params

  def setParams(self,**p):
    """Update parameters. Only given """
    # Check, if any key was passed, which is not known
    for k in p.keys():
      if k not in self.validKeys():
        raise ValueError("Only the following keys are supported: "+self.validKeys().__str__())

    # When an interaction is newly activated, all required keys must be given
    if not self.isActive():
      for k in self.requiredKeys():
        if k not in p:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
    
    # If this instance refers to an interaction defined in the espresso core,
    # load the parameters from there

    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._params=self._getParamsFromEsCore()
    
    # Put in values given by the user
    self._params.update(p)
    
    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._setParamsInEsCore()

  def validateParams(self):
    return True

  def _getParamsFromEsCore(self):
    raise Exception("Subclasses of NonBondedInteraction must define the _getParamsFromEsCore() method.")
  
  def _setParamsInEsCore(self):
    raise Exception("Subclasses of NonBondedInteraction must define the setParamsFromEsCore() method.")
  
  def defaultParams(self):
    raise Exception("Subclasses of NonBondedInteraction must define the defaultParams() method.")
  
  def isActive(self):
    # If this instance refers to an actual interaction defined in the es core, load
    # current parameters from there
    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._params=self._getParamsFromEsCore()
    raise Exception("Subclasses of NonBondedInteraction must define the isActive() method.")

  
  def typeName(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the typeName() method.")
  
  def validKeys(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the validKeys() method.")
  
  def requiredKeys(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the requiredKeys() method.")

# Lennard Jones

cdef class LennardJonesInteraction(NonBondedInteraction):
  if LENNARD_JONES == 1:
      def validateParams(self):
        if self._params["epsilon"]<0:
          raise ValueError("Lennard-Jones eps has to be >=0")
        if self._params["sigma"]<0:
          raise ValueError("Lennard-Jones sigma has to be >=0")
        if self._params["cutoff"]<0:
          raise ValueError("Lennard-Jones cutoff has to be >=0")
        return True
    
      def _getParamsFromEsCore(self):
        cdef IA_parameters* iaParams
        iaParams =  get_ia_param(self._partTypes[0],self._partTypes[1]) 
        return { \
          "epsilon": iaParams.LJ_eps, \
          "sigma": iaParams.LJ_sig, \
          "cutoff": iaParams.LJ_cut, \
          "shift": iaParams.LJ_shift, \
          "offset": iaParams.LJ_offset, \
          "min": iaParams.LJ_min }
           
    
      def isActive(self):
        return (self._params["epsilon"] >0)
      
      def _setParamsInEsCore(self):
        # Handle the case of shift="auto"
        if self._params["shift"]=="auto": 
          # Calc shift
          self._params["shift"]= -( (self._params["sigma"]/self._params["cutoff"])**12 - (self._params["sigma"]/self._params["cutoff"])**6 )
        
        if lennard_jones_set_params(self._partTypes[0],self._partTypes[1],\
                                            self._params["epsilon"], \
                                            self._params["sigma"], \
                                            self._params["cutoff"], \
                                            self._params["shift"], \
                                            self._params["offset"], \
    					0.0, \
                                            self._params["min"]):
          raise Exception("Could not set Lennard Jones parameters")					
      
      def defaultParams(self):
        self._params={\
          "epsilon":0.,\
          "sigma":0.,\
          "cutoff":0.,\
          "shift":0.,\
          "offset":0.,\
          "min":0.}
    
      def typeName(self): 
        return "LennardJones" 
      
      def validKeys(self): 
        return "epsilon","sigma","cutoff","shift","offset","min"
      
      def requiredKeys(self): 
        return "epsilon","sigma","cutoff","shift" 


class NonBondedInteractionHandle(object):
  """Provides access to all Non-bonded interactions between 
  two particle types.""" 

  type1=-1
  type2=-1

  # Here, one line per non-bonded ia
  lennardJones=None


  def __init__(self, _type1, _type2):
    """Takes two particle types as argument"""
    if not (isinstance(_type1,int) and isinstance(_type2,int)):
      raise TypeError("The particle types have to be of type integer.")
    self.type1=_type1
    self.type2=_type2
    
    
    # Here, add one line for each nonbonded ia
    self.lennardJones =LennardJonesInteraction(_type1,_type2)

  
  


cdef class NonBondedInteractions:
  """Access to non-bonded interaction parameters via [i,j], where i,j are particle 
  types. Returns NonBondedInteractionHandle.
  Also: access to force capping
  """
  def __getitem__(self,key):
    if not isinstance(key,tuple):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    if len(key) != 2 or (not isinstance(key[0],int)) or (not isinstance(key[1],int)):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    return NonBondedInteractionHandle(key[0],key[1])
    
  def setForceCap(self,cap):
   if forcecap_set_params(cap):
     raise Exception("Could not set forcecap")

  def getForceCap(self):
    return force_cap






cdef class BondedInteraction(object):
  def __init__(self, *args, **kwargs):
    """Either called with an interaction id, in which case, the interaction will represent
       the bonded interaction as it is defined in Espresso core
       Or called with keyword arguments describing a new interaction."""
    # Interaction id as argument
    if len(args)==1 and isinstance(args[0],int):
      bondId=args[0]
      # Check, if the bond in Espresso core is really defined as a FENE bond
      if bonded_ia_params[bondId].type != self.typeNumber():
        raise Exception("The bond with this id is not defined as a "+self.typeName()+" bond in the Espresso core.")
 
      self._bondId=bondId

      # Load the parameters currently set in the Espresso core
      self._params=self._getParamsFromEsCore()
      self._bondId=bondId
    
    # Or have we been called with keyword args describing the interaction
    elif len(args)==0:
      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
      
      self.params = kwargs
      
      # Validation of parameters
      self.validateParams()
      
    else: 
      raise Exception("The constructor has to be called either with a bond id (as interger), or with a set of keyword arguments describing a new interaction")

   
      
   


  def isValid(self):
    """Check, if the data stored in the instance still matches what is in Espresso"""
    # Check if the bond type in Espresso still matches the bond type saved in this class
    if bonded_ia_params[self._bondId].type != self.typeNumber():
      return False

    # check, if the bond parameters saved in the class still match those saved in Espresso
    tempParams =self._getParamsFromEsCore()
    if self._params != tempParams:
      return False
    
    # If we're still here, the instance is valid
    return True
 
 
  property params:
    def __get__(self):
      return self._params

    def __set__(self,p):
      # Check, if any key was passed, which is not known
      for k in p.keys():
        if k not in self.validKeys():
          raise ValueError("Only the following keys are supported: "+self.validKeys().__str__)
      
      # Initialize default values
      self.setDefaultParams()
      # Put in values given by the user
      self._params.update(p)

  def validateParams(self):
    return True

  def _getParamsFromEsCore(self):
    raise Exception("Subclasses of BondedInteraction must define the _getParamsFromEsCore() method.")
  
  def _setParamsInEsCore(self):
    raise Exception("Subclasses of BondedInteraction must define the setParamsFromEsCore() method.")
  
  def setDefaultParams(self):
    raise Exception("Subclasses of BondedInteraction must define the setDefaultParams() method.")

  def typeNumber(self): 
    raise Exception("Subclasses of BondedInteraction must define the typeNumber() method.")
   
  
  def typeName(self): 
    raise Exception("Subclasses of BondedInteraction must define the typeName() method.")
  
  def validKeys(self): 
    raise Exception("Subclasses of BondedInteraction must define the validKeys() method.")
  
  def requiredKeys(self): 
    raise Exception("Subclasses of BondedInteraction must define the requiredKeys() method.")



class BondedInteractionNotDefined(BondedInteraction):
  def typeNumber(self):
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

  def typeName(self): 
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

  def validKeys(self):
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

  def requiredKeys(self):
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

  def setDefaultParams(self):
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

  def _getParamsFromEsCore(self):
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

  def _setParamsInEsCore(self):
    raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)



class FeneBond(BondedInteraction):

  def typeNumber(self):
    return 0

  def typeName(self): 
    return "FENE"

  def validKeys(self):
    return "k","d_r_max","r_0"

  def requiredKeys(self): 
    return "k","d_r_max"

  def setDefaultParams(self):
    self._params = {"r_0":0.} 
    # Everything else has to be supplied by the user, anyway

  def _getParamsFromEsCore(self):
    return \
      {"k":bonded_ia_params[self._bondId].p.fene.k,\
       "d_r_max":bonded_ia_params[self._bondId].p.fene.drmax,\
       "r_0":bonded_ia_params[self._bondId].p.fene.r0}

  def _setParamsInEsCore(self):
    fene_set_params(self._bondId,self._params["k"],self._params["d_r_max"],self._params["r_0"])


class HarmonicBond(BondedInteraction):
  def typeNumber(self):
    return 1

  def typeName(self): 
    return "HARMONIC"

  def validKeys(self):
    return "k","r_0","r_cut"

  def requiredKeys(self): 
    return "k","r_0"

  def setDefaultParams(self):
    self._params = {"k'":0.,"r_0":0.,"r_cut":0.} 

  def _getParamsFromEsCore(self):
    return \
      {"k":bonded_ia_params[self._bondId].p.harmonic.k,\
       "r_0":bonded_ia_params[self._bondId].p.harmonic.r,\
       "r_cut":bonded_ia_params[self._bondId].p.harmonic.r_cut}

  def _setParamsInEsCore(self):
    harmonic_set_params(self._bondId,self._params["k"],self._params["r_0"],self._params["r_cut"])
   

class Dihedral(BondedInteraction):
  def typeNumber(self):
    return 5

  def typeName(self): 
    return "DIHEDRAL"

  def validKeys(self):
    return "mult","bend","phase"

  def requiredKeys(self): 
    return "mult","bend", "phase"

  def setDefaultParams(self):
    self._params = {"mult'":1.,"bend":0.,"phase":0.} 

  def _getParamsFromEsCore(self):
    return \
      {"mult":bonded_ia_params[self._bondId].p.dihedral.mult,\
       "bend":bonded_ia_params[self._bondId].p.dihedral.bend,\
       "phase":bonded_ia_params[self._bondId].p.dihedral.phase}

  def _setParamsInEsCore(self):
    dihedral_set_params(self._bondId,self._params["mult"],self._params["bend"],self._params["phase"])
    

IF TABULATED == 1:
  class Tabulated(BondedInteraction):
    def typeNumber(self):
      return 6

    def typeName(self): 
      return "TABULATED"

    def validKeys(self):
      return "type", "filename", "npoints", "minval", "maxval", "invstepsize"

    def requiredKeys(self): 
      return "type", "filename", "npoints", "minval", "maxval", "invstepsize"

    def setDefaultParams(self):
      self._params = {"type":1, "filename":"", "npoints":0, "minval":0, "maxval":1, \
        "invstepsize":1} 

    def _getParamsFromEsCore(self):
      return \
        {"type":bonded_ia_params[self._bondId].p.tab.type,\
         "filename":bonded_ia_params[self.bondID].p.tab.filename,\
         "npoints":bonded_ia_params[self._bondId].p.tab.npoints,\
         "minval":bonded_ia_params[self._bondId].p.tab.minval,\
         "maxval":bonded_ia_params[self._bondId].p.tab.maxval,\
         "invstepsize":bonded_ia_params[self._bondId].p.tab.invstepsize}

    def _setParamsInEsCore(self):
      tabulated_bonded_set_params(self._bondId, self._params["type"], self._params["filename"])

   
IF TABULATED != 1:
  class Tabulated(BondedInteraction):
    def typeNumber(self):
      raise Exception("TABULATED has to be defined in myconfig.hpp.")

    def typeName(self): 
      raise Exception("TABULATED has to be defined in myconfig.hpp.")

    def validKeys(self):
      raise Exception("TABULATED has to be defined in myconfig.hpp.")

    def requiredKeys(self):
      raise Exception("TABULATED has to be defined in myconfig.hpp.")

    def setDefaultParams(self):
      raise Exception("TABULATED has to be defined in myconfig.hpp.")

    def _getParamsFromEsCore(self):
      raise Exception("TABULATED has to be defined in myconfig.hpp.")

    def _setParamsInEsCore(self):
      raise Exception("TABULATED has to be defined in myconfig.hpp.")
   
   
class Subt_Lj(BondedInteraction):
  IF LENNARD_JONES == 1:
    def typeNumber(self):
      return 7

    def typeName(self): 
      return "SUBT_LJ"

    def validKeys(self):
      return "r", "k"

    def requiredKeys(self): 
      return "r", "k"

    def setDefaultParams(self):
      self._params = {"k":0, "r":0} 

    def _getParamsFromEsCore(self):
      return \
        {"k":bonded_ia_params[self._bondId].p.subt_lj.k,\
         "r":bonded_ia_params[self._bondId].p.subt_lj.r}

    def _setParamsInEsCore(self):
      subt_lj_set_params(self._bondId,self._params["k"],self._params["r"])

IF BOND_VIRTUAL == 1:
  class Virtual(BondedInteraction):
    def typeNumber(self):
      return 9

    def typeName(self): 
      return "VIRTUAL"

    def validKeys(self):
      return

    def requiredKeys(self):
      return

    def setDefaultParams(self):
      pass

    def _getParamsFromEsCore(self):
      pass

    def _setParamsInEsCore(self):
      virtual_set_params(self._bondId)

ELSE:
  class Virtual(BondedInteractionNotDefined):
    name="BOND_VIRTUAL"      
   
IF BOND_ENDANGLEDIST == 1:
  class Endangledist(BondedInteraction):
    def typeNumber(self):
      return 11

    def typeName(self): 
      return "ENDANGLEDIST"

    def validKeys(self):
      return "bend", "phi0", "distmin", "distmax"

    def requiredKeys(self):
      return "bend", "phi0", "distmin", "distmax"

    def setDefaultParams(self):
      self._params = {"bend":0, "phi0":0, "distmin":0, "distmax":1}

    def _getParamsFromEsCore(self):
      return \
        {"bend":bonded_ia_params[self._bondId].p.endangledist.bend,\
         "phi0":bonded_ia_params[self._bondId].p.endangledist.phi0,\
         "distmin":bonded_ia_params[self._bondId].p.endangledist.distmin,\
         "distmax":bonded_ia_params[self._bondId].p.endangledist.distmax}

    def _setParamsInEsCore(self):
      endangledist_set_params(self._bondId, self._params["bend"], self._params["phi0"], self._params["distmin"],\
      self._params["distmax"])

ELSE:
  class Endangledist(BondedInteractionNotDefined):
    name="BOND_ENDANGLEDIST"
         
IF OVERLAPPED == 1:
  class Overlapped(BondedInteraction):
    def typeNumber(self):
      return 12

    def typeName(self): 
      return "OVERLAPPED"

    def validKeys(self):
      return "overlap_type", "filename"

    def requiredKeys(self):
      return "overlap_type", "filename"

    def setDefaultParams(self):
      self._params = {"overlap_type":0, "filename":""}

    def _getParamsFromEsCore(self):
      return \
        {"bend":bonded_ia_params[self._bondId].p.overlap.type,\
         "phi0":bonded_ia_params[self._bondId].p.overlap.filename}

    def _setParamsInEsCore(self):
      overlapped_bonded_set_params(self._bondId, self._params["overlap_type"], self._params["filename"])
      
ELSE:
  class Overlapped(BondedInteractionNotDefined):
    name="OVERLAPPED"
   
IF BOND_ANGLE == 1:   
  class Angle_Harmonic(BondedInteraction):
    def typeNumber(self):
      return 13

  def typeName(self): 
    return "ANGLE_HARMONIC"

  def validKeys(self):
    return "bend", "phi0"

  def requiredKeys(self): 
    return "bend", "phi0"

  def setDefaultParams(self):
    self._params = {"bend":0, "phi0":0} 

  def _getParamsFromEsCore(self):
    return \
      {"bend":bonded_ia_params[self._bondId].p.angle_harmonic.bend,\
       "phi0":bonded_ia_params[self._bondId].p.angle_harmonic.phi0}

  def _setParamsInEsCore(self):
    angle_harmonic_set_params(self._bondId,self._params["bend"],self._params["phi0"])
ELSE:
  class Angle_Harmonic(BondedInteractionNotDefined):
    name="BOND_ANGLE"
 
IF BOND_ANGLE == 1:   
  class Angle_Cosine(BondedInteraction):
    def typeNumber(self):
      return 14

  def typeName(self): 
    return "ANGLE_COSINE"

  def validKeys(self):
    return "bend", "phi0"

  def requiredKeys(self): 
    return "bend", "phi0"

  def setDefaultParams(self):
    self._params = {"bend":0, "phi0":0} 

  def _getParamsFromEsCore(self):
    return \
      {"bend":bonded_ia_params[self._bondId].p.angle_cosine.bend,\
       "phi0":bonded_ia_params[self._bondId].p.angle_cosine.phi0}

  def _setParamsInEsCore(self):
    angle_cosine_set_params(self._bondId,self._params["bend"],self._params["phi0"])
ELSE:
  class Angle_Cosine(BondedInteractionNotDefined):
    name="BOND_ANGLE"
   
IF BOND_ANGLE == 1:      
  class Angle_Cossquare(BondedInteraction):
    def typeNumber(self):
      return 15

  def typeName(self): 
    return "ANGLE_COSSQUARE"

  def validKeys(self):
    return "bend", "phi0"

  def requiredKeys(self): 
    return "bend", "phi0"

  def setDefaultParams(self):
    self._params = {"bend":0, "phi0":0} 

  def _getParamsFromEsCore(self):
    return \
      {"bend":bonded_ia_params[self._bondId].p.angle_cossquare.bend,\
       "phi0":bonded_ia_params[self._bondId].p.angle_cossquare.phi0}

  def _setParamsInEsCore(self):
    angle_cossquare_set_params(self._bondId,self._params["bend"],self._params["phi0"])
ELSE:
  class Angle_Cossquare(BondedInteractionNotDefined):
    name="BOND_ANGLE"
      
class Stretching_Force(BondedInteraction):
  def typeNumber(self):
    return 16

  def typeName(self): 
    return "STRETCHING_FORCE"

  def validKeys(self):
    return "r0", "ks"

  def requiredKeys(self): 
    return "r0", "ks"

  def setDefaultParams(self):
    self._params = {"r0":1., "ks":0} 

  def _getParamsFromEsCore(self):
    return \
      {"r0":bonded_ia_params[self._bondId].p.stretching_force.r0,\
       "ks":bonded_ia_params[self._bondId].p.stretching_force.ks}

  def _setParamsInEsCore(self):
    stretching_force_set_params(self._bondId,self._params["r0"],self._params["ks"])
   
   
class Area_Force_Local(BondedInteraction):
  def typeNumber(self):
    return 17

  def typeName(self): 
    return "AREA_FORCE_LOCAL"

  def validKeys(self):
    return "A0_l", "ka_l"

  def requiredKeys(self): 
    return "A0_l", "ka_l"

  def setDefaultParams(self):
    self._params = {"A0_l":1., "ka_l":0} 

  def _getParamsFromEsCore(self):
    return \
      {"A0_l":bonded_ia_params[self._bondId].p.area_force_local.A0_l,\
       "ka_l":bonded_ia_params[self._bondId].p.area_force_local.ka_l}

  def _setParamsInEsCore(self):
    area_force_local_set_params(self._bondId,self._params["A0_l"],self._params["ka_l"])
   
   
class Bending_Force(BondedInteraction):
  def typeNumber(self):
    return 18

  def typeName(self): 
    return "BENDING_FORCE"

  def validKeys(self):
    return "phi0", "kb"

  def requiredKeys(self): 
    return "phi0", "kb"

  def setDefaultParams(self):
    self._params = {"phi0":1., "kb":0} 

  def _getParamsFromEsCore(self):
    return \
      {"phi0":bonded_ia_params[self._bondId].p.bending_force.phi0,\
       "kb":bonded_ia_params[self._bondId].p.bending_force.kb}

  def _setParamsInEsCore(self):
    bending_force_set_params(self._bondId,self._params["phi0"],self._params["kb"])
   
   
class Volume_Force(BondedInteraction):
  def typeNumber(self):
    return 19

  def typeName(self): 
    return "VOLUME_FORCE"

  def validKeys(self):
    return "V0", "kv"

  def requiredKeys(self): 
    return "V0", "kv"

  def setDefaultParams(self):
    self._params = {"V0":1., "kv":0} 

  def _getParamsFromEsCore(self):
    return \
      {"V0":bonded_ia_params[self._bondId].p.volume_force.V0,\
       "kv":bonded_ia_params[self._bondId].p.volume_force.kv}

  def _setParamsInEsCore(self):
    volume_force_set_params(self._bondId,self._params["V0"],self._params["kv"])
   
   
class Area_Force_Global(BondedInteraction):
  def typeNumber(self):
    return 20

  def typeName(self): 
    return "AREA_FORCE_GLOBAL"

  def validKeys(self):
    return "A0_g", "ka_g"

  def requiredKeys(self): 
    return "A0_g", "ka_g"

  def setDefaultParams(self):
    self._params = {"A0_g":1., "ka_g":0} 

  def _getParamsFromEsCore(self):
    return \
      {"A0_g":bonded_ia_params[self._bondId].p.area_force_global.A0_g,\
       "ka_g":bonded_ia_params[self._bondId].p.area_force_global.ka_g}

  def _setParamsInEsCore(self):
    area_force_global_set_params(self._bondId,self._params["A0_g"],self._params["ka_g"])
   
   
class Stretchlin_Force(BondedInteraction):
  def typeNumber(self):
    return 21

  def typeName(self): 
    return "STRETCHLIN_FORCE"

  def validKeys(self):
    return "r0", "kslin"

  def requiredKeys(self): 
    return "r0", "kslin"

  def setDefaultParams(self):
    self._params = {"r0":1., "kslin":0} 

  def _getParamsFromEsCore(self):
    return \
      {"r0":bonded_ia_params[self._bondId].p.stretchlin_force.r0,\
       "kslin":bonded_ia_params[self._bondId].p.stretchlin_force.kslin}

  def _setParamsInEsCore(self):
    stretchlin_force_set_params(self._bondId,self._params["r0"],self._params["kslin"])

    


bondedInteractionClasses = {0:FeneBond, 1:HarmonicBond, 5:Dihedral, 6:Tabulated, 7:Subt_Lj,\
    9:Virtual, 11:Endangledist, 12:Overlapped,\
    13:Angle_Harmonic, 14:Angle_Cosine, 15:Angle_Cossquare, 16:Stretching_Force, 17:Area_Force_Local,\
    18:Bending_Force, 19:Volume_Force, 20:Area_Force_Global, 21:Stretchlin_Force}






class BondedInteractions:
  """Represents the non-bonded interactions. Individual interactions can be accessed using
  NonBondedInteractions[i], where i is the bond id. Will return an instance o
  BondedInteractionHandle"""
  def __getitem__(self, key):
    if not isinstance(key,int):
      raise ValueError("Index to BondedInteractions[] hast to ba an integer referring to a bond id")

    # Find out the type of the interaction from Espresso
    bondType = bonded_ia_params[key].type

    # Check if the bonded interaction exists in Espresso core
    if bondType == -1:
      raise ValueError("The bonded interaction with the id "+str(key)+" is not yet defined.")

    # Find the appropriate class representing such a bond
    bondClass = bondedInteractionClasses[bondType]
    print bondType
    print "  "

    # And return an instance of it, which refers to the bonded interaction id in Espresso
    return bondClass(key)
     
  def __setitem__(self,key,value):
    # Validate arguments
   
    # type of key must be int
    if not isinstance(key,int):
      raise ValueError("Index to BondedInteractions[] has to ba an integer referring to a bond id")

    # Value must be subclass off BondedInteraction
    if not isinstance(value,BondedInteraction):
      raise ValueError("Only subclasses of BondedInteraction can be assigned.")

    # Save the bond id in the BondedInteraction instance
    value._bondId=key

    # Set the parameters of the BondedInteraction instance in the Es core
    value._setParamsInEsCore()






