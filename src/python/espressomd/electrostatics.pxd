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
# Handling of interactions

include "myconfig.pxi"
from _system cimport *
cimport numpy as np
from utils cimport *


IF ELECTROSTATICS == 1:
  cdef extern from "interaction_data.hpp":
    int coulomb_set_bjerrum(double bjerrum)

  IF P3M == 1:
    cdef extern from "p3m-common.hpp":
      ctypedef struct p3m_parameter_struct:
        double alpha_L
        double r_cut_iL
        int    mesh[3]
        double mesh_off[3]
        int    cao
        int    inter
        double accuracy
        double epsilon
        double cao_cut[3]
        double a[3]
        double ai[3]
        double alpha
        double r_cut
        int    inter2
        int    cao3
        double additional_mesh[3]

      ctypedef struct p3m_send_mesh:
        int s_dim[6][3]
        int s_ld[6][3]
        int s_ur[6][3]
        int s_size[6]
        int r_dim[6][3]
        int r_ld[6][3]
        int r_ur[6][3]
        int r_size[6]
        int max

      ctypedef struct p3m_local_mesh:
        int dim[3]
        int size
        int ld_ind[3]
        double ld_pos[3]
        int inner[3]
        int in_ld[3]
        int in_ur[3]
        int margin[6]
        int r_margin[6]
        int q_2_off
        int q_21_off


    cdef extern from "p3m.hpp":
      int p3m_set_params(double r_cut, int *mesh, int cao, double alpha, double accuracy)
      void p3m_set_tune_params(double r_cut, int mesh[3], int cao, double alpha, double accuracy, int n_interpol)
      int p3m_set_mesh_offset(double x, double y, double z)
      int p3m_set_eps(double eps)
      int p3m_set_ninterpol(int n)
      int p3m_adaptive_tune(char** log)

 
      ctypedef struct p3m_data_struct:
        p3m_parameter_struct params
      

      p3m_data_struct p3m_get_params()

    # Convert C arguments into numpy array
    cdef inline python_p3m_set_mesh_offset(mesh_off):
      cdef double mesh_offset[3]
      mesh_offset[0] = mesh_off[0]
      mesh_offset[1] = mesh_off[1]
      mesh_offset[2] = mesh_off[2]
      return p3m_set_mesh_offset(mesh_offset[0],mesh_offset[1],mesh_offset[2])

      
    cdef inline python_p3m_adaptive_tune():
      cdef char* log = NULL
      cdef int response
      response = p3m_adaptive_tune(&log)
      return response, log


    cdef inline python_p3m_set_params(p_r_cut, p_mesh, p_cao, p_alpha, p_accuracy, p_n_interpol):
      cdef int mesh[3]
      cdef double r_cut
      cdef int cao
      cdef double alpha
      cdef double accuracy
      cdef int n_interpol
      r_cut = p_r_cut
      cao = p_cao
      alpha = p_alpha
      accuracy = p_accuracy
      n_interpol = p_n_interpol
      if isinstance(p_mesh,int):
        mesh[0]=p_mesh
        mesh[1]=p_mesh
        mesh[2]=p_mesh
      else:
        mesh=p_mesh

      p3m_set_params(r_cut, mesh, cao, alpha, accuracy)

    cdef inline python_p3m_set_tune_params(p_r_cut, p_mesh, p_cao, p_alpha, p_accuracy, p_n_interpol):
      # cdef inline python_p3m_set_tune_params():
      cdef int mesh[3]
      cdef double r_cut
      cdef int cao
      cdef double alpha
      cdef double accuracy
      cdef int n_interpol
      r_cut = p_r_cut
      cao = p_cao
      alpha = p_alpha
      accuracy = p_accuracy
      n_interpol = p_n_interpol
      if isinstance(p_mesh,int):
        mesh[0]=p_mesh
        mesh[1]=p_mesh
        mesh[2]=p_mesh
      else:
        mesh=p_mesh

      p3m_set_tune_params(r_cut, mesh, cao, alpha, accuracy, n_interpol)
