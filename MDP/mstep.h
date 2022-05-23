/*  ---------------------------------------------------------------------
    Copyright 2012 Marc Toussaint
    email: mtoussai@cs.tu-berlin.de
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a COPYING file of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>
    -----------------------------------------------------------------  */
#ifndef RAI_mstep_h
#define RAI_mstep_h

#include <Core/array.h>

//enum MstepType{ MstepNoisyMax, MstepExact, MstepNone, MstepCopyExpectations };

void standardMstep(arr& param, const arr& post, uint left, double noise=0);

void noisyMaxMstep_old(arr& param, const arr& post, uint left);
void noisyMaxMstep(arr& param, const arr& post, uint left, double rate, double noise=0);
void noisyMaxMstep2(arr& param, const arr& likelihood, uint left, double noise=0);
void mstep_11rule(arr& param, const arr& post, uint left, double rate, double noise=0);

#ifdef  RAI_IMPLEMENTATION
#  include "mstep.cpp"
#endif

#endif
