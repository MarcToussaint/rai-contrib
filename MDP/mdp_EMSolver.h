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
#ifndef RAI_solver_h
#define RAI_solver_h

#include <Core/util.h>
#include "mdp.h"

namespace mdp {

enum FscType { FscPlain, FscHierarchical, FscReactive };

struct EMSolver {
  //----- public options of the solver
  //general
  rai::String problemFile;    ///< file name of the problem
  rai::String fscFile;        ///< file name for the controller
  rai::String outputPrefix;   ///< path-prefix for the outputfile
  uint seed;                 ///< random seed
  uint EMiterations;         ///< #iterations when loop() is called
  uint evaluationCheckHorizon;  ///< use a (flat) policy evaluation in each iteration to check the reward
  //controller options
  FscType fscType;           ///< fsc type [[todo make an enum]]
  uintA levels;              ///< array with #nodes in each level
  //Mstep options
  MstepType  mstepType;      ///< type of Mstep
  double mstepRate;          ///< Mstep ``convergence rate''
  double mstepNoise;         ///< Mstep noise level
  //Estep options
  bool estepStructured;      ///< use structured inference instead of collapsing everything to a single variable
  bool estepIncremental;     ///< reuse the alpha and beta from the previous iteration (-> can use much smaller horizon)
  uint estepHorizon;         ///< number of time slices propagated during the E-step
  //obsolete (still here for testing)
  bool obsolete_forceLevel1;
  
  //----- internal variables
  rai::String outfilename;
  ofstream outfile;
  arr alpha, beta;
  double tic;
  MDP_structured mdps;
  FSC_structured fsc;
  uint k;
  
  //----- methods
  void getParameters();
  void reportParameters(std::ostream& os);
  
  void clear(){ alpha.clear(); beta.clear(); clearMDP(mdps); clearFSC(fsc); }
  void initProblem();
  void initFsc();
  void resetTimer();
  void step();
  void loop();
  void gnuplot(bool byTime);
  
  void obsolete_loop_lev12(); ///< old routines for 1 or 2 levels (very useful for testing!)
  
};

}//end of namespace mdp

#ifdef  RAI_IMPLEMENTATION
#include "mdp_EMSolver.cpp"
#endif

#endif
