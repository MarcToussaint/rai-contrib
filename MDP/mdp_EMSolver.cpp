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
#include "mdp_EMSolver.h"
#include <iomanip>

void mdp::EMSolver::getParameters(){
  rai::getParameter(problemFile, "problemFile", rai::String("paint.95.POMDP.arr"));
  rai::getParameter(fscFile, "fscFile", rai::String("z.fsc"));
  rai::getParameter(outputPrefix, "outputPrefix", rai::String("data/test"));
  rai::getParameter(seed, "seed", (uint)0);
  rai::getParameter(EMiterations, "EMiterations", (uint)200);
  rai::getParameter(evaluationCheckHorizon, "evaluationCheckHorizon", (uint)0);
  
  rai::getParameter((int&)fscType, "fscType", 0);
  rai::getParameter(levels, "levels", TUP(10));
  
  rai::getParameter((int&)mstepType, "mstepType", 0);
  rai::getParameter(mstepRate, "mstepRate", .3);
  rai::getParameter(mstepNoise, "mstepNoise", 1e-5);
  
  rai::getParameter(estepHorizon, "estepHorizon", (uint)200);
  rai::getParameter(estepStructured, "estepStructured", true);
  rai::getParameter(estepIncremental, "estepIncremental", true);
  
  rai::getParameter(obsolete_forceLevel1, "forceLevel1", false);;
}

void mdp::EMSolver::reportParameters(std::ostream& os){
  cout <<"\n *** EMSolver parameter setting"
       <<"\ngeneral:"
       <<"\n  problemFile  = " <<problemFile
       <<"\n  fscFile      = " <<fscFile
       <<"\n  outputPrefix = " <<outputPrefix
       <<"\n  seed         = " <<seed
       <<"\n  EMiterations = " <<EMiterations
       <<"\n  evaluationCheckHorizon = " <<evaluationCheckHorizon
       <<"\ncontroller parameters:"
       <<"\n  fscType = " <<fscType
       <<"\n  levels  = " <<levels
       <<"\nM-step parameters:"
       <<"\n  mstepType  = " <<mstepType
       <<"\n  mstepRate  = " <<mstepRate
       <<"\n  mstepNoise = " <<mstepNoise
       <<"\nE-step parameters:"
       <<"\n  estepIncremental = " <<estepIncremental
       <<"\n  estepStructured  = " <<estepStructured
       <<"\n  estepHorizon     = " <<estepHorizon
       <<"\nobseolete:"
       <<"\n  forceLevel1 = " <<obsolete_forceLevel1
       <<endl;
}

void mdp::EMSolver::initProblem(){
  outfilename.clear() <<outputPrefix;
  //rai::IOraw=true;
  //levels.write(outfilename, "-");
  //rai::IOraw=false;
  outfilename <<"." <<seed;
  cout <<"output filename = " <<outfilename <<endl;
  outfile.close();
  rai::open(outfile, outfilename);
  
  readMDP(mdps, problemFile);
}

void mdp::EMSolver::initFsc(){
  rnd.seed(seed);
  switch(fscType){
    case FscPlain:
      //standardInitFsc_structured_levels(fsc, mdp, levels);
      standardInitFsc_structured_levels(fsc, mdps, levels);
      break;
    case FscHierarchical:
      //standardInitFsc_structured_hierarchical(fsc, mdp, levels);
      NIY;
      break;
    case FscReactive:
      //standardInitFsc_structured_reactive(fsc, mdp, levels);
      NIY;
      break;
  }
  //ofstream z(outfilename+".init-fsc");
  //writeFSC_fg(fsc, z, false);
  alpha.clear();
  beta.clear();
  outfile <<endl;
  k=0;
}

void mdp::EMSolver::resetTimer(){
  tic=rai::cpuTime();
}

void mdp::EMSolver::step(){
  cout <<k <<' ';
  cout <<std::setprecision(8);
  double R;
  R=pomdpEM_structured(mdps, fsc,
                       estepHorizon, estepStructured, alpha.N>0 && estepIncremental,
                       mstepType, mstepRate, mstepNoise,
                       false, &alpha, &beta, &cout);
  outfile <<k <<' ' <<rai::timerRead(false, tic) <<' ' <<R <<endl;
  if(false && evaluationCheckHorizon){ //NIY...
    uint horizon=evaluationCheckHorizon;
    if(horizon==1) horizon=2*estepHorizon;
    FSC_lev1 fsc1;
    MDP mdp;
    collapseFSC(fsc1, fsc);
    collapseToFlat(mdp, mdps);
    cout <<"evaluation = " <<evaluateFsc1(fsc1, mdp, horizon) <<endl;
  }
  k++;
}

void mdp::EMSolver::loop(){
  tic=rai::cpuTime();
  for(; k<EMiterations;) step();
  //cout <<"final R = " <<R <<"\ntotal time = " <<rai::timerRead(false, tic) <<endl;
  cout <<"total time = " <<rai::timerRead(false, tic) <<endl;
}

void mdp::EMSolver::gnuplot(bool byTime){
  rai::String cmd;
  if(byTime) cmd <<"plot '" <<outfilename <<"' us 2:3";
  else       cmd <<"plot '" <<outfilename <<"' us 1:3";
  cout <<"gnuplot command: " <<cmd <<endl;
  ::gnuplot(cmd.p);
}

void mdp::EMSolver::obsolete_loop_lev12(){
  MDP mdp;
  mdp::collapseToFlat(mdp, mdps);
  
  //init policy
  FSC_lev1 fsc1;
  FSC_lev2 fsc2;
  rnd.seed(seed);
  CHECK(levels.N==1 || levels.N==2, "old solver can only handle 1 or 2 levels!");
  if(levels.N==1){
    standardInitFsc1(fsc1, mdp, levels(0));
    writeFSC_lev1(STRING(outfilename <<".init-fsc1"), fsc1, false);
  }
  if(levels.N==2){
    standardInitFsc2(fsc2, mdp, levels(0), levels(1), fscType==FscHierarchical);
    writeFSC_lev2(STRING(outfilename <<".init-fsc2"), fsc2, false);
  }
  
  //iterate EM
  double R=0.;
  tic=rai::cpuTime();
  if(levels.N==1 || obsolete_forceLevel1){ //use specialized flat controller optimization
    if(levels.N==2) collapse2levelFSC(fsc1, fsc2);
    for(uint k=0; k<EMiterations; k++){
      if(evaluationCheckHorizon){
        uint horizon=evaluationCheckHorizon;
        if(horizon==1) horizon=2*estepHorizon;
        cout <<"evaluation = " <<evaluateFsc1(fsc1, mdp, horizon) <<endl;
      }
      cout <<k <<' ';
      //cout <<std::setprecision(8);
      R=pomdpEM_lev1(mdp, fsc1, estepHorizon, estepStructured, alpha.N>0 && estepIncremental,
                     (MstepType)mstepType, mstepRate, mstepNoise,
                     false, &alpha, &beta, &cout);
      outfile <<k <<' ' <<rai::timerRead(false, tic) <<' ' <<R <<endl;
    }
  }else{
    for(uint k=0; k<EMiterations; k++){
      if(evaluationCheckHorizon){
        uint horizon=evaluationCheckHorizon;
        if(horizon==1) horizon=2*estepHorizon;
        collapse2levelFSC(fsc1, fsc2);
        cout <<"evaluation = " <<evaluateFsc1(fsc1, mdp, horizon) <<endl;
      }
      cout <<k <<' ';
      R=pomdpEM_lev2(mdp, fsc2, estepHorizon, estepStructured, (MstepType)mstepType, false, &cout);
      outfile <<k <<' ' <<rai::timerRead(false, tic) <<' ' <<R <<endl;
    }
  }
  cout <<"final R = " <<R <<"\ntotal time = " <<rai::timerRead(false, tic) <<endl;
}
