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
#include "mdp.h"
#include "mstep.h"
#include <Infer/infer.h>



double mdp::pomdpEM_structured(
  const MDP& mdp,
  FSC_structured& fsc,
  uint estepHorizon,
  bool estepStructured,
  bool estepIncremental,
  MstepType mstepType,
  double mstepRate,
  double mstepNoise,
  bool adaptP0,
  arr* alpha, arr* beta,
  ostream *os){
  
  rai::timerStart();
  
  uint i;
  //----- rescale rewards if necessary
  arr mdp_Rax = mdp.Rax;
#if rescaleRewards
  double Rmin=mdp_Rax.min(), Rmax=mdp_Rax.max();
  if(rescaleRewards || (mstepType!=MstepNoisyMax && Rmin<0.)){
    //if(!rescaleRewards) RAI_MSG("can't handle neg rewards in case of exact M-step -- I'm enforcing rescaling of rewards!");
    for(i=0; i<mdp_Rax.N; i++) mdp_Rax.elem(i) = (mdp_Rax.elem(i)-Rmin)/(Rmax-Rmin);
  }else{
    Rmin=0.; Rmax=1.;
  }
#endif
  
  
  //----- find the variable ids for the mdp world:
  infer::Variable *x=NULL, *y=NULL, *a=NULL, *x_=NULL, *y_=NULL;
  for(i=0; i<fsc.vars.N; i++) if(fsc.vars(i)->name=="state(t)"){         x =fsc.vars(i); break; }  if(i==fsc.vars.N) HALT("something's really wrong!");
  for(i=0; i<fsc.vars.N; i++) if(fsc.vars(i)->name=="observation(t)"){   y =fsc.vars(i); break; }  if(i==fsc.vars.N) HALT("something's really wrong!");
  for(i=0; i<fsc.vars.N; i++) if(fsc.vars(i)->name=="action(t)"){        a =fsc.vars(i); break; }  if(i==fsc.vars.N) HALT("something's really wrong!");
  for(i=0; i<fsc.vars.N; i++) if(fsc.vars(i)->name=="state(t+1)"){       x_=fsc.vars(i); break; }  if(i==fsc.vars.N) HALT("something's really wrong!");
  for(i=0; i<fsc.vars.N; i++) if(fsc.vars(i)->name=="observation(t+1)"){ y_=fsc.vars(i); break; }  if(i==fsc.vars.N) HALT("something's really wrong!");
  
  infer::VariableList mdp_leftVars={x};
  infer::VariableList mdp_rightVars={x_};
  
  //----- define factors for the mdp components
  //start
  infer::Factor Fx({x});       Fx.setP(mdp.Px);
  infer::Factor Fy({y});       Fy.setUniform();
  //transition
  infer::Factor Fxax({x_, a, x});  Fxax.setP(mdp.Pxax);
  infer::Factor Fyxa({y_, x_, a}); Fyxa.setP(mdp.Pyxa);
  //reward
  infer::Factor FRax({a, x});     FRax.setP(mdp_Rax);
  
  infer::FactorList mdp_transitions = {&Fyxa, &Fxax};
  infer::FactorList mdp_inits       = {&Fy, &Fx};
  infer::FactorList mdp_rewards     = {&FRax};
  
  infer::VariableList leftVars=(fsc.leftVars, mdp_leftVars);
  infer::VariableList rightVars=(fsc.rightVars, mdp_rightVars);
  infer::VariableList tail_headVars=(rightVars, leftVars);
  
  infer::FactorList allTransitions = (fsc.transFacs, mdp_transitions);
  infer::FactorList allRewards = (mdp_rewards, allTransitions);
  infer::FactorList allInits = (fsc.initFacs, mdp_inits);
  
  infer::Factor Falpha(leftVars);
  infer::Factor Fbeta(rightVars);
  arr PT;
  double PR, ET;
  if(!estepStructured){
    //----- collapse to unstructured model for generic inference
    uint dz = 1;  for(i=0; i<leftVars.N; i++) dz *= leftVars(i)->dim;
    //get transition matrix
    infer::Factor Fzz;
    eliminationAlgorithm(Fzz, allTransitions, tail_headVars);
    Fzz.P.reshape(dz, dz);
    //get reward vector
    infer::Factor FRz;
    eliminationAlgorithm(FRz, allRewards, leftVars);
    FRz.P.reshape(dz);
    //get start vector
    infer::Factor Fz;
    eliminationAlgorithm(Fz, allInits, leftVars);
    Fz.P.reshape(dz);
    
    //Fz >>FILE("z.gFz");
    //FRz >>FILE("z.gFRz");
    //Fzz >>FILE("z.gFzz");
    
    //----- E-STEP
    arr _alpha, _beta;
    if(estepIncremental){
      _alpha.referTo(*alpha);  _alpha.reshape(dz);
      _beta .referTo(*beta);   _beta .reshape(dz);
    }
    
    infer::inferMixLengthUnstructured(_alpha, _beta, PT, PR, ET,
				      Fz.P, FRz.P, Fzz.P, mdp.gamma, estepHorizon,
				      estepIncremental);
    Falpha.setP(_alpha);
    Fbeta .setP(_beta);
    
    if(!estepIncremental){
      if(alpha)(*alpha)=Falpha.P;
      if(beta)(*beta) =Fbeta.P;
    }
  }else{
    //----- use factor lists for generic inference
    infer::FactorList temporary;
    //eliminateVariable(allTransitions, temporary, a);
    //eliminateVariable(allTransitions, temporary, y_);
    //eliminateVariable(allRewards, temporary, a);
    
    if(estepIncremental){
      Falpha.setP(*alpha);
      Fbeta .setP(*beta);
    }
    
    inferMixLengthStructured(Falpha, Fbeta, PT, PR, ET,
                             leftVars, rightVars,
                             allInits,
                             allRewards,
                             allTransitions,
                             mdp.gamma, estepHorizon,
                             estepIncremental);
                             
    if(alpha)(*alpha)=Falpha.P;
    if(beta)(*beta) =Fbeta.P;
    
    for(uint i=0; i<temporary.N; i++) delete temporary(i);
  }
  
  if(os)(*os) <<"E: " <<rai::timerRead(true) <<"sec, M: " <<std::flush;
  
  //----- M-STEP
  //term2: derived from the full two-time-slice model (beta*P_(x'|x)*alpha)
  infer::FactorList twotimeslice; NIY// = ({&Fbeta}, fsc.transFacs, mdp_transitions, {&Falpha});
  
  //term1: derived from the immediate reward model
  infer::FactorList immediateR; NIY// = (mdp_rewards, fsc.transFacs, mdp_transitions, {&Falpha});
  
  //loop through all transition factors of the controller
  for(i=0; i<fsc.transFacs.N; i++){
    //term2: terms from the two-time-slice model
    infer::Factor X_term2;
    eliminationAlgorithm(X_term2, twotimeslice, fsc.transFacs(i)->variables);
    
    //term1: terms from immediate reward
    infer::Factor X_term1;
    eliminationAlgorithm(X_term1, immediateR  , fsc.transFacs(i)->variables);
    
    //get the expectations by adding both terms
    arr X;
    X = X_term2.P*::exp(X_term2.logP) + X_term1.P*::exp(X_term1.logP);
    
    //do the Mstep
    fsc.transFacs(i)->logP=0.;
    switch(mstepType){
      case MstepNoisyMax:
        noisyMaxMstep(fsc.transFacs(i)->P, X, 1, mstepRate, mstepNoise);
        tensorCheckCondNormalization(fsc.transFacs(i)->P, 1);
        break;
      case MstepExact:
        standardMstep(fsc.transFacs(i)->P, X, 1, mstepNoise);
        tensorCheckCondNormalization(fsc.transFacs(i)->P, 1);
        break;
      case MstepNone:
        break;
      case MstepCopyExpectations:
        fsc.transFacs(i)->P = X;
        break;
      default:
        NIY;
    }
  }
  
  if(adaptP0) NIY;
  
  //in case we rescaled, reset
  double expR;
#if rescaleRewards
  expR=(PR*(Rmax-Rmin)+Rmin)/(1.-mdp.gamma);
#else
  expR=PR/(1.-mdp.gamma);
#endif
  
  //----- rest is cosmetics
  //report
  if(os)(*os) <<rai::timerRead() <<"sec, " <<std::flush;
  if(os)
    (*os) <<" P(r=1)=" <<PR
    <<", Exp(T)=" <<ET <<"/" <<::log(PR)/::log(mdp.gamma)
    <<", Exp(R)=" <<expR
    <<endl;
    
  return expR;
}
