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



double mdp::pomdpEM_lev1(
  const MDP& mdp,
  FSC_lev1& fsc,
  uint estepHorizon,
  bool estepStructured,
  bool estepIncremental,
  MstepType mstepType,
  double mstepRate,
  double mstepNoise,
  bool adaptP0,
  arr* alpha, arr* beta,
  ostream *os){
  
//   cout <<"estepHorizon=" <<estepHorizon <<endl;
//   cout <<"estepStructured=" <<estepStructured <<endl;
//   cout <<"estepIncremental=" <<estepIncremental <<endl;
//   cout <<"mstepType=" <<mstepType <<endl;
//   cout <<"mstepRate=" <<mstepRate <<endl;
//   cout <<"mstepNoise=" <<mstepNoise <<endl;
//   cout <<"adaptP0=" <<adaptP0 <<endl;

  rai::timerStart();
  
  CHECK(mdp.Px.nd==1 && mdp.Pxax.nd==3 &&
        mdp.Pyxa.nd==3 && mdp.Rax.nd==2 &&
        fsc.Pa0.nd==2 && fsc.P0y0.nd==3, "");
        
  uint
  dx=mdp.Pxax.d0,
     da=mdp.Pxax.d1,
        dy=mdp.Pyxa.d0,
           d0=fsc.P0y0.d0;
           
  //----- rescale rewards if necessary
  arr mdp_Rax = mdp.Rax;
#if rescaleRewards
  double Rmin=mdp_Rax.min(), Rmax=mdp_Rax.max();
  if(rescaleRewards || (mstepType!=MstepNoisyMax && Rmin<0.)){
    //if(!rescaleRewards) RAI_MSG("can't handle neg rewards in case of exact M-step -- I'm enforcing rescaling of rewards!");
    for(uint i=0; i<mdp_Rax.N; i++) mdp_Rax.elem(i) = (mdp_Rax.elem(i)-Rmin)/(Rmax-Rmin);
  }else{
    Rmin=0.; Rmax=1.;
  }
#endif
  
  //----- define the factor model
  infer::Variable x(dx , "state(t)");
  infer::Variable y(dy , "observation(t)");
  infer::Variable n0(d0 , "node0(t)");
  infer::Variable a(da , "action(t)");
  infer::Variable x_(dx , "state(t+1)");
  infer::Variable y_(dy , "observation(t+1)");
  infer::Variable n0_(d0 , "node0(t+1)");
  //start
  infer::Factor Fx({&x}        , mdp.Px);
  infer::Factor Fy({&y});      Fy.setUniform();
  infer::Factor F0({&n0}       , fsc.P0);
  //transition
  infer::Factor Fa0({&a, &n0}     , fsc.Pa0);
  infer::Factor Fxax({&x_, &a, &x}   , mdp.Pxax);
  infer::Factor Fyxa({&y_, &x_, &a}  , mdp.Pyxa);
  infer::Factor F0y0({&n0_, &y_, &n0}, fsc.P0y0);
  //reward
  infer::Factor FRax({&a, &x}      , mdp_Rax);
  
  infer::VariableList leftVars={&n0 , &x};
  infer::VariableList rightVars={&n0_, &x_};
  infer::VariableList tail_headVars=cat(rightVars, leftVars);
  
  //infer::FactorList allTransitions={&Fa0, &Fxax, &Fyxa, &F0y0};
  infer::FactorList allTransitions = {&F0y0, &Fyxa, &Fxax, &Fa0};
  infer::FactorList allRewards = {&FRax, &Fa0};
  infer::FactorList allInits = {&F0, &Fx, &Fy};
  
  infer::Factor Falpha(leftVars);
  infer::Factor Fbeta(rightVars);
  arr PT;
  double PR, ET;
  if(!estepStructured){
    //----- collapse to unstructured model for generic inference
    uint dz = d0*dx;
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
  infer::FactorList twotimeslice = {&Falpha, &Fa0, &Fxax, &Fyxa, &F0y0, &Fbeta};
  
  infer::Factor X0y0_term2;
  eliminationAlgorithm(X0y0_term2, twotimeslice, {&n0_, &y_, &n0});
  infer::Factor Xa0_term2;
  eliminationAlgorithm(Xa0_term2, twotimeslice, {&a, &n0});
  
  //consider the 1st term (alpha*R_x)
  infer::FactorList immediateR   = {&Falpha, &Fa0, &Fxax, &Fyxa, &F0y0, &FRax};
  
  infer::Factor X0y0_term1;
  eliminationAlgorithm(X0y0_term1 , immediateR, {&n0_, &y_, &n0});
  infer::Factor Xa0_term1;
  eliminationAlgorithm(Xa0_term1, immediateR, {&a, &n0});
  
  arr X0y0;
#if 0 //this parameter does not depend on immediate reward (only term2 is relevant)
  X0y0 = X0y0_term2.P*::exp(X0y0_term2.logP)
         + X0y0_term1.P*::exp(X0y0_term1.logP);
#else
  X0y0 = X0y0_term2.P*::exp(X0y0_term2.logP);
#endif
         
  arr Xa0;
  Xa0 = Xa0_term2.P*::exp(Xa0_term2.logP)
        + Xa0_term1.P*::exp(Xa0_term1.logP);
        
  //do the M-step!
  switch(mstepType){
    case MstepNoisyMax:
      noisyMaxMstep(fsc.P0y0, X0y0, 1, mstepRate, mstepNoise);
      noisyMaxMstep(fsc.Pa0 , Xa0 , 1, mstepRate, mstepNoise);
      //if(adaptP0) noisyMaxMstep(fsc.P0   , X0.P   , 1);
      tensorCheckCondNormalization(fsc.Pa0 , 1);
      tensorCheckCondNormalization(fsc.P0y0, 1);
      break;
    case MstepExact:
      standardMstep(fsc.P0y0, X0y0, 1, mstepNoise);
      standardMstep(fsc.Pa0 , Xa0 , 1, mstepNoise);
      //if(adaptP0) standardMstep(fsc.P0   , X0.P   , 1);
      tensorCheckCondNormalization(fsc.Pa0, 1);
      tensorCheckCondNormalization(fsc.P0y0, 1);
      break;
    case Mstep11Rule:
      mstep_11rule(fsc.P0y0, X0y0, 1, mstepRate, mstepNoise);
      mstep_11rule(fsc.Pa0 , Xa0 , 1, mstepRate, mstepNoise);
      //if(adaptP0) mstep_11rule(fsc.P0   , X0   , 1, 1.1);
      tensorCheckCondNormalization(fsc.P0 , 1);
      tensorCheckCondNormalization(fsc.Pa0 , 1);
      tensorCheckCondNormalization(fsc.P0y0, 1);
      break;
    case MstepNone:
      break;
    case MstepCopyExpectations:
      fsc.P0y0= X0y0;
      fsc.Pa0 = Xa0 ;
      break;
    default:
      NIY;
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
