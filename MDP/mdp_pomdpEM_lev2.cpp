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



#ifndef rescaleRewards
#  define rescaleRewards false
#endif

#undef useStructure
//#define useStructure

double mdp::pomdpEM_lev2(
  const MDP& mdp,
  FSC_lev2& fsc,
  uint T,
  bool structuredEstep,
  bool maxMstep,
  bool adaptP0,
  std::ostream *os){
  
  CHECK(mdp.Px.nd==1 && mdp.Pxax.nd==3 &&
        mdp.Pyxa.nd==3 && mdp.Rax.nd==2 &&
        fsc.Pa0.nd==2 && fsc.P01y0.nd==4, "");
        
  uint
  dx=mdp.Pxax.d0,
     da=mdp.Pxax.d1,
        dy=mdp.Pyxa.d0,
           d0=fsc.P0.N,
              d1=fsc.P1.N;
              
  //----- rescale rewards if necessary
  arr mdp_Rax = mdp.Rax;
  double Rmin=mdp_Rax.min(), Rmax=mdp_Rax.max();
  if(rescaleRewards || (!maxMstep && Rmin<0.)){
    //if(!rescaleRewards) RAI_MSG("can't handle neg rewards in case of exact M-step -- I'm enforcing rescaling of rewards!");
    for(uint i=0; i<mdp_Rax.N; i++) mdp_Rax.elem(i) = (mdp_Rax.elem(i)-Rmin)/(Rmax-Rmin);
  }else{
    Rmin=0.; Rmax=1.;
  }
  
  rai::timerStart();
  
  //----- define the factor model
  infer::Variable x(dx , "state(t)");
  infer::Variable y(dy , "observation(t)");
  infer::Variable n1(d1 , "node1(t)");
  infer::Variable n0(d0 , "node0(t)");
  infer::Variable a(da , "action(t)");
  infer::Variable x_(dx , "state(t+1)");
  infer::Variable y_(dy , "observation(t+1)");
  infer::Variable n1_(d1 , "node1(t+1)");
  infer::Variable n0_(d0 , "node0(t+1)");
  //start
  infer::Factor Fx({&x}        , mdp.Px);
  infer::Factor F1({&n1}       , fsc.P1);
  infer::Factor F0({&n0}       , fsc.P0);
  //transition
  infer::Factor Fa0({&a, &n0}         , fsc.Pa0);
  infer::Factor Fxax({&x_, &a, &x}       , mdp.Pxax);
  infer::Factor Fyxa({&y_, &x_, &a}      , mdp.Pyxa);
  infer::Factor F1y01({&n1_, &y_ , &n0, &n1}, fsc.P1y01);
  infer::Factor F01y0({&n0_, &n1_, &y_, &n0}, fsc.P01y0);
  //reward
  infer::Factor FRax({&a, &x}      , mdp_Rax);
  
  infer::Factor Falpha({&n0 , &n1 , &x});
  infer::Factor Fbeta({&n0_, &n1_, &x_});
  arr PT;
  double PR, ET;
  if(!structuredEstep){
    //----- collapse to unstructured model for generic inference
    //get transition matrix
    infer::Factor F01x01x;
    eliminationAlgorithm(F01x01x, {&Fa0, &Fxax, &Fyxa, &F1y01, &F01y0}, {&n0_, &n1_, &x_, &n0, &n1, &x});
    F01x01x.P.reshape(d0*d1*dx, d0*d1*dx);
    //get reward vector
    infer::Factor FR01x;
    infer::Factor tmp({&n1}); tmp.setOne();
    eliminationAlgorithm(FR01x, {&tmp, &Fa0, &FRax}, {&n0, &n1, &x});
    FR01x.P.reshape(d0*d1*dx);
    //get start vector
    infer::Factor F01x;
    eliminationAlgorithm(F01x, {&F0, &F1, &Fx}, {&n0, &n1, &x});
    F01x.P.reshape(d0*d1*dx);
    
    //----- E-STEP
    arr alpha, beta;
    infer::inferMixLengthUnstructured(alpha, beta, PT, PR, ET,
				      F01x.P, FR01x.P, F01x01x.P, mdp.gamma, T);
    Falpha.setP(alpha);
    Fbeta .setP(beta);
  }else{
    //----- use factor lists for generic inference
    infer::FactorList trans = {&Fa0, &Fxax, &Fyxa, &F1y01, &F01y0};
    infer::FactorList newed;
    eliminateVariable(trans, newed, &a);
    //eliminateVariable(trans, newed, y_);
    
    infer::Factor tmp({&n1}); tmp.setOne();
    infer::FactorList rewards = {&FRax, &Fa0, &tmp};
    eliminateVariable(rewards, newed, &a);
    
    infer::inferMixLengthStructured(Falpha, Fbeta, PT, PR, ET,
				    {&n0 , &n1, &x}, {&n0_, &n1_, &x_},
				    {&F0, &F1, &Fx},
				    rewards,
				    trans,
				    mdp.gamma, T);
                             
    for(uint i=0; i<newed.N; i++) delete newed(i);
  }
  
  if(os)(*os) <<"E: " <<rai::timerRead(true) <<"sec, M: " <<std::flush;
  
  //----- M-STEP
  //consider the 2nd term (alpha*P_(x'|x)*beta)
  infer::FactorList twotimeslice = {&Falpha, &Fa0, &Fxax, &Fyxa, &F1y01, &F01y0, &Fbeta};
  
  infer::Factor X1y01_term2;
  eliminationAlgorithm(X1y01_term2, twotimeslice, {&n1_, &y_ , &n0, &n1});
  infer::Factor X01y0_term2;
  eliminationAlgorithm(X01y0_term2, twotimeslice, {&n0_, &n1_, &y_, &n0});
  infer::Factor Xa0_term2;
  eliminationAlgorithm(Xa0_term2, twotimeslice, {&a, &n0});
  
  //consider the 1st term (alpha*R_x)
  infer::FactorList immediateR   = {&Falpha, &Fa0, &Fxax, &Fyxa, &F1y01, &F01y0, &FRax};
  
  infer::Factor X1y01_term1;
  eliminationAlgorithm(X1y01_term1 , immediateR, {&n1_, &y_ , &n0, &n1});
  infer::Factor X01y0_term1;
  eliminationAlgorithm(X01y0_term1 , immediateR, {&n0_, &n1_, &y_, &n0});
  infer::Factor Xa0_term1;
  eliminationAlgorithm(Xa0_term1, immediateR, {&a, &n0});
  
  arr X1y01;
  X1y01 = X1y01_term2.P*::exp(X1y01_term2.logP);
  
  arr X01y0;
  X01y0 = X01y0_term2.P*::exp(X01y0_term2.logP);
  
  arr Xa0;
  Xa0 = Xa0_term2.P*::exp(Xa0_term2.logP)
        + Xa0_term1.P*::exp(Xa0_term1.logP);
        
  //do the M-step!
  if(maxMstep){
    noisyMaxMstep_old(fsc.P1y01, X1y01, 1);
    noisyMaxMstep_old(fsc.P01y0, X01y0, 1);
    noisyMaxMstep_old(fsc.Pa0  , Xa0  , 1);
    //if(adaptP0) noisyMaxMstep(fsc.P0   , X0.P   , 1);
  }else{
    standardMstep(fsc.P1y01, X1y01, 1);
    standardMstep(fsc.P01y0, X01y0, 1);
    standardMstep(fsc.Pa0  , Xa0  , 1);
    //if(adaptP0) standardMstep(fsc.P0   , X0.P   , 1);
  }
  
  tensorCheckCondNormalization(fsc.Pa0  , 1, 1e-10);
  tensorCheckCondNormalization(fsc.P1y01, 1, 1e-10);
  tensorCheckCondNormalization(fsc.P01y0, 1, 1e-10);
  
  //----- rest is cosmetics
  //report
  if(os)(*os) <<rai::timerRead() <<"sec, " <<std::flush;
  if(os)(*os) <<" P(r=1)=" <<PR
    <<", Exp(T)=" <<ET <<"/" <<::log(PR)/::log(mdp.gamma)
    <<", Exp(R)=" <<(PR*(Rmax-Rmin)+Rmin)/(1.-mdp.gamma)
    <<endl;
    
  return (PR*(Rmax-Rmin)+Rmin)/(1.-mdp.gamma);
}
