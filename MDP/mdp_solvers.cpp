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
#include <Core/util.h>
#include "mdp.h"
#include "mstep.h"
#include "math.h"

#define Pt(t) ::pow(mdp.gamma, t)
#ifndef rescaleRewards
#  define rescaleRewards 0
#endif

arr ALPHA, BETA, GAMMA, LIKE;

//===========================================================================
//
// policy manipulation
//

void mdp::randomPolicy(arr& pi, const MDP& mdp){
  uint x, a, X=mdp.Px.N, A=mdp.Pxax.d1;
  pi.resize(A, X);
  for(x=0; x<X; x++) for(a=0; a<A; a++) pi(a, x)=1./A;
}

void mdp::maxPolicyMap(uintA& piMap, const arr& pi){
  uint x, X=pi.d1;
  arr tpi;
  op_transpose(tpi, pi);
  piMap.resize(X);
  for(x=0; x<X; x++) piMap(x) = argmax(tpi[x]);
}

void mdp::getPxx(arr& Pxx, const arr& Pxax, const arr& pi){
  uint x, y, a, X=pi.d1, A=pi.d0;
  Pxx.resize(X, X);
  Pxx.setZero();
  for(y=0; y<X; y++) for(a=0; a<A; a++) for(x=0; x<X; x++)
        Pxx(y, x) += Pxax(y, a, x) * pi(a, x);
}

void mdp::getMaxPxxMap(uintA& PxxMap, const arr& Pxax, const arr& pi){
  arr tmp;
  uint x, X=pi.d1;
  PxxMap.resize(X);
  getPxx(tmp, Pxax, pi);
  transpose(tmp);
  for(x=0; x<X; x++) PxxMap(x) = argmax(tmp[x]);
}

//===========================================================================
//
// Dynamic Programming
//

double mdp::MDP::Qvalue(uint a, uint x, const arr& V) const{
  uint y, j;
  uintA& neigh = neighbors.elem(x);
  double Qax=0.;
  for(j=0; j<neigh.N; j++){
    y=neigh(j);
    Qax+=Pxax(y, a, x) * V(y);
  }
  Qax *= gamma;
  Qax += Rax(a, x);
  return Qax;
}

void mdp::MDP::Qfunction(arr& Q, const arr& V) const{
  uint x, a, X=Px.N, A=Pxax.d1;
  Q.resize(A, X).setZero();
  for(x=0; x<X; x++) for(a=0; a<A; a++) Q(a, x) = Qvalue(a, x, V);
}

void mdp::MDP::valueIteration(arr& V) const{
  uint x, a, X=Px.N, A=Pxax.d1;
  if(V.N!=X) V.resize(X).setZero();
  double Qax;
  arr Vnew(X);
  for(x=0; x<X; x++) for(a=0; a<A; a++){
      Qax = Qvalue(a, x, V);
      if(!a || Qax>Vnew(x)) Vnew(x)=Qax;
  }
  V=Vnew;
}

void mdp::MDP::policyEvaluation(arr& V, const arr& pi) const{
  uint x, a, X=Px.N, A=Pxax.d1;
  if(V.N!=X){ V.resize(X); V.setZero(); }
  arr Vnew;
  Vnew.resize(X);
  Vnew.setZero();
  for(x=0; x<X; x++) for(a=0; a<A; a++){
      Vnew(x) += pi(a, x) * Qvalue(a, x, V);
  }
  V=Vnew;
}

void mdp::MDP::maxPolicy(arr& pi, const arr& V) const{
  uint x, a, X=Px.N, A=Pxax.d1;
  uint amax=0;
  double Qax, Qmax=0.;
  pi.resize(A, X);
  pi.setZero();
  for(x=0; x<X; x++){
    for(a=0; a<A; a++){
      Qax=Qvalue(a, x, V);
      if(!a || Qax>Qmax){ amax=a; Qmax=Qax; }
    }
    pi(amax, x)=1.;
  }
}

#undef ERR
arr *ERR;
static bool PQcompare(const uint& a, const uint& b){
  if((*ERR)(a)<(*ERR)(b)) return true;
  if((*ERR)(a)>(*ERR)(b)) return false;
  return a<=b;
}

void mdp::MDP::prioritizedSweeping(arr& V, double VerrThreshold) const{
  uint x, a, X=Px.N, A=Pxax.d1;
  
  if(V.N!=X){ V.resize(X); V.setZero(); }
  
  arr Vnew(X);
  Vnew.setZero();
  
  uintA *neigh;
  uint j, y;
  double Qay;
  
  uintA queue;
  queue.memMove=true;
  arr   ERR(X);
  ::ERR = &ERR;
  boolA inQueue(X);
  inQueue=false;
  
  //add max reward to the que
  arr Rx;
  tensorMaxMarginal(Rx, Rax, uintA{1});
  x=argmax(Rx);
  Vnew(x)=Rx(x);
  ERR(x) =Vnew(x)-V(x);
  queue.insertInSorted(x, PQcompare);
  inQueue(x)=true;
  
  uint t;
  for(t=0; queue.N; t++){
    x=queue.popLast();
    inQueue(x)=false;
    V(x) += 1. * (Vnew(x) - V(x)); //update towards correct value;
    
    neigh=&neighbors(x);
    for(j=0; j<neigh->N; j++){
      y=(*neigh)(j);
      for(a=0; a<A; a++){
        Qay=Qvalue(a, y, V);
        if(!a || Qay>Vnew(y)) Vnew(y)=Qay;
      }
      if(inQueue(y)){ queue.removeValueInSorted(y, PQcompare); inQueue(y)=false; }
      ERR(y) = Vnew(y)-V(y);
      if(ERR(y) > VerrThreshold){ queue.insertInSorted(y, PQcompare); inQueue(y)=true; }
    }
    //if(!(t%X)) glDisplayRedBlue(V, global_maze.d0, global_maze.d1, false);
  }
  //report to console:
  std::cout <<"PS: stopped after " <<t <<" updates" <<std::endl;
  // <<" cost=" <<cost
  // <<" value of start state=" <<from->V <<std::endl;
}

void weightedAddLog(arr& x, double& xLog, double w, const arr& y, double yLog){
  w *= ::exp(yLog-xLog);
  for(uint i=0; i<x.N; i++) x.elem(i) = x.elem(i) + w*y.elem(i);
  //lognormScale(x, xLog);
}

void mdp::mdpEM(const MDP& mdp, arr& pi, arr& hatBeta, uint Tmax, float cutoffTimeFactor, MstepType mstepType, arr *Pvisited){
  uint x, y, a, X=mdp.Px.N, A=mdp.Pxax.d1;
  
  uint cutoffTime=0;
  
  arr mdp_Rax = mdp.Rax;
  double Rmin=min(mdp_Rax), Rmax=max(mdp_Rax);
  if(rescaleRewards || (mstepType!=MstepNoisyMax && Rmin<0.)){
    if(!rescaleRewards) RAI_MSG("can't handle neg rewards in case of exact M-step -- I'm enforcing rescaling of rewards!");
    for(uint i=0; i<mdp_Rax.N; i++) mdp_Rax.elem(i) = (mdp_Rax.elem(i)-Rmin)/(Rmax-Rmin);
  }else{
    Rmin=0.; Rmax=1.;
  }
  
  arr Rx(X);
  Rx.setZero();
  for(x=0; x<X; x++) for(a=0; a<A; a++) Rx(x) += mdp_Rax(a, x)*pi(a, x);
  
  arr Pxx(X, X);
  Pxx.setZero();
  for(x=0; x<X; x++) for(y=0; y<X; y++)
      for(a=0; a<A; a++) Pxx(y, x) += mdp.Pxax(y, a, x) * pi(a, x);
      
  arr Pxx_back(X, X);
  arr Px_norm(X);
  if(Pvisited && Pvisited->N==X){
    Pxx_back.setZero();
    Px_norm.setZero();
    for(x=0; x<X; x++) for(y=0; y<X; y++) Pxx_back(x, y) += Pxx(y, x) * Pvisited->elem(x);
    for(x=0; x<X; x++) for(y=0; y<X; y++) Px_norm(y) += Pxx(y, x) * Pvisited->elem(x);
    tensorCondNormalize(Pxx_back, 1);
  }else{
    op_transpose(Pxx_back, Pxx);
  }
  
  arr alpha(X), beta(X), beta_neutral(X);
  double alphaLog=0., betaLog=0.;
  alpha = mdp.Px;
  beta  = Rx;
  beta_neutral = 1.;
  
  arr alpha_(X), beta_(X), beta_neutral_(X);
  ALPHA.resize(Tmax, X); ALPHA.setZero(); ALPHA[0] = alpha;
  BETA .resize(Tmax, X); BETA .setZero(); BETA[0]  = beta;
  
  
  arr hat_alpha(X), hat_beta(X);
  double hat_alphaLog=0., hat_betaLog=0.;
  hat_alpha= alpha;
  hat_beta = beta;
  
  //if(timePriorType()=='w') cutoffTime=timeWindowH+1;
  //if(!prune) forNodes(n, G){ addA(n); addB(n); } addB(to);
  
  uintA *neigh;
  
  uint t, j;
  arr L(2*Tmax);  L.setZero();
  L(0) = sum(alpha % beta)*::exp(alphaLog+betaLog);
  
  ofstream os("z.EM");
  
  for(t=1; t<Tmax; t++){
    //prune in second phase: states for which gamma will be zero
    /*
    if(prune && cutoffTime && t>cutoffTime/2){
    forNodes_save(nn, nns, Aset) if(nn->x->beta [cutoffTime-t]==0) delA(nn);
    forNodes_save(nn, nns, Bset) if(nn->x->alpha[cutoffTime-t]==0) delB(nn);
    }
    */
    
    // BETA bwd propagation
    beta_.setZero();
    for(x=0; x<X; x++){
      neigh=&mdp.neighbors(x);
      for(j=0; j<neigh->N; j++){
        y=(*neigh)(j);
        beta_(y) += Pxx_back(y, x) * beta(x);
        //beta_(y) += Pxx(x, y) * beta(x);
        //for(a=0;a<A;a++)  beta_(y) += mdp.Pxax(x, a, y) * pi(a, y) * beta(x);
      }
    }
    beta = beta_;
    //for(x=0;x<X;x++) if(beta(x)) addB(x);
    //lognormScale(beta, betaLog);
    
    // BETA_neutral bwd propagation
    beta_neutral_.setZero();
    for(x=0; x<X; x++){
      neigh=&mdp.neighbors(x);
      for(j=0; j<neigh->N; j++){
        y=(*neigh)(j);
        beta_neutral_(y) += Pxx_back(y, x) * beta_neutral(x);
      }
    }
    beta_neutral = beta_neutral_;
    
    L(t+t-1)   = sum(alpha % beta/beta_neutral)*::exp(alphaLog+betaLog);
    
    // ALPHA fwd propagation
    alpha_.setZero();
    for(x=0; x<X; x++){
      neigh=&mdp.neighbors(x);
      for(j=0; j<neigh->N; j++){
        y=(*neigh)(j);
        alpha_(y) += Pxx(y, x) * alpha(x);
        //for(a=0;a<A;a++)  alpha_(y) += mdp.Pxax(y, a, x) * pi(a, x) * alpha(x);
      }
    }
    alpha = alpha_;
    //for(x=0;x<X;x++) if(alpha(x)) addA(x);
    //lognormScale(alpha, alphaLog);
    
    L(t+t)   = sum(alpha % beta/beta_neutral)*::exp(alphaLog+betaLog);
    
    weightedAddLog(hat_alpha, hat_alphaLog, Pt(t), alpha, alphaLog);
    weightedAddLog(hat_beta , hat_betaLog , Pt(t), beta /beta_neutral, betaLog);
    
    ALPHA[t] = alpha;
    BETA[t]  = beta;
    //prune in all phases: too small alphas or betas
    if(false){
      /*maxb=0.;
      forNodes(nn, Bset) if(!maxb || nn->x->beta[t]>maxb) maxb=nn->x->beta[t];
      forNodes_save(nn, nns, Bset) if(nn->x->beta [t]<1e-3*maxb) delB(nn);
        */
    }
    
    //choose cutoff time
    if(!cutoffTime && t>10 && L(t+t)>0.){
      cutoffTime=(uint)ceil(2.*cutoffTimeFactor*t);
      if(cutoffTime>Tmax) cutoffTime=Tmax;
      std::cout <<"fronts meet at time " <<t <<" - cutoff time set to " <<cutoffTime <<std::endl;
    }
    
    //stopping criterion
    if(cutoffTime && t>cutoffTime) break;
    //CHECK(Aset && Bset.N, "no active sites");
    
    //glDisplayRedBlue(beta-alpha, global_maze.d0, global_maze.d1, false);
    //glDisplayRedBlue(alpha, global_maze.d0, global_maze.d1, false);
    showAB(alpha, beta);
    os <<t <<' ' <<L(t+t-1) <<' ' <<L(t+t) <<std::endl;
    //Z-file: alpha-norm beta-norm time Z(t) beta_START(t) alpha_GOAL(t) Z(t)-alternative
  }
  cout <<"E: " <<rai::timerRead() <<"sec, M: " <<std::flush;
  
  /*if(Pvisited && Pvisited->N==X){
    for(x=0;x<X;x++) hat_beta(x) = rai::DIV(hat_beta(x), Px_norm(x));
  }*/
  
  hatBeta=hat_beta * ::exp(hat_betaLog);
  
  
  //***** M-STEP, part 1: compute the expectations
  
  arr Xax(A, X);
  Xax.setZero();
  
  for(x=0; x<X; x++) for(a=0; a<A; a++){
      neigh=&mdp.neighbors(x);
      for(j=0; j<neigh->N; j++){
        y=(*neigh)(j);
        Xax(a, x) += mdp.Pxax(y, a, x) * hat_beta(y);
      }
    }
  Xax *= ::exp(hat_betaLog);
  Xax *= mdp.gamma;
  Xax += mdp_Rax;
  
  if(mstepType==MstepNoisyMax){
    noisyMaxMstep2(pi, Xax, 1, .0);
  }
  if(mstepType==MstepExact){
    standardMstep(pi, Xax, 1, .0);
  }
  
  LIKE = L;
  
  //report
  cout <<rai::timerRead() <<"sec, " <<std::flush;
  double Ltot=0., PtNorm=0., Texp=0.;
  arr PtL(2*Tmax);
  for(t=0; t<2*Tmax; t++){ PtL(t)=Pt(t)*L(t); Ltot+=PtL(t); Texp+=t*PtL(t); PtNorm+=Pt(t); }
  PtL  /= PtNorm;
  Texp /= Ltot;
  Ltot /= PtNorm; //(in case Pt was not normalized)
  cout <<" P(r=1)=" <<Ltot
       <<", Exp(T)=" <<Texp <<"/" <<::log(Ltot*PtNorm)/::log(mdp.gamma)
       <<", Exp(R)=" <<(Ltot*(Rmax-Rmin)+Rmin)*PtNorm
       <<endl;
       
  //glDisplayRedBlue(hat_beta, global_maze.d0, global_maze.d1, true);
  //glDisplayRedBlue(hat_alpha, global_maze.d0, global_maze.d1, true);
}

/// calculates Gvis correctly with the posterior over states
void mdp::calcPvisited(arr& Pvisited, const MDP& mdp){
  uint t, T, x, Tmax=ALPHA.d0, X=mdp.Px.N;
  double Z, p;
  
  GAMMA.resizeAs(ALPHA);
  GAMMA.setZero();
  Pvisited.resize(X);
  Pvisited.setZero();
  for(x=0; x<X; x++){
    for(T=0; T<Tmax; T++) if(LIKE(T)){
        p=0.;
        for(t=0; t<=T; t++) p+=(1.-p) * 1./LIKE(T) * ALPHA(t, x)*BETA(T-t, x);
        GAMMA(T, x)=p;
        Pvisited(x)+=Pt(T)*LIKE(T)*GAMMA(T, x);
      }
  }
  Z=0.;
  for(T=0; T<Tmax; T++) Z+=Pt(T)*LIKE(T);
  Pvisited/=Z;
}
