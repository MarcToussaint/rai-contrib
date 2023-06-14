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
#include "BinaryBP.h"
#include <Core/util.h>


//===========================================================================
//
// basic equations
//

/*double ATANH(double x){
  //if(fabs(x)>.999) cout <<"high atanh(x), x=" <<x ;
  return atanh(x);
}
double TANH(double x){
  //if(fabs(x)>10.) cout <<"high tanh(x), x=" <<x ;
  return tanh(x);
}*/
#define TANH tanh
#define ATANH atanh
#undef EXP
double RATIO(double a, double b){
  double c=(1.-rai::sqr(TANH(a)));
  double d=(1.-rai::sqr(TANH(b)));
  //if(fabs(d)<1e-5) cout <<"low quotient! 1/x, x=" <<d;
  return c/d;
}
double myLog(double x){ //return log(x);
  if(x<1e-10) return -23.02585093;
  if(x>1e+10) return +23.02585093;
  return log(x);
}
double EXP(double x){ //return exp(x);
  if(x<-23.02585093) return 1e-10;
  if(x>+23.02585093) return 1e+10;
  return exp(x);
}

double p_to_ratio(double p){
  if(p<.5){ if(p<1e-10) p=1e-10; return .5*log(p/(1.-p)); } else    { p=1.-p; if(p<1e-10) p=1e-10; return -.5*log(p/(1.-p)); }
}
double ratio_to_p(double a){ double p1=exp(a); return p1/(p1+1./p1); }
arr    ratios_to_ps(const arr& a){ arr p(a.N);  for(uint i=0; i<a.N; i++) p(i)=ratio_to_p(a(i));  return p; }
double table_to_nodeExp(const arr& p){ return p_to_ratio(p(1)/(p(0)+p(1))); }
void   nodeExp_to_Table(arr& p, double a){ p.resize(2); p(1)=exp(a); p(0)=1./p(1); } //p/=p(0)+p(1); }
void   pairExp_to_Table(arr& p, double J, double a, double b){
  p.resize(2, 2);
  p(0, 0) = exp(J - a - b);
  p(0, 1) = exp(-J - a + b);
  p(1, 0) = exp(-J + a - b);
  p(1, 1) = exp(J + a + b);
  //p /= sum(p);
}
double table_to_pairExp(const arr& p){
  return .25 * log((p(0, 0)*p(1, 1)) / (p(0, 1)*p(1, 0)));
}
void   table_to_pairExp(double& J, double& a, double& b, const arr& p){
  J=.25 * log((p(0, 0)*p(1, 1)) / (p(0, 1)*p(1, 0)));
  a=.25 * log((p(1, 0)*p(1, 1)) / (p(0, 0)*p(0, 1)));
  b=.25 * log((p(0, 1)*p(1, 1)) / (p(0, 0)*p(1, 0)));
}
double KLD_node_probVsProb(double p, double q){
  return log(cosh(q)/cosh(p)) + (p-q)*tanh(p);
}
double KLD_node_probVsFactor(double b, double f){
  return -log(2.*cosh(b)) + (b-f)*tanh(b);
}
double KLD_pair_probVsFactor(double J, double a, double b, double thi, double thj){
  double tJ=tanh(J), ta=tanh(a), tb=tanh(b);
  double cJ=cosh(J), ca=cosh(a), cb=cosh(b);
  double down=1+tJ*ta*tb;
#if 1
  return -log(4.*(cJ*ca*cb*down)) + (a-thi)*(ta+tJ*tb)/down + (b-thj)*(tb+tJ*ta)/down;
#else //different ways to compute the same thing...
  arr B, F;
  pairExp_to_Table(B, J, a, b);
  double sumB=sum(B);
  B/=sumB;
  pairExp_to_Table(F, J, thi, thj);
  double D=0.;
  for(uint x=0; x<2; x++) for(uint y=0; y<2; y++)  D += B(x, y) * log(B(x, y)/F(x, y));
  cout <<"D1=" <<D;
  D=0.;
  for(uint x=0; x<2; x++) for(uint y=0; y<2; y++)  D += B(x, y) * ((2.*x-1.)*(a-thi)+ (2.*y-1.)*(b-thj));
  D -= log(sumB);
  cout <<" D3=" <<D;
  D = (B(1, 0)+B(1, 1)-B(0, 0)-B(0, 1))*(a-thi)
      + (B(0, 1)+B(1, 1)-B(0, 0)-B(1, 0))*(b-thj)
      - log(sumB);
  cout <<" D2=" <<D;
  D = -log(4.*(cJ*ca*cb*down)) + (a-thi)*(ta+tJ*tb)/down + (b-thj)*(tb+tJ*ta)/down;
  cout <<" D4=" <<D <<endl;
  return D;
#endif
}
double KLD(const BinaryPairFG &B, const BinaryPairFG &F){
  CHECK_EQ(B.f_i.N ,F.f_i.N, "");
  CHECK_EQ(B.f_ij.N,F.f_ij.N, "");
  uint i, N=B.f_i.d0, E=B.edges.d0, x, y;
  uintA degree(N); degree.setZero();
  arr b, f;
  double D = 0.;
  for(i=0; i<E; i++){  degree(B.edges(i, 0))++;  degree(B.edges(i, 1))++;  }
  for(i=0; i<E; i++){
    b.referToDim(B.f_ij, i);
    f.referToDim(F.f_ij, i);
    for(x=0; x<2; x++) for(y=0; y<2; y++)  D += b(x, y) * log(b(x, y)/f(x, y));
  }
  for(i=0; i<N; i++){
    b.referToDim(B.f_i, i);
    f.referToDim(F.f_i, i);
    for(x=0; x<2; x++) D += (1.-degree(i)) * b(x) * log(b(x)/f(x));
  }
  return D;
}

//===========================================================================
//
// BinaryPairFG
//

void BinaryPairFG::write(std::ostream& os){
  os <<"logZ=" <<logZ <<"\nedges=\n" <<edges <<"\nf_ij=\n" <<f_ij <<"\nf_i=\n" <<f_i <<endl;
}

//===========================================================================
//
// BinaryBPNet
//

void BinaryBPNet::randomizeWeightsUniform(double theta_range, double J_range, bool attractive){
  for(node *n: nodes) n->theta=rnd.uni(-theta_range, theta_range);
  if(attractive){
    for(edge *e: edges) e->J = rnd.uni(0., J_range);
  }else{
    for(edge *e: edges) e->J = rnd.uni(-J_range, J_range);
  }
  for(edge *e: edges) e->tanhJ = TANH(e->J);
  zeroMessages();
}

void BinaryBPNet::randomizeWeightsGauss(double theta_sig, double J_sig, bool attractive){
  for(node *n: nodes) n->theta = theta_sig*rnd.gauss();
  for(edge *e: edges) e->J = J_sig*rnd.gauss();
  if(attractive) for(edge *e: edges) e->J = fabs(e->J);
  for(edge *e: edges) e->tanhJ = TANH(e->J);
  zeroMessages();
}

void BinaryBPNet::zeroMessages(){
  NIY//graphConnectUndirected(nodes, edges);
  for(node *n: nodes){
    n->b = n->theta;
    n->conditioned=false;
  }
  for(edge *e: edges){
    e->mf=e->mb=0.;
  }
  steps=0;
}

double BinaryBPNet::addInputEvidence(const arr& input, bool condition){
  double logZ = 0.;
  double th;
  for(uint i=0; i<input.N; i++){
    if(!condition){
      th=p_to_ratio(input.elem(i));
      nodes(i)->theta += th;
      logZ += -log(2.*cosh(th));
    }else{
      nodes(i)->theta = p_to_ratio(input.elem(i));
      nodes(i)->conditioned = true;
    }
  }
  return logZ;
}

double BinaryBPNet::addOutputEvidence(const arr& output, bool condition){
  double logZ = 0.;
  double th;
  for(uint i=0; i<output.N; i++){
    if(!condition){
      th=p_to_ratio(output.elem(i));
      nodes(nodes.N-output.N+i)->theta += th;
      logZ += -log(2.*cosh(th));
    }else{
      nodes(nodes.N-output.N+i)->theta = p_to_ratio(output.elem(i));
      nodes(nodes.N-output.N+i)->conditioned=true;
    }
  }
  return logZ;
}
/*void unconditionOutput(const arr& output){
for(uint i=0;i<output.N;i++){
bp.nodes(bp.nodes.N-output.N+i)->theta = 0.;
bp.nodes(bp.nodes.N-output.N+i)->conditioned=false;
}
}*/

void BinaryBPNet::zeroDeltas(uint dim){
#ifdef GradTypeArr
#  define ZERO(z) z.resize(dim); z.setZero();
#else
#  define ZERO(z) z=0.;
#endif
  for(node *n: nodes){ ZERO(n->delT) ZERO(n->G) ZERO(n->G_sum) }
  for(edge *e: edges){ ZERO(e->delJ) ZERO(e->delf) ZERO(e->delb) ZERO(e->Gf) ZERO(e->Gb) ZERO(e->G_sum) }
#undef ZERO
  firstGrad=true;
}

void BinaryBPNet::addBeliefDeltas(const arr& delta){
#ifndef GradTypeArr
  CHECK_EQ(delta.N,nodes.N, "");
  for(auto n:enumerated(nodes)){
    n->delT += delta(n.count);
    for(edge *e: n->edges){
      if(e->to==n()) e->delf += delta(n.count);
      else         e->delb += delta(n.count);
    }
  }
#else
  HALT("this feature requires NO GradTypeArr compile flag");
#endif
}

void BinaryBPNet::addThetaDeltas(const arr& delta){
#ifndef GradTypeArr
  CHECK_EQ(delta.N,nodes.N, "");
  for(auto n: enumerated(nodes)) n->delT += delta(n.count);
#else
  HALT("this feature requires NO GradTypeArr compile flag");
#endif
}

void BinaryBPNet::addJDeltas(const arr& delta){
#ifndef GradTypeArr
  CHECK_EQ(delta.N,edges.N, "");
  for(auto e: enumerated(edges)) e->delJ += delta(e.count);
#else
  HALT("this feature requires NO GradTypeArr compile flag");
#endif
}

void BinaryBPNet::setNodeBeliefDeltas(){
#ifdef GradTypeArr
  zeroDeltas(nodes.N);
  for(node *n: nodes){
    n->delT(i) += 1.;
    for_list(Type,  e,  n->edges){
      if(e->to==n) e->delf(i) += 1.;
      else         e->delb(i) += 1.;
    }
  }
#else
  HALT("this feature requires GradTypeArr compile flag");
#endif
}

void BinaryBPNet::setPairBeliefDeltas(){
#ifdef GradTypeArr
  zeroDeltas(edges.N);
  for(edge *e: edges){ //pair
    double tanh_b = TANH(e->mf) * TANH(e->mb) / e->tanhJ;
    //e->from->delT(i) = (e->tanhJ * (1.-rai::sqr(TANH(collect(e->from, e)))) * TANH(collect(e->to, e))) / (1.-rai::sqr(tanh_b));
    //e->to  ->delT(i) = (e->tanhJ * TANH(collect(e->from, e)) * (1.-rai::sqr(TANH(collect(e->to, e))))) / (1.-rai::sqr(tanh_b));
    e->delJ(i)       = - tanh_b / e->tanhJ * RATIO(e->J, tanh_b);
    //e->delJ(i)       = tanh_b / e->tanhJ * RATIO(e->J, tanh_b);
    e->delf(i)       = TANH(e->mb) / e->tanhJ * RATIO(e->mf, tanh_b);
    e->delb(i)       = TANH(e->mf) / e->tanhJ * RATIO(e->mb, tanh_b);
  }
#else
  HALT("this feature requires GradTypeArr compile flag");
#endif
}

void BinaryBPNet::stepBP(){
#if 0
  for(edge *e: edges){ e->mf_old=e->mf;  e->mb_old=e->mb; }
  for(node *n: nodes){ n->b=collect(n); }
  for(edge *e: edges){
    e->mf = ATANH(e->tanhJ * TANH(e->from->b-e->mb_old));
    e->mb = ATANH(e->tanhJ * TANH(e->to  ->b-e->mf_old));
  }
#else
  for(node *n: nodes){ n->b_old[steps]=n->b; n->b=n->theta; }
  for(edge *e: edges){
    e->mf_old[steps] = e->mf;
    e->mb_old[steps] = e->mb;
    e->to->b   += e->mf = ATANH(e->tanhJ * TANH(e->from->b_old[steps]-e->mb_old[steps]));
    e->from->b += e->mb = ATANH(e->tanhJ * TANH(e->to  ->b_old[steps]-e->mf_old[steps]));
  }
#endif
  steps++;
}

void BinaryBPNet::stepGradBP(){
#if 0
  GradType e__Gf;
  for(node *n: nodes){ n->G_old[steps-1]=n->G; n->G=0.; }
  for(edge *e: edges){ //G_{jil} update equation
    e__Gf = e->Gf;
    e->Gf_old[steps-1] = e->Gf;
    e->Gb_old[steps-1] = e->Gb;
    e->from->G += e->Gf = e->tanhJ * RATIO(e->from->b_old[steps-1]-e->mb_old[steps-1], e->mf) * (e->to  ->G_old[steps-1] - e->Gb + e->delf);
    e->to->G   += e->Gb = e->tanhJ * RATIO(e->to  ->b_old[steps-1]-e->mf_old[steps-1], e->mb) * (e->from->G_old[steps-1] - e__Gf + e->delb);
  }
#else
  steps--;
  for(node *n: nodes){ n->G_old=n->G; n->G=0.; }
  for(edge *e: edges){ //G_{jil} update equation
    e->Gf_old = e->Gf;
    e->Gb_old = e->Gb;
    if(firstGrad){
      e->from->G += e->Gf = e->tanhJ * RATIO(e->from->b_old[steps]-e->mb_old[steps], e->mf) * e->delf;
      e->to->G   += e->Gb = e->tanhJ * RATIO(e->to  ->b_old[steps]-e->mf_old[steps], e->mb) * e->delb;
    }else{
      e->from->G += e->Gf = e->tanhJ * RATIO(e->from->b_old[steps]-e->mb_old[steps], e->mf_old[steps+1]) * (e->to  ->G_old - e->Gb_old);
      e->to->G   += e->Gb = e->tanhJ * RATIO(e->to  ->b_old[steps]-e->mf_old[steps], e->mb_old[steps+1]) * (e->from->G_old - e->Gf_old);
    }
  }
  for(node *n: nodes){  n->G_sum += n->G;  }
  for(edge *e: edges){
    if(firstGrad){
      e->G_sum += TANH(e->from->b_old[steps]-e->mb_old[steps]) * RATIO(e->J, e->mf) * e->delf;
      e->G_sum += TANH(e->to  ->b_old[steps]-e->mf_old[steps]) * RATIO(e->J, e->mb) * e->delb;
    }else{
      e->G_sum += TANH(e->from->b_old[steps]-e->mb_old[steps]) * RATIO(e->J, e->mf_old[steps+1]) * (e->to  ->G_old - e->Gb_old);
      e->G_sum += TANH(e->to  ->b_old[steps]-e->mf_old[steps]) * RATIO(e->J, e->mb_old[steps+1]) * (e->from->G_old - e->Gf_old);
    }
  }
  firstGrad=false;
#endif
}

double BinaryBPNet::nodeBelief(const node *n){
  return n->b;
}

double BinaryBPNet::edgeBelief(const edge *e){
  return ATANH(TANH(e->mf) * TANH(e->mb) / e->tanhJ);
}

void BinaryBPNet::getNodeBeliefs(arr &b){
  b.resize(nodes.N);
  for(auto n: enumerated(nodes)) b(n.count)=nodeBelief(n());
}

void BinaryBPNet::getNodeBeliefTables(arr &b, bool addOn){
  if(!addOn || !b.N){  b.resize(nodes.N, 2);  b.setZero();  }
  CHECK(b.nd==2 && b.d0==nodes.N && b.d1==2, "");
  arr P;
  for(auto n: enumerated(nodes)){
    nodeExp_to_Table(P, nodeBelief(n()));
    P/=sum(P);
    b[n.count] += P;
  }
}

void BinaryBPNet::getPairBeliefs(arr &b){
  b.resize(edges.N);
  for(auto e: enumerated(edges)) b(e.count)=edgeBelief(e());
}

void BinaryBPNet::getPairBeliefTables(arr &b, bool addOn){
  if(!addOn || !b.N){  b.resize(edges.N, 2, 2);  b.setZero();  }
  CHECK(b.nd==3 && b.d0==edges.N && b.d1==2 && b.d2==2, "");
  arr P;
  for(auto e: enumerated(edges)){
    pairExp_to_Table(P, e->J, e->from->b-e->mb, e->to->b-e->mf);
    P/=sum(P);
    b[e.count] += P;
  }
}

//===========================================================================
//
// converts
//

void BinaryBPNet::getFG(BinaryPairFG &FG, bool betheForm){
  uint N=nodes.N, E=edges.N;
  
  FG.logZ=0.;
  FG.f_i .resize(N, 2);
  FG.f_ij.resize(E, 2, 2);
  FG.edges.resize(E, 2);
  
  for(auto e: enumerated(edges)){  FG.edges(e.count, 0)=e->ifrom;  FG.edges(e.count, 1)=e->ito;  }
  if(betheForm){
    for(auto e: enumerated(edges)) pairExp_to_Table(FG.f_ij[e.count].noconst(), e->J, e->from->theta, e->to->theta);
  } else {
    for(auto e: enumerated(edges)) pairExp_to_Table(FG.f_ij[e.count].noconst(), e->J, 0., 0.);
  }
  for(auto n: enumerated(nodes)){  nodeExp_to_Table(FG.f_i[n.count].noconst(), n->theta);  }
}

void BinaryBPNet::getBeliefFG(BinaryPairFG &FG, bool addOn, bool normalized){
  uint N=nodes.N, E=edges.N;
  arr P;
  
  if(addOn){
    CHECK_EQ(FG.f_i .d0,N, "");
    CHECK_EQ(FG.f_ij.d0,E, "");
    CHECK_EQ(FG.edges.N,E, "");
  }else{
    FG.logZ=0.;
    FG.f_i .resize(N, 2);    FG.f_i .setZero();
    FG.f_ij.resize(E, 2, 2);  FG.f_ij.setZero();
    FG.edges.resize(E, 2);
    for(auto e: enumerated(edges)){  FG.edges(e.count, 0)=e->ifrom;  FG.edges(e.count, 1)=e->ito;  }
  }
  
  for(auto e: enumerated(edges)){
    pairExp_to_Table(P, e->J, e->from->b-e->mb, e->to->b-e->mf);
    if(normalized) P/=sum(P);
    FG.f_ij[e.count] += P;
  }
  
  for(auto n: enumerated(nodes)){
    nodeExp_to_Table(P, n->b);
    if(normalized) P/=sum(P);
    FG.f_i[n.count] += P;
  }
  
}

double BinaryBPNet::Bethe(){
  double logZ=0.;
  for(edge *e: edges){
    double a=e->from->b-e->mb, b=e->to->b-e->mf, thi=e->from->theta, thj=e->to->theta;
    double tJ=e->tanhJ, ta=tanh(a), tb=tanh(b);
    double cJ=cosh(e->J), ca=cosh(a), cb=cosh(b);
    double match=1.+tJ*ta*tb;
    double tmji = tJ*tb, tmij = tJ*ta;
    double bi=a+atanh(tmji), bj=b+atanh(tmij);
    double tbi=tanh(bi), tbj=tanh(bj);
    logZ +=
      -log(4.*(cJ*ca*cb*match))
      + (a-thi)*tbi
      + (b-thj)*tbj;
  }
  for(node *n: nodes){
    double b=n->b, th=n->theta;
    logZ += (1.-n->edges.N) * (-log(2.*cosh(b)) + (b-th)*tanh(b));
  }
  return logZ;
}

void BinaryBPNet::addBetheGradientDeltas(){
#ifdef GradTypeArr
  HALT("this feature requires NO GradTypeArr compile flag");
#endif
  //double logZ=0.;
  for(edge *e: edges){
    double a=e->from->b-e->mb, b=e->to->b-e->mf, thi=e->from->theta, thj=e->to->theta;
    double tJ=e->tanhJ, ta=tanh(a), tb=tanh(b);
    //double cJ=cosh(e->J), ca=cosh(a), cb=cosh(b);
    double match=1.+tJ*ta*tb;
    double tmji=tJ*tb, tmij=tJ*ta;
    double bi=a+atanh(tmji), bj=b+atanh(tmij);
    double tbi=tanh(bi), tbj=tanh(bj);
    //logZ += -log(4.*(cJ*ca*cb*match)) + (a-thi)*tbi + (b-thj)*tbj;
    double one_tmji2 = 1.-tmji*tmji, one_tmij2 = 1.-tmij*tmij;
    double one_tbi2 = 1.-tbi*tbi, one_tbj2 = 1.-tbj*tbj;
    double one_tJ2 = 1.-tJ*tJ;
    double one_ta2 = 1.-ta*ta, one_tb2 = 1.-tb*tb;
    double delta_J
    = -tJ - (ta*tb*one_tJ2)/match
      + (a-thi)*one_tbi2 * (tb * one_tJ2/one_tmji2)
      + (b-thj)*one_tbj2 * (ta * one_tJ2/one_tmij2);
    double delta_a
    = -ta - (tmji*one_ta2)/match
      + (a-thi)*one_tbi2
      + (b-thj)*one_tbj2 * (tJ*one_ta2/one_tmij2);
    double delta_b
    = -tb - (tmij*one_tb2)/match
      + (a-thi)*one_tbi2 * (tJ*one_tb2/one_tmji2)
      + (b-thj)*one_tbj2;
      
    e->delJ += delta_J;
    
    node *n=e->from;
    n->delT += delta_a;
    for(edge *ee: n->edges){
      if(ee==e) continue; //exclude edge j->i
      if(ee->to==n) ee->delf += delta_a + tbi; //the tbi stems from the (a-thi)*tbi term...
      else          ee->delb += delta_a + tbi;
    }
    
    n=e->to;
    n->delT += delta_b;
    for(edge *ee: n->edges){
      if(ee==e) continue; //exclude edge i->j
      if(ee->to==n) ee->delf += delta_b + tbj;
      else          ee->delb += delta_b + tbj;
    }
  }
  for(node *n: nodes){
    double b=n->b, th=n->theta;
    double tb=tanh(b);
    //logZ += (1.-n->edges.N) * (-log(2.*cosh(b)) + (b-th)*tanh(b));
    double delta_th = (1.-n->edges.N) * (-tb + (b-th)*(1.-tb*tb));
    double delta_mu = (1.-n->edges.N) * (b-th)*(1.-tb*tb);
    n->delT += delta_th;
    for(edge *e: n->edges){
      if(e->to==n) e->delf += delta_mu;
      else         e->delb += delta_mu;
    }
  }
}

void BinaryBPNet::setT(const arr &w){
  CHECK_EQ(w.N,nodes.N, "");
  for(auto n: enumerated(nodes)) n->theta = w(n.count);
}

void BinaryBPNet::setJ(const arr &w){
  CHECK_EQ(w.N,edges.N, "");
  for(auto e: enumerated(edges)){ e->J = w(e.count);  e->tanhJ=TANH(e->J);  }
}

void BinaryBPNet::getT(arr &w){
  w.resize(nodes.N);
  for(auto n: enumerated(nodes)) w(n.count) = n->theta;
}

void BinaryBPNet::getJ(arr &w){
  w.resize(edges.N);
  for(auto e: enumerated(edges)) w(e.count) = e->J;
}

void BinaryBPNet::getGradT(arr &g){
  CHECK(!steps, "");
#ifdef GradTypeArr
  g.resize(nodes.N, nodes(0)->delT.N);
  for(node *n: nodes){
    g[n.count]() = n->delT;
    g[n.count]() += n->G_sum;
    //g[i]() *= n->discount;
  }
#else
  g.resize(nodes.N);
  for(auto n: enumerated(nodes)){
    g(n.count) = n->delT;
    g(n.count) += n->G_sum;
    //g(i) *= n->discount;
  }
#endif
}

void BinaryBPNet::getGradJ(arr &g){
#ifdef GradTypeArr
  g.resize(edges.N, nodes(0)->delT.N);
  for(auto e: enumerated(edges)){
    g[e.count]() = e->delJ;
    g[e.count]() += e->G_sum;
  }
#else
  g.resize(edges.N);
  for(auto e: enumerated(edges)){
    g(e.count) = e->delJ;
    g(e.count) += e->G_sum;
  }
#endif
}


void BinaryBPNet::getPerturbationGradT(arr &g){
#ifdef GradTypeArr
  uint i;  node *n;
  CHECK_EQ(nodes(0)->delT.N,nodes.N, "need to set node belief deltas to compute perturbation grad");
  g.resize(nodes.N);
  for(node *n: nodes){
    g(i) = n->delT(i);
    g(i) += n->G_sum(i);
  }
#else
  HALT("this feature requires GradTypeArr compile flag");
#endif
}

void BinaryBPNet::getNeighborPerturbationError(arr &g, const arr &J0, bool moduloNodeError){
#ifdef GradTypeArr
  uint i;  edge *e;
  CHECK_EQ(nodes(0)->delT.N,nodes.N, "need to set node belief deltas to compute perturbation grad");
  g.resize(edges.N, 2);
  //g.setZero();
  for(edge *e: edges){
    //-- these are the true perturbation gradients
    g(i, 0) = e->to  ->delT(e->ifrom) + e->to  ->G_sum(e->ifrom); //how does 'to' vary when I perturb at 'from'
    g(i, 1) = e->from->delT(e->ito) + e->from->G_sum(e->ito);
    if(moduloNodeError){
      //-- neglect effects from the node gradient:
      g(i, 0) *= e->from->delT(e->ifrom) + e->from->G_sum(e->ifrom); //devide the effect of "how does 'from' vary when I perturb at 'from'"
      g(i, 1) *= e->to  ->delT(e->ito) + e->to  ->G_sum(e->ito);
    }
    //-- now substract the desired (theoretical) gradients
    g(i, 0) -= TANH(J0(i)) * RATIO(e->from->b-e->mb, e->mf) ; //J0 determines the desired effect
    g(i, 1) -= TANH(J0(i)) * RATIO(e->to  ->b-e->mf, e->mb) ;
  }
#else
  HALT("this feature requires GradTypeArr compile flag");
#endif
}

void BinaryBPNet::getPerturbationGradJ(arr &g){
#ifdef GradTypeArr
  uint i;  edge *e;
  CHECK_EQ(nodes(0)->delT.N,nodes.N, "need to set node belief deltas to compute perturbation grad");
  g.resize(edges.N, nodes.N);
  arr dFdm;
  for(edge *e: edges){
    g[i]() = 0.;
    g(i, e->ito)   += (1.-rai::sqr(TANH(e->from->b-e->mb))) * RATIO(e->J, e->mf) * (e->to  ->G-e->Gb)(e->ito);
    g(i, e->ifrom) += (1.-rai::sqr(TANH(e->to  ->b-e->mf))) * RATIO(e->J, e->mb) * (e->from->G-e->Gf)(e->ifrom);
  }
#else
  HALT("this feature requires GradTypeArr compile flag");
#endif
}


void BinaryBPNet::report(std::ostream &os){
  arr theta, J, b, M;
  getT(theta);
  getNodeBeliefs(b);
  J.resize(nodes.N, nodes.N);  J.setZero();
  M.resize(nodes.N, nodes.N);  M.setZero();
  for(edge *e: edges){
    J(e->ifrom, e->ito) = e->J;
    J(e->ito, e->ifrom) = e->J;
    M(e->ifrom, e->ito) = e->mf;
    M(e->ito, e->ifrom) = e->mb;
  }
  cout <<"** BinaryBPNet:"
       <<"\n  theta=  " <<theta
       <<"\n  beliefs=" <<b
       <<"\n  J=\n" <<J
       <<"\n  messages=\n" <<M <<endl;
}

void BinaryBPNet::getSamples(uintA &samples, uint S){
  arr phi(nodes.N);
  samples.resize(S, nodes.N);
  samples.setZero();
  {for(auto n: enumerated(nodes)) phi(n.count)=n->theta;}
  for(uint s=0; s<S; s++){
    {for(auto n: enumerated(nodes)) n->theta=phi(n.count);} //reset all evidences
    zeroMessages();                                 //reset all messages
    for(auto n: enumerated(nodes)){
      for(uint t=0; t<100; t++) stepBP();  //do inference for some time
      //get the belief at node i
      //translate to a probability
      if(rnd.uni()<ratio_to_p(nodeBelief(n()))) samples(s, n.count)=1; //sample a state of node i
      if(samples(s, n.count)) n->theta=1000.; else n->theta=-1000.; //set hard evidence at node i
    }
  }
  for(auto n: enumerated(nodes)) n->theta=phi(n.count); //reset all evidences
}

/////////////////////////////////////////////////////////////////////

void BinaryBPGrid::discount(float gamma){
  u*=gamma; d*=gamma; l*=gamma; r*=gamma;
}

void BinaryBPGrid::step(uint iter){
  CHECK_EQ(phi.nd,2, "");
  uint Y=phi.d0, X=phi.d1;
  phi.reshape(Y*X);
  
  //-- fwd & bwd messages
  if(u.N==Y*X){
    u.reshape(Y*X);
    d.reshape(Y*X);
    l.reshape(Y*X);
    r.reshape(Y*X);
  }else{
    u.resize(Y*X);
    d.resize(Y*X);
    l.resize(Y*X);
    r.resize(Y*X);
    u=0.; d=0.; l=0.; r=0.;
  }
  floatA ud(Y*X), lr(Y*X);
  
//#define MOD if((x&1)^(k&1))
#define MOD

  uint k, y, x, i, j;
  for(k=0; k<iter; k++){
    //cout <<'.' <<std::flush;
    //recompute all messages
    for(y=0; y<Y; y++){
      for(x=0; x<X; x++) MOD { i=y*X+x; ud(i)=phi(i)+u(i)+d(i); }
        for(x=1; x<X; x++) MOD { i=y*X+x; j=i-1;  msg_eq(r(i), r(j)+ud(j));  }
          for(x=X-1; x--;)  MOD { i=y*X+x; j=i+1;  msg_eq(l(i), l(j)+ud(j));  }
          }
    for(x=0; x<X; x++){
      for(y=0; y<Y; y++) MOD { i=y*X+x; lr(i)=phi(i)+l(i)+r(i); }
        for(y=1; y<Y; y++) MOD { i=y*X+x; j=i-X;  msg_eq(d(i), d(j)+lr(j));  }
          for(y=Y-1; y--;)  MOD { i=y*X+x; j=i+X;  msg_eq(u(i), u(j)+lr(j));  }
          }
  }
  //cout <<"done" <<endl;
  
  //-- posterior beliefs
  b.resize(Y*X);
  for(i=0; i<X*Y; i++){
    b(i) = phi(i)+u(i)+d(i)+l(i)+r(i);
  }
  phi.reshape(Y, X);
  b.reshape(Y, X);
}
