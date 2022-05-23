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
#ifndef RAI_BinaryBP_h
#define RAI_BinaryBP_h

#include <math.h>
#include <Core/array.h>

#ifdef GradTypeArr
#define GradType arr
#else
#  define GradType double
#endif

#define HIST 200

struct BinaryPairFG;  ///<minimalistic representation of a binary pair-wise FG
struct BinaryBPNet;
struct BinaryBPGrid;

double ratio_to_p(double a);
double p_to_ratio(double p);

double table_to_nodeExp(const arr& p);
void   nodeExp_to_Table(arr& p, double a);

double table_to_pairExp(const arr& p);
void   table_to_pairExp(double& J, double& a, double& b, const arr& p);
void   pairExp_to_Table(arr& p, double J, double a, double b);

double KLD_ratio(double a, double b);

double KLD(const BinaryPairFG &b, const BinaryPairFG &f);

struct BinaryPairFG { //minimalistic representation of a binary pair-wise FG
  double logZ;
  arr f_i, f_ij;
  uintA edges;
  void write(std::ostream& os);
};

struct BinaryBPNet {
  struct node;
  struct edge;
  rai::Array<node*> nodes;
  rai::Array<edge*> edges;
  uint steps;
  bool firstGrad;
  
  BinaryBPNet(){ firstGrad=false; }
  
  //-- initialization
  void randomizeWeightsUniform(double theta_range, double J_range, bool attractive);
  void randomizeWeightsGauss(double theta_sig, double J_sig, bool attractive);
  void zeroMessages();
  
  //-- adding evidence (add to thetas!!)
  double addInputEvidence(const arr& input , bool condition=false);
  double addOutputEvidence(const arr& output, bool condition=false);
  
  //-- set gradient deltas
  void zeroDeltas(uint dim);
  void addBeliefDeltas(const arr& delta);
  void addThetaDeltas(const arr& delta);
  void addJDeltas(const arr& delta);
  void setNodeBeliefDeltas();
  void setPairBeliefDeltas();
  
  //-- make one inference or grad prop step
  void stepBP();
  void stepGradBP();
  
  //set and get model parameters
  void setT(const arr &w);
  void setJ(const arr &w);
  void setTandJ(const arr &w){ setT(w.sub(0, nodes.N-1)); setJ(w.sub(nodes.N, -1)); }
  void getT(arr &w);
  void getJ(arr &w);
  void getTandJ(arr &w){ arr wJ; getT(w); getJ(wJ); w.append(wJ); }
  void getGradT(arr &g);
  void getGradJ(arr &g);
  void getGradTandJ(arr &g){ arr gJ; getGradT(g); getGradJ(gJ); g.append(gJ); g=~g; }
  
  //-- read out results
  double nodeBelief(const node *n);
  double edgeBelief(const edge *e);
  void getNodeBeliefs(arr &b);
  void getPairBeliefs(arr &b);
  void getNodeBeliefTables(arr &b, bool addOn);
  void getPairBeliefTables(arr &b, bool addOn);
  void getPerturbationGradT(arr &g);
  void getPerturbationGradJ(arr &g);
  void getNeighborPerturbationError(arr& g, const arr& J0, bool moduloNodeError);
  
  //-- interface to FG
  void getFG(BinaryPairFG &FG, bool betheForm);
  void getBeliefFG(BinaryPairFG &FG, bool addOn, bool normalized);
  
  //-- high level
  double Bethe();
  void addBetheGradientDeltas();
  
  //-- I/O
  void report(std::ostream &os);
  
  //-- generate samples
  void getSamples(uintA &samples, uint S);
};

struct BinaryBPNet::node {
  //uint id;
  rai::Array<edge*> edges;
  double theta;
  double b, b_old[HIST];
  bool conditioned;
  GradType G, G_old, G_sum;
  GradType delT;
};
struct BinaryBPNet::edge {
  uint index;
  node *from, *to;
  uint ifrom, ito;
  double J, tanhJ;          //parameter J
  double mf, mb;            //fwd & bwd messages
  double mf_old[HIST], mb_old[HIST];    //fwd & bwd messages
  GradType delJ, delf, delb; //partials del F/del J, for J, m_fwd, m_bwd
  GradType Gf, Gb;          //gradient messages fwd & bwd
  GradType G_sum;          //gradient messages fwd & bwd
  GradType Gf_old, Gb_old;          //gradient messages fwd & bwd
};


struct BinaryBPGrid {
  floatA phi, b, u, d, l, r;
  
  float tanh_J;
  
  void msg_eq(float &msg, float in){ msg=atanh(tanh_J*tanh(in)); }
  
  void init(){ phi.clear(); b.clear(); u.clear(); d.clear(); l.clear(); r.clear(); }
  void step(uint iter=1);
  void discount(float gamma=.8);
};

#ifdef  RAI_IMPLEMENTATION
#  include "BinaryBP.cpp"
#endif

#endif
