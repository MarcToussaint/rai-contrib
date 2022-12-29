#pragma once

#include "computeNode.h"

#include <Core/array.h>

//===========================================================================
// problem

struct CB_Node {
  CB_Node *parent=0;
  std::shared_ptr<rai::ComputeNode> comp;

  bool childrenComplete=false; //all possible children are complete
  bool branchComplete=false; //all possible decsendants are complete
  uint R=0;        //#children
  uint c_children=0; //#compute invested in childrent \sum_i:ch ch->c
  uint n_children=0; //#complete children
  double comp_n=0.;
//  arr D;           //data at leafs

  double y_tot=0., y_num=0., y_ucb=0.;
  double mean_Y=0., mean_n=0., mean_ucb=0.;
  double eff=0.;
  double c_tot=0., c_soFar=0.;
  double score=0.;
  bool isSelected=false, isBest=false;

  rai::Array<shared_ptr<CB_Node>> children;

  CB_Node(CB_Node* _parent, shared_ptr<rai::ComputeNode> _comp);

  void write(ostream& os) const;
};
stdOutPipe(CB_Node)


//===========================================================================
// solver

struct CBSolverOptions {
  enum SolverMethod { noMethod=0, SCE_Thresholded, SCE_RoundRobin, SCE_IterativeLimited };

  RAI_PARAM_ENUM("CB/", SolverMethod, method1, SCE_Thresholded)
  RAI_PARAM_ENUM("CB/", SolverMethod, method2, SCE_RoundRobin)
  RAI_PARAM("CB/", int, verbose, 1)
  RAI_PARAM("CB/", double, gamma, 1.)
  RAI_PARAM("CB/", double, beta, 1.)
  RAI_PARAM("CB/", double, epsilon, .1)
  RAI_PARAM("CB/", double, theta, .1)
  RAI_PARAM("CB/", double, rr_sampleFreq, 10.)
  RAI_PARAM("CB/", double, rr_computeFreq, 3.)
};

struct UILE_Solver{
  CB_Node root;

  rai::Array<CB_Node*> all;
  rai::Array<CB_Node*> terminals;
  rai::Array<CB_Node*> nonTerminals;
  rai::Array<CB_Node*> solutions;

  //parameter
  CBSolverOptions opt;

  //variables
  uint steps=0;
  double y_baseline = -1., y_now=-1, c_now=0.;
  double regret=0.;

  //reporting
  shared_ptr<ofstream> fil;
  uint filLast=0;

  UILE_Solver(const shared_ptr<rai::ComputeNode>& _root);

  void step();
  void run(double costLimit){ costLimit += totalCost(); while(totalCost()<costLimit) step(); }
  void runTrivial(uint k, double maxEffortPerCompute=10.);
  void report();
  double totalCost(){ return root.c_tot + opt.epsilon*root.y_num; }

private:

  void query(CB_Node *n);
  CB_Node* select_Thresholded();
  CB_Node* select_RoundRobin();
  CB_Node* selectBestCompute_IterativeLimited();
  CB_Node* selectBestCompute_RoundRobin();

  void clearScores();
  CB_Node* getBestCompute();
  CB_Node* getBestExpand();
  CB_Node* getBestSample_Flat();
  CB_Node* getBestSample_UCT();

  CB_Node* getCheapestIncompleteChild(CB_Node *r);

  //--
  uint rr_sample=0, rr_compute=0, rr_compComp=0, rr_compExp=0;
  rai::Array<CB_Node*> lifo;
  uint limit_R=0;
  double limit_c=1.;
  rai::Array<CB_Node*> rr_computeFifo;


};

//===========================================================================
// helpers

void printTree(ostream& os, CB_Node& root);
template<class T> uint getDepth(T* n){
  int i=0;
  while(n->parent){ n=n->parent; i++; }
  return i;
}
