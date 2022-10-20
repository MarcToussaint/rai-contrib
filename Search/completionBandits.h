#pragma once

#include "computeNode.h"

#include <Core/array.h>

//===========================================================================
// problem

struct CB_Node {
  uint ID=0;
  CB_Node *parent=0;
  std::shared_ptr<rai::ComputeNode> comp;

  bool childrenComplete=false; //all possible children are complete
  bool branchComplete=false; //all possible decsendants are complete
  uint R=0;        //#children
  uint c_children=0; //#compute invested in childrent \sum_i:ch ch->c
  uint n_children=0; //#complete children
  double thetaImprovement=.1; //improvement threshold
  double comp_n=0.;
//  arr D;           //data at leafs

  double y_tot=0., y_num=0., y_ucb=0.;
  double mean_Y=0., mean_n=0., mean_ucb=0.;
  double eff=0.;
  double c_tot=0.;
  double score=0.;
  bool isSelected=false, isBest=false;

  rai::Array<shared_ptr<CB_Node>> children;

  CB_Node(CB_Node* _parent, shared_ptr<rai::ComputeNode> _comp);

  void write(ostream& os) const;
};
stdOutPipe(CB_Node)


//===========================================================================
// solver

struct UILE_Solver{
  CB_Node root;

  rai::Array<CB_Node*> all;
  rai::Array<CB_Node*> terminals;
  rai::Array<CB_Node*> nonTerminals;

  //parameter
  int verbose=1;
  double gamma;
  double beta;
  double epsilon;
  bool useUCBasData;
  double costCoeff;

  //variables
  uint steps=0;
  double y_baseline = -1., y_now=-1;

  //reporting
  shared_ptr<ofstream> fil;

  UILE_Solver(const shared_ptr<rai::ComputeNode>& _root);

  void query(CB_Node *n);
  void report();

  void run(uint k){ for(uint i=0;i<k;i++) step(); }
  void step();
  double totalCost(){ return root.c_tot + epsilon*root.y_num; }

  CB_Node* select();
  CB_Node* selectNew();

  CB_Node* select_doubleMcts();

  CB_Node* select_terminalFlat();
  CB_Node* select_terminalUCT();

  CB_Node* select_compParent_IE();
  CB_Node* select_compParent_CR();
  CB_Node* select_compParent_Tree();

  CB_Node* select_computeChild(CB_Node *r);
  CB_Node* getCheapestIncompleteChild(CB_Node *r);
  double get_expandThreshold(CB_Node *r);

  void clearScores();
  CB_Node* getBestCompute();
  CB_Node* getBestExpand();
  CB_Node* getBestSample();


  void runTrivial(uint k, double maxEffortPerCompute=10.);

  public:
  void backupMeans();


  private:

};

//===========================================================================
// helpers

double ExpectedImprovement(const arr& y);
void printTree(ostream& os, CB_Node& root);
template<class T> uint getDepth(T* n){
  int i=0;
  while(n->parent){ n=n->parent; i++; }
  return i;
}
