#pragma once

#include "computeNode.h"

#include <Core/array.h>

//===========================================================================
// problem

struct CB_Node {
  uint ID=0;
  CB_Node *parent=0;
  shared_ptr<rai::ComputeNode> comp;

  bool isClosed=false; //all possible children are complete
  double c=0.;     //cost invested into completion of THIS node
  double l=-1.;    //lower bound computed at completion
  uint R=0;        //#children
  uint c_children=0; //#compute invested in childrent \sum_i:ch ch->c
  uint n_children=0; //#complete children
  double thetaImprovement=.1; //improvement threshold
//  arr D;           //data at leafs

  double data_Y=0., data_n=0., data_ucb=0.;
  double mean_Y=0., mean_n=0., mean_ucb=0.;
  double eff=0.;
  double comp_C=0., comp_n=0.;

  rai::Array<shared_ptr<CB_Node>> children;

  CB_Node(CB_Node* _parent, shared_ptr<rai::ComputeNode> _comp);

  void write(ostream& os) const;
};
stdOutPipe(CB_Node)


//===========================================================================
// solver

struct UILE_Solver{
  CB_Node& root;

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
  double c_total=0.;

  //reporting
  shared_ptr<ofstream> fil;

  UILE_Solver(CB_Node& _root);

  void query(CB_Node *n);
  void report();

  void run(uint k){ for(uint i=0;i<k;i++) step(); }
  void step();
  double totalCost(){ return c_total + epsilon*root.data_n; }

  CB_Node* select();

  CB_Node* select_doubleMcts();

  CB_Node* select_terminalFlat();
  CB_Node* select_terminalUCT();

  CB_Node* select_compParent_IE();
  CB_Node* select_compParent_CR();
  CB_Node* select_compParent_Tree();

  CB_Node* select_computeChild(CB_Node *r);
  CB_Node* getCheapestIncompleteChild(CB_Node *r);
  double get_novelThreshold(CB_Node *r);

  void runTrivial(uint k, double maxEffortPerCompute=10.);

  public:
  void backupMeans();


  private:

};

//===========================================================================
// helpers

double ExpectedImprovement(const arr& y);
void printTree(ostream& os, CB_Node& root);
uint getDepth(CB_Node* n);
