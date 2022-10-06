#pragma once

#include "ConfigurationProblem.h"
#include "methods.h"

#include <Algo/ann.h>

/// just a data structure, no algorithms
struct RRT_SingleTree : GLDrawer {
  ANN ann;         //ann stores all points added to the tree in ann.X
  uintA parent;    //for each point we also store the index of the parent node
  rai::Array<ptr<QueryResult>> queries;
  arr disp3d;
  Mutex drawMutex;
  uint nearestID = UINT_MAX; //nearest node from the last 'getProposalToward' call!

  RRT_SingleTree(const arr& q0, const ptr<QueryResult>& q0_qr);

  //core method
  double getNearest(const arr& target);
  arr getProposalTowards(const arr& target, double stepsize);

  arr getNewSample(const arr& target, double stepsize, double p_sideStep, bool& isSideStep, const uint recursionDepth);

  //trivial
  uint add(const arr& q, uint parentID, const ptr<QueryResult>& _qr);

  //trivial access routines
  uint getParent(uint i){ return parent(i); }
  uint getNumberNodes(){ return ann.X.d0; }
  uint getDim(){ return ann.X.d1; }
  arr getNode(uint i){ return ann.X[i].copy(); }
  arr getLast(){ return ann.X[ann.X.d0-1].copy(); }
  arr getRandomNode(){ return ann.X[rnd(ann.X.d0)].copy(); }
  arr getPathFromNode(uint fromID);

  void glDraw(OpenGL &gl);

public:
  arr getSideStep(std::shared_ptr<QueryResult> qr);
};

//===========================================================================

///algorithms
struct PathFinder_RRT : PathFinder {
  double stepsize;
  uint maxIters=5000;
  uint verbose;
  bool intermediateCheck;

  uint n_backStep=0, n_backStepGood=0, n_sideStep=0, n_sideStepGood=0, n_forwardStep=0, n_forwardStepGood=0, n_rndStep=0, n_rndStepGood=0;


  PathFinder_RRT(ConfigurationProblem& _P, const arr& starts, const arr& goals, double _stepsize = .2, uint _verbose=0, bool _intermediateCheck=false);

  void planForward(const arr& q0, const arr& qT);
  arr planConnect(const arr& q0, const arr& qT, double p_forwardStep=.5, double p_sideStep=.0, double p_backwardStep=.0); //default numbers: equivalent to standard bidirect

  bool growTreeTowardsRandom(RRT_SingleTree& rrt);
  bool growTreeToTree(RRT_SingleTree& rrt_A, RRT_SingleTree& rrt_B, double p_forwardStep, double p_sideStep, double p_backwardStep);

  virtual ptr<PathResult> run(double timeBudget=1.);

private:
  rai::Configuration DISP;
};

void revertPath(arr& path);
