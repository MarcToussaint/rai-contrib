#pragma once

#include <Kin/kin.h>
#include <KOMO/objective.h>

struct QueryResult{
  //goal features
  arr goal_y, goal_J;

  //collision features
  uintA collisions;
  arr coll_y, coll_J;
  arr normal_y, normal_J;
  arr side_J;

  bool isGoal=true;
  bool isFeasible=true;

  //optional a 3D coordinate for display
  arr disp3d;

  arr getSideStep();
  arr getForwardStep();
  arr getBackwardStep(double relativeStepLength=1.1); //1.1 means overstepping the constraint IK a bit

  void write(ostream& os) const;
  void writeDetails(ostream& os) const;
};
stdOutPipe(QueryResult)

struct PathProblem {
  rai::Configuration C;
  arr q0, limits, max_step;
  rai::Array<Objective*> goalObjectives;
  rai::Array<Objective*> pathObjectives;

  int display=0;
  uint evals=0;
  bool computeCollisions;

  PathProblem(const rai::Configuration& _C, bool _computeCollisions=true);

  void resetStart(const arr& _q0);

  Objective* addObjective(const FeatureSymbol& feat, const StringA& frames, ObjectiveType type, const arr& scale=NoArr, const arr& target=NoArr);

  ptr<QueryResult> query(const arr& x);

  //helpers
  void cutStepToMax(arr& d);
};


struct GoalStateProblem : ConstrainedProblem {
  PathProblem &P;
  double scaleCollisions = 1e1;

  GoalStateProblem(PathProblem& _P) : P(_P) {}

  virtual void phi(arr& phi, arr& J, arr& H, ObjectiveTypeA& tt, const arr& x);
};

struct BallPathProblem : PathProblem {
  BallPathProblem(PathProblem& org, double d);
};

struct PathProblem_Simplification : PathProblem {
  PathProblem_Simplification(PathProblem& org, uint numCutLinks);
};

arr goalOptim(PathProblem& P);

bool makePoseFeasible(arr& x, PathProblem& P, uint trials=3, double IKstepSize=1.1);

