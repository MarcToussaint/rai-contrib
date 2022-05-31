#pragma once

#include <Kin/kin.h>
#include <KOMO/objective.h>
#include <Optim/NLP.h>

#include <unordered_map>

#include "Animation.h"

struct ConfigurationProblem;

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

  void getViolatedContacts(arr& y, arr& J, double margin=0.);
  arr getSideStep();
  arr getForwardStep();
  arr getBackwardStep(double relativeStepLength=1.1, double margin=0., const arr& nullRef=NoArr); //1.1 means overstepping the constraint IK a bit

  void write(ostream& os) const;
  void writeDetails(ostream& os, const ConfigurationProblem& P, double margin=0.) const;
};
stdOutPipe(QueryResult)

struct ConfigurationProblem {
  rai::Configuration C;
  arr q0, limits, max_step;
  rai::Array<shared_ptr<GroundedObjective>> objectives;

  int display=0;
  uint evals=0;
  bool computeCollisions;
  double collisionTolerance;

  // ignore collision-pairs
  bool activeOnly = true;
  bool isActuated(const rai::Frame *f);
  std::unordered_map<uint, bool> actuated;

  ConfigurationProblem(const rai::Configuration& _C, bool _computeCollisions=true, double _collisionTolerance=1e-3);

  shared_ptr<GroundedObjective> addObjective(const FeatureSymbol& feat, const StringA& frames, ObjectiveType type, const arr& scale=NoArr, const arr& target=NoArr);

  shared_ptr<QueryResult> query(const arr& x);

  bool checkConnection(const arr &start, const arr &end, const uint disc=3, const bool binary=false);
  arr sample(const arr &start={}, const arr &goal={}, const double c_max=0);
};

struct TimedConfigurationProblem : ConfigurationProblem{
  rai::Animation A;

  TimedConfigurationProblem(const rai::Configuration& _C, const rai::Animation &_A);
  ConfigurationProblem getConfigurationProblemAtTime(const double t);

  ptr<QueryResult> query(const arr& x, const double t, const double tMin=0);
  bool checkEdge(const arr& x0, const double t0, const arr &x1, const double t1, const uint discretization=3);
};

struct GoalStateProblem : NLP {
  ConfigurationProblem &P;
  double scaleCollisions = 1e1;
  int nEq=-1, nIneq=-1;

  GoalStateProblem(ConfigurationProblem& _P) : P(_P) {}

  virtual uint getDimension(){ return P.q0.N; }
  virtual void getFeatureTypes(ObjectiveTypeA& tt);
  virtual void evaluate(arr& phi, arr& J, const arr& x);
};

bool makePoseFeasible(arr& x, ConfigurationProblem& P, double IKstepSize=1.1, double maxQStepSize=.1, uint trials=3);


