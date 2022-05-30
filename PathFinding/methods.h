#pragma once

#include "ConfigurationProblem.h"

#include <KOMO/komo.h>

namespace rai{
  struct Spline;
}


void makeSpline(rai::Spline& S, const arr& path, double initialDuration);

struct PathResult{
  arr path;
  arr Xpath;
  int feasible=-1;
  double duration=-1.;
  double cost=-1., ineq=-1., eq=-1.;
  bool ineqFeasible(double eps=1e-3){ return ineq<=eps && eq<=eps; }

  PathResult() {}
  PathResult(bool infeasible) : feasible(0) { CHECK(infeasible==false, ""); }
  PathResult(const arr& _path) : feasible(1), path(_path) { if(!path.N) feasible=0; }

  void write(ostream& os) const{
    if(path.N) os <<" path-dim:" <<path.dim();
    if(feasible>=0) os <<" feasible:" <<feasible;
    if(duration>=0) os <<" duration:" <<duration;
    if(cost>=0) os <<" cost:" <<cost <<" ineq:" <<ineq <<" eq:" <<eq;
  }
};
stdOutPipe(PathResult)

struct PathMethod {
  PathMethod(ConfigurationProblem& _P, const ptr<PathResult>& _initialPath) : P(_P), initialPath(_initialPath) {}

  virtual ptr<PathResult> run() = 0;

protected:
  ConfigurationProblem& P;
  ptr<PathResult> initialPath;
};

struct GoalSampler {
  GoalSampler(ConfigurationProblem& p) : p(p) {}

  struct Result{
    bool feasible;
    arr goal;
    ptr<QueryResult> qr;
  };
  ptr<Result> run();

  ConfigurationProblem& p;
};

struct PathInitializer {
  PathInitializer(ConfigurationProblem& p, const arr& goal);

  arr run();
};

struct PathFinder {
  PathFinder(ConfigurationProblem& _P, const arr& _starts, const arr& _goals) : P(_P), starts(_starts), goals(_goals) {}
  virtual ~PathFinder(){}

  virtual void resetStartAndGoals(const arr& _starts, const arr& _goals){ NIY }

  virtual ptr<PathResult> run(double timeBudget=1.) = 0;

protected:
  ConfigurationProblem& P;
  arr starts, goals;
};

struct PathSmoother {
  PathSmoother(ConfigurationProblem& _P, const arr& path) : P(_P), initialPath(path) {}

  arr run();

protected:
  ConfigurationProblem& P;
  const arr initialPath;
};

struct TimeOptimizer {
  TimeOptimizer(const arr& _path) : path(_path) {}

  ptr<PathResult> run(double initialDuration);


public:

private:
  arr path;
};

struct PathOptimizer : KOMO {
  PathOptimizer(ConfigurationProblem& p, const arr& initialPath, double duration);

  struct Return{
    arr path;
    rai::Graph report;
    double sos_sumOfSqr;
    double eq_sumOfAbs;
    double ineq_sumOfPos;
  };
  ptr<Return> run(double timeBudget=1.);

private:
  ConfigurationProblem& P;
};

struct MultiAgentKOMO : KOMO {
  uintA agentBusyTil;
  rai::Array<ptr<Objective>> allAgentObjectives;

  MultiAgentKOMO(rai::Configuration& C, uint na, bool computeCollisions=true);

  void runForSingleAgent(uint agent, ConfigurationProblem& p, const arr& initialPath, double duration);
  void pastePickAndPlaceAndRun(uint agent, uint object, const arr& goalFrameState, ptr<PathResult> initPickPath, ptr<PathResult> initPlacePath);
  void runForAllAgents();
};
