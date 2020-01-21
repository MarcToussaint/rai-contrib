#pragma once

#include "pathProblem.h"

#include <KOMO/komo.h>

struct PathResult{
  arr path;
  int feasible=-1;
  double duration=-1.;
  double cost=-1., ineq=-1., eq=-1.;
  bool ineqFeasible(double eps=1e-3){ return ineq<=eps && eq<=eps; }

  PathResult() {}
  PathResult(bool infeasible) : feasible(0) { CHECK(infeasible==false, ""); }
  PathResult(const arr& _path) : path(_path) {}

  void write(ostream& os) const{
    if(path.N) os <<" path-dim:" <<path.dim();
    if(feasible>=0) os <<" feasible:" <<feasible;
    if(duration>=0) os <<" duration:" <<duration;
    if(cost>=0) os <<" cost:" <<cost <<" ineq:" <<ineq <<" eq:" <<eq;
  }
};
stdOutPipe(PathResult)

struct PathMethod {
  PathMethod(PathProblem& _P, const ptr<PathResult>& _initialPath) : P(_P), initialPath(_initialPath) {}

  virtual ptr<PathResult> run() = 0;

protected:
  PathProblem& P;
  ptr<PathResult> initialPath;
};

struct GoalSampler {
  GoalSampler(PathProblem& p) : p(p) {}

  struct Result{
    bool feasible;
    arr x;
    ptr<QueryResult> qr;
  };
  ptr<Result> run();

protected:
  PathProblem& p;
};

struct PathInitializer {
  PathInitializer(PathProblem& p, const arr& goal);

  arr run();
};

struct PathFinder {
  PathFinder(PathProblem& _P, const arr& _goals) : P(_P), goals(_goals) {}
  virtual ~PathFinder(){}

  virtual void resetStartAndGoals(const arr& q0, const arr& goals){ NIY }

  virtual ptr<PathResult> run(double timeBudget=1.) = 0;

protected:
  PathProblem& P;
  arr goals;
};

struct PathSmoother {
  PathSmoother(PathProblem& _P, const arr& path) : P(_P), initialPath(path) {}

  arr run();

protected:
  PathProblem& P;
  const arr initialPath;
};

struct TimeOptimizer {
  TimeOptimizer(const arr& _path) : path(_path) {}

  ptr<PathResult> run(double initialDuration);

private:
  arr path;
};

struct PathOptimizer : KOMO {
  PathOptimizer(PathProblem& p, const arr& initialPath, double duration);

  struct Return{
    arr path;
    Graph report;
    double sos_sumOfSqr;
    double eq_sumOfAbs;
    double ineq_sumOfPos;
  };
  ptr<Return> run(double timeBudget=1.);

private:
  PathProblem& P;
};

struct MultiAgentKOMO : KOMO {
  uintA agentBusyTil;
  rai::Array<ptr<Objective>> allAgentObjectives;

  MultiAgentKOMO(rai::Configuration& C, uint na, bool computeCollisions=true);

  void runForSingleAgent(uint agent, PathProblem& p, const arr& initialPath, double duration);
  void pastePickAndPlaceAndRun(uint agent, uint object, const arr& goalFrameState, KOMO& seqkomo, ptr<PathResult> initPickPath, ptr<PathResult> initPlacePath);
  void runForAllAgents();
};
