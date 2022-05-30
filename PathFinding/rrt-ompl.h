#include "methods.h"

struct PathFinder_OMPL : PathFinder {
  double stepsize;
  uint verbose;

  uint n_backStep=0, n_backStepGood=0, n_sideStep=0, n_sideStepGood=0, n_forwardStep=0, n_forwardStepGood=0, n_rndStep=0, n_rndStepGood=0;


  PathFinder_OMPL(ConfigurationProblem& _P, const arr& goals, double _stepsize = .2, uint _verbose=0) : PathFinder(_P, {}, goals), stepsize(_stepsize), verbose(_verbose) {}

  virtual ptr<PathResult> run(double timeBudget=1.);

private:
  rai::Configuration DISP;
};

