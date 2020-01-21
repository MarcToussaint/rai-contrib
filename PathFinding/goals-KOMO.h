#include "methods.h"

struct GoalSampler_KOMO : GoalSampler {
  GoalSampler_KOMO(PathProblem& p) : GoalSampler(p) {}

  ptr<Result> run(double noise=.1);
};
