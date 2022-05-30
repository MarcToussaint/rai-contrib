#include "methods.h"

struct GoalSampler_KOMO : GoalSampler {
  GoalSampler_KOMO(ConfigurationProblem& p) : GoalSampler(p) {}

  ptr<Result> run(double noise=.1);
};
