#pragma once

#include <ompl/base/StateValidityChecker.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/State.h>
#include <ompl/base/ScopedState.h>

#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/PlannerTerminationCondition.h>

#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/planners/rrt/InformedRRTstar.h>

#include <ompl/geometric/planners/bitstar/BITstar.h>

#include <ompl/geometric/planners/prm/PRM.h>
#include <ompl/geometric/planners/prm/PRMstar.h>

#include <ompl/geometric/PathGeometric.h>
#include <ompl/geometric/PathSimplifier.h>

namespace ob = ompl::base;
namespace og = ompl::geometric;

#include "pathProblem.h"

struct PathProblem_OMPL : ob::StateValidityChecker {
  PathProblem& P;
  arr x;
  uint numQueries = 0;

  PathProblem_OMPL(PathProblem& _P, const ob::SpaceInformationPtr& si)
    : ob::StateValidityChecker(si),
      P(_P) {}

  virtual ~PathProblem_OMPL() = default;

  virtual bool isValid(const ompl::base::State *state) const{
    const double* _state = state->as<ob::RealVectorStateSpace::StateType>()->values;

    ((PathProblem_OMPL*)this)->numQueries++;
    if(!(numQueries%100)) cout <<"#queries: " <<numQueries <<endl;

    arr x(_state, P.q0.N, true);
    auto qr = P.query(x);

    return qr->isFeasible;
  }
};
