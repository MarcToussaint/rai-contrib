#pragma once

#include "methods.h"
#include "ConfigurationProblem.h"

#include <Kin/kin.h>
#include <Kin/frame.h>
#include <Algo/astar.h>
#include <Core/array.h>

#include <utility>
#include <stdlib.h>


struct PathSmoother_RH : PathSmoother {
  PathSmoother_RH(ConfigurationProblem& _P, const arr& initialPath, double _duration=20, uint _horizon=10)
    :PathSmoother(_P, initialPath),
    options(NOOPT),
    horizon(_horizon),
    totalDuration(_duration){
      if(initialPath.d0 < horizon){
        horizon = initialPath.d0;
      }
    };

  arr run(int verbose=0);

  rai::OptOptions options;
  uint horizon;
  double totalDuration;
};


struct PathSmoother_PartialShortcutter: public PathSmoother{
  PathSmoother_PartialShortcutter(ConfigurationProblem& _P, const arr& path): PathSmoother(_P, path){};
  
  arr run();
};

struct PathSmoother_Shortcutter: public PathSmoother{
  PathSmoother_Shortcutter(ConfigurationProblem& _P, const arr& path): PathSmoother(_P, path){};
  arr run();
};

struct PathSmoother_Perturber: public PathSmoother{
  PathSmoother_Perturber(ConfigurationProblem& _P, const arr& path): PathSmoother(_P, path){};
  arr run();
};

//struct PathSmoother_Kalman: public PathSmoother{
//};
