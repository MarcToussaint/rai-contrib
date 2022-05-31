
#include <KOMO/objective.h>
#include <KOMO/komo.h>
#include <Kin/viewer.h>

#include "smoother.h"

const double pathLength(const arr& path){
  double cost = 0.;
  for (uint i=0; i<path.d0-1; ++i){
    cost += length(path[i] - path[i+1]);
  }

  return cost;
}

/*
arr PathSmoother_FullOptimizer::run(){
  KOMO komo(C, true);

  komo.setTiming(1, path.d0, totalDuration, 2);
  komo.setSquaredQAccVelHoming(0, -1, 1, 0, 0);
  
  komo.add_collision(true);

  // add objectives
  komo.addObjective({1, -1},
                FS_qItself,
                {},
                OT_eq,
                {1e1},
                path[path.d0-1]);
  komo.addObjective({1, -1}, FS_qItself, {}, OT_eq, {}, {}, 1);

  komo.reset();
  
  komo.setConfiguration(-2, path[0]);
  komo.setConfiguration(-1, path[0]);
  for(uint i=0; i<path.d0; ++i){
    komo.setConfiguration(i, path[i]);
  }

  std::cout << "test" << std::endl;
  options.stopIters = 10000;
  komo.run(options);

  //komo.getReport(true, 1);
  //for(uint i=0;i<2;i++) {
  //  komo.displayTrajectory(.5, true);
  //}

  std::cout << "test" << std::endl;
  const arrA optPath = komo.getPath_q();
  arr smoothed(optPath.d0, path.d1);
  for(uint i=0; i<optPath.d0; ++i){
    smoothed[i] = optPath(i);
  }

  return smoothed;
}
*/
arr PathSmoother_RH::run(int verbose){
  bool disp = false;
  std::unique_ptr<rai::ConfigurationViewer> V;

  arr smoothed = initialPath;

  const double horizonDuration = totalDuration / smoothed.d0 * horizon;
  
  // TODO: Check better iteration termination criteria
  //Marc: perhaps have different parameters for the first optimization (conservative), and rest (with low damping)
  options.stopIters = 10;
//  options.damping = 1e-1;

  // set up KOMO problem for part of the path
  KOMO komo;
  komo.setModel(P.C, true);
  komo.setTiming(1., horizon, horizonDuration, 2);
  komo.opt.verbose = verbose;

  CHECK_EQ(komo.T, horizon, "");

  if(disp){
    komo.opt.verbose = 6;
  }

  for(uint i=1; i<=initialPath.d0 - horizon; ++i){
    if(verbose>1) LOG(0) << "Smoother Iteration " << i;

    komo.add_collision(true);
    komo.add_qControlObjective({}, 2, 1.);
    //-- set prefix configurations
    if(i<=1){
      komo.setConfiguration_qOrg(-2, smoothed[0]);
      komo.setConfiguration_qOrg(-1, smoothed[0]);
    }else{
      komo.setConfiguration_qOrg(-2, smoothed[i-2]);
      komo.setConfiguration_qOrg(-1, smoothed[i-1]);
    }
    //-- initialize other configurations with the previous horizon
    for(uint j=0; j<horizon; ++j){
      komo.setConfiguration_qOrg(j, smoothed[i+j]);
    }
    komo.run_prepare(0.); //mt: only this makes komo actually adopt the set configurations as decision variable komo.x!!! TODO: setConfiguration should komo.x.clear()...!

    // final obj
    komo.addObjective({1}, FS_qItself, {}, OT_eq, {1e1}, smoothed[i+komo.T-1]);

    // for velocity=0 at the end
    komo.addObjective({1}, FS_qItself, {}, OT_eq, {1e1}, NoArr, 1);
    

    komo.run(options);
    
    // get results from komo
    for(uint j=0; j<horizon; ++j){
      smoothed[i+j] = komo.getConfiguration_qOrg(j);
    }
    
    if(disp){
      std::cout << komo.getReport(true, 0) << std::endl;
      if(!V) V = make_unique<rai::ConfigurationViewer>();
      V->setConfiguration(komo.pathConfig, "smoothing...", true);
//      rai::wait();
    }

    komo.clearObjectives();
  }

  //copy last configuration, to ensure exact goal
  smoothed[-1] = initialPath[-1];

  return smoothed;
}

arr PathSmoother_Perturber::run(){
  arr smoothedPath = initialPath;

  uint no_improvement_counter = 0;
  const uint max_iter = 10000;
  for(uint i=0; i<max_iter; ++i){
    arr perturbation(smoothedPath.d1);
    //rndUniform(perturbation, -0.2, .2);
    rndGauss(perturbation, 0.1);

    //const uint k = rand() % (path.d0-2) + 1;
    const uint k = i % (initialPath.d0-2) + 1;
    arr p = smoothedPath[k] + perturbation;
    
    // check if perturbation is feasible
    auto qr = P.query(p);
    if(!qr->isFeasible){
      continue;
    }

    // check if cost decreased
    arr perturbedPath(3, smoothedPath.d1);
    perturbedPath.append(smoothedPath[k-1]);
    perturbedPath.append(p);
    perturbedPath.append(smoothedPath[k+1]);

    arr partialPath(3, smoothedPath.d1);
    partialPath.append(smoothedPath[k-1]);
    partialPath.append(smoothedPath[k]);
    partialPath.append(smoothedPath[k+1]);

    double perturbedCost = pathLength(perturbedPath);
    double cost = pathLength(partialPath);

    if(cost > perturbedCost){
      smoothedPath[k] = p;
      no_improvement_counter = 0;
    }
    else{
      no_improvement_counter++;
    }

    // proxy measure for convergence
    if(no_improvement_counter >= 100){
      std::cout << "converged after " << k << " iterations" << std::endl;
      break;
    }

  }

  return smoothedPath;
}

arr constructShortcutPath(const arr& path, 
    const uint i, 
    const uint j, 
    const std::vector<uint> short_ind){
  arr p(j-i, path.d1);

  for(uint l=0; l<j-i; ++l){
    for(uint k=0; k<path.d1; ++k){
      if(std::find(short_ind.begin(), short_ind.end(), k) != short_ind.end()){
        const double a = 1. * l / (j-i);
        p(l, k) = path(i, k) + a * (path(j, k) - path(i, k));
      }
      else{
        p(l, k) = path(l+i, k);
      }
    }
  }

  return p;
}

arr PathSmoother_PartialShortcutter::run(){
  arr smoothedPath = initialPath;

  std::vector<double> costs;
  costs.push_back(pathLength(initialPath));

  const uint max_iter = 10000;
  const uint numPoints = 5;
  for(uint k=0; k<max_iter; ++k){
    // choose random indices
    int i, j;
    do{
      i = rand() % initialPath.d0;
      j = rand() % initialPath.d0;
    } while(abs(j - i) <= 1);

    if(i > j){
      std::swap(i, j);
    }

    // choose, which indices to shortcut
    std::vector<uint> ind;
    for(uint q=0; q<smoothedPath.d1; ++q){
      const double r = 1. * rand() / RAND_MAX;
      if(r > 1./smoothedPath.d1){
        ind.push_back(q);
      }
    }
    
    // construct the new path
    auto p = constructShortcutPath(smoothedPath, i, j, {});
    auto ps = constructShortcutPath(smoothedPath, i, j, ind);

    bool shortcutFeasible = true;

    // check if the new path is feasible (interpolate)
    for(uint n=0; (int)n<j-i-1; ++n){
      const arr dir = ps[n+1] - ps[n];
      for(uint l=0; l<numPoints; ++l){
        arr point = ps[n] + dir * static_cast<double>(1.0*l/numPoints);
        auto qr = P.query(point);

        if(!qr->isFeasible){
          shortcutFeasible = false;
          break;
        }
      }
    }
   
    // check if the new path is shorter
    if(shortcutFeasible && pathLength(p) > pathLength(ps)){
      for(uint n=1; (int)n<j-i; ++n){
        smoothedPath[i+n] = ps[n];
      }
    }
    
    auto c = pathLength(smoothedPath);
    costs.push_back(c);

    // proxy measure for convergence
    const uint conv = 100;
    if(k > conv && abs(costs.back() - costs[costs.size()-conv-1]) < 1e-6){
      std::cout << "converged after " << k << " iterations" << std::endl;
      break;
    }
  }

  return smoothedPath;
}


arr PathSmoother_Shortcutter::run(){
  // TODO: This implementation is really inefficient
  arr smoothedPath = initialPath;

  std::vector<double> costs;
  costs.push_back(pathLength(initialPath));

  const uint max_iter = 10000;
  const uint numPoints = 2;
  for(uint k=0; k<max_iter; ++k){
    // choose random indices
    int i, j;
    do{
      i = rand() % initialPath.d0;
      j = rand() % initialPath.d0;
    } while(abs(j - i) <= 1);

    if(i > j){
      std::swap(i, j);
    }

    const arr p1 = smoothedPath[i];
    const arr p2 = smoothedPath[j];

    arr p = constructShortcutPath(smoothedPath, i, j, {});

    // check if the direct connection is feasible by sampling the path
    bool shortcutFeasible = true;
    const arr dir = p2 - p1;

    for(uint n=1; n<=numPoints * (j-i); ++n){
      arr point = p1 + dir * static_cast<double>(1.0*n/(numPoints * (j-i)));
      auto qr = P.query(point);

      if(!qr->isFeasible){
        shortcutFeasible = false;
        break;
      }
    }

    if(shortcutFeasible && pathLength(p) > length(p2 - p1)){
      for(uint n=1; (int)n<j-i; ++n){
        smoothedPath[i+n] = p1 + dir * static_cast<double>(1.0*n/(j-i));
      }
    }
    auto c = pathLength(smoothedPath);
    costs.push_back(c);

    // proxy measure for convergence
    const uint conv = 100;
    if(k > conv && abs(costs.back() - costs[costs.size()-conv-1]) < 1e-6){
      std::cout << "converged after " << k << " iterations" << std::endl;
      break;
    }
  }

  return smoothedPath;
}
