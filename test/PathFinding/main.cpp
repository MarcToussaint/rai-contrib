#include <contrib/PathFinding/rrt-ompl.h>
#include <contrib/PathFinding/goals-KOMO.h>
#include <Kin/viewer.h>

void test(){
  rai::Configuration C("scene.g");

  rai::ConfigurationViewer V;

  PathProblem P(C);
  P.addObjective(FS_positionDiff, {"endeff", "target"}, OT_eq, {1e1});

  P.display = 0;

  GoalSampler_KOMO GS(P);

  arr goals;
  for(uint k=0;k<10;k++){
    auto goal = GS.run(1.);
    C.setJointState(goal->x);
    if(goal->feasible){
      goals.append(goal->x);
      V.setConfiguration(C, "found goal", true);
    }else{
      V.setConfiguration(C, "infeasible...", true);
    }
  }
  goals.reshape(-1, P.q0.N);

  P.limits = repmat(arr(1,2,{-4,4}), P.q0.N, 1); //overwrite limits... not nice...

  PathFinder_OMPL rrt(P, goals, .2);
  auto path = rrt.run();
  V.setPath(C, path->path, "PathFinder_RRT", true);
  while(V.playVideo(true, 3.));
}


int main(int argc,char **argv){
  rai::initCmdLine(argc, argv);

  test();

  return 0;
}
