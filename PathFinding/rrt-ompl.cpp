#include "rrt-ompl.h"
#include "pathProblem_OMPL.h"

struct RobotSpace : public ob::RealVectorStateSpace{
  rai::Configuration &C;
  RobotSpace(rai::Configuration &_C, const int _dim) : RealVectorStateSpace(_dim), C(_C) {}

  // implement check for all vertices on how much they moved
  double distance(const ob::State *s1, const ob::State *s2) const override{
#if 0
    // set joint position and get position
    std::vector<std::string> frame_names{"KUKA1_gripperCenter",
                                         "KUKA1_mobileXY",
                                         "KUKA1_iiwa_joint_1",
                                         "KUKA1_iiwa_joint_2",
                                         "KUKA1_iiwa_joint_3",
                                         "KUKA1_iiwa_joint_4",
                                         "KUKA1_iiwa_joint_5",
                                         "KUKA1_iiwa_joint_6"};

    double max_dist = 0;
    for(auto& name: frame_names){
      auto p1 = getPos(s1, name);
      auto p2 = getPos(s2, name);

      double dist = rai::sqrDistance(p1, p2);
      if(dist > max_dist) max_dist = dist;
    }

    return max_dist;
#else
    return RealVectorStateSpace::distance(s1, s2);
#endif
  }

  arr getPos(const ob::State *s, std::string frame_name) const{
    arr c1(9);

    for(unsigned int i=0; i<9; ++i){
      c1(i) = s->as<ob::RealVectorStateSpace::StateType>()->values[i];
    }

    C.setJointState(c1);

    auto f = C.getFrame(frame_name.c_str());
    auto p1 = f->getPosition();

    return p1;
  }
};

ptr<PathResult> PathFinder_OMPL::run(double timeBudget) {
  //---------- OMPL setup
  ompl::RNG::setSeed(42);

  uint dim = P.q0.N;
  auto space = std::make_shared<RobotSpace>(P.C, dim);

  ob::RealVectorBounds bounds(dim);

  for(uint i=0; i < P.limits.d0; ++i){
    if(P.limits(i, 0) == P.limits(i, 1)){
      bounds.low[i] = -10.;
      bounds.high[i] = 10.;
    }
    else{
      bounds.low[i] = P.limits(i, 0);
      bounds.high[i] = P.limits(i, 1);
    }
  }

  //TODO
//  bounds.setLow(-4);
//  bounds.setHigh(4);
  LOG(0) <<"BOUNDS:" <<P.limits <<endl;

  space->setBounds(bounds);

  ob::SpaceInformationPtr si = std::make_shared<ob::SpaceInformation>(space);

  si->setStateValidityChecker(std::make_shared<PathProblem_OMPL>(P, si));
  si->setStateValidityCheckingResolution(0.03); // 3%
  si->setup();

  ob::ScopedState<> start(space);
  start.operator=(P.q0.vec());

  ob::ScopedState<> goal(space);
  if(goals.nd==1) goals.reshape(1,-1);
  goal.operator=(goals[0].vec());

  LOG(0) <<"GOAL:" <<goals[0] <<endl;
  LOG(0) <<"GOAL query:" <<*P.query(goals[0]) <<endl;

  auto pdef(std::make_shared<ob::ProblemDefinition>(si));
  pdef->setGoalState(start);
  for(uint i=0;i<goals.d0;i++){
    goal.operator=(goals[0].vec());
    pdef->addStartState(goal);
  }
//  pdef->setStartAndGoalStates(start, goal);

  std::vector<ob::SpaceInformationPtr> siVec;
  siVec.push_back(si);
  auto planner(std::make_shared<og::RRTConnect>(si));
  //auto planner(std::make_shared<ompl::multilevel::QRRT>(siVec));
  //auto planner(std::make_shared<og::RRTstar>(si));
  //auto planner(std::make_shared<og::InformedRRTstar>(si));
  //auto planner(std::make_shared<og::PRM>(si));
  //auto planner(std::make_shared<og::PRMstar>(si));
  //auto planner(std::make_shared<og::BITstar>(si));

  //planner->setIntermediateStates(true);
  planner->setRange(.5);
  //planner->setGoalBias(0.5);

  planner->setProblemDefinition(pdef);
  planner->setup();

  // std::cout << "Range: " << planner->getRange() << std::endl;
  // planner->setRange(stepsize);


  //auto termcond = ob::exactSolnPlannerTerminationCondition(pdef);
  auto termcond = ob::timedPlannerTerminationCondition(10.0);
  ob::PlannerStatus solved = planner->solve(termcond);

  if (solved){
    // get the goal representation from the problem definition (not the same as the goal state)
    // and inquire about the found path
    ob::PathPtr path = pdef->getSolutionPath();

    auto &path2 = static_cast<og::PathGeometric&>(*path);

    /*auto ps = og::PathSimplifier(si);
      ps.smoothBSpline(path2, 3);

        bool simplify = ps.simplifyMax(path2);
        if(!simplify){
          std::cout << "simplification failed" << std::endl;
        }*/

    arr patharr(path2.getStateCount(), P.q0.N);
    for(uint i=0;i<patharr.d0;i++){
      for(uint j=0;j<patharr.d1;j++){
        patharr(i,j) = path2.getState(i)->as<ob::RealVectorStateSpace::StateType>()->values[j];

      }
    }

    return make_shared<PathResult>(patharr);
  }

  return make_shared<PathResult>(false);
}
