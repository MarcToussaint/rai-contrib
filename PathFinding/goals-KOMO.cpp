#include "goals-KOMO.h"

ptr<GoalSampler::Result> GoalSampler_KOMO::run(double noise){

  GoalStateProblem p_goal(p);

  arr x=p.q0;
  if(noise>0.){
    rndGauss(x, noise, true);
  }
#if 0
  p_goal.scaleCollisions=0.;
  {
    arr dual;
    OptConstrained opt(x, dual, p_goal, 2);
    opt.run();
  }
  p_goal.scaleCollisions=1e1;
#endif
  {
    arr dual;
    OptConstrained opt(x, dual, p_goal.ptr());
    opt.run();
  }

  //final query
  auto r = make_shared<Result>();
  r->qr = p.query(x);
  r->feasible = r->qr->isGoal && r->qr->isFeasible;
  r->goal = x;
  LOG(0) <<*r->qr;

  return r;
}
