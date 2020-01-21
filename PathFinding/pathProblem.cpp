#include <Kin/proxy.h>
#include <Optim/constrained.h>

#include "pathProblem.h"

#define PENETRATION_TOLERANCE 1e-3 //1e-3

PathProblem::PathProblem(const rai::Configuration& _C, bool _computeCollisions)
  : C(_C), computeCollisions(_computeCollisions) {

  q0 = C.getJointState();
  limits = C.getLimits();
  max_step = zeros(limits.d0);

  for(rai::Joint *j: C.activeJoints) {
    uint i=j->qIndex;
    uint d=j->qDim();
    if(d){
      switch(j->type) {
        case rai::JT_transXY:
        case rai::JT_transXYPhi:
        case rai::JT_free:
        case rai::JT_transX:
        case rai::JT_transZ:
        case rai::JT_trans3:
        case rai::JT_quatBall: 
        case rai::JT_hingeX: 
        case rai::JT_hingeY: 
        case rai::JT_hingeZ: 
        case rai::JT_rigid: {
          for(uint k=0; k<d; k++) max_step(i+k) = 1.;
        } break;

        default: NIY
      };
    }
  }
  // cout <<max_step <<endl;
}

Objective* PathProblem::addObjective(const FeatureSymbol& feat, const StringA& frames, ObjectiveType type, const arr& scale, const arr& target){
  ptr<Feature> f = symbols2feature(feat, frames, C, scale, target, 0);

  Objective *ob = new Objective(f, type, f->shortTag(C));

  goalObjectives.append(ob);
  return ob;
}

ptr<QueryResult> PathProblem::query(const arr& x){
  C.setJointState(x);
  if(computeCollisions){
    C.stepSwift();
    //  C.stepFcl();
    C.totalCollisionPenetration();
  }
//  C.reportProxies(std::cout, .001);
  evals++;
  if(display) C.watch(display>1);


  ptr<QueryResult> qr = make_shared<QueryResult>();

  //collision features
  uint N = C.proxies.N;
  qr->collisions.resize(N, 2);
  qr->coll_y.resize(N, 1);
  qr->coll_J.resize(N, 1, x.N);
  qr->normal_y.resize(N, 3);
  qr->normal_J.resize(N, 3, x.N);
  qr->side_J.resize(N, 3, x.N);

  uint i=0;
  for(const rai::Proxy& p:C.proxies){
    qr->collisions[i] =  TUP(p.a->ID, p.b->ID);
    arr Jp1, Jp2, Jx1, Jx2;
    {
      C.jacobian_pos(Jp1, C(p.a->ID), p.coll->p1);
      C.jacobian_pos(Jp2, C(p.b->ID), p.coll->p2);
      C.jacobian_angular(Jx1, C(p.a->ID));
      C.jacobian_angular(Jx2, C(p.b->ID));
    }
    p.coll->kinDistance(qr->coll_y[i](), qr->coll_J[i](), Jp1, Jp2);
    p.coll->kinNormal(qr->normal_y[i](), qr->normal_J[i](), Jp1, Jp2, Jx1, Jx2);

    arr a, b, Ja, Jb;
    {
      C.kinematicsPos(a, Ja, C(p.a->ID));
      C.kinematicsPos(b, Jb, C(p.b->ID));
    }
#if 0
    arr z = a-b;
    z /= length(z);
#else
    arr z = qr->normal_y[i];
#endif
    qr->side_J[i] = (eye(3) - (z^z))* (Ja - Jb);

    i++;
  }
  CHECK_EQ(i, N, "");
  qr->coll_J.reshape(qr->coll_y.N, x.N);

  //is feasible?
  qr->isFeasible = (!qr->coll_y.N || min(qr->coll_y)>=-PENETRATION_TOLERANCE);

  //goal features
  N=0;
  for(Objective *ob:goalObjectives) N += ob->map->__dim_phi(C);
  qr->goal_y.resize(N);
  qr->goal_J.resize(N, x.N);

  i=0;
  arr z, Jz;
  for(Objective *ob : goalObjectives){
    ob->map->__phi(z, Jz, C);
    for(uint j=0;j<z.N;j++){
      qr->goal_y(i+j) = z(j);
      qr->goal_J[i+j] = Jz[j];
    }
    i += z.N;
  }
  CHECK_EQ(i, N, "");

  //is goal?
  qr->isGoal= (absMax(qr->goal_y)<1e-2);

  //display (link of last joint)
  qr->disp3d = C.activeJoints.last()->frame->getPosition();

  return qr;
}

void PathProblem::cutStepToMax(arr& d){
  double alpha=1.;
  for(uint i=0;i<d.N;i++){
    double a = max_step.elem(i) / fabs(d.elem(i));
    if(a<alpha) alpha=a;
  }
  d *= alpha;
}

void GoalStateProblem::phi(arr& phi, arr& J, arr& H, ObjectiveTypeA& tt, const arr& x){
  auto qr = P.query(x);
  phi = qr->goal_y;
  if(!!J) J = qr->goal_J;
  if(!!tt) tt = consts<ObjectiveType>(OT_eq, qr->goal_y.N);

  if(scaleCollisions>0.){
    arr _y, _J;
    accumulateInequalities(_y, _J, -qr->coll_y, -qr->coll_J);
    phi.append( scaleCollisions*_y );
    if(!!J) J.append( scaleCollisions*_J );
    if(!!tt) tt.append( OT_ineq ); //consts<ObjectiveType>(OT_ineq, qr->coll_y.N) );
  }

  CHECK_EQ(phi.N, J.d0, "");
}

arr goalOptim(PathProblem& P){
  GoalStateProblem G(P);

  G.scaleCollisions=0.;
  arr x=P.q0;
  {
    arr dual;
    OptConstrained opt(x, dual, G, 2);
    opt.run();
  }

  G.scaleCollisions=1e0;
  {
    arr dual;
    OptConstrained opt(x, dual, G, 2);
    opt.run();
  }

  return x;
}

BallPathProblem::BallPathProblem(PathProblem& org, double d)
  : PathProblem(org.C){

  //remove actuated frames
  FrameL robots;
  for(rai::Frame *f: C.frames) if(f->joint){ robots.append(f); f->getRigidSubFrames(robots); }
  arr cen = zeros(3);
  cout <<"removing robot: ";
  for(rai::Frame *f: robots){
    cout <<f->name <<' ';
    cen += f->getPosition();
  }
  cout <<endl;
  cen/=(double)robots.N;
  while(robots.N) delete robots.popLast();

  //add a ball robot
  if(!C.getFrameByName("world", false)) C.addFrame("world");
  rai::Frame *b = C.addFrame("ball", "world");
  b->setShape(rai::ST_sphere, {d/2});
  b->setColor({1.,1.,0.});
  b->setPosition(cen);
  b->setJoint(rai::JT_trans3);
  b->setContact(1);

  q0 = C.getJointState();
  cout <<C.getJointNames() <<endl;

  arr lo = zeros(3), hi = zeros(3);
  for(rai::Frame* f: C.frames){
    arr pos = f->getPosition();
    lo = elemWiseMin(lo, pos);
    hi = elemWiseMax(hi, pos);
  }
  lo -= 2*d;
  hi += 2*d;

  limits = catCol(lo, hi);
  max_step = consts<double>(10.*d, 3);
}

PathProblem_Simplification::PathProblem_Simplification(PathProblem& org, uint numCutLinks)
  : PathProblem(org.C){

  //get frames in topologically sorted order
  FrameL sortedFrames = C.calc_topSort();

  //go reversly to collect frames to be removed
  FrameL toBeDeleted;
  for(int i=sortedFrames.N;i--;){
    if(toBeDeleted.N>=numCutLinks) break;
    rai::Frame *f = sortedFrames(i);
    if(f->joint && f->joint->dim>0) toBeDeleted.append(f);
  }

  //delete these frames and all children
  for(rai::Frame *f:toBeDeleted){
    FrameL children;
    children.append(f);
    f->getSubtree(children);
    for(rai::Frame *c:children){
      cout <<"removing frame '" <<c->name <<"'" <<endl;
      delete c;
    }
  }

  C.checkConsistency();
  cout <<C <<endl;

  q0 = C.getJointState();
  limits = C.getLimits();
  max_step = zeros(limits.d0);

  CHECK(q0.N > 0, "this simplification removed all DOFs!");

  for(rai::Joint *j: C.activeJoints) {
    uint i=j->qIndex;
    uint d=j->qDim();
    if(d){
      switch(j->type) {
        case rai::JT_transXY:
        case rai::JT_trans3:
        case rai::JT_hingeX: {
          for(uint k=0; k<d; k++) max_step(i+k) = 1.;
        } break;

        default: NIY
      };
    }
  }
  cout <<max_step <<endl;
}

arr QueryResult::getSideStep(){
  arr s = randn(3);
  s /=length(s);

  arr S(side_J.d0, 3);
  for(uint i=0;i<S.d0;i++) S[i] = s;

  arr J = side_J;

  S.reshape(-1);
  J.reshape(S.N, -1);

#if 0
  arr U, sig ,V;
  svd(U, sig, V, J);
  arr d = ~V * sig % V * randn(V.d1); //random step in input space of J!
#else
  arr JI = ~J; //pseudoInverse(J);
  arr d = JI * S;
#endif

  if(length(d)<1e-10) HALT("???");

  return d;
}

arr QueryResult::getForwardStep(){
  arr goal_JI = pseudoInverse(goal_J);
  arr d = goal_JI * (-goal_y);
  return d;
}

arr QueryResult::getBackwardStep(double relativeStepLength){
  CHECK(!isFeasible, "");
  CHECK(coll_y.N>0, "");

  uintA violated;
  for(uint i=0;i<coll_y.N;i++) if(coll_y.elem(i)<0.) violated.append(i);

  arr y = coll_y.sub(violated);
  arr J = coll_J.sub(violated);

  y.reshape(-1);
  J.reshape(y.N, -1);

  arr JI = pseudoInverse(J);
  arr d = JI * (-relativeStepLength * y); // 10% overstepping to feasibility

  return d;
}

void QueryResult::write(std::ostream& os) const{
  os <<"query: h_goal: " <<sumOfAbs(goal_y)
    <<" g_coll: " <<sum(elemWiseHinge(-coll_y))
   <<" isGoal: " <<isGoal
  <<" isFeasible: " <<isFeasible;
}

void QueryResult::writeDetails(std::ostream& os) const{
  write(os);
  for(uint i=0;i<coll_y.N;i++){
    if(coll_y.elem(i)<0.){
      os <<"\ncoll " <<i <<':' <<collisions[i] <<" y:" <<coll_y.elem(i) <<" normal:" <<normal_y[i];
    }
  }
  os <<endl;
}

bool makePoseFeasible(arr& x, PathProblem& P, uint trials, double IKstepSize){
  ptr<QueryResult> qr = P.query(x);
  for(uint k=0;k<trials;k++){
    if(qr->isFeasible) break;
    arr delta = qr->getBackwardStep(IKstepSize);
    x += delta;
    qr = P.query(x);
  }
  return qr->isFeasible;
}
