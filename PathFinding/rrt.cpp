#include "rrt.h"

#include <GL/gl.h>
#include <Kin/viewer.h>

#include <chrono>

RRT_SingleTree::RRT_SingleTree(const arr& q0, const ptr<QueryResult>& q0_qr){
  if(!q0_qr->isFeasible) LOG(0) <<"rooting RRT with infeasible start configuration -- that's likely to fail: query is:\n" <<*q0_qr;
  add(q0, 0, q0_qr);
}

uint RRT_SingleTree::add(const arr& q, uint parentID, const ptr<QueryResult>& _qr){
  drawMutex.lock(RAI_HERE);
  ann.append(q);
  parent.append(parentID);
  queries.append(_qr);
  disp3d.append(_qr->disp3d);
  disp3d.reshape(-1,3);

  CHECK_EQ(parent.N, ann.X.d0, "");
  CHECK_EQ(queries.N, ann.X.d0, "");
  CHECK_EQ(disp3d.d0, ann.X.d0, "");
  drawMutex.unlock();
  return parent.N-1;
}

double RRT_SingleTree::getNearest(const arr& target){
  //find NN
  nearestID = ann.getNN(target);
  return length(target - ann.X[nearestID]);
}

arr RRT_SingleTree::getProposalTowards(const arr& target, double stepsize){
  //find NN
  nearestID = ann.getNN(target);

  //compute default step
  arr delta = target - ann.X[nearestID]; //difference vector between q and nearest neighbor
  double dist = length(delta);
  if(dist>stepsize)  delta *= stepsize/dist;

  return getNode(nearestID) + delta;
}

arr RRT_SingleTree::getNewSample(const arr& target, double stepsize, double p_sideStep, bool& isSideStep, const uint recursionDepth){
  //find NN
  nearestID = ann.getNN(target);
  std::shared_ptr<QueryResult> qr = queries(nearestID);

  //compute default step
  arr delta = target - getNode(nearestID);
  double dist = length(delta);
  if(dist>stepsize) delta *= stepsize/dist;

  //without side stepping, we're done
  isSideStep = false;
  if(p_sideStep<=0. || recursionDepth >= 3) return getNode(nearestID) + delta;

  //check whether this is a predicted collision
  bool predictedCollision=false;
  if(qr->coll_y.N){
    arr y = qr->coll_y + qr->coll_J * delta;
    if(min(y)<0.) predictedCollision = true;
  }

  if(predictedCollision && p_sideStep>0. && rnd.uni()<p_sideStep){
    isSideStep=true;

    //compute new target
    arr d = qr->getSideStep();
    d *= rnd.uni(stepsize,2.) / length(d);
    arr targ = getNode(nearestID) + d;
    bool tmp;
    return getNewSample(targ, stepsize, p_sideStep, tmp, recursionDepth + 1);
  }else{
    return getNode(nearestID) + delta;
  }

  HALT("shouldn't be here");
  return NoArr;
}

arr RRT_SingleTree::getPathFromNode(uint fromID){
  arr path;
  uint node = fromID;
  for(;;){
    path.append(ann.X[node]);
    if(!node) break;
    node = getParent(node);
  }
  path.reshape(-1, ann.X.d1);
  return path;
}

void RRT_SingleTree::glDraw(OpenGL& gl){
  glColor(.0, .0, .0);
  glLineWidth(2.f);
  glBegin(GL_LINES);
  drawMutex.lock(RAI_HERE);
  for(uint i=1;i<getNumberNodes();i++){
    glVertex3dv(&disp3d(parent(i),0));
    glVertex3dv(&disp3d(i,0));
  }
  drawMutex.unlock();
  glEnd();
  glLineWidth(1.f);
}

//===========================================================================

bool PathFinder_RRT::growTreeTowardsRandom(RRT_SingleTree& rrt){
  const arr start = rrt.ann.X[0];
  arr t(rrt.getNode(0).N);
  rndUniform(t,-RAI_2PI,RAI_2PI,false);

  arr q = rrt.getProposalTowards(t, stepsize);

  auto qr = P.query(q);
  if(qr->isFeasible){
    if (intermediateCheck && !P.checkConnection(start, q, 20, true)){
      return false;
    }

    rrt.add(q, rrt.nearestID, qr);
    return true;
  }
  return false;
}

bool PathFinder_RRT::growTreeToTree(RRT_SingleTree& rrt_A, RRT_SingleTree& rrt_B, double p_forwardStep, double p_sideStep, double p_backwardStep){
  bool isSideStep, isForwardStep;
  //decide on a target: forward or random
  arr t;
  if(rnd.uni()<p_forwardStep){
    t = rrt_B.getRandomNode();
    isForwardStep = true;
  }else{
    t.resize(rrt_A.getNode(0).N);
    rndUniform(t,-RAI_2PI,RAI_2PI,false);
    isForwardStep = false;
  }

  //sample configuration towards target, possibly sideStepping
  arr q = rrt_A.getNewSample(t, stepsize, p_sideStep, isSideStep, 0);

  //evaluate the sample
  auto qr = P.query(q);
  if(isForwardStep){  n_forwardStep++; if(qr->isFeasible) n_forwardStepGood++; }
  if(!isForwardStep){  n_rndStep++; if(qr->isFeasible) n_rndStepGood++; }
  if(isSideStep){  n_sideStep++; if(qr->isFeasible) n_sideStepGood++; }

  //if infeasible, make a backward step from the sample configuration
  if(!qr->isFeasible && p_backwardStep>0. && rnd.uni()<p_backwardStep){
    t = q + qr->getBackwardStep();
    q = rrt_A.getNewSample(t, stepsize, p_sideStep, isSideStep, 0);
    qr = P.query(q);
    n_backStep++; if(qr->isFeasible) n_backStepGood++;
    if(isSideStep){  n_sideStep++; if(qr->isFeasible) n_sideStepGood++; }
  };

  // TODO: add checking motion
  if(qr->isFeasible){
    const arr start = rrt_A.ann.X[rrt_A.nearestID];
    if (intermediateCheck && !P.checkConnection(start, q, 20, true)){
      return false;
    }

    rrt_A.add(q, rrt_A.nearestID, qr);
    double dist = rrt_B.getNearest(q);
    if(dist<stepsize) return true;
  }
  return false;
}

//===========================================================================

PathFinder_RRT::PathFinder_RRT(ConfigurationProblem& _P, const arr& starts, const arr& goals, double _stepsize, uint _verbose, bool _intermediateCheck)
  : PathFinder(_P, starts, goals),
    stepsize(_stepsize),
    verbose(_verbose),
    intermediateCheck(_intermediateCheck) {
}

void PathFinder_RRT::planForward(const arr& q0, const arr& qT){
  DISP.clear();
  DISP.copy(P.C);

  arr q = q0;

  RRT_SingleTree rrt(q0, P.query(q0));
  DISP.gl()->add(rrt);

  bool success=false;

  uint i;
  for(i=0;i<100000;i++){
    //let rrt0 grow
    bool added = growTreeTowardsRandom(rrt);
    if(added){
      if(length(rrt.getLast() - qT)<stepsize) success = true;
    }
    if(success) break;

    //some output
    if (verbose > 2){
      if(!(i%100)){
        DISP.setJointState(q); //updates display (makes it slow)
        DISP.watch(false);
        cout <<"RRT samples=" <<i <<" tree size = " <<rrt.getNumberNodes() <<std::endl;
      }
    }
  }

  if(!success) return;

  if (verbose > 0){
    std::cout <<"SUCCESS!"
              <<"\n  tested samples=" <<2*i
              <<"\n  #tree-size=" <<rrt.getNumberNodes()
     << std::endl;
  }

  arr path = rrt.getPathFromNode(rrt.nearestID);
  revertPath(path);

  //display
  if(verbose > 1){
    std::cout << "path-length= " << path.d0 <<endl;
    DISP.proxies.clear();

    for(uint t=0;t<path.d0;t++){
      DISP.setJointState(path[t]);
      //DISP.watch();
      DISP.watch(false);
      rai::wait(.1);
    }
  }

  path >>FILE("z.path");
}

arr PathFinder_RRT::planConnect(const arr& q0, const arr& qT, double p_forwardStep, double p_sideStep, double p_backwardStep){
  auto q0ret = P.query(q0);
  auto qTret = P.query(qT);
  if(!q0ret->isFeasible){ if(verbose>0) LOG(0) <<"initializing with infeasible q0"; if(verbose>1) q0ret->writeDetails(cout, P); return {}; }
  if(!qTret->isFeasible){ if(verbose>0) LOG(0) <<"initializing with infeasible qT"; if(verbose>1) qTret->writeDetails(cout, P); return {}; }
  RRT_SingleTree rrt0(q0, q0ret);
  RRT_SingleTree rrtT(qT, qTret);

  if(verbose>2){
    DISP.clear();
    DISP.copy(P.C);
    DISP.gl()->add(rrt0);
    DISP.gl()->add(rrtT);
  }

  bool success=false;

  uint i;
  for(i=0;i<maxIters;i++){
    success = growTreeToTree(rrt0, rrtT, p_forwardStep, p_sideStep, p_backwardStep);
    if(success) break;

    success = growTreeToTree(rrtT, rrt0, p_forwardStep, p_sideStep, p_backwardStep);
    if(success) break;

    //some output
    if(verbose>2){
      if(!(i%100)){
        DISP.setJointState(rrt0.getLast());
        DISP.watch(verbose>4, STRING("planConnect it " <<i));
        std::cout <<"RRT queries=" <<P.evals <<" tree sizes = " <<rrt0.getNumberNodes()  <<' ' <<rrtT.getNumberNodes() <<std::endl;
      }
    }
  }

  if(verbose>0){
    if(success)
      cout <<"\nSUCCESS!" <<endl;
    else
      cout <<"\nFAIL!" <<endl;
    cout <<"  RRT queries=" <<P.evals <<" tree sizes = " <<rrt0.getNumberNodes()  <<' ' <<rrtT.getNumberNodes() <<std::endl;
    cout <<"  forwardSteps: " <<(100.*n_forwardStepGood/n_forwardStep) <<"%/" <<n_forwardStep;
    cout <<"  backSteps: " <<(100.*n_backStepGood/n_backStep) <<"%/" <<n_backStep;
    cout <<"  rndSteps: " <<(100.*n_rndStepGood/n_rndStep) <<"%/" <<n_rndStep;
    cout <<"  sideSteps: " <<(100.*n_sideStepGood/n_sideStep) <<"%/" <<n_sideStep;
    cout <<endl;
  }

  if(!success) return NoArr;

  arr path = rrt0.getPathFromNode(rrt0.nearestID);
  arr pathT = rrtT.getPathFromNode(rrtT.nearestID);

  revertPath(path);
  path.append(pathT);

  //display
  if(verbose>1){
    cout <<"  path-length=" <<path.d0 <<endl;
    if(verbose>2){
      DISP.proxies.clear();
      for(uint t=0;t<path.d0;t++){
        DISP.setJointState(path[t]);
        DISP.watch(false, STRING("rrt result "<<t));
        rai::wait(.1);
      }
      DISP.watch(true);
      DISP.clear();
    }
  }

  path >>FILE("z.path");
  return path;
}

void revertPath(arr& path){
  uint N = path.d0;
  arr x;
  for(uint i=0;i<N/2;i++){
    x = path[i];
    path[i] = path[N-1-i];
    path[N-1-i] = x;
  }
}

ptr<PathResult> PathFinder_RRT::run(double timeBudget){
  CHECK(starts.N, "no goal given")
  CHECK(goals.N, "no goal given")
  arr path = planConnect(starts, goals, .5, .0, .0); //.9, .9);

  return make_shared<PathResult>(path);
}

