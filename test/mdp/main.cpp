/*
    Copyright (C) 2007, 2008  Marc Toussaint.
    email: mtoussai@cs.tu-berlin.de

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    A copy of the GNU Lesser General Public License can usually be found
    at http://www.gnu.org/copyleft/lesser.html; if not, write to the Free
    Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
    02110-1301 USA
*/

#include <MDP/mdp.h>
#include <Core/util.h>
#include <iomanip>

enum { ValueIteration, PrioritizedSweeping, PolicyIteration, EM, POMDP_EM, POMDP_EM_hierarchical };

struct MDP_Parameters{
  RAI_PARAMt("", int, mode, double, 0)
  RAI_PARAM("", rai::String, problem, "15x20holes.ppm")
  RAI_PARAM("", double, gamma, .999)
  RAI_PARAMt("", uint, method, double, 3)
  RAI_PARAMt("", uint, d0, double, 0)
  RAI_PARAMt("", uint, d1, double, 0)
  RAI_PARAMt("", uint, seed, double, 1)
  RAI_PARAMt("", uint, iterations, double, 200)

  RAI_PARAMt("", uint, estepHorizon, double, 100)
  RAI_PARAM("", bool, estepStructured, true)
  RAI_PARAM("", bool, estepIncremental, false)
  RAI_PARAMt("", mdp::MstepType, mstepType, double, 0)
  RAI_PARAM("", double, mstepRate, .3)
  RAI_PARAM("", double, mstepNoise, 1e-5)

  RAI_PARAM("", bool, forceLevel1, false)
};



//rai::getParameter(PARAMS.mode,"mode",0);
//rai::getParameter(PARAMS.seed,"seed",(uint)0);
//rai::getParameter(PARAMS.method,"method",(uint)3);
//rai::getParameter(PARAMS.d0,"d0",(uint)0);
//rai::getParameter(PARAMS.d1,"d1",(uint)0);
//rai::getParameter(PARAMS.problem,"problem",rai::String("15x20holes.ppm"));
//rai::getParameter(PARAMS.iterations,"iterations",(uint)200);
//rai::getParameter(PARAMS.estepHorizon,"estepHorizon",(uint)100);
//rai::getParameter(PARAMS.estepStructured,"estepStructured",true);
//rai::getParameter(PARAMS.estepIncremental,"estepIncremental",false);
//rai::getParameter((int&)PARAMS.mstepType,"mstepType",0);
//rai::getParameter(PARAMS.mstepRate,"mstepRate",.3);
//rai::getParameter(PARAMS.mstepNoise,"mstepNoise",1e-5);
//rai::getParameter(PARAMS.forceLevel1,"forceLevel1",false);
//rai::getParameter(PARAMS.gamma,"gamma",.999);


void loadProblem(mdp::MDP& mdp, MDP_Parameters& PARAMS){
  //load the problem
  uint n=strlen(PARAMS.problem);
  if(!strcmp(PARAMS.problem,"simpleMaze")){
    cout <<"loading standard problem: simpleMaze" <<endl;
    mdp::generateStandardProblem(mdp,mdp::miniMaze); //tinyMaze //heavenAndHell
    mdp.gamma=PARAMS.gamma;
    mdp::createNeighorList(mdp);
  }else if(!strcmp(PARAMS.problem+n-3,"ppm")){
    cout <<"loading ppm image maze file " <<PARAMS.problem <<endl;
    mdp::readImageMaze(mdp, PARAMS.problem);
    //tunnelRaxToPx(mdp.Pxax,mdp.Rax,mdp.Px);
    mdp::tunnelRaxTo(mdp.Pxax,mdp.Rax,-1);
    mdp::addActionNoise(mdp.Pxax,rai::getParameter<double>("mazeNoise"));
    mdp::showMaze();
    mdp.gamma=PARAMS.gamma;
    mdp::createNeighorList(mdp);
  }else if(!strcmp(PARAMS.problem+n-9,"POMDP.arr")){
    cout <<"reading POMDP file: " <<PARAMS.problem <<endl;
    mdp::readMDP_arr(mdp, PARAMS.problem);
  }else HALT("don't know how to load problem " <<PARAMS.problem);
}

void run_pomdpEM_lev2(const mdp::MDP& mdp, MDP_Parameters& PARAMS, bool hierarchy=false){
  //init policy
  uint d0=PARAMS.d0,d1=PARAMS.d1,da=mdp.Pxax.d1,dy=mdp.Pyxa.d0;
  mdp::FSC_lev2 fsc;
  fsc.P01y0.resize(TUP(d0,d1,dy,d0));
  fsc.P1y01.resize(TUP(d1,dy,d0,d1));
  mdp::zeroNodeStart(fsc.P0,d0);
  mdp::zeroNodeStart(fsc.P1,d1);
  mdp::oneNodeOneAction(fsc.Pa0,da,d0,1.,1.,100.);
  mdp::generalNode0Transition(fsc.P01y0,d0,d1,dy,.1,.1,0.);
  mdp::generalNode1Transition(fsc.P1y01,d1,dy,d0,.1,.1,1.);

  if(hierarchy){
    initHierarchical(TUP(da,d1,d0/2).min(),fsc);
  }else{
    fsc.hierarchical=false;
  }

  //prepare output file
  rai::String filename;
  filename <<"data/OUT";
  if(!hierarchy) filename<<"2"; else filename<<"H";
  filename<<"n";
  filename <<"-"<<PARAMS.d0<<"-"<<PARAMS.d1<<"-"<<PARAMS.problem<<"-"<<PARAMS.estepHorizon<<"-"<<PARAMS.seed;
  cout <<"output filename = " <<filename <<endl;
  ofstream out(filename);
  
  //iterate EM
  double Like=0.;
  double tic=rai::cpuTime();
  if(PARAMS.forceLevel1 || d1==1){ //use specialized flat controller optimization
    mdp::FSC_lev1 fsc1;
    collapse2levelFSC(fsc1,fsc);
    for(uint k=0;k<PARAMS.iterations;k++){
      cout <<k <<' ';
      Like=pomdpEM_lev1(mdp,fsc1,
			PARAMS.estepHorizon,PARAMS.estepStructured,PARAMS.estepIncremental,
			PARAMS.mstepType,PARAMS.mstepRate,PARAMS.mstepNoise,
			false,
			NULL,NULL,NULL);
      out <<k <<' ' <<rai::timerRead(false,tic) <<' ' <<Like <<endl;
    }
  }else{
    for(uint k=0;k<PARAMS.iterations;k++){
      cout <<k <<' ';
      Like=pomdpEM_lev2(mdp,fsc,PARAMS.estepHorizon,PARAMS.estepStructured,PARAMS.mstepType,false);
      out <<k <<' ' <<rai::timerRead(false,tic) <<' ' <<Like <<endl;
    }
  }
  cout <<"\n\ntotal time = " <<rai::timerRead(false,tic) <<endl;
}

void testSolvers(uint method, const mdp::MDP& mdp, MDP_Parameters& PARAMS){
  arr V,pi,Pvisited;
  uint i,j;
  
  switch(method){
  case ValueIteration:
    cout <<"*** value iteration" <<endl;
    for(i=0;i<PARAMS.iterations;i++){
      mdp.valueIteration(V);
      mdp.maxPolicy(pi,V);
      cout <<"\r #" <<std::setw(3) <<i <<" V(x_0)=" <<scalarProduct(mdp.Px,V) <<std::flush;
      if(i<50 || !(i%10)) plotPolicyAndValue(pi,V,mdp,i<100);
    }
    cout <<endl;
    plotPolicyAndValue(pi,V,mdp,true);
    break;
  case PrioritizedSweeping:
    cout <<"*** prioritized sweeping" <<endl;
    V.setZero();
    mdp.prioritizedSweeping(V);
    mdp.maxPolicy(pi,V);
    plotPolicyAndValue(pi,V,mdp,true);
    break;
  case PolicyIteration:
    cout <<"*** policy iteration" <<endl;
    randomPolicy(pi,mdp);
    V.setZero();
    for(i=0;i<PARAMS.iterations/10;i++){
      for(j=0;j<10;j++) mdp.policyEvaluation(V,pi);
      mdp.maxPolicy(pi,V);
      cout <<"\r #" <<std::setw(3) <<i <<" V(x_0)=" <<scalarProduct(mdp.Px,V) <<std::flush;
      if(i<10) plotPolicyAndValue(pi,V,mdp,i<100);
    }
    cout <<endl;
    plotPolicyAndValue(pi,V,mdp,true);
    break;
  case EM:
    cout <<"*** EM policy optimization" <<endl;
    randomPolicy(pi,mdp);
    V.setZero();
    for(i=0;i<PARAMS.iterations;i++){
      mdpEM(mdp,pi,V,PARAMS.estepHorizon,5,mdp::MstepNone,&Pvisited);
      mdpEM(mdp,pi,V,PARAMS.estepHorizon,5,PARAMS.mstepType,NULL); //&Pvisited);
      //maxPolicy(pi,V,mdp);
      calcPvisited(Pvisited,mdp);
      cout <<"\r #" <<std::setw(3) <<i <<" V(x_0)=" <<scalarProduct(mdp.Px,V) <<std::flush;
      if(i<100){
        plotPolicyAndValue(pi,Pvisited,mdp,i<100);
        plotPolicyAndValue(pi,V,mdp,i<100);
      }
    }
    cout <<endl;
    plotPolicyAndValue(pi,V,mdp,true);
    break;
  case POMDP_EM:
    CHECK(PARAMS.d0 && PARAMS.d1,"need to specify parameters d0 and d1 (number of nodes on base and abstract level) greater than zero");
    run_pomdpEM_lev2(mdp, PARAMS, false);
    break;
  case POMDP_EM_hierarchical:
    CHECK(PARAMS.d0 && PARAMS.d1,"need to specify parameters d0 and d1 (number of nodes on base and abstract level) greater than zero");
    run_pomdpEM_lev2(mdp, PARAMS, true);
    break;
  default: HALT("don't know this method: " <<method);
  }
}

int main(int argc,char** argv){
  rai::initCmdLine(argc,argv);

  if(!rai::checkParameter<double>("mode")){
    char buf[256];
    ifstream is;
    rai::open(is,"README");
    while(is.good()){ is.getline(buf,256); cout <<buf <<endl; }
    return 1;
  }

  MDP_Parameters PARAMS;
  rnd.seed(PARAMS.seed);

  cout <<std::setprecision(5);
  
  mdp::MDP mdp;

  switch(PARAMS.mode){
  case 1: mdp::convertTonyToMDP(PARAMS.problem); break;
  case 2: scanArrFile(PARAMS.problem); break;
  case 3:
    loadProblem(mdp,PARAMS);
    writeMDP_arr(mdp,"z.mdp");
    testSolvers(PARAMS.method, mdp, PARAMS);
    break;
  default:
    NIY;
  }
  
  return 0;
}


