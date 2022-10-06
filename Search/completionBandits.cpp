#include "completionBandits.h"

#include <math.h>
#include <Core/util.h>
#include <Core/graph.h>

static uint CB_Node_ID=0;

CB_Node::CB_Node(CB_Node* _parent, shared_ptr<rai::ComputeNode> _comp)
  : ID(CB_Node_ID++), parent(_parent), comp(_comp) {
  thetaImprovement = rai::getParameter<double>("CB/thetaImprovement", .1);
}

void CB_Node::write(std::ostream& os) const {
  os <<'#' <<ID <<'_' <<comp->name <<" n=" <<c_children <<" R=" <<R;
//  if(comp->isTerminal) os <<" data: " <<D;
}

//===========================================================================

UILE_Solver::UILE_Solver(CB_Node& _root) : root(_root) {
  gamma = rai::getParameter<double>("CB/gamma", 1.);
  beta = rai::getParameter<double>("CB/beta", 1.);
  epsilon = rai::getParameter<double>("CB/epsilon", .1);
  useUCBasData = rai::getParameter<bool>("CB/useUCBasData", false);
  costCoeff = rai::getParameter<double>("CB/costCoeff", .01);

  root.comp->isComplete=true;
  root.comp->isTerminal=false;
  nonTerminals.append(&root);
}

void UILE_Solver::step(){
  y_now=-1.;
  query(select());
  report();
  steps++;
}

void UILE_Solver::query(CB_Node *n){
  if(!n) return;

  if(verbose>0) LOG(0) <<"querying " <<n->ID <<'_' <<n->comp->name;

  if(!n->comp->isComplete){ //incomplete

    double c = n->comp->compute();
//    n->comp->isComplete = n->isComplete();
    n->c += c;
    c_total += c;
    n->parent->c_children++;
//    n->comp_n += c;
//    n->parent->comp_n += c;

    if(n->comp->isComplete){
      if(n->l<0.) n->l = n->comp->costHeuristic();
      if(verbose>0) LOG(0) <<"computed " <<n->ID <<'_' <<n->comp->name <<" -> complete with c=" <<n->c <<" l=" <<n->l;
      CHECK_GE(n->l, 0., "lower bound was not computed");
      if(n->comp->isTerminal) terminals.append(n);
      else nonTerminals.append(n);
    }else{
      if(verbose>0) LOG(0) <<"computed " <<n->ID <<'_' <<n->comp->name <<" -> still incomplete with c=" <<n->c;
    }

    //backup compute costs
    CB_Node* p=n->parent;
    while(p){
      p->comp_C += c;
      if(n->comp->isComplete && n->comp->isTerminal) p->comp_n += 1.;
      p=p->parent;
    }

    //backup closedness
    if(n->comp->isComplete && n->comp->isTerminal){
      n->isClosed=true;
      CB_Node* p=n->parent;
      while(p){
        if(p->comp->getNumDecisions()<0 || (int)p->R<p->comp->getNumDecisions()) break; //not all possible children expanded
        for(auto& ch:p->children) if(!ch->isClosed){ p=0; break; }//not all children closed
        if(!p) break;
        p->isClosed=true;
        p=p->parent;
      }
    }

//  }else{ //complete
  }

//  n->isTerminal = n->isTerminal();
  if(n->comp->isComplete && n->comp->isTerminal){
    CHECK(n->comp->isTerminal, "");

    double y = n->comp->sample();

    if(y_baseline<0. || y>y_baseline) y_baseline=y;
    y_now=y;
    if(verbose>0) LOG(0) <<"sampled " <<n->ID <<'_' <<n->comp->name <<" -> return " <<y;

    while(n){
      n->data_Y += y;
      n->data_n += 1.;
      n=n->parent;
    }
  }
}

void UILE_Solver::runTrivial(uint k, double maxEffortPerCompute){
  CB_Node *n = &root;

  for(uint i=0;i<k;i++){
    query(n);
    report();
    steps++;

    if(n->comp->isComplete){
      if(n->comp->isTerminal) n=&root;
      else n=select_computeChild(n);
    }else{
      if(n->c > maxEffortPerCompute){
        if(verbose>0) LOG(0) <<"compute " <<n->ID <<'_' <<n->comp->name <<" -> *** aborted with c=" <<n->c;
        n=&root;
      }
    }

    printTree(cout, root); rai::wait();
  }
}

void UILE_Solver::backupMeans(){
  //clear all data
  for(CB_Node* n:nonTerminals){
    n->mean_Y=0.; n->mean_n=0.;
    //n->mean_Y=y_baseline; n->mean_n=1.;
  }
  //standard backup from all terminals
  for(CB_Node* n:terminals){
    double y = n->data_Y/n->data_n; //->ucb; OPTION!
//    if(useUCBasData) y = n->ucb;
    n=n->parent;
    while(n){
      n->mean_Y += y;
      n->mean_n += 1.;
      n=n->parent;
    }
  }
}

CB_Node* UILE_Solver::select_terminalFlat(){
  if(!terminals.N) return 0;

  //-- compute UCB1 score for all terminals
  arr score(terminals.N);
  uint i=0;
  for(CB_Node* n:terminals){
    if(n->data_n>0.){
//      double parent_num=0;
//      for(auto& ch: n->parent->children) parent_num += ch->data_n;
      n->data_ucb = n->data_Y/n->data_n + beta * ::sqrt(2.*::log(root.data_n) / n->data_n);
    }else{
      n->data_ucb = double(1<<10);
    }
    score(i++) = n->data_ucb;
  }
  if(verbose>0) LOG(0) <<"terminal's data scores: " <<score;

  return terminals(argmax(score));
}

CB_Node* UILE_Solver::select_terminalUCT(){
  if(!terminals.N) return 0;

  //-- select terminal node using tree policy
  CB_Node* n = &root;
  while(!n->comp->isTerminal){
    //for all children with data compute UCB1 score
    arr alpha(n->children.N);
    uint i=0;
    for(auto& n:n->children){
      if(n->data_n>0.){
        n->data_ucb = n->data_Y/n->data_n + beta * ::sqrt(2.*::log(n->parent->data_n) / n->data_n);
      }else{
        n->data_ucb = -1.;
      }
      alpha(i++) = n->data_ucb;
    }

    // pick the child with highest
    n = n->children(argmax(alpha)).get();
  }
  return n;
}

CB_Node* UILE_Solver::select_compParent_IE(){
  if(!nonTerminals.N){
    if(verbose>0) LOG(0) <<"no nonTerminals -- select none";
    return 0;
  }

  //-- count 'samples'(=children) per nonTerminal, and total
  uint nTotal=0;
  for(CB_Node* n:nonTerminals){
    uint ncompl=0;
    for(auto& ch: n->children) if(ch->comp->isComplete) ncompl++;
    n->n_children = ncompl;
    nTotal += ncompl;
  }

  //-- compute UCB1 score for all non-terminals
  arr alpha(nonTerminals.N);
  uint i=0;
  //pick non-terminal node with highest UCB1 (could be in the middle!!), or with highest UI/LE
  for(CB_Node* n:nonTerminals){
    if(n->n_children && n->data_n>0.){
#if 0
      double parent_num=0;
      if(n->parent){
        for(auto& ch: n->parent->children) parent_num += ch->c_children;
      }else{
        parent_num = n->c_children;
      }
      n->mean_ucb = n->mean_Y/n->mean_n + beta * ::sqrt(2.*::log(parent_num) / n->c_children);
#else
      n->mean_ucb = n->data_Y/n->data_n + beta * ::sqrt(2.*::log(nTotal) / n->n_children);
#endif
    }else{
      n->mean_ucb = double(1<<10);
    }
    //lowest compute of children
    double LE=get_novelThreshold(n);
    for(auto& ch:n->children) LE += 1. + costCoeff*ch->c;// + ch->effortHeuristic();
    LE /= double(1+n->children.N);
    n->eff = LE;
    alpha(i++) = (n->mean_ucb - y_baseline) / LE;
  }

  if(verbose>0) LOG(0) <<"nonTerminal's mean alphas: " <<alpha;

  return nonTerminals(argmax(alpha));
}

CB_Node* UILE_Solver::select_compParent_CR(){
  if(!nonTerminals.N){
    if(verbose>0) LOG(0) <<"no nonTerminals -- select none";
    return 0;
  }

  CB_Node* m=0;
  double mScore=0.;
  for(CB_Node* n:nonTerminals){
//    double score = sqrt(n->comp_C)*sqrt(n->R);
    double score = n->R;
    if(!m || score<mScore){ m=n; mScore=score; }
  }

  return m;
}

CB_Node* UILE_Solver::select_compParent_Tree(){
  if(!nonTerminals.N){
    if(verbose>0) LOG(0) <<"no nonTerminals -- select none";
    return 0;
  }

  //-- select terminal node using tree policy
  CB_Node* n = &root;
  CB_Node* best = 0;
  if(nonTerminals.contains(n)) best=n;
  while(!n->comp->isTerminal){
    //for all children with data compute UCB1 score
    CB_Node* m=0;
    double mScore=0.;
    for(auto& ch:n->children){
      if(ch->isClosed) continue;
      double score = ch->comp_C;
      if(!m || score<mScore){ m=ch.get(); mScore=score; }
    }

    if(!m) break; //all children closed
    if(!m->comp->isComplete) break; //m itself (as compParent) needs to be complete

    if(nonTerminals.contains(m)){
      if(!best || best->R>m->R) best=m;
    }

    n = m; //iterate down
  }
  return best;
}

CB_Node* UILE_Solver::getCheapestIncompleteChild(CB_Node *r){
  CB_Node* m=0;
  for(auto& ch:r->children) if(!ch->comp->isComplete && (!m || ch->c < m->c)) m=ch.get();
  return m;
}

double UILE_Solver::get_novelThreshold(CB_Node *r){
//  return gamma*::sqrt(r->c_children);
  return gamma*(r->R);
}


CB_Node* UILE_Solver::select_computeChild(CB_Node *r){
  //-- get best INCOMPLETE child
  CB_Node* j = getCheapestIncompleteChild(r);

  if(r->comp->getNumDecisions()<0 || (int)r->R < r->comp->getNumDecisions()){
    if(!j || j->c>get_novelThreshold(r)){ //not good enough -> create a new child
      r->children.append(make_shared<CB_Node>(r, r->comp->getNewChild(r->R)));
      j = r->children(-1).get();
      r->R++;
      if(verbose>0) LOG(0) <<"created new child ID:" <<j->ID <<" of type '" <<rai::niceTypeidName(typeid(*j->comp)) <<"'";
    }
  }

  if(!j){ //all possible children are complete
    if(verbose>0) LOG(0) <<"all children of r=" <<r->comp->name <<" complete -- removing it from nonTerminals -- selected none";
    nonTerminals.removeValue(r);
  }else{
    if(verbose>0) LOG(0) <<"chose j=" <<j->comp->name <<" as cheepest compute child of r=" <<r->comp->name;
  }

  return j;
}

CB_Node* UILE_Solver::select(){

  //-- baseline
#if 0 //recompute y_max as mean...
  y_max = 0.;
  for(CB_Node* n:terminals) if(n->data_n){
    double mean = n->data_Y/n->data_n;
    if(mean>y_max) y_max = mean;
  }
#else
  y_baseline = 1.;
#endif

  //-- SELECT TERMINAL
//  CB_Node *n = select_terminalFlat();
  CB_Node *n = select_terminalUCT();

  //select and threshold
  if(n && (n->data_ucb - y_baseline)>root.thetaImprovement) return n;

  //-- SELECT NON-TERMINAL
  //backup all these scores down the tree (as in MCTS)
  CB_Node *r = select_compParent_IE();
//  CB_Node *r = select_compParent_CR();
//  CB_Node *r = select_compParent_Tree();
  if(!r) return n;

  //-- PICK CHILD OF NON-TERMINAL (OR CREATE NOVEL)
  return select_computeChild(r);
}

void UILE_Solver::report(){
  double groundTruthMean=0.;
  double groundTruthBest=0.;
  for(CB_Node* n:terminals){
//    double y = n->groundTruthMean();
    double y = n->data_Y/n->data_n;
    groundTruthMean += y;
    if(y>groundTruthBest) groundTruthBest=y;
  }
  if(terminals.N) groundTruthMean/= double(terminals.N);
  if(fil) (*fil) <<steps <<' ' <<totalCost() <<' ' <<y_now <<' ' <<(root.data_n>0.?root.data_Y/root.data_n:0.) <<' ' <<y_baseline <<' ' <<groundTruthMean <<' ' <<groundTruthBest <<' ' <<terminals.N <<' ' <<nonTerminals.N <<endl;
}


rai::Array<CB_Node*> getAllNodes(CB_Node& root){
  rai::Array<CB_Node*> queue;
  queue.append(&root);

  uint i=0;
  while(i<queue.N){
    CB_Node *n = queue(i);
    n->ID=i;
    i++;

    if(n->parent) CHECK_EQ(queue(n->parent->ID), n->parent, "");

    for(uint j=0;j<n->children.N;j++){
      queue.append(n->children(j).get());
    }
  }
  return queue;
}

void printTree(std::ostream& os, CB_Node& root){

  rai::Array<CB_Node*> T = getAllNodes(root);

  rai::Graph G;
  for(uint i=0;i<T.N;i++){
    CB_Node *n = T(i);
    rai::NodeL par;
    if(n->parent) par.append(G.elem(n->parent->ID));
    rai::Graph& sub = G.newSubgraph(n->comp->name, par, {});

    sub.newNode<double>("c", {}, n->c);
//    sub.newNode<double>("c_children", {}, n->c_children);
    sub.newNode<double>("l", {}, n->l);
//    sub.newNode<double>("R", {}, n->R);
    sub.newNode<double>("D_n", {}, n->data_n);
    sub.newNode<double>("D_mean", {}, n->data_Y/n->data_n);
    sub.newNode<double>("D_ucb", {}, n->data_ucb);
    sub.newNode<double>("C_ucb", {}, n->mean_ucb);
    sub.newNode<double>("C_eff", {}, n->eff);
    sub.newNode<double>("C_tot", {}, n->comp_C);
    sub.newNode<double>("C_n", {}, n->comp_n);
    sub.newNode<bool>("closed", {}, n->isClosed);

//    if(n->D.N) sub.newNode<arr>("data", {}, n->D);
    if(n->comp->isTerminal){
//      sub.newNode<rai::String>("dotstyle", {}, ", color=red");
      G.getRenderingInfo(sub.isNodeOfGraph).dotstyle <<", shape=box";
    }
  }

  G.checkConsistency();
  cout <<G <<endl;
  G.writeDot(FILE("z.dot"));
  rai::system("dot -Tpdf z.dot > z.pdf");
}

double ExpectedImprovement(const arr& y){
  uint n=y.N;
  arr ysort = y;
  for(uint i=0;i<y.N;i++) if(ysort(i)<0.){ ysort(i) = 100.; n--; } //discard nil returns
  std::sort(ysort.p, ysort.p+y.N);

  if(n<3){
    if(n<2) return 1.;
    return ysort(1) - ysort(0);
  }
  double EI=0.;
  EI += (1./n)*(1./(n-1)) * (ysort(2) - ysort(0));
  EI += (1./n)*(double(n-2)/(n-1)) * (ysort(1) - ysort(0));
  EI += (1./n)*(1./(n-1)) * (ysort(2) - ysort(1));

  if(EI<0.) HALT(ysort);
  //  CHECK_GE(EI, 0., "");

  return EI;
}


uint getDepth(CB_Node* n){
  int i=0;
  while(n->parent){ n=n->parent; i++; }
  return i;
}
