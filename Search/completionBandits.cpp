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

UILE_Solver::UILE_Solver(const shared_ptr<rai::ComputeNode>& _root) : root(0, _root) {
  gamma = rai::getParameter<double>("CB/gamma", 1.);
  beta = rai::getParameter<double>("CB/beta", 1.);
  epsilon = rai::getParameter<double>("CB/epsilon", .1);
  useUCBasData = rai::getParameter<bool>("CB/useUCBasData", false);
  costCoeff = rai::getParameter<double>("CB/costCoeff", .01);

  root.comp->isComplete=true;
  root.comp->isTerminal=false;
  all.append(&root);
  nonTerminals.append(&root);
}

void UILE_Solver::query(CB_Node *n){
  if(!n) return;

  if(verbose>0) LOG(0) <<"querying #" <<n->ID <<'_' <<n->comp->name;

  //== complete & non-terminal -> create new child, then query it directly
  if(n->comp->isComplete && !n->comp->isTerminal){
    shared_ptr<CB_Node> child = make_shared<CB_Node>(n, n->comp->getNewChild(n->R));
    n->R++;
    n->children.append(child);
    all.append(child.get());
    if(verbose>0) LOG(0) <<"created new child ID:" <<child->ID <<" of type '" <<rai::niceTypeidName(typeid(*child->comp)) <<"'";
    query(child.get());
    return;
  }

  //== incomplete -> compute
  if(!n->comp->isComplete){

    double time = n->comp->timedCompute();
    n->comp->c += time;
    n->comp_n += 1.;
    n->parent->c_children++;

    if(n->comp->isComplete){
      if(n->comp->l<0.) n->comp->l = n->comp->valueHeuristic();
      if(verbose>0) LOG(0) <<"computed #" <<n->ID <<'_' <<n->comp->name <<" -> complete with c=" <<n->comp->c <<" l=" <<n->comp->l;
      CHECK_GE(n->comp->l, 0., "lower bound was not computed");
      if(n->comp->isTerminal){
        terminals.append(n);
      }else{
        if(n->comp->l<1e10){
          nonTerminals.append(n);
        }else{
          n->childrenComplete = true;
          n->branchComplete = true;
          CB_Node* p=n;
          while(p){
            p->y_tot += -1.;
            p->y_num += 1.;
            p=p->parent;
          }
        }
      }
    }else{
      if(verbose>0) LOG(0) <<"computed #" <<n->ID <<'_' <<n->comp->name <<" -> still incomplete with c=" <<n->comp->c;
    }

    //backup compute costs
    CB_Node* p=n->parent;
    while(p){
      p->c_tot += time;
      p=p->parent;
    }

    //backup completeness
    if(n->comp->isComplete){
      CB_Node* p=n->parent;
      bool allComplete=true;
      if(p->comp->getNumDecisions()<0 || (int)p->R<p->comp->getNumDecisions()){
        allComplete=false; //not all possible children expanded
      }else{
        for(auto& ch:p->children) if(!ch->comp->isComplete){ allComplete=false; break; }//not all children closed
      }
      if(allComplete) p->childrenComplete=true;
    }

    //backup closedness
    if(n->comp->isComplete && n->comp->isTerminal) n->branchComplete=true;
    if(n->branchComplete){
      CB_Node* p=n->parent;
      while(p){
        if(p->comp->getNumDecisions()<0 || (int)p->R<p->comp->getNumDecisions()) break; //not all possible children expanded
        for(auto& ch:p->children) if(!ch->branchComplete){ p=0; break; }//not all children closed
        if(!p) break;
        p->branchComplete=true;
        p=p->parent;
      }
    }
  }

  //== complete & terminal -> sample
  if(n->comp->isComplete && n->comp->isTerminal){
    CHECK(n->comp->isTerminal, "");

    double y = n->comp->sample();

    if(y_baseline<0. || y>y_baseline) y_baseline=y;
    y_now=y;
    if(verbose>0) LOG(0) <<"sampled #" <<n->ID <<'_' <<n->comp->name <<" -> return " <<y;

    while(n){
      n->y_tot += y;
      n->y_num += 1.;
      n=n->parent;
    }
  }
}

void UILE_Solver::step(){
  y_now=-1.;
  CB_Node *n = selectNew();
  n->isSelected = true;
  query(n);
  report();
  steps++;
}

void UILE_Solver::runTrivial(uint k, double maxEffortPerCompute){
  CB_Node *n = &root;

  for(uint i=0;i<k;i++){
    query(n);
    report();
    steps++;

    if(n->comp->isComplete){
      if(n->comp->l>=1e10){
        if(verbose>0) LOG(0) <<"compute #" <<n->ID <<'_' <<n->comp->name <<" -> *** infeasible with c=" <<n->comp->c;
        n=&root;
      }else{
        if(n->comp->isTerminal) n=&root;
        else n=select_computeChild(n);
      }
    }else{
      if(n->comp->c > maxEffortPerCompute){
        if(verbose>0) LOG(0) <<"compute #" <<n->ID <<'_' <<n->comp->name <<" -> *** aborted with c=" <<n->comp->c;
        n=&root;
      }
    }

//    if(n==&root){ printTree(cout, root); rai::wait(); }
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
    double y = n->y_tot/n->y_num; //->ucb; OPTION!
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
    if(n->y_num>0.){
//      double parent_num=0;
//      for(auto& ch: n->parent->children) parent_num += ch->y_num;
      n->y_ucb = n->y_tot/n->y_num + beta * ::sqrt(2.*::log(root.y_num) / n->y_num);
    }else{
      n->y_ucb = double(1<<10);
    }
    score(i++) = n->y_ucb;
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
    CB_Node* best=0;
    for(auto& n:n->children){
      if(n->y_num>0.){
        n->y_ucb = n->y_tot/n->y_num + beta * ::sqrt(2.*::log(n->parent->y_num) / n->y_num);
      }else{
        n->y_ucb = -1.;
      }
      if(!best || n->y_ucb>=best->y_ucb) best=n.get();
      if(n->comp->isTerminal) n->score = n->y_ucb - y_baseline;
    }

    // pick the child with highest
    n = best;
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
    if(n->n_children && n->y_num>0.){
#if 0
      double parent_num=0;
      if(n->parent){
        for(auto& ch: n->parent->children) parent_num += ch->c_children;
      }else{
        parent_num = n->c_children;
      }
      n->mean_ucb = n->mean_Y/n->mean_n + beta * ::sqrt(2.*::log(parent_num) / n->c_children);
#else
      n->mean_ucb = n->y_tot/n->y_num + beta * ::sqrt(2.*::log(nTotal) / n->n_children);
#endif
    }else{
      n->mean_ucb = double(1<<10);
    }
    //lowest compute of children
    double LE=get_expandThreshold(n);
    for(auto& ch:n->children) LE += 1. + costCoeff*ch->comp->c;// + ch->effortHeuristic();
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
//    double score = sqrt(n->c_tot)*sqrt(n->R);
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
  if(!n->childrenComplete/*nonTerminals.contains(n)*/) best=n;
  while(!n->comp->isTerminal){
    //for all non-closed children compute score
    CB_Node* m=0;
    double mScore=0.;
    for(auto& ch:n->children){
      if(ch->branchComplete) continue;
      double score = ch->c_tot;
      if(!m || score<mScore){ m=ch.get(); mScore=score; }
    }

    if(!m) break; //all children closed
    if(!m->comp->isComplete) break; //m itself (as compParent) needs to be complete

    if(!m->childrenComplete/*nonTerminals.contains(m)*/){
      if(!best || best->R>m->R) best=m;
    }

    n = m; //iterate down
  }
  return best;
}

CB_Node* UILE_Solver::getCheapestIncompleteChild(CB_Node *r){
  CB_Node* m=0;
  for(auto& ch:r->children) if(!ch->comp->isComplete && (!m || ch->comp->c < m->comp->c)) m=ch.get();
  return m;
}

double UILE_Solver::get_expandThreshold(CB_Node *r){
//  return gamma*::sqrt(r->c_children);
  return gamma*(r->R);
}

void UILE_Solver::clearScores() {
  for(CB_Node *n:all){
    n->score = -1.;
    n->isSelected = false;
    n->isBest = false;
  }
}

CB_Node* UILE_Solver::getBestCompute(){
  CB_Node* best=0;
  for(CB_Node *n:all){
    if(!n->comp->isComplete){
      n->score = 1./( n->comp->c + n->comp->effortHeuristic() );
      if(!best || n->score>=best->score) best=n;
    }
  }
  if(best) best->isBest=true;
  return best;
}

CB_Node* UILE_Solver::getBestExpand(){
  CB_Node* best=0;
  for(CB_Node *n:all){
    int Rmax = n->comp->getNumDecisions();
    if(n->comp->isComplete && !n->comp->isTerminal && n->comp->l<1e9 && (Rmax<0 || (int)n->R<Rmax)){
      n->score = 1./( sqrt(n->y_num+1.) * (n->R+1.) * n->comp->correlationHeuristic() ); // * (n->comp->effortHeuristic());
      if(!best || n->score>=best->score) best=n;
    }
  }
  if(best) best->isBest=true;
  return best;
}


CB_Node* UILE_Solver::select_computeChild(CB_Node *r){
  //-- get best INCOMPLETE child
  CB_Node* j = getCheapestIncompleteChild(r);

  if(r->comp->getNumDecisions()<0 || (int)r->R < r->comp->getNumDecisions()){
    if(!j || j->comp->c>get_expandThreshold(r)){ //not good enough -> create a new child
      return r; //r-decision: create new child
//      r->children.append(make_shared<CB_Node>(r, r->comp->getNewChild(r->R)));
//      j = r->children(-1).get();
//      r->R++;
//      if(verbose>0) LOG(0) <<"created new child ID:" <<j->ID <<" of type '" <<rai::niceTypeidName(typeid(*j->comp)) <<"'";
    }
  }

  if(!j){ //all possible children are complete
    if(verbose>0) LOG(0) <<"all children of r=" <<r->comp->name <<" complete -- removing it from nonTerminals -- selected none";
    nonTerminals.removeValue(r);
    r->childrenComplete=true;
  }else{
    if(verbose>0) LOG(0) <<"chose j=" <<j->comp->name <<" as cheepest compute child of r=" <<r->comp->name;
  }

  return j;
}

CB_Node* UILE_Solver::select(){

  //-- baseline
#if 0 //recompute y_max as mean...
  y_max = 0.;
  for(CB_Node* n:terminals) if(n->y_num){
    double mean = n->y_tot/n->y_num;
    if(mean>y_max) y_max = mean;
  }
#else
  y_baseline = 1.;
#endif

  //-- SELECT TERMINAL
//  CB_Node *n = select_terminalFlat();
  CB_Node *n = select_terminalUCT();

  //select and threshold
  if(n && (n->y_ucb - y_baseline)>root.thetaImprovement) return n;

  //-- SELECT NON-TERMINAL
  //backup all these scores down the tree (as in MCTS)
//  CB_Node *r = select_compParent_IE();
//  CB_Node *r = select_compParent_CR();
  CB_Node *r = select_compParent_Tree();
  if(!r) return n;

  //-- PICK CHILD OF NON-TERMINAL (OR CREATE NOVEL)
  return select_computeChild(r);
}

CB_Node* UILE_Solver::selectNew(){
  clearScores();

  //sample?
  CB_Node *s = select_terminalUCT();
  //select and threshold
  if(s && s->score>root.thetaImprovement) return s;

  //compute?
  CB_Node *c = getBestCompute();
  if(c && c->comp->c < gamma*sqrt(root.c_tot)) return c;

  //expand?
  CB_Node *e = getBestExpand();
  return e;
}

void UILE_Solver::report(){
  double groundTruthMean=0.;
  double groundTruthBest=0.;
  for(CB_Node* n:terminals){
//    double y = n->groundTruthMean();
    double y = n->y_tot/n->y_num;
    groundTruthMean += y;
    if(y>groundTruthBest) groundTruthBest=y;
  }
  if(terminals.N) groundTruthMean/= double(terminals.N);
  if(fil) (*fil) <<steps <<' ' <<totalCost() <<' ' <<y_now <<' ' <<(root.y_num>0.?root.y_tot/root.y_num:0.) <<' ' <<y_baseline <<' ' <<groundTruthMean <<' ' <<groundTruthBest <<' ' <<terminals.N <<' ' <<nonTerminals.N <<endl;
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

    sub.newNode<double>("score", {}, n->score);
    sub.newNode<double>("c", {}, n->comp->c);
    sub.newNode<double>("comp_n", {}, n->comp_n);
//    sub.newNode<double>("c_children", {}, n->c_children);
    if(n->comp->l>=0.){
      sub.newNode<double>("l", {}, n->comp->l);
    }
//    sub.newNode<double>("R", {}, n->R);
    if(n->comp->l<1e9){
      if(n->y_num){
        sub.newNode<double>("y_mean", {}, n->y_tot/n->y_num);
        sub.newNode<double>("y_num", {}, n->y_num);
      }
        //    sub.newNode<double>("y_ucb", {}, n->y_ucb);
        //    sub.newNode<double>("C_ucb", {}, n->mean_ucb);
        //    sub.newNode<double>("C_eff", {}, n->eff);
      if(n->c_tot){
        sub.newNode<double>("c_tot", {}, n->c_tot);
      }
      if(n->childrenComplete) sub.newNode<bool>("childrenCpl");
      if(n->branchComplete) sub.newNode<bool>("branchCpl");
    }

//    if(n->D.N) sub.newNode<arr>("data", {}, n->D);
    if(!n->comp->isComplete) G.getRenderingInfo(sub.isNodeOfGraph).dotstyle <<", shape=box, style=dashed";
    else if(n->comp->isTerminal) G.getRenderingInfo(sub.isNodeOfGraph).dotstyle <<", shape=box, style=rounded";
    if(n->isSelected) G.getRenderingInfo(sub.isNodeOfGraph).dotstyle <<", color=red";
    else if(n->isBest) G.getRenderingInfo(sub.isNodeOfGraph).dotstyle <<", color=orange";
  }

  G.checkConsistency();
  G.write(FILE("z.tree"));
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
