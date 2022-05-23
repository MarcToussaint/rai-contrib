/*  ---------------------------------------------------------------------
    Copyright 2012 Marc Toussaint
    email: mtoussai@cs.tu-berlin.de
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a COPYING file of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>
    -----------------------------------------------------------------  */

#include "infer.h"
#include <algorithm>
#include <iomanip>


#define LOG_NORM_SCALE
int DEBUG_INFER_LEVEL=0;
#define DEBUG_INFER(l, x) if(l<=DEBUG_INFER_LEVEL){ x; }
uint infer::VarCount=0;

#define FOR1D(x,i) for(i=0;i<x.N;i++)
#define FOR2D(x,i,j) for(i=0;i<x.d0;i++) for(j=0;j<x.d1;j++)

//dummy to encorde pre-main initialization
struct infer_Init {
  infer_Init(){ infer::VariableList::memMove=1; infer::FactorList::memMove=1; infer::MessagePairList::memMove=1; }
} infer_init;


//===========================================================================
//
// Declarations of various minor functions
//


//void get_vars(infer::VariableList& V, const infer::FactorList& factors);
void getNeighbors(infer::Factor* f, infer::FactorList& neighbors);
void convert_id2config(uintA& config, uint id, const uintA& varDimensions);
uint convert_config2id(uintA& config, const uintA& varDimensions);
uint convert_config2id(boolA truths);
int get_list_index(infer::FactorGraph& fg, infer::Factor* f);
uint get_list_index_unsigned(infer::FactorGraph& fg, infer::Factor* f);


void get_vars(infer::VariableList& V, const infer::FactorList& factors){
  V.clear();
  uint i, k;
  FOR1D(factors, i){
    FOR1D(factors(i)->variables, k){
      V.setAppend(factors(i)->variables(k));
    }
  }
}

void getNeighbors(infer::Factor* f, infer::FactorList& neighbors){
  neighbors.clear();
  uint i;
  FOR1D(f->messages, i){
    if(f->messages(i)->f1 == f)
      neighbors.append(f->messages(i)->f2);
    else
      neighbors.append(f->messages(i)->f1);
  }
}

infer::MessagePair* getMessagePair(const infer::Factor* f1, const infer::Factor* f2){
  uint i;
  FOR1D(f1->messages, i){
    if(f1->messages(i)->f1 == f2  ||  f1->messages(i)->f2 == f2)
      return f1->messages(i);
  }
  return NULL;
}

infer::Factor* getMessage(const infer::Factor* f_from, const infer::Factor* f_to){
  infer::MessagePair* msg_pair = getMessagePair(f_from, f_to);
  CHECK(msg_pair != NULL,  "Message pair not found");
  if(msg_pair->f1 == f_from)
    return &msg_pair->m12;
  else
    return &msg_pair->m21;
}

void writeMessage(const infer::Factor* f_from, const infer::Factor* f_to, ostream& os){
  os <<*getMessage(f_from, f_to);
}


int get_list_index__orig(infer::FactorGraph& fg, infer::Factor* f_orig){
  return fg.F.findValue(f_orig);
}

uint get_list_index_unsigned__orig(infer::FactorGraph& fg, infer::Factor* f_orig){
  uint i=0;
  while(i < fg.F.N){
    if(fg.F(i) == f_orig)
      break;
    i++;
  }
  CHECK(i<fg.F.N, "infer::Factor over " <<f_orig->varIds <<" not found in infer::FactorGraph!" <<endl);
  return i;
}





/* maybe useful...
void checkSpaceRegistry(iSpace *S){
  if(!S) S=&globalInferSpace;
  infer::Variable *v;
  infer::Factor *f;
  uint i, j, k;
  for(i=0;i<S->variables.N;i++){
    v=S->variables(i);
    CHECK_EQ(v->space,S, "variale points to another space");
    CHECK_EQ(v->id,i, "variable id is not the position globalVariableList registry");
    for(j=0;j<v->factors.N;j++){
      f=v->factors(j);
      k=f->varIds.findValue(i);
      CHECK(k<f->varIds.N, "variable's factor doesn't include the variable's id");
      CHECK_EQ(v->dim,f->dim(k), "variable's factor thinks variable has different dim");
    }
  }
  for(i=0;i<S->factors.N;i++){
    f=S->factors(i);
    CHECK_EQ(f->space,S, "factor points to another space");
    for(j=0;j<f->varIds.N;j++){
      v=S->variables(f->varIds(j));
      k=v->factors.findValue(f);
      //CHECK(k<v->factors.N, "factor's variable doesn't include the factor in its list"); [[REDO!]]
      CHECK_EQ(f->dim(j),v->dim, "factor's variable thinks it is of different dim");
    }
  }
}

void iSpace::deleteAll(){
  infer::Factor *f;
  while(factors.N){
    f=factors(0);
    while(f->messages.N) delete f->messages(0);
    delete f;
  }
  while(variables.N) delete variables(0);
}

void iSpace::write(std::ostream& os) const{
  os <<"<iSpace>"
    <<" #variable="   <<variables.N <<endl;
  // <<" #msg_pairs=" <<g.messages.N <<endl;
  uint i, j;
  for(i=0;i<variables.N;i++){
    os <<'[' <<i <<"] " <<variables(i)->name <<' ' <<variables(i)->dim;
    for(j=0;j<variables(i)->factors.N;j++) os <<',' <<variables(i)->factors(j)->varIds;
    os <<endl;
    CHECK_EQ(variables(i)->id,i, "");
  }
  os <<" #factors="    <<factors.N <<endl;
  for(i=0;i<factors.N;i++){
    os <<'[' <<i <<"] " <<*factors(i);
    //for(j=0;j<factors(i)->varIds.N;j++) os <<',' <<variables(factors(i)->varIds(j))->name;
    os <<endl;
  }
  os <<endl;
  / *
  for(i=0;i<g.messages.N;i++){
  os <<'[' <<g.messages(i).from <<'-' <<g.messages(i).to <<"] " <<g.messages.elem(i) <<endl;
  }
  * /
  os <<"\n</iSpace>" <<endl;
}


void iSpace::writeVariables(std::ostream& os) const {
  uint i;
  for(i=0; i<variables.N; i++){
    os <<"[" <<variables(i)->id <<"]  " <<variables(i)->name <<endl;
  }
}

  */




//===========================================================================
//
// infer::Variable
//

infer::Variable::Variable(){
  factors.memMove=true;
  messages.memMove=true;
  id=-1;
}

infer::Variable::Variable(uint _dim, const char *_name){
  id=VarCount++;
  dim=_dim;
  name=_name;
}

infer::Variable::Variable(uint _dim, const char *_name, uint _id){
  id=_id;
  dim=_dim;
  name=_name;
}

infer::Variable::~Variable(){
  if(factors.N) RAI_MSG("you shouldn't destroy variables that are still linked to factors");
}

void infer::Variable::write(ostream& os) const {
  os <<"variable " <<name <<" { dim=" <<dim;
  if(ats.N){ os <<", "; listWrite(ats, os, ","); }
  os <<" }";
}







//===========================================================================
//
// infer::Factor
//

void disconnectFactor(infer::Factor &f){
  for(infer::Variable *v:  f.variables) v->factors.removeValue(&f);
  f.variables.clear();
  f.varIds.clear();
  f.dim.clear();
}

void infer::Factor::relinkTo(const infer::VariableList& vars){
  for(Variable *v:  variables) v->factors.removeValue(this);
  variables=vars;
  for(Variable *v: variables) v->factors.append(this);
  for(uint i=0; i<vars.N; i++) varIds(i)=variables(i)->id;
  for(uint i=0; i<vars.N; i++) CHECK_EQ(dim(i),variables(i)->dim, "relinking to variables with different dimension!");
}

void infer::Factor::init(const infer::VariableList& vars){
  disconnectFactor(*this);
  variables=vars;
  varIds.resize(vars.N);
  dim.resize(vars.N);
  for(uint i=0; i<vars.N; i++) varIds(i)=vars(i)->id;
  for(uint i=0; i<vars.N; i++) dim(i)   =vars(i)->dim;
  for(Variable *v:  variables) v->factors.append(this);
}

void infer::checkConsistent(const infer::Factor &f){
  CHECK_EQ(f.variables.N,f.varIds.N, "");
  CHECK_EQ(f.variables.N,f.dim.N, "");
  for(auto v: f.variables.itEnumerated()){
    CHECK_EQ(v->id ,f.varIds(v.count), "");
    CHECK_EQ(v->dim,f.dim(v.count), "");
  }
}

void infer::checkConsistent(const infer::FactorList& F){
  for(Factor *f: F) checkConsistent(*f);
}


infer::Factor::Factor(){
  specialType=NONE;
  messages.memMove=true;
  logP=0.;
}

infer::Factor::Factor(const infer::VariableList& variables, const char *_name){
  specialType=NONE;
  init(variables);
  setOne();
  if(_name) name=_name;
}

infer::Factor::Factor(const infer::VariableList& variables, const arr& q, const char *_name){
  specialType=NONE;
  init(variables);
  setP(q);
  if(_name) name=_name;
}

infer::Factor::~Factor(){
  disconnectFactor(*this);
}

void infer::Factor::operator=(const infer::Factor& q){
  init(q.variables);
  P=q.P; logP=q.logP;
}

void infer::Factor::setP(const arr& p){
  P = p;
  CHECK_EQ(P.N,product(dim), "infer::Factor set with ill-dimensioned array");
  P.reshape(dim);
  logP = 0.;
#ifdef LOG_NORM_SCALE
  lognormScale(P, logP);
#endif
}

void infer::Factor::setText(const char* str){
  rai::String(str) >>P;
  CHECK_EQ(P.N,product(dim), "infer::Factor set with ill-dimensioned array");
  P.reshape(dim);
  logP = 0.;
#ifdef LOG_NORM_SCALE
  lognormScale(P, logP);
#endif
}

void infer::Factor::setOne(){
  P.resize(dim);
  P.setUni(1.);
  logP=0.;
#ifdef LOG_NORM_SCALE
  lognormScale(P, logP);
#endif
}

void infer::Factor::setUniform(){
  double normalizedValue = 1./product(dim);
  //uint i;
  //FOR1D(dim, i)
  //  normalizedValue /= dim(i);
  P.setUni(normalizedValue);
  logP=0.;
#ifdef LOG_NORM_SCALE
  lognormScale(P, logP);
  CHECK_EQ(logP , 0., "Normalizing in lognormScale failed.");
#endif
}

void infer::Factor::setRandom(){
  rndUniform(P, .1, 1., false);
  logP=rnd.uni();
#ifdef LOG_NORM_SCALE
  lognormScale(P, logP);
#endif
}

void infer::Factor::setEvidence(uint e){
  P=0.;
  P(e)=1.;
  logP=0.;
}

bool infer::Factor::operator==(infer::Factor& q){
  //uint i=P.maxIndex();
  //double ratio=P(i)/q.P(i);
  lognormScale(P, logP, true);
  lognormScale(q.P, q.logP, true);
  bool c1 = varIds==q.varIds;
  bool c2 = dim==q.dim;
  bool c3 = maxDiff(P, q.P, 0)<1e-10;
  bool c4 = fabs(logP-q.logP)/(1+fabs(logP)+fabs(q.logP))<1e-10;
  bool r = c1  &&  c2  &&  c3  &&  c4;
  //if(!r){
  //RAI_MSG("unequal table factors: " <<id <<q.id <<dim <<q.dim <<maxDiff(P, q.P) <<' ' <<logP <<' ' <<q.logP);
  //}
  return r;
}

void infer::Factor::getP(arr& p) const {
  p = ::exp(logP)*P;
}

void infer::Factor::write(std::ostream& os, bool brief) const {
  os <<"factor " <<name <<" (";
  for(uint v=0; v<variables.N; v++){ if(v) os <<' ';  os <<variables(v)->name; }
  os <<") {";
  if(specialType){
    os <<" specialType=" <<specialType;
    if(ats.N){ os <<", "; listWrite(ats, os, ","); }
    os <<' ';
  }else{
    arr p;
    getP(p);
    char SEP=(p.nd<=1)?' ':'\n';
    os <<SEP <<"P=" <<SEP;
    p.write(os, " ", "\n ", "[]", false);
    if(ats.N){ os <<',' <<SEP; listWrite(ats, os, ","); }
    os <<SEP;
  }
  os <<'}' <<std::endl;
}

void infer::Factor::writeNice(std::ostream& os) const {
  uint i;
  os <<"[ ";
  FOR1D(variables, i)
  os<<variables(i)->name <<"  ";
  os <<"]";
}

void infer::Factor::writeExtremelyNice(std::ostream& os) const {
  writeNice(os); os<<endl;
  uint i, k;
  for(i=0; i<P.N; i++){
    os<<i <<": " <<P.elem(i) <<" (*" <<exp(logP) <<")     ";
    uintA config;
    P.getIndexTuple(config, i);
    FOR1D(variables, k){
      os<<variables(k)->name <<"=" <<config(k) <<"   ";
    }
    os<<endl;
  }
}

void infer::Factor::checkCondNormalization(uint left, double tol){
  tensorCheckCondNormalization_with_logP(P, left, logP, tol);
}

uint infer::Factor::numberNonZeroEntries(){
  uint num = 0;
  uint i;
  FOR1D(P, i){
    if(fabs(P.elem(i) > 10e-10))
      num++;
  }
  return num;
}



//===========================================================================
//
// MessagePair
//

infer::MessagePair::MessagePair(){
  f1=f2=NULL;
}

infer::MessagePair::MessagePair(infer::Factor *_f1, infer::Factor *_f2){
  init(_f1, _f2);
}

infer::MessagePair::MessagePair(infer::Variable *_v1, infer::Variable *_v2, infer::Factor *_v_to_v_fac){
  init(_v1, _v2, _v_to_v_fac);
}

infer::MessagePair::MessagePair(infer::Factor *_f1, infer::Variable *_v2){
  init(_f1, _v2);
}

infer::MessagePair::~MessagePair(){
  if(f1) f1->messages.removeValue(this);
  if(f2) f2->messages.removeValue(this);
  if(v1) v1->messages.removeValue(this);
  if(v2) v2->messages.removeValue(this);
}

void infer::MessagePair::init(infer::Factor *_f1, infer::Factor *_f2){
  f1 = _f1;
  f2 = _f2;
  v1 = v2 = NULL;
  v_to_v_fac = NULL;
  setSection(variables, f1->variables, f2->variables);
  m12.init(variables);  m12.setOne();
  m21.init(variables);  m21.setOne();
  f1->messages.append(this);
  f2->messages.append(this);
}

void infer::MessagePair::init(infer::Variable *_v1, infer::Variable *_v2, infer::Factor *_v_to_v_fac){
  f1 = f2 = NULL;
  v1 = _v1;
  v2 = _v2;
  v_to_v_fac=_v_to_v_fac;
  variables.clear();
  m12.init({v2});  m12.setOne();
  m21.init({v1});  m21.setOne();
  v1->messages.append(this);
  v2->messages.append(this);
}

void infer::MessagePair::init(infer::Factor *_f1, infer::Variable *_v2){
  f1 = _f1;
  f2 = NULL;
  v1 = NULL;
  v2 = _v2;
  v_to_v_fac = NULL;
  variables={v2};
  m12.init(variables);  m12.setOne();
  m21.init(variables);  m21.setOne();
  f1->messages.append(this);
  v2->messages.append(this);
}

void infer::MessagePair::write(ostream& os) const {
  os <<"MessagePair ";
  if(f1) os <<" fac1 " <<f1->varIds;
  if(v1) os <<" var1 [ " <<v1->id <<" ]";
  if(f2) os <<" fac2 " <<f2->varIds;
  if(v2) os <<" var2 [ " <<v2->id <<" ]";
}

void infer::MessagePair::writeIds(ostream& os) const {
  NIY;
  //rai::IOraw=true;
  //os <<"[" <<f1->varIds <<" |" <<f2->varIds <<" ]";
  //rai::IOraw=false;
}


void checkMessagePairConsistency(rai::Array<infer::MessagePair*> messages){
  uint i, k;
  int idx;
  infer::Factor* f;
  FOR1D(messages, i){
    f = messages(i)->f1;
    FOR1D(f->messages, k){
      idx = messages.findValue(f->messages(k));
      CHECK(idx>=0, "infer::Factor references an unknown message pair!");
      CHECK(f->messages(k)->f1 == f  ||  f->messages(k)->f2 == f, "infer::Factor references a foreign message pair!");
    }
    f = messages(i)->f2;
    FOR1D(f->messages, k){
      idx = messages.findValue(f->messages(k));
      CHECK(idx>=0, "infer::Factor references an unknown message pair!");
      CHECK(f->messages(k)->f1 == f  ||  f->messages(k)->f2 == f, "infer::Factor references a foreign message pair!");
    }
  }
}


//===========================================================================
//
// infer::FactorGraph
//


void infer::FactorGraph::write(std::ostream& os, bool writeBeliefs, bool writeMessages) const {
  uint i;
  os <<"<infer::FactorGraph>" <<endl;
  // <<" #variable="   <<g.variables.dims.N
  //        <<" #factors="    <<this->factors.N
  //             <<" B_v="    <<this->B_v.N
  //        <<" #msg_pairs=" <<this->messages.N <<endl;
  //os <<g.variables <<endl;
  os <<"Original factors (" <<this->F.N <<"):" <<endl;
  for(i=0; i<this->F.N; i++){
    os <<"[" <<i <<"] " <<this->F(i)->varIds <<"  " <<this->F(i) <<endl;
//     os <<*F(i);
    F(i)->writeExtremelyNice(os);
//     os <<"msg: " <<F(i)->messages <<endl;
  }
  os <<endl;
  if(writeBeliefs){
    os <<"Beliefs-cliques (" <<this->B_c.N <<"):" <<endl;
    for(i=0; i<this->B_c.N; i++){
      os <<"[" <<i <<"] " <<this->B_c(i)->varIds <<endl;
      os <<*B_c(i);
    }
    os <<endl;
    os <<"Beliefs-vars (" <<this->B_v.N <<"):" <<endl;
    for(i=0; i<this->B_v.N; i++){
      os <<"[" <<i <<"] " <<this->B_v(i)->varIds <<endl;
      os <<*B_v(i);
    }
    os <<endl;
  }
  if(writeMessages){
    os <<"messages (" <<this->messages.N <<"):" <<endl;
    //     uint i1, i2;
    for(i=0; i<this->messages.N; i++){
      os <<"[" <<i <<"] " <<this->messages(i) <<"  ";
      //       FOR1D(factors, i1){
      //           if(factors(i1)==messages(i)->f1) break;
      //       }
      //       FOR1D(factors, i2){
      //           if(factors(i2)==messages(i)->f2) break;
      //       }
      //       os<<i1 <<"-" <<i2 <<"  ";
      os<<this->messages(i)->variables  <<" [ " <<this->messages(i)->f1->varIds <<"---" <<this->messages(i)->f2->varIds <<"]  " <<messages(i)->f1 <<"  " <<messages(i)->f2 <<endl;
      os<<"m12  " <<messages(i)->m12;
      os<<"m21  " <<messages(i)->m21;
    }
  }
  os <<"</infer::FactorGraph>" <<endl;
}


void infer::FactorGraph::writeNice(std::ostream& os) const {
  uint i, j;
  os <<"Factors (specified by their vars):" <<endl;
  for(i=0; i<this->B_c.N; i++){
    os <<"[" <<i <<"] " <<endl;
    os<<this->B_c(i)->varIds.N <<"  ";
    os <<this->B_c(i)->varIds <<endl;
    FOR1D(B_c(i)->variables, j){os <<B_c(i)->variables(j)->name <<"  ";}
    os<<endl;
  }
}


void infer::FactorGraph::writeMessagePairs(std::ostream& os) const {
  uint s;
  FOR1D(this->messages, s){
    os <<*this->messages(s);
  }
}


void infer::FactorGraph::writeVariableBeliefs(std::ostream& os) const {
  uint i;
  os <<"infer::Variable beliefs:" <<endl;
  FOR1D(B_v, i){
    os <<*B_v(i) <<endl;
  }
}


void infer::FactorGraph::checkCondNormalization_B_c(double tol){
  uint i;
  FOR1D(B_c, i){
    B_c(i)->checkCondNormalization(1, tol);
  }
}


void infer::FactorGraph::setCliqueBeliefs(const infer::FactorList& fs_orig){
  uint i;
  // build belief factors if needed
  infer::Factor* b;
  if(B_c.N == 0  &&  fs_orig.N > 0){
    FOR1D(fs_orig, i){
      b = new infer::Factor(fs_orig(i)->variables);
      B_c.append(b);
    }
//     RAI_MSG("Rebuilding clique beliefs");
  }
  CHECK_EQ(fs_orig.N , B_c.N, "Number of original factors does not fit number of clique beliefs.");
  FOR1D(B_c, i){
    B_c(i)->P = fs_orig(i)->P;
    B_c(i)->logP = fs_orig(i)->logP;
  }
}


void infer::FactorGraph::resetMessages(){
  uint i;
  FOR1D(messages, i){
    messages(i)->m12.setOne();
    messages(i)->m21.setOne();
  }
}


void infer::FactorGraph::resetVariableFactors(){
  uint i;
  FOR1D(F_v, i){
    F_v(i)->setOne();
  }
}

void infer::FactorGraph::resetVariableBeliefs(){
  uint i;
  FOR1D(B_v, i){
    B_v(i)->setOne();
  }
}


void infer::FactorGraph::deleteAll(){
  uint i;
  FOR1D(messages, i) delete messages(i);
  //FOR1D(F, i)   delete F(i);
  FOR1D(F_v, i) delete F_v(i);
  FOR1D(B_c, i) delete B_c(i);
  FOR1D(B_v, i) delete B_v(i);
}


double infer::FactorGraph::computeBeliefs(){
  double change, maxChange = 0.0;
  arr P_old;
  uint i;
  FOR1D(F, i){
    P_old = B_c(i)->P;
    collectBelief(*B_c(i), *F(i), 0);
    change = absMax(B_c(i)->P - P_old);
    if(change > maxChange) maxChange = change;
  }
  FOR1D(F_v, i){
    P_old = B_v(i)->P;
    collectBelief(*B_v(i), *F_v(i), 0);
    change = absMax(B_v(i)->P - P_old);
    if(change > maxChange) maxChange = change;
  }
  return maxChange;
}


infer::Factor* infer::FactorGraph::getBelief(infer::Factor* f_orig){
  int idx;
  idx = F.findValue(f_orig); // TODO teuer, da sehr haeufig aufgerufen!
  if(idx >= 0)
    return B_c(idx);
  else {
    idx = F_v.findValue(f_orig);
    CHECK(idx>=0, "unknown factor " <<f_orig);
    return B_v(idx);
  }
}

void infer::FactorGraph::setF2Bmap(){
  if(!F2B.empty()) F2B.clear();
  uint i;
  FOR1D(F, i){F2B[F(i)] = B_c(i);}
  FOR1D(F_v, i){F2B[F_v(i)] = B_v(i);}
}

void infer::FactorGraph::addF2Bmap(infer::Factor* f, infer::Factor* b){
  F2B[f] = b;
}

infer::Factor* infer::FactorGraph::getBelief_fast(infer::Factor* f_orig){
  return F2B[f_orig];
}

void infer::FactorGraph::setV2Fv(){
  if(!V2Fv.empty()) V2Fv.clear();
  uint i, k;
  FOR1D(V, i){
    FOR1D(F_v, k){
      if(F_v(k)->varIds(0)==V(i)->id){
        V2Fv[V(i)->id] = F_v(k);
        break;
      }
    }
  }
}

void infer::FactorGraph::setV2F(){
  if(!V2F.empty()) V2F.clear();
  uint i, k;
  FOR1D(V, i){
    infer::FactorList dummy;
    V2F[V(i)->id] = dummy;
    FOR1D(F, k){
      if(F(k)->varIds(0)==V(i)->id){
        V2F[V(i)->id].append(F(k));
      }
    }
  }
}

void infer::FactorGraph::addV2Fmap(infer::Factor* f){
  V2F[f->varIds(0)].append(f);
}


// void infer::FactorGraph::checkFaithfulness(){
//   uint DEBUG=0;
//   if(DEBUG>0){cout <<"checkFaithfulness [START]" <<endl;}
//   uint i, k;
//   // for clique factors
//   if(DEBUG>0){cout <<"Checking normal factors:" <<endl;}
//   FOR1D(B_c, i){
//     if(DEBUG>0){cout <<"Belief is:" <<endl <<*B_c(i);}
//     infer::Factor product(B_c(i)->varIds);
//     product.setOne();
//     if(DEBUG>0){cout <<"Multiplying:" <<endl;}
//     if(DEBUG>0){cout <<"* orig factor:" <<endl <<*F(i);}
//     tensorMultiply(product, *F(i));
//     FOR1D(B_c(i)->messages, k){
//       if(B_c(i)->messages(k)->f1 == B_c(i)){
//         if(DEBUG>0){cout <<"* msg:" <<endl <<B_c(i)->messages(k)->m21;}
//         tensorMultiply(product, B_c(i)->messages(k)->m21);
//       }
//       else if(B_c(i)->messages(k)->f2 == B_c(i)){
//         if(DEBUG>0){cout <<"* msg:" <<endl <<B_c(i)->messages(k)->m12;}
//         tensorMultiply(product, B_c(i)->messages(k)->m12);
//       }
//       else HALT("Impossible!");
//     }
//     if(DEBUG>0){cout <<"Product is:" <<endl <<product;}
//     if(!(product == *B_c(i))){
//       cerr <<"Product:" <<endl <<product <<endl;
//       cerr <<"*B_c(i):" <<endl <<*B_c(i) <<endl;
//       CHECK_EQ(product , *B_c(i), "Faithfulness failed for factor over " <<B_c(i)->varIds);
//     }
//   }
//   if(DEBUG>0){cout <<"Checking variable factors:" <<endl;}
//   // for variable factors
//   FOR1D(B_v, i){
//     if(DEBUG>0){cout <<"Belief is:" <<endl <<*B_v(i);}
//     infer::Factor product(B_v(i)->varIds);
//     product.setOne();
//     if(DEBUG>0){cout <<"Multiplying:" <<endl;}
//     FOR1D(B_v(i)->messages, k){
//       if(B_v(i)->messages(k)->f1 == B_v(i)){
//         if(DEBUG>0){cout <<"* msg:" <<endl <<B_v(i)->messages(k)->m21;}
//         tensorMultiply(product, B_v(i)->messages(k)->m21);
//       }
//       else if(B_v(i)->messages(k)->f2 == B_v(i)){
//         if(DEBUG>0){cout <<"* msg:" <<endl <<B_v(i)->messages(k)->m12;}
//         tensorMultiply(product, B_v(i)->messages(k)->m12);
//       }
//       else HALT("Impossible!");
//     }
//     if(DEBUG>0){cout <<"Product is:" <<endl <<product;}
//     if(!(product == *B_v(i))){
//       cerr <<"Product:" <<endl <<product <<endl;
//       FOR1D(product.P, k){cout <<product.P.elem(k) <<endl;}
//       cerr <<"*B_v(i):" <<endl <<*B_v(i) <<endl;
//       FOR1D(B_v(i)->P, k){cout <<B_v(i)->P.elem(k) <<endl;}
//       CHECK_EQ(product , *B_v(i), "Faithfulness failed for factor over " <<B_v(i)->varIds);
//     }
//   }
//   if(DEBUG>0){cout <<"checkFaithfulness [END]" <<endl;}
// }


void get_vars(infer::FactorGraph& fg, uintA& varIds){
  varIds.clear();
  uint i;
  FOR1D(fg.B_c, i) varIds.setAppend(fg.B_c(i)->varIds);
}


// Edge directions:
// edges between clique factors (cf) and variable factors (vf)
// cf-->vf  iff  variable in vf is the one that is conditioned on in cf (the first one there)
// cf<--vf  otherwise
#if 0
void infer::LoopyBP::constructBipartiteinfer::FactorGraph(infer::FactorGraph& fg, const infer::FactorList& factors){
  uint i, j;
  //copy factors
  fg.F = factors;
  //generate new variable factors
  get_vars(fg.V, factors);
  infer::Factor* f;
  fg.F_v.clear();
  FOR1D(fg.V, i){
    f = new infer::Factor(TUP(*fg.V(i)));
    f->setOne();
    fg.F_v.append(f);
  }
  //connect them with messages
  infer::Variable *v;
  infer::MessagePair* s;
  FOR1D(fg.F, i){
    f=fg.F(i);
    FOR1D(f->varIds, j){
      v = getvar(f->varIds(j));
      s = new infer::MessagePair(f, v->factors.last());
      fg.messages.append(s);
    }
  }
  // beliefs
  fg.B_c.clear();
  fg.setCliqueBeliefs(factors);
  fg.B_v.clear();
  FOR1D(fg.V, i){
    f = new infer::Factor(TUP(*fg.V(i)));
    f->setOne();
    fg.B_v.append(f);
  }
  
  fg.resetMessages();
  fg.setF2Bmap();
  fg.setV2F();
  fg.setV2Fv();
}
#else
void infer::LoopyBP_obsolete::constructBipartiteFactorGraph(infer::FactorGraph& fg, const infer::FactorList& factors){
  uint i, k;
  // (1) variables
  get_vars(fg.V, factors);
  // (2) clique factors
  // original factors
  fg.F = factors;
  // beliefs
  fg.B_c.clear();
  fg.setCliqueBeliefs(factors);
  // (3) variable factors & (4) message pairs
  fg.F_v.clear();
  fg.B_v.clear();
  infer::Factor* f;
  FOR1D(fg.V, i){
    f = new infer::Factor({fg.V(i)});
    f->setOne();
    fg.F_v.append(f);
    FOR1D(fg.F, k){
      if(fg.F(k)->variables.findValue(fg.V(i)) >= 0){
        infer::MessagePair* s;
        if(fg.F(k)->variables(0) == fg.V(i))
          s = new infer::MessagePair(fg.F(k), f);
        else
          s = new infer::MessagePair(f, fg.F(k));
        fg.messages.append(s);
      }
    }
    // beliefs
    f = new infer::Factor({fg.V(i)});
    f->setOne();
    fg.B_v.append(f);
  }
  fg.resetMessages();
  fg.setF2Bmap();
  fg.setV2F();
  fg.setV2Fv();
}
#endif


//===========================================================================
//
// Algorithms on Factors and MessagePairs
//


/// collects a belief at a factor, optionally exluding one incoming message
void infer::collectBelief(infer::Factor& belief, const infer::Factor& f, const infer::MessagePair *exclude){
  belief = f;
  for(MessagePair *s: f.messages){
    if(s==exclude) continue;
    if(s->f1==&f){
      tensorMultiply(belief, s->m21);
    }else{
      CHECK(s->f2==&f, "");
      tensorMultiply(belief, s->m12);
    }
  }
}

/// collects a belief at a factor, optionally exluding one incoming message
void infer::collectBelief(infer::Factor& belief, infer::Variable *v, const MessagePair *exclude){
  belief.init({v});
  belief.setOne();
  for(MessagePair *s: v->messages){
    if(s==exclude) continue;
    if(s->v1==v){
      tensorMultiply(belief, s->m21);
    }else{
      CHECK_EQ(s->v2,v, "");
      tensorMultiply(belief, s->m12);
    }
  }
}

/// collects a belief at f1 and assigned m12 to its marginal
void infer::recomputeMessage_12(MessagePair& sep){
  CHECK((sep.f1 && !sep.v1) || (!sep.f1 && sep.v1), "");
  CHECK((sep.f2 && !sep.v2) || (!sep.f2 && sep.v2), "");
  infer::Factor belief;
  if(sep.f1){ //factor-to-variable or factor-to-factor
    collectBelief(belief, *sep.f1, &sep);
    tensorMarginal(sep.m12, belief, sep.variables);
  } else if(sep.v1 && sep.v2){ //variable-to-variable
    collectBelief(belief, sep.v1, &sep);
    tensorProductMarginal(sep.m12, *sep.v_to_v_fac, belief, {sep.v1});
  } else if(sep.v1 && sep.f2){ //variable-to-factor
    collectBelief(sep.m12, sep.v1, &sep);
  }else{ NIY; }
}

/// collects a belief at f2 and assigned m21 to its marginal
void infer::recomputeMessage_21(MessagePair& sep){
  CHECK((sep.f1 && !sep.v1) || (!sep.f1 && sep.v1), "");
  CHECK((sep.f2 && !sep.v2) || (!sep.f2 && sep.v2), "");
  infer::Factor belief;
  if(sep.f2){ //factor-to-variable or factor-to-factor
    collectBelief(belief, *sep.f2, &sep);
    tensorMarginal(sep.m21, belief, sep.variables);
  } else if(sep.v2 && sep.v1){ //variable-to-variable
    collectBelief(belief, sep.v2, &sep);
    tensorProductMarginal(sep.m21, *sep.v_to_v_fac, belief, {sep.v2});
  } else if(sep.v2 && sep.f1){ //variable-to-factor
    collectBelief(sep.m21, sep.v2, &sep);
  }else{ NIY; }
}

/// checks that marginals of connected factors are equal and the same as m12*m21
bool infer::checkConsistency(const MessagePair& sep){
  infer::Factor fmu1, fmu2, f1_marg, f2_marg, m12_m21;
  
  collectBelief(fmu1, *sep.f1, 0);
  tensorMarginal(f1_marg, fmu1, sep.variables);
  
  collectBelief(fmu2, *sep.f2, 0);
  tensorMarginal(f2_marg, fmu2, sep.variables);
  
  m12_m21=sep.m12;
  tensorMultiply(m12_m21, sep.m21);
  
  CHECK_EQ(f1_marg,f2_marg, "marginals inconsistent");
  CHECK_EQ(f1_marg,m12_m21, "marginals inconsistent");
  CHECK_EQ(m12_m21,f2_marg, "marginals inconsistent");
  return true;
}

/// checks a whole list of messages
bool infer::checkConsistencyBatch(const MessagePairList& msgs){
  uint i;
  FOR1D(msgs, i) checkConsistency(*msgs(i));
  return true;
}

#if 1

// Marc's eq (5) and (7)
void infer::computeMessage_noDiv(infer::Factor& f_from, infer::Factor& f_to){
  uint DEBUG = 0;
  if(DEBUG>0){cout <<"computeMessage_noDivision [START]" <<endl;}
  uint i;
  FOR1D(f_from.messages, i){
    if(f_from.messages(i)->f1 == &f_to){
      recomputeMessage_21(*f_from.messages(i));
      break;
    } else if(f_from.messages(i)->f2 == &f_to){
      recomputeMessage_12(*f_from.messages(i));
      break;
    }
  }
  CHECK(f_from.messages.N > i, "Factors don't have a MessagePair!");
  if(DEBUG>0){cout <<"computeMessage_noDivision [END]" <<endl;}
}


void infer::computeMessage_withDiv(infer::Factor& f_from, infer::Factor& f_to){
  uint DEBUG = 0;
  if(DEBUG>0){cout <<"calcMessage [START]" <<endl;}
  // find MessagePair
  infer::MessagePair* s = getMessagePair(&f_from, &f_to);
  if(&f_from == s->f2){
    if(DEBUG > 0){cout <<f_from.varIds <<" --> " <<s->f1->varIds<<endl;}
    if(DEBUG > 2){
      cout <<"Message 2 --> 1" <<endl;
      cout <<"f_from: " <<endl <<f_from;
    }
    infer::Factor belief;
    collectBelief(belief, f_from, 0);
    tensorMarginal(s->m21, belief, s->variables);
    if(DEBUG>2){
      cout <<"Belief over f_from: " <<endl <<belief;
      cout <<"Old msg in other direction (m12): " <<endl <<s->m12;
    }
    tensorDivide(s->m21, s->m12);
    if(DEBUG > 1) cout <<"New msg (m21): " <<endl <<s->m21;
  }
  // message m12
  else {
    if(DEBUG > 0){cout <<f_from.varIds <<" --> " <<s->f2->varIds<<endl;}
    if(DEBUG > 2){
      cout <<"Message 1 --> 2" <<endl;
      cout <<"f_from: " <<endl <<f_from <<endl;
    }
    infer::Factor belief;
    collectBelief(belief, f_from, 0);
    tensorMarginal(s->m12, belief, s->variables);
    if(DEBUG>2){
      cout <<"Belief over f_from: " <<endl <<belief;
      cout <<"Old msg in other direction (m21): " <<endl <<s->m21;
    }
    tensorDivide(s->m12, s->m21);
    if(DEBUG > 1) cout <<"New msg (m12): " <<endl <<s->m12;
  }
  if(DEBUG>0){cout <<"lbp_calcMessage [END]" <<endl;}
}


void infer::computeMessage_noDiv(infer::Factor& f_from, infer::Factor& b_from, infer::Factor& f_to){
  HALT("this is actually with div...");
  uint DEBUG = 0;
  if(DEBUG>0){cout <<"calcMessage [START]" <<endl;}
  // find MessagePair
  infer::MessagePair* s = getMessagePair(&f_from, &f_to);
  if(&f_from == s->f2){
    if(DEBUG > 0){cout <<f_from.varIds <<" --> " <<s->f1->varIds<<endl;}
    if(DEBUG > 2){
      cout <<"Message 2 --> 1" <<endl;
      cout <<"f_from: " <<endl <<f_from;
    }
    tensorMarginal(s->m21, b_from, s->variables);
    if(DEBUG>2){
      cout <<"Marginalised f_from: " <<endl <<s->m21;
      cout <<"Old msg in other direction (m12): " <<endl <<s->m12;
    }
    tensorDivide(s->m21, s->m12);
    if(DEBUG > 1) cout <<"New msg (m21): " <<endl <<s->m21;
  }
  // message m12
  else {
    if(DEBUG > 0){cout <<f_from.varIds <<" --> " <<s->f2->varIds<<endl;}
    if(DEBUG > 2){
      cout <<"Message 1 --> 2" <<endl;
      cout <<"f_from: " <<endl <<f_from <<endl;
    }
    tensorMarginal(s->m12, b_from, s->variables);
    if(DEBUG>2){
      cout <<"Marginalised f_from: " <<endl <<s->m12;
      cout <<"Old msg in other direction (m21): " <<endl <<s->m21;
    }
    tensorDivide(s->m12, s->m21);
    if(DEBUG > 1) cout <<"New msg (m12): " <<endl <<s->m12;
  }
  if(DEBUG>0){cout <<"lbp_calcMessage [END]" <<endl;}
}






double infer::passMessage(infer::Factor& f_from, infer::Factor& f_to, infer::Factor& b_to, MsgCalc calcMsgType){
  uint DEBUG = 0;
  
  // store old
  // delete this later (for efficiency reasons)
  infer::Factor belief_to_old;
  collectBelief(belief_to_old, f_to, 0);
  arr p_old = belief_to_old.P;
  
  if(DEBUG>0){
    cout <<"passMessage: " <<endl;
    cout <<"f_from: " <<f_from <<endl;
    cout <<"f_to: " <<f_to <<endl;
  }
  
  // (1) calc new message
  if(calcMsgType == WITH_DIV)
    computeMessage_withDiv(f_from, f_to);
  else
    computeMessage_noDiv(f_from, f_to);
  // (2) recalculate belief accordingly
  collectBelief(b_to, f_to, 0);
  
  
  if(DEBUG>0){
    cout <<"Message:" <<endl;
    writeMessage(&f_from, &f_to, cout); cout <<endl;
    cout <<"Updated b_to:" <<endl <<b_to <<endl;
    cout <<" --> " <<absMax(b_to.P - p_old) <<endl;
  }
  
  // calc change (and ignore log_P !)
  return absMax(b_to.P - p_old);
}

#endif

// double distributeMessages(infer::FactorGraph& fg, infer::Factor& f, MsgCalc calcMsgType){
//   uint DEBUG = 0;
//   uint i, index;
//   double change, max_change=0.0;
//   if(DEBUG>0){cout <<"Distributing messages from factor over " <<f.varIds<<endl;}
//   if(DEBUG>1){cout <<f <<endl;}
//   FOR1D(f.messages, i){
//     infer::Factor* neighbor;
//     infer::Factor* neighbor_belief;
//     if(f.messages(i)->f1 == &f){
//       index = get_list_index_unsigned(fg, f.messages(i)->f2);
//       neighbor = fg.
//       if(DEBUG>1){cout <<"to (before) "; f.messages(i)->f2->writeExtremelyNice(cout);cout <<endl;}
//       change = passMessage_1(f, *f.messages(i)->f2, calcMsgType); // based on Marc's eq (2) / (3)
//       if(max_change<change) max_change=change;
//       if(DEBUG>1){cout <<"to (after) "; f.messages(i)->f2->writeExtremelyNice(cout);cout <<endl;}
//     }
//     else {
//       if(DEBUG>1){cout <<"to (before) "; f.messages(i)->f1->writeExtremelyNice(cout);cout <<endl;}
//       change = fg.passMessage_1(f, *f.messages(i)->f1, calcMsgType); // based on Marc's eq (2) / (3)
//       if(max_change<change) max_change=change;
//       if(DEBUG>1){cout <<"to (after) "; f.messages(i)->f1->writeExtremelyNice(cout);cout <<endl;}
//     }
//   }
//   return max_change;
// }


void getMarginal(infer::Factor& marginal, const infer::VariableList& marginalVars, infer::FactorGraph& fg){
  // search for factors containing marginalVars
  uint f, v;
  FOR1D(fg.B_c, f){
    FOR1D(marginalVars, v){
      if(fg.B_c(f)->variables.findValue(marginalVars(v))<0)
        break;
    }
    if(marginalVars.N > v)
      continue;
    else {
      tensorMarginal(marginal, *fg.B_c(f), marginalVars);
      return;
    }
  }
  HALT("No factor with all marginalVars " <<marginalVars <<" -- "); //getvar(marginalVars(0))->name);
}


/* return a list of posteriors for each variable.
  We step through all variables, get the (last) factor over
  this variable, collect a belief for it, and marginalize it */
void infer::getVariableBeliefs(rai::Array<arr>& post, const infer::VariableList& vars){
  uint i, N=vars.N;
  post.resize(N);
  infer::Factor *f, belief, marg;
  for(i=0; i<N; i++){
    if(!vars(i)->factors.N){ post(i).resize(vars(i)->dim); post(i) = 1.; continue; }
    f=vars(i)->factors.last(); //simply take last factor in factor list!
    collectBelief(belief, *f, NULL);
    tensorMarginal(marg, belief, {vars(i)});
    post(i) = marg.P;
  }
}






void printFactorArrayViaVariableSets(infer::FactorList factors){
  uint i;
  FOR1D(factors, i){
    cout <<factors(i)->varIds <<" ";
  }
  cout <<endl;
}







//===========================================================================
//
//  infer::Factor Algebra

double robustDivide(double a, double b){
  if(!a){ if(!b) return 1.; else return 0.; }
  //if(a==b) return 1.;
  if(b) return a/b;
  HALT("division by zero: " <<a <<"/" <<b);
  return 0;
}

/** @brief If id is a tuple of k variable identifiers, and mid is a _subset_ of l of those,
then pick will contain l indices, shouting which slot of id corresponds to each slot of mid. For instance
if we have a product of two sensors \f$A_{ijklm} B_{kli}\f$ then pick will be <2 3 0> shouting that
the 0th slot of B is the 2nd of A, that the 1st slot of B is the 3rd of A, and that the 3rd slot of B is the 0th of A. */
void getPick(uintA& pick, const infer::VariableList& base_vars, const infer::VariableList& multiplier_vars){
  uint i=0, k=0;
  pick.resize(multiplier_vars.N);
  for(k=0; k<multiplier_vars.N; k++){
    for(i=0; i<base_vars.N; i++) if(base_vars(i)==multiplier_vars(k)) break;
    CHECK(i<base_vars.N, "getPick error: base_vars=" <<base_vars <<" multiplier_vars=" <<multiplier_vars);
    pick(k)=i;
  }
  DEBUG_INFER(3, cout <<"** getPick: base_vars=" <<base_vars <<" multiplier_vars=" <<multiplier_vars <<" pick=" <<pick <<endl);
}

void infer::tensorProduct(infer::Factor& f, const infer::Factor& a, const infer::Factor& b){
  infer::VariableList fvars;
  setUnion(fvars, a.variables, b.variables);
  f.init(fvars);
  f.P.resize(f.dim);
  uintA pickA, pickB;
  getPick(pickA, f.variables, a.variables);
  getPick(pickB, f.variables, b.variables);
  tensorEquation(f.P, a.P, pickA, b.P, pickB, 0);
  f.logP = a.logP + b.logP;
  lognormScale(f.P, f.logP);
}

void infer::tensorProductMarginal(infer::Factor& f, const infer::Factor& a, const infer::Factor& b, const infer::VariableList& s){
  infer::VariableList fvars;
  setUnion(fvars, a.variables, b.variables);
  setMinus(fvars, s);
  f.init(fvars);
  f.P.resize(f.dim);
  uintA pickA, pickB;
  infer::VariableList all;
  all=f.variables;
  for(uint i=0; i<s.N; i++) all.append(s(i)); //all.append(s);
  getPick(pickA, all, a.variables);
  getPick(pickB, all, b.variables);
  tensorEquation(f.P, a.P, pickA, b.P, pickB, s.N);
  f.logP = a.logP + b.logP;
  lognormScale(f.P, f.logP);
}

void infer::tensorMarginal(infer::Factor& m, const infer::Factor& f, const infer::VariableList& marginalVars){
  m.init(marginalVars);
  uintA pick;
  getPick(pick, f.variables, m.variables);
  ::tensorMarginal(m.P, f.P, pick);
  m.logP=f.logP;
  lognormScale(m.P, m.logP);
}

void infer::tensorMaxMarginal(infer::Factor& m, const infer::Factor& f, const infer::VariableList& marginalVars){
  m.init(marginalVars);
  uintA pick;
  getPick(pick, f.variables, m.variables);
  ::tensorMaxMarginal(m.P, f.P, pick);
  //do we need to normalize in the max-marginal case?
  //normalizeDist(Y);
  m.logP=f.logP;
  lognormScale(m.P, m.logP);
}

void infer::tensorMultiply(infer::Factor& f, const infer::Factor& m){
  if(m.variables==f.variables){
    f.P *= m.P;
  }else{
    uintA pick;
    getPick(pick, f.variables, m.variables);
    ::tensorMultiply(f.P, m.P, pick);
  }
  f.logP += m.logP;
  lognormScale(f.P, f.logP);
}

void infer::tensorDivide(infer::Factor& f, const infer::Factor& m){
  if(m.variables==f.variables){
    f.P /= m.P;
  }else{
    uintA pick;
    getPick(pick, f.variables, m.variables);
    ::tensorDivide(f.P, m.P, pick);
  }
  f.logP -= m.logP;
  //   cout <<"infer.cpp.tensorDivide vormScalen: " <<f.P <<endl;
  lognormScale(f.P, f.logP);
  //   cout <<"infer.cpp.tensorDivide nachmScalen: " <<f.P <<endl;
}

void infer::tensorAdd(infer::Factor& f, const infer::Factor& m){
  arr mP=m.P;
  mP *= ::exp(m.logP-f.logP); //get m.P on the same log scale as f!
  if(m.variables==f.variables){
    f.P += mP;
  }else{
    uintA pick;
    getPick(pick, f.variables, m.variables);
    ::tensorAdd(f.P, mP, pick);
  }
  lognormScale(f.P, f.logP);
}

void infer::tensorInvertMultiply(infer::Factor& f, const infer::Factor& m){
  CHECK_EQ(m.variables,f.variables, "infer::Factor invMultiply needs identical variables");
  for(uint i=0; i<f.P.N; i++) f.P.elem(i)=robustDivide(m.P.elem(i), f.P.elem(i));
  f.logP = -f.logP + m.logP;
  lognormScale(f.P, f.logP);
}

void infer::tensorWeightedAdd(infer::Factor& f, double w, const infer::Factor& m){
  CHECK_EQ(m.variables,f.variables, "infer::Factor weightedAdd needs identical variables");
  w *= ::exp(m.logP-f.logP);
  for(uint i=0; i<f.P.N; i++) f.P.elem(i) = f.P.elem(i) + w*m.P.elem(i);
  lognormScale(f.P, f.logP);
}

/*void writeInfo(ostream& os, const infer::Factor& f){
  uint i;
  os <<"vars=(";
  for(i=0;i<f.varIds.N;i++){
    if(i) os <<' ';
    os <<f.varIds(i) <<'.' <<globalVariableList(f.varIds(i))->name;
    CHECK_EQ(f.varIds(i),globalVariableList(f.varIds(i))->id, "identity mismatch!!");
  }
  os <<") dims=<";
  for(i=0;i<f.varIds.N;i++){
    if(i) os <<' ';
    os <<f.dim(i);
    CHECK_EQ(f.dim(i),globalVariableList(f.varIds(i))->dim, "dimensionality mismatch!!");
  }
  os <<"> P.dims=";
  f.P.writeDim(os);
  os <<endl;
}*/




//===========================================================================
//
// ELIMINATION ALGORITHM methods
//
void infer::getJoint(infer::Factor& joint, const infer::FactorList& factors){
  DEBUG_INFER(1, cout <<RAI_HERE <<endl);
  uint i;
  //get tuple of vars
  infer::VariableList jointVars;
  for(i=0; i<factors.N; i++) jointVars.setAppend(factors(i)->variables);
  DEBUG_INFER(2, cout <<"  jointVars=" <<jointVars <<endl);
  //compute joint
  joint.init(jointVars);
  joint.setOne();
  for(i=0; i<factors.N; i++) tensorMultiply(joint, *(factors(i)));
}

/** Computes an order for the elimination of variables.
The order is such that in each iteration among all variables that have not been eliminated yet
the one which has the fewest links to other not-deleted variables is chosen.
This is equivalent to saying that per iteration the variable chosen for elimination is the one
that would create the smallest clique if all factors that involve this variable were multiplied.
*/
void infer::computeEliminationOrder(infer::VariableList& elimOrder, const infer::FactorList& factors, const infer::VariableList& elimVars){
  int DEBUG_INFER_LEVEL = 0;
  DEBUG_INFER(1, cout <<RAI_HERE <<endl);
  DEBUG_INFER(2, cout <<"  input factors=\n" <<factors <<endl);
  DEBUG_INFER(1, cout <<"variables to eliminate=" <<elimVars <<endl);
  
  infer::VariableList vars;
  get_vars(vars, factors);
  
  elimOrder.resize(elimVars.N);
  uint f, v, e;
  
  // Determine for each variable the set of variables it is linked to.
  rai::Array<infer::VariableList> connectedVarSets(elimVars.N);
  for(v=0; v<elimVars.N; v++){
    for(f=0; f<factors.N; f++){
      if(factors(f)->variables.findValue(elimVars(v))>-1){
        connectedVarSets(v).setAppend(factors(f)->variables);
      }
    }
    DEBUG_INFER(2, cout <<"  neighbors of v" <<elimVars(v) <<"=" <<connectedVarSets(v) <<endl);
  }
  // Calculate elimination order
  boolA used(elimOrder.N);
  used=false;
  DEBUG_INFER(1, cout <<"calculation of elimination order: " <<endl);
  uint cost, minCost=-1;
  int minIndex;
  for(e=0; e<elimOrder.N; e++){
    DEBUG_INFER(2, cout <<"  (" <<e <<")" <<endl);
    // Determine variable with smallest set of connected remaining variables
    minIndex=-1;
    for(v=0; v<elimVars.N; v++){
      if(used(v)) continue; //don't consider already used vars
#if 1
      //determine cost as the #variables in the generated clique
      cost=connectedVarSets(v).N;
#else
      //determine cost as the size (=size of probability table) of the generated clique
      cost=1;
      for(uint i=0; i<connectedVarSets(v).N; i++)
        if(connectedVarSets(v)(i)!=v)  cost*=globalVariableList(connectedVarSets(v)(i))->dim;
#endif
      DEBUG_INFER(2, cout <<"  testing v" <<elimVars(v) <<": cost=" <<cost <<endl);
      if(minIndex==-1 || minCost>cost){
        minCost = cost;
        minIndex = v;
      }
    }
    elimOrder(e) = elimVars(minIndex);
    used(minIndex) = true;
    DEBUG_INFER(1, cout <<"eliminating v" <<elimOrder(e) <<" '" <<elimOrder(e)->name <<"' (cost=" <<minCost <<")" <<endl);
    // Delete eliminated variable out of sets of connected variables of other vars
    for(v=0; v<elimVars.N; v++){
      connectedVarSets(v).removeValue(elimOrder(e, false));
      DEBUG_INFER(2, cout <<"  neighbors of v" <<elimVars(v) <<"=" <<connectedVarSets(v) <<endl);
    }
  }
  DEBUG_INFER(1, cout <<"elimination order: " <<elimOrder <<endl);
}

/** eliminates a single variable from the factor list.

At input, factors describes the full model; at output, factors contains the reduced model which include unchanged old factors and some newed factors.
The newed factors are additionally appended to the newed_factors (to allow for external cleanup) */
void infer::eliminateVariable(infer::FactorList& factors, infer::FactorList& newed_factors, infer::Variable *var){
  uint f;
  infer::FactorList referencedFactors;
  factors.memMove=true;
  referencedFactors.memMove=true;
  //collect referenced factors
  for(f=0; f<factors.N; f++)
    if(factors(f)->variables.findValue(var)!=-1) referencedFactors.append(factors(f));
    
  //compute joint variable tuple
  infer::VariableList jointVars;
  for(f=0; f<referencedFactors.N; f++)
    jointVars.setAppend(referencedFactors(f)->variables);
  //setUnion(jointVars, jointVars, referencedFactors(f)->varIds);//[mt]
  
  checkConsistent(factors);
  
  //compute joint clique and remove! the used factors from the list
  infer::Factor jointFactor(jointVars);
  //jointFactor.setOne(); is done in the constructor
  for(f=0; f<referencedFactors.d0; f++){
    tensorMultiply(jointFactor, *(referencedFactors(f)));
    factors.removeValue(referencedFactors(f));
  }
  checkConsistent(jointFactor);
  
  //eliminate the single variable from the joint factor
  jointVars.removeValue(var);
  infer::Factor *newFactor = new infer::Factor(jointVars);
  tensorMarginal(*newFactor, jointFactor, jointVars);
  factors.append(newFactor);
  newed_factors.append(newFactor);
  checkConsistent(*newFactor);
  
  //DEBUG_INFER(1, checkSpaceRegistry());
}

/** marginalizes a factor list over all variables except the "remaining_vars". The output is a
single factor over the remaining_vars with the marginal. The factors list remains unchanged. */
void infer::eliminationAlgorithm(infer::Factor& posterior, const infer::FactorList& factors, const infer::VariableList& remaining_vars){
  DEBUG_INFER(1, cout <<RAI_HERE <<endl);
  uint i, f;
  checkConsistent(factors);
  
  // determine which variables need to be eliminated: ALL \ post_vars
  infer::VariableList facVars;
  for(f=0; f<factors.N; f++) facVars.setAppend(factors(f)->variables);
  
  DEBUG_INFER(3, cout <<"  all facs=" <<factors <<endl);
  DEBUG_INFER(2, cout <<"  factor vars=" <<facVars <<"\n  remaining vars=" <<remaining_vars <<endl);
  
  infer::VariableList elimVars=facVars; elimVars.memMove=true;
  elimVars.setAppend(remaining_vars); //in case the posterior wants more variables that the factors are defined over
  for(i=0; i<remaining_vars.N; i++) elimVars.removeValue(remaining_vars(i));
  
  DEBUG_INFER(2, cout <<"  elim vars=" <<elimVars <<endl);
  
  // determine order in which variables are eliminated.
  infer::VariableList elimOrder;
  computeEliminationOrder(elimOrder, factors, elimVars);
  
  DEBUG_INFER(2, cout <<"  elimination order=" <<elimOrder <<endl);
  
  // eliminate in this order
  infer::FactorList factors_copy(factors);
  infer::FactorList newedFactors;
  for(i=0; i<elimOrder.N; i++) eliminateVariable(factors_copy, newedFactors, elimOrder(i));
  
  // calculate posterior
  posterior.init(remaining_vars);
  posterior.setOne();
  /*if(remaining_vars.N==0){ //DON'T HANDLE SCALAR OUTPUT SPECIAL...
    for(i=0;i<factors_copy.N;i++){
      CHECK_EQ(factors_copy.N , 1, "too many factors created");
    posterior.P.scalar()=sum(factors_copy(0)->P) * exp(factors_copy(0)->logP);
  }else{*/
  for(i=0; i<factors_copy.N; i++){
    DEBUG_INFER(3, cout <<"  remaining factor " <<i <<"=" <<*factors_copy(i) <<endl);
    CHECK(&posterior!=factors_copy(i), "one of the input factors is also the output factor of the elimination algorithm - that doesn't work!!");
    tensorMultiply(posterior, *factors_copy(i));
  }
  //}
  
  DEBUG_INFER(3, cout <<"  posterior=" <<posterior <<endl);
  
  //cleanup newed factores
  for(i=0; i<newedFactors.N; i++) delete newedFactors(i);
}


void moralize(infer::FactorList& factorsOfDirectedGraph, infer::FactorList& factorsOfMoralizedGraph){
  // first variable is child, remaining variables are parents
  NIY;
}


//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//===========================================================================
//
// JUNCTION TREE methods
//

void checkJunctionTreeProperty_dfs(infer::Factor* node, infer::Factor* parent, infer::FactorGraph& junctionTree, boolA& containsId, infer::Variable *id){
  int DEBUG = 0;
  uint i;
  containsId(get_list_index_unsigned__orig(junctionTree, node)) = 0;
  infer::FactorList neighbors;
  getNeighbors(node, neighbors);
  FOR1D(neighbors, i){
    if(neighbors(i) == parent)
      continue;
    if(DEBUG > 0)
      cout <<"Checking direction " <<node->varIds <<" --> " <<neighbors(i)->varIds <<endl;
    // check whether neighbor contains varible
    if(neighbors(i)->variables.findValue(id) != -1){
      if(DEBUG > 0)
        cout <<"   --> going on that way" <<endl;
      // check that variable contained in MessagePair
#ifndef RAI_NOCHECK
      infer::MessagePair* s = getMessagePair(node, neighbors(i));
#endif
      CHECK(s->variables.findValue(id) != -1, "infer::Variable " <<id <<" is not contained in MessagePair " <<s->variables <<" for factors over " <<node->varIds <<" and " <<neighbors(i)->varIds <<"!" <<endl);
      // continue dfs here
      checkJunctionTreeProperty_dfs(neighbors(i), node, junctionTree, containsId, id);
    }else{
      if(DEBUG > 0)
        cout <<"   --> dead end!" <<endl;
      continue;
    }
  }
}

// Jordan, Chapter 17, p.22: For every pair of cliques V and W, all cliques on (unique) path between V and W contain section(V, W).
// --> checks that each variable induces a single subtree in junction tree
void infer::JunctionTree::checkJunctionTreeProperty(FactorGraph& junctionTree){
  int DEBUG = 0;
  if(DEBUG > 0)
    cout <<"checkJunctionTreeProperty [START]" <<endl;
    
  uint v, f;
  infer::Variable *id;
  boolA containsId(junctionTree.B_c.N);
  FOR1D(junctionTree.V, v){
    id = junctionTree.V(v);
    if(DEBUG > 0)
      cout <<" id=" <<id <<endl;
    FOR1D(junctionTree.B_c, f){
      if(junctionTree.B_c(f)->variables.findValue(id) != -1)
        containsId(f) = 1;
      else
        containsId(f) = 0;
    }
    if(DEBUG > 0)
      cout <<" containsId_before=" <<containsId <<endl;
    f=0;
    while(f<junctionTree.B_c.N && !containsId(f))
      f++;
    CHECK(f<junctionTree.B_c.N, " infer::Variable ID " <<id <<" has not been found in graph." <<endl);
    checkJunctionTreeProperty_dfs(junctionTree.B_c(f), 0, junctionTree, containsId, id);
    if(DEBUG > 0)
      cout <<" containsId_after=" <<containsId <<endl;
    FOR1D(containsId, f)
    CHECK_EQ(containsId(f) , 0, " junction tree property violated for variable " <<id <<endl);
  }
  if(DEBUG > 0)
    cout <<"checkJunctionTreeProperty [END]" <<endl;
}



void infer::JunctionTree::addEvidence(FactorGraph& junctionTree, infer::Factor& evid){
  uint f;
  infer::VariableList varSection;
  FOR1D(junctionTree.F, f){
    setSection(varSection, evid.variables, junctionTree.F(f)->variables);
    if(varSection.N > 0){
      infer::MessagePair* s = new infer::MessagePair(junctionTree.F(f), &evid);
      junctionTree.F.append(&evid);
      // update orig factor
      junctionTree.B_c.append(new infer::Factor());
      *junctionTree.B_c.last() = evid;
      junctionTree.messages.append(s);
      break;
    }
  }
}





void infer::JunctionTree::buildTriangulatedCliques(const infer::FactorList& factors, infer::FactorList& triangulatedCliques){
  bool DEBUG = false;
  bool DEBUG_VERBOSE = false;
  
  uint v, v2, v3, f;
  infer::FactorList intermediateFactors = factors;
  boolA intermediateFactors_removed(intermediateFactors.N);
  FOR1D(intermediateFactors_removed, v)
  intermediateFactors_removed(v) = 0;
  
  // determine existing variables
  infer::VariableList vars;
  for(f=0; f < factors.N; f++){
    setUnion(vars, vars, factors(f)->variables);
  }
  std::sort(vars.p, vars.p + vars.N);
  if(DEBUG)
    cout <<"Existing vars: " <<vars <<endl;
    
    
  // map varID -> order over participating vars
  std::map<infer::Variable*, uint> var_id2order;
  FOR1D(vars, v){
    var_id2order[vars(v)] = v;
  }
  if(DEBUG){
    std::map<infer::Variable*, uint>::iterator itVars;
    for(itVars = var_id2order.begin() ; itVars != var_id2order.end(); itVars++){
      cout <<"var_id2order[" <<itVars->first <<"] = " <<itVars->second <<endl;
    }
  }
  
  // determine elimination ordering
  //  uintA elimOrder(6);
  //  for(v=0; v<6; v++){
  //    elimOrder(v) = 5-v;
  //  }
  infer::VariableList elimOrder;
  computeEliminationOrder(elimOrder, intermediateFactors, vars);
  //     elimOrder <<"[   9   14   7   8   10   11   12   13   0    1    2    15   16   3   4   5   6]";
  //     elimOrder <<"[   12   17   10   11   13   14   15   16  7   0   1   2   8   9   3   4   5   6 ]";
  //     elimOrder <<"[   19   29   16   17   18   20   21   22   23   24   25   26   27   28   0   12   1   2   13   14   15   3   4   5   6   7   8   9   10   11]";
  //     elimOrder = global_elimOrder;
  
  if(DEBUG){
    cout <<"Elimination order: " <<elimOrder <<endl;
    if(DEBUG_VERBOSE){
      FOR1D(elimOrder, v){
        cout <<elimOrder(v)->name <<"  ";
      }
      cout <<endl;
    }
  }
  
  // edges matrix
  // indexed by var ids
  rai::Array<byte> edges(elimOrder.N, elimOrder.N);
  FOR2D(edges, v, v2){
    edges(v, v2) = 0;
  }
  // init edges matrix with existing edges in original graph
  FOR1D(intermediateFactors, f){
    FOR1D(intermediateFactors(f)->variables, v){
      for(v2 = v+1; v2 < intermediateFactors(f)->variables.N; v2++){
        edges(var_id2order[intermediateFactors(f)->variables(v)], var_id2order[intermediateFactors(f)->variables(v2)]) = 1;
        edges(var_id2order[intermediateFactors(f)->variables(v2)], var_id2order[intermediateFactors(f)->variables(v)]) = 1;
      }
    }
  }
  //  if(DEBUG)
  //    cout <<"Edges matrix:" <<endl <<edges <<endl;
  
  
  // UNDIRECTED_GRAPH_ELIMINATE
  
  FOR1D(elimOrder, v){
    if(DEBUG){
      cout <<"--------------------------------------------------" <<endl;
      cout <<"(" <<v <<")" <<" Eliminating " <<elimOrder(v) <<endl;
      if(DEBUG_VERBOSE){cout <<"(" <<v <<")" <<" Eliminating " <<elimOrder(v)->name <<endl;}
    }
    // determine remaining neighbors
    infer::VariableList remainingNeighbors;
    FOR1D(edges, v2){
      if(edges(var_id2order[elimOrder(v)], var_id2order[elimOrder(v2)]))
        remainingNeighbors.append(elimOrder(v2));
    }
    if(DEBUG){
      cout <<"Remaining neighbors:" <<remainingNeighbors <<endl;
      //if(DEBUG_VERBOSE){cout <<"Remaining neighbors: "; printvars(remainingNeighbors); cout <<endl;}
    }
    
    // connect remaining neighbors
    FOR1D(remainingNeighbors, v2){
      for(v3 = v2+1; v3<remainingNeighbors.N; v3++){
        edges(var_id2order[remainingNeighbors(v2)], var_id2order[remainingNeighbors(v3)]) = 1;
        edges(var_id2order[remainingNeighbors(v3)], var_id2order[remainingNeighbors(v2)]) = 1;
      }
    }
    //    if(DEBUG)
    //      cout <<"Edges matrix (remaining neighbors connected):" <<endl <<edges <<endl;
    
    // set up new factor with all neighbors
    if(DEBUG)
      cout <<"Setting up new factor [START]" <<endl;
    infer::VariableList vars;
    vars.append(elimOrder(v));
    setUnion(vars, vars, remainingNeighbors);
    if(DEBUG)
      cout <<"Variables [" <<vars.N <<"]: " <<vars <<endl;
    //if(DEBUG_VERBOSE){cout <<"Variables: "; printvars(vars);cout <<endl;}
    // check whether "vars" of new factors are subset of existing factor f
    bool containedInOtherFactor = 0;
    FOR1D(triangulatedCliques, f){
      if(numberSharedElements(triangulatedCliques(f)->variables, vars) == vars.N){
        containedInOtherFactor = true;
        break;
      }
    }
    infer::Factor* newFactor;
    // if so, include into f
    if(containedInOtherFactor){
      newFactor = triangulatedCliques(f);
      if(DEBUG)
        cout <<"--> included in " <<newFactor->varIds <<endl;
      //if(DEBUG_VERBOSE){cout <<"--> included in "; printvars(newFactor->varIds);cout <<endl;}
    }
    // else construct new factor
    else {
      if(DEBUG)
        cout <<"Spaeter: Variables [" <<vars.N <<"]: " <<vars <<endl;
      // TODO weg
      //             if(vars.N == 28) continue;
      
      newFactor = new infer::Factor(vars);
      newFactor->setOne();
      // add new factor to triangulated set
      triangulatedCliques.append(newFactor);
      if(DEBUG)
        cout <<"--> true new factor" <<endl;
    }
    
    if(DEBUG)
      cout <<"Calculating new potential:" <<endl;
    // calc new potential
    for(f=0; f < intermediateFactors.N; f++){
      if(intermediateFactors_removed(f) != 1
          && intermediateFactors(f)->variables.findValue(elimOrder(v)) > -1){
        if(DEBUG)
          cout <<" ** Using factor: " <<*intermediateFactors(f);
        tensorMultiply(*newFactor, *(intermediateFactors(f)));
        if(DEBUG)
          cout <<"  --> updated newFactor: " <<*newFactor;
        intermediateFactors_removed(f) = 1;
      }
      if(DEBUG){
        cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" <<endl;
        cout <<"All factors thus far:" <<endl;
        FOR1D(triangulatedCliques, v2){
          cout <<v2 <<":  " <<triangulatedCliques(v2)->varIds.N <<endl;
        }
      }
    }
    if(DEBUG){
      cout <<"Final newFactor: " <<endl <<*newFactor;
      //if(DEBUG_VERBOSE){cout <<" over vars "; printvars(newFactor->varIds);cout <<endl;}
      cout <<"Setting up new factor [END]" <<endl;
    }
    
    // remove node from graph
    for(v2=0; v2 < elimOrder.N; v2++){
      edges(var_id2order[elimOrder(v)], var_id2order[elimOrder(v2)]) = 0;
      edges(var_id2order[elimOrder(v2)], var_id2order[elimOrder(v)]) = 0;
    }
    
    if(DEBUG){
      //      cout <<"Edges matrix (node removed):" <<endl <<edges <<endl;
      cout <<"Removed factors: " <<intermediateFactors_removed <<endl;
    }
  }
}


// Kruskal algorithm
void infer::JunctionTree::buildMaxSpanningTree(infer::FactorList& factors, const infer::VariableList& vars, FactorGraph& cliqueTree){
  uint DEBUG = 0;
  if(DEBUG >= 1){
    cout <<"========================================" <<endl;
    cout <<"buildMaxSpanningTree" <<endl;
    cout <<"--------------------" <<endl;
  }
  
  uint f, f2, maxId;
  uint max;
  infer::Factor *fac1, *fac2;
  infer::VariableList tempA;
  
  if(DEBUG >= 1){
    cout <<"input factors [START]" <<endl;
    FOR1D(factors, f){
      cout <<"(" <<f <<") " <<factors(f)->varIds <<endl;
    }
    cout <<"input factors [END]" <<endl;
  }
  
  // check that each factor has a non-empty MessagePair for some other factor
  // = shares variables with some other factor
  // = not disconnected from tree
  boolA sharesVariables(factors.N);
  FOR1D(sharesVariables, f)
  sharesVariables(f) = false;
  FOR1D(factors, f){
    if(sharesVariables(f))
      continue;
    FOR1D(factors, f2){
      if(f == f2)
        continue;
      setSection(tempA, factors(f)->variables, factors(f2)->variables);
      if(tempA.N > 0){
        sharesVariables(f) = true;
        sharesVariables(f2) = true;
        break;
      }
    }
    // might be a single factor that contains all vars
    if(!sharesVariables(f)){
      FOR1D(vars, f2){
        if(factors(f)->variables.findValue(vars(f2)) < 0)
          break;
      }
      if(f2 == factors(f)->variables.N)
        sharesVariables(f) = true;
    }
    if(!sharesVariables(f))
      HALT("While building max spanning tree: factor node over " <<factors(f)->varIds <<" is disconnected from junction tree!");
  }
  
  
  rai::Array<MessagePair*> candidate_msg_pairs;
  boolA candidate_msg_pairs_contained;
  FOR1D(factors, f){
    for(f2 = f+1; f2 < factors.N; f2++){
      setSection(tempA, factors(f)->variables, factors(f2)->variables);
      if(tempA.N > 0){
        MessagePair* s = new MessagePair(factors(f), factors(f2));
        candidate_msg_pairs.append(s);
        candidate_msg_pairs_contained.append(1);
      }
    }
  }
  
  if(DEBUG >= 1){
    cout <<"candidate_msg_pairs [START]" <<endl;
    FOR1D(candidate_msg_pairs, f){
      cout <<*(candidate_msg_pairs(f)) <<endl;
    }
    cout <<"candidate_msg_pairs [END]" <<endl;
  }
  
  // define elementary clusters
  // factorID --> Cluster (made of factor IDs)
  rai::Array< uintA > clusters(factors.N);
  FOR1D(factors, f){clusters(f).append(f);}
  if(DEBUG >= 1){
    cout <<"Elementary clusters [START]" <<endl;
    FOR1D(factors, f){
      cout <<"(" <<f <<") infer::Factor over " <<factors(f)->varIds <<" connected to:" <<endl;
      FOR1D(clusters(f), f2){
        cout <<factors((clusters(f))(f2))->varIds <<"  ";
      }
      cout <<endl;
      // ids only:
      cout <<"In ids: #" <<f <<" -> " <<clusters(f) <<endl;
    }
    cout <<"Elementary clusters [END]" <<endl;
  }
  
  // THE ALGORITHM
  // big loop
  if(DEBUG >= 1){
    cout <<endl <<"---------------------------------------" <<endl;
    cout <<"THE ALGORITHM [START]" <<endl <<endl;
  }
  while(cliqueTree.messages.N < factors.N -1){
    // determine S with max weight
    maxId = 1000;
    max = 0;
    FOR1D(candidate_msg_pairs, f){
      if(candidate_msg_pairs_contained(f)
          && candidate_msg_pairs(f)->variables.N > max){
        maxId = f;
        max = candidate_msg_pairs(f)->variables.N;
      }
    }
    if(maxId == 1000){
      HALT("Constructing max spanning tree failed: No more MessagePair candidates with |MessagePair set| > 0 found!");
    }
    
    if(DEBUG >= 1){
      cout <<"++++++++++++++" <<endl;
      cout <<"Next MessagePair candidate:" <<endl;
      cout <<*(candidate_msg_pairs(maxId));
    }
    fac1 = candidate_msg_pairs(maxId)->f1;
    fac2 = candidate_msg_pairs(maxId)->f2;
    // index of fac1 in factors
    f = factors.findValue(fac1);
    // index of fac2 in factors
    f2 = factors.findValue(fac2);
    // Combine both trees
    if(numberSharedElements(clusters(f), clusters(f2)) == 0){
      if(DEBUG >= 1)
        cout <<" --> included" <<endl;
      cliqueTree.messages.append(candidate_msg_pairs(maxId));
      uint k;
      FOR1D(clusters(f), k){clusters(clusters(f)(k)).setAppend(clusters(f2));}
      FOR1D(clusters(f2), k){clusters(clusters(f2)(k)).setAppend(clusters(f));}
    }else{
      // delete if not needed
      delete candidate_msg_pairs(maxId);
      
      if(DEBUG >= 1)
        cout <<" --> discarded" <<endl;
    }
    
    if(DEBUG >= 1){
      cout <<"Elementary clusters [START]" <<endl;
      FOR1D(factors, f){
        cout <<"(" <<f <<") infer::Factor over " <<factors(f)->varIds <<" connected to:" <<endl;
        FOR1D(clusters(f), f2){
          cout <<factors((clusters(f))(f2))->varIds <<"  ";
        }
        cout <<endl;
        // ids only:
        cout <<"In ids: #" <<f <<" -> " <<clusters(f) <<endl;
      }
      cout <<"Elementary clusters [END]" <<endl;
    }
    
    candidate_msg_pairs_contained(maxId) = false;
  }
  
  if(DEBUG >= 1)
    cout <<"THE ALGORITHM [END]" <<endl <<endl;
    
  // delete remaining (unused) messages
  FOR1D(candidate_msg_pairs, f){
    if(candidate_msg_pairs_contained(f)){
      delete candidate_msg_pairs(f);
    }
  }
  
  // factors in clique tree uebertragen
  cliqueTree.F = factors;
  cliqueTree.setCliqueBeliefs(factors);
}




/**
Jordan, Chap. 17, p.22: "In a junction tree, local consistency implies global consistency."
*/

// bool JunctionTree::checkLocalConsistency(MessagePair& s, bool crashIfInconsistent){
//   uint DEBUG = 1;
//   if(DEBUG > 0)
//     cout <<"Consistency between " <<s.f1->varIds <<" and " <<s.f2->varIds <<" with MessagePair " <<s.varIds;
//   // test marginals over MessagePair variables on both factors
//   infer::Factor f1_marg(s.varIds), f2_marg(s.varIds);
//   tensorMarginal( f1_marg, *(s.f1), s.varIds);
//   tensorMarginal( f2_marg, *(s.f2), s.varIds);
//   if(DEBUG > 0){
//     cout <<endl <<"f1_marg:" <<endl <<f1_marg;
//     cout <<"f2_marg:" <<endl <<f2_marg;
//   }
//   if(f1_marg == f2_marg){
//     if(DEBUG > 0)
//       cout <<" -> consistent" <<endl;
//     // test that messages product is also the same
//     infer::Factor messagesProduct;
//     tensorProduct( messagesProduct, s.m12, s.m21);
//     if(f1_marg == messagesProduct){
//       if(DEBUG > 0)
//         cout <<" Also consistent with m12*m21 = " <<messagesProduct <<"." <<endl;
//       return true;
//     }
//     else {
//       if(DEBUG > 0)
//         cout <<" INCONSISTENT with m12*m21 = " <<messagesProduct <<"." <<endl;
//       HALT("Local inconsistency for "  <<s.f1->varIds <<" and " <<s.f2->varIds);
//       return false;
//     }
//   }
//   else {
//     if(DEBUG > 0)
//       cout <<" -> INCONSISTENT !!!" <<endl;
//     if(crashIfInconsistent){
//       HALT("Local inconsistency for "  <<s.f1->varIds <<" and " <<s.f2->varIds);
//     }
//     return false;
//   }
// }




// uses Marc's eq (12)
// void JunctionTree::jt_passMessage(MessagePair& s, infer::Factor& f_from, FactorGraph& jt){
//  uint DEBUG = 0;
//     if(DEBUG > 0){
//         cout <<"Pass message [START]" <<endl;
//         cout <<&s<<endl;
//     }
//
//  // message m21
//  if(&f_from == s.f2){
//         if(DEBUG > 0){cout <<f_from.varIds <<" --> " <<s.f1->varIds<<endl;}
//    if(DEBUG > 2){
//      cout <<"f_from: " <<endl <<f_from;
//             cout <<"f_to_old: " <<endl <<*s.f1;
//             cout <<"Message 2 --> 1" <<endl;
//    }
//
//         // Calculate new msg
//         jt.calcMessage_withDivision(f_from, *s.f1);
//    // Update target potential (~ incorporate message ~ rescaling)
//    tensorMultiply(*s.f1, s.m21);
//
//         if(DEBUG > 1) cout <<"New msg (m21): " <<endl <<s.m21;
//         if(DEBUG > 2) cout <<"f_to: " <<endl <<*s.f1;
//  }
//  // message m12
//  else {
//         if(DEBUG > 0){cout <<f_from.varIds <<" --> " <<s.f2->varIds<<endl;}
//    if(DEBUG > 2){
//             cout <<"f_from: " <<endl <<f_from <<endl;
//             cout <<"f_to_old: " <<endl <<*s.f2 <<endl;
//             cout <<"Message 1 --> 2" <<endl;
//    }
//
//         // Calculate new msg
//         jt.calcMessage_withDivision(f_from, *s.f2);
//    // Update target potential (~ incorporate message ~ rescaling)
//         tensorMultiply(*s.f2, s.m12);
//
//         if(DEBUG > 1) cout <<"New msg (m12): " <<endl <<s.m12;
//         if(DEBUG > 2) cout <<"f_to: " <<endl <<*s.f2;
//  }
//
//  if(DEBUG > 0){
//    cout <<"Pass message [END]" <<endl;
//  }
// }



/** Start message collecting from all children-subtrees (i.e., except its parent subtree). !*/
void recursiveCollectEvidence(infer::Factor* node, infer::Factor* parent, infer::FactorGraph& junctionTree){
  uint DEBUG = 0;
  if(DEBUG > 0)
    cout <<"Collect evidence for " <<node->varIds <<" [START]" <<endl;
  uint i;
  infer::FactorList neighbors;
  getNeighbors(node, neighbors);
  if(DEBUG > 0){
    cout <<" Neighbors: ";
    printFactorArrayViaVariableSets(neighbors);
  }
  //   cout <<junctionTree <<endl;
  FOR1D(neighbors, i){
    if(neighbors(i) == parent)
      continue;
    recursiveCollectEvidence(neighbors(i), node, junctionTree);
    // update with respect to this direction
    infer::MessagePair* s = getMessagePair(node, neighbors(i));
    //    cout <<"got MessagePair: " <<std::flush <<s <<endl;
    // message m21
    if(node == s->f1){
      if(DEBUG > 0) cout <<"Pass message " <<s->f2->varIds <<" --> " <<node->varIds <<endl;
      infer::Factor* b_to = junctionTree.B_c(get_list_index_unsigned__orig(junctionTree, s->f1));
      infer::passMessage(*(s->f2), *(s->f1), *b_to, NO_DIV);
    }
    // message 12
    else {
      if(DEBUG > 0) cout <<"Pass message " <<s->f1->varIds <<" --> " <<node->varIds <<endl;
      infer::Factor* b_to = junctionTree.B_c(get_list_index_unsigned__orig(junctionTree, s->f2));
      infer::passMessage(*(s->f1), *(s->f2), *b_to, NO_DIV);
    }
  }
  if(DEBUG > 0)
    cout <<"Collect evidence for " <<node->varIds <<" [END]" <<endl;
}


/** Starts message distributing on node into the direction of all its children (i.e., except its parent). */
void recursiveDistributeEvidence(infer::Factor* node, infer::Factor* parent, infer::FactorGraph& junctionTree){
  uint DEBUG = 0;
  if(DEBUG > 0)
    cout <<"Distribute evidence for " <<node->varIds <<" [START]" <<endl;
  uint i;
  infer::FactorList neighbors;
  getNeighbors(node, neighbors);
  FOR1D(neighbors, i){
    if(neighbors(i) == parent)
      continue;
    // update with respect to this direction
    infer::MessagePair* s = getMessagePair(node, neighbors(i));
    // message m12
    if(node == s->f1){
      // neighbors(i) == s->f2
      if(DEBUG > 0) cout <<"Pass message " <<node->varIds <<" --> " <<s->f2->varIds <<endl;
      infer::Factor* b_to = junctionTree.B_c(get_list_index_unsigned__orig(junctionTree, s->f2));
      passMessage(*(s->f1), *(s->f2), *b_to, NO_DIV);
    }
    // message 21
    else {
      // neighbors(i) == s->f1
      if(DEBUG > 0) cout <<"Pass message " <<node->varIds <<" --> " <<s->f1->varIds <<endl;
      infer::Factor* b_to = junctionTree.B_c(get_list_index_unsigned__orig(junctionTree, s->f1));
      infer::passMessage(*(s->f2), *(s->f1), *b_to, NO_DIV);
    }
    recursiveDistributeEvidence(neighbors(i), node, junctionTree);
  }
  if(DEBUG > 0)
    cout <<"Distribute evidence for " <<node->varIds <<" [END]" <<endl;
}



void infer::JunctionTree::collectAndDistributeInference(FactorGraph& junctionTree){
  // reset junction tree factors to original values
  junctionTree.resetMessages();
  junctionTree.resetCliqueBeliefs();
  
  uint DEBUG = 0;
  if(DEBUG > 0){
    cout <<"---------------------------" <<endl;
    cout <<"Propagating probabilities:" <<endl;
  }
  // determine root
  infer::Factor* p_root = junctionTree.F(0);
  if(DEBUG > 0){
    cout <<"Root: " <<p_root->varIds <<endl;
  }
  // collect evidence
  if(DEBUG > 0){
    cout <<"\n\n*******************************\n******************************\n**********************\n";
    cout <<"Collecting evidence:" <<endl;
  }
  recursiveCollectEvidence(p_root, 0, junctionTree);
  // distribute evidence
  if(DEBUG > 0){
    cout <<"\n\n*******************************\n******************************\n**********************\n";
    cout <<"Distributing evidence:" <<endl;
  }
  recursiveDistributeEvidence(p_root, 0, junctionTree);
  
  
  //   if(!checkConsistency(junctionTree, true)){
  //     HALT("Propagation of probabilities failed. (Local inconsistencies.)");
  //   }
}


void infer::JunctionTree::constructJunctionTree(FactorGraph& junctionTree, const infer::FactorList& factors, const infer::VariableList& vars){
  uint DEBUG = 0;
  if(DEBUG>0){cout <<"constructJunctionTree [START]" <<endl;}
  
  //construct cliques of triangulated tree
  infer::FactorList triangulatedCliques;
  buildTriangulatedCliques(factors, triangulatedCliques);
  
  if(DEBUG > 0){
    cout <<"Input factors:" <<endl;
    listWrite(factors, cout);
    uint f;
    cout <<"Factors built from triangulated graph:" <<endl;
    FOR1D(triangulatedCliques, f){
      cout <<"(" <<f <<") " <<triangulatedCliques(f)->varIds <<endl;
    }
    cout <<endl;
  }
  
  //build junction tree
  buildMaxSpanningTree(triangulatedCliques, vars, junctionTree);
  checkMessagePairConsistency(junctionTree.messages);
  if(DEBUG>0){cout <<"Message pairs consistent." <<endl;}
  if(DEBUG>0) junctionTree.write(cout);
  junctionTree.resetMessages();
  if(DEBUG > 0){
    cout <<"Junction tree before probability propagation:" <<endl;
    cout <<junctionTree <<endl;
    if(DEBUG > 1)
      junctionTree.writeMessagePairs(cout);
  }
  checkJunctionTreeProperty(junctionTree);
}


void infer::JunctionTree::junctionTreeInference(FactorGraph& junctionTree, const infer::FactorList& factors, const infer::VariableList& vars){
  constructJunctionTree(junctionTree, factors, vars);
  collectAndDistributeInference(junctionTree);
  
  uint DEBUG = 0;
  
  if(DEBUG > 0){
    cout <<"Junction tree after probability propagation:" <<endl;
    cout <<junctionTree <<endl;
    if(DEBUG > 1)
      junctionTree.writeMessagePairs(cout);
  }
  if(DEBUG>0){cout <<"constructJunctionTree [END]" <<endl;}
}








//===========================================================================
//
//     Loopy BP methods
//

void infer::LoopyBP_obsolete::loopy_belief_propagation(FactorGraph& fg, const infer::FactorList& factors){
  uint DEBUG = 0;
  constructBipartiteFactorGraph(fg, factors);
  uint MAX_STEPS = 20;
  double change;
  double STOP_TOL = 0.001;
  uint i;
  for(i=0; i<MAX_STEPS; i++){
    if(DEBUG>0){cout <<"+++ Round " <<(i+1) <<" +++" <<endl;}
    change = passAllEdges(fg, PARALLEL);
    if(change < STOP_TOL) break;
    if(DEBUG>0){fg.writeVariableBeliefs();}
    if(DEBUG>1) fg.write(cout);
  }
}


double infer::LoopyBP_obsolete::passAllEdges(FactorGraph& fg, PassType type){
  if(type == PARALLEL){
    return passAllEdges_parallel(fg);
  } else
    NIY;
  return 0;
}


double infer::LoopyBP_obsolete::passAllEdges_parallel(FactorGraph& fg){
  uint DEBUG =0;
  // (1) cliques --> vars
  uint i;
  if(DEBUG>0){cout <<"***** Propagating forth (factor -> variable) *****" <<endl;}
  FOR1D(fg.F, i){
    shoutMessages(*fg.F(i), NO_DIV);
  }
  if(DEBUG>0){cout <<"***** Propagating back (variable -> factor) *****" <<endl;}
  // (2) vars --> cliques
  FOR1D(fg.F_v, i){
    shoutMessages(*fg.F_v(i), NO_DIV);
  }
  double maxChange = fg.computeBeliefs();
  if(DEBUG>0){cout <<"Belief change:  maxChange=" <<maxChange <<endl;}
  return maxChange;
}

void infer::LoopyBP_obsolete::shoutMessages(infer::Factor& f, MsgCalc calcMsgType){
  uint i;
  FOR1D(f.messages, i){
    infer::Factor* f_to;
    if(f.messages(i)->f1 == &f)
      f_to = f.messages(i)->f2;
    else
      f_to = f.messages(i)->f1;
    if(calcMsgType == WITH_DIV)
      computeMessage_withDiv(f, *f_to);
    else
      computeMessage_noDiv(f, *f_to);
  }
}



void check_exactlyOneConditional(infer::VariableList& vars, infer::FactorList& facs){
  uint i, k;
  CHECK_EQ(vars.N,facs.N, "#vars != #facs");
  FOR1D(vars, i){
    FOR1D(facs, k){
      if(facs(k)->variables(0) == vars(i)) break;
    }
    CHECK(k<facs.N, "No factor for variable " <<vars(i)->id);
  }
}

void check_atLeastOneConditional(infer::VariableList& vars, infer::FactorList& facs){
  uint i, k;
  FOR1D(vars, i){
    FOR1D(facs, k){
      if(facs(k)->variables(0) == vars(i)) break;
    }
    CHECK(k<facs.N, "No factor for variable " <<vars(i)->id);
  }
}



// =======================================================================
//
//  Loopy BP - other approach (MT)
//

void infer::connectThemUp(infer::VariableList& V, infer::FactorList& F){
  RAI_MSG("you shouldn't use this anymore!!");
  for(Factor *f: F) checkConsistent(*f);
#if 0
  NIY;
  infer::Variable *v;
  infer::Factor *f;
  uint i, j;
  for_list(Variable, v,  V) v->factors.clear();
  for_list(Factor, f,  F) f->variables.clear();
  for_list(Factor, f,  F) for(j=0; j<f->varIds.N; j++){
    v=V(f->varIds(j));
    v->factors.append(f);
    f->variables.append(v);
  }
#endif
}

void infer::LoopyBP::clear(){
  vars.clear();
  facs.clear();
  listDelete(msgs);
}

void infer::LoopyBP::initBipartite(const infer::VariableList& _vars, const infer::FactorList& _facs){
  CHECK(!msgs.N, "delete list before");
  vars=_vars;
  facs=_facs;
  for(Factor *f: facs){
    for(uint j=0; j<f->variables.N; j++){
      msgs.append(new MessagePair(f, f->variables(j)));
    }
  }
}

void infer::LoopyBP::initPairwise(const infer::VariableList& _vars, const infer::FactorList& _facs){
  CHECK(!msgs.N, "delete list before");
  vars=_vars;
  facs=_facs;
  for(Factor *f: facs){
    CHECK(f->variables.N<=2, "only for pair-wise networks!");
    if(f->variables.N==1) msgs.append(new MessagePair(f, f->variables(0)));
    if(f->variables.N==2) msgs.append(new MessagePair(f->variables(0), f->variables(1), f));
  }
}

void infer::LoopyBP::getVarBeliefs(rai::Array<infer::Factor>& beliefs, bool normalized){
  uint i;
  beliefs.resize(vars.N);
  for(i=0; i<vars.N; i++){
    collectBelief(beliefs(i), vars(i), NULL);
    if(normalized) beliefs(i).normalize();
  }
}

void infer::LoopyBP::getVarBelief(infer::Factor& belief, infer::Variable *v, bool normalized){
  collectBelief(belief, v, NULL);
  if(normalized) belief.normalize();
}

void infer::LoopyBP::step(){
  recomputeBatchOfMessages(msgs, false);
  recomputeBatchOfMessages(msgs, true);
}

infer::LoopyBP::~LoopyBP(){
  listDelete(msgs);
}

void loopyBP_bipartite(const infer::VariableList& vars, const infer::FactorList& facs, uint T){
  infer::LoopyBP lbp;
  lbp.initBipartite(vars, facs);
  
  rai::Array<infer::Factor> beliefs(vars.N);
  uint t;
  for(t=0; t<T; t++){
    lbp.getVarBeliefs(beliefs);
    beliefs.write(cout, "\n");
    lbp.step();
  }
}

void loopyBP_pairwise(const infer::VariableList& vars, const infer::FactorList& facs, uint T){
  infer::LoopyBP lbp;
  lbp.initPairwise(vars, facs);
  
  rai::Array<infer::Factor> beliefs(vars.N);
  uint t;
  for(t=0; t<T; t++){
    lbp.getVarBeliefs(beliefs);
    beliefs.write(cout, "\n");
    lbp.step();
  }
}

// =======================================================================
//
//  mean field
//

void meanField_collectBeliefs(arr& beliefs, const infer::VariableList& vars){
  HALT("that's broke");
#if 0
  if(beliefs.N!= vars.N){
    beliefs.resize(vars.N);
    beliefs = .5;
  }
  
  uint i, j;
  infer::Variable *v;
  infer::Factor *f;
  arr b(2);
  for_list(infer::Variable, v,  vars){
    b.setZero();
    for_list(infer::Factor, f,  v->factors){
      if(f->variables(0)==v){
        b += f->P * ARR(1.-beliefs(f->variables(1)), beliefs(f->variables(1)));
      }else{
        b +=ARR(1.-beliefs(f->variables(0)), beliefs(f->variables(0))) * f->P;
      }
    }
    beliefs(i) = 1./(1.+exp(-b(1)));
  }
#endif
}

void infer::LoopyBP::step_meanfield(){
  arr beliefs;
  meanField_collectBeliefs(beliefs, vars);
  cout <<beliefs <<endl;
}

//=======================================================================
//
//  inference on trees
//

std::ostream& operator<<(std::ostream& os, const TreeNode& t){
  return os <<"par=" <<t.parent <<" dim=" <<t.dim <<" P=" <<t.P <<endl;
}

void randomTree(Tree& tree, uint N, uint K, uint roots){
  tree.resize(N);
  uint n;
  for(n=0; n<N; n++){
    tree(n).dim=K;                 //K-dim variable
    if(n<roots) tree(n).parent=-1;
    else        tree(n).parent=rnd(n);  //random parent
    if(n<roots) tree(n).P.resize(K); //K-times-K joint factor (this, parent)
    else        tree(n).P.resize(K, K);
    RAI_MSG("fully random roots!");
    /*if(n<roots){ // deterministic roots for testing purposes
      tree(n).P(0) = 0.0;
      tree(n).P(1) = 1.0;
    } else*/
    rndUniform(tree(n).P, .1, 1., false);
  }
}


/* conversion: every tree node corresponds to
- a variable
- a pairwise factor connecting it to its parent
- and a Message connecting its factor to the parent's factor

everything is ordered precisely */

void tree2FactorGraph(infer::FactorGraph& fg, const rai::Array<TreeNode>& tree){
  uint i, N=tree.N;
  fg.V.resize(N);
  fg.F.resize(N);
  for(i=0; i<N; i++){  //variables
    fg.V(i) = new infer::Variable(tree(i).dim, STRING("tree_node_" <<std::setfill('0') <<std::setw(3) <<i), i);
  }
  for(i=0; i<N; i++){  //factors
    if(tree(i).parent<0) fg.F(i) = new infer::Factor({fg.V(i)});
    else                 fg.F(i) = new infer::Factor({fg.V(i), fg.V(tree(i).parent)});
    fg.F(i)->setP(tree(i).P);
  }
  for(i=0; i<N; i++){  //messages
    if(tree(i).parent>=0) fg.messages.append(new infer::MessagePair(fg.F(i), fg.F(tree(i).parent)));
  }
  fg.resetMessages();
  fg.resetCliqueBeliefs();
}

// void getPosteriors(rai::Array<arr>& posteriors, FactorGraph& model){
//   uint i;
//   FOR1D(model.V, i){
//     uintA var_ids;
//     var_ids.append(model.V(i)->id);
//     infer::Factor marginal;
//     getMarginal(marginal, var_ids, model);
//     posteriors.append(marginal.P);
//   }
// }

void treeInference(rai::Array<arr>& posteriors, const Tree& tree){
  infer::FactorGraph model;
  tree2FactorGraph(model, tree);
  treeInference(model.F(0), true);
  infer::getVariableBeliefs(posteriors, model.V);
}

void treeInference(rai::Array<arr>& posteriors, const Tree& forest, uintA& roots){
  infer::FactorGraph model;
  tree2FactorGraph(model, forest);
  for(uint i=0; i<roots.N; i++){
    treeInference(model.F(roots(i)), true);
  }
//   model.write(cout, true, true);
//   model.writeMessagePairs(cout);
  infer::getVariableBeliefs(posteriors, model.V);
}


void write(Tree& tree){
  uint i;
  cout <<"tree: ";
  FOR1D(tree, i){
    cout <<tree(i).parent <<" ";
  }
  cout <<endl;
}


// =======================================================================
//
//  Exact Tree Inference (by MT for NP)
//

/* message orders allows us to compute a series of messages
   in a specified order and direction. A message order is stored
   as a list of messages and a directionality with each of them. */


/// compute a series of messages in a given order, optionally inversely (also inverting each message's directionality)
void infer::recomputeBatchOfMessages(MessagePairList &msgs, bool invert_order){
  uint i, N=msgs.N;
  if(!invert_order) for(i=0; i<N; i++) recomputeMessage_12(*msgs(i));
  else              for(i=N; i--;)    recomputeMessage_21(*msgs(i));
}

/// compute a series of messages in a given order, optionally inversely (also inverting each message's directionality)
void infer::recomputeBatchOfMessages(MessagePairList &msgs, const boolA &msgFlips, bool invert_order){
  uint i, N=msgs.N;
  if(!invert_order) for(i=0; i<N; i++){
      if(!msgFlips(i)) recomputeMessage_12(*msgs(i));
      else             recomputeMessage_21(*msgs(i));
    } else for(i=N; i--;){
      if(!msgFlips(i)) recomputeMessage_21(*msgs(i));
      else             recomputeMessage_12(*msgs(i));
    }
}

/// for inference in a tree, we first construct an ordering descending from the root
void infer::constructTreeMessageOrder(MessagePairList& msgs, boolA &msgFlips, const infer::Factor *root){
  uint i, j;
  infer::Factor *f;
  MessagePair *m;
  
  //--- construct a tree oder
  //first ad all neighbor links of root to the order
  for(j=0; j<root->messages.N; j++){
    m=root->messages(j);
    msgs.append(m);
    if(m->f1==root) msgFlips.append(false);
    else            msgFlips.append(true);
  }
  //step through the existing links and append sublinks dynamically
  for(i=0; i<msgs.N; i++){
    if(!msgFlips(i)) f=msgs(i)->f2; else f=msgs(i)->f1; //f is the sub-factor of order(i).m
    for(j=0; j<f->messages.N; j++){
      m=f->messages(j);                                     //m is one sub-message of f
      if(m==msgs(i)) continue;                     //discard it when it is going up
      msgs.append(m);
      if(m->f1==f) msgFlips.append(false);
      else         msgFlips.append(true);
    }
  }
}

/// inference on a tree is then trivial: compute a descending-from-root order,
/// pass all messages backward (collecting towards the root)
/// pass all messages forward  (distributing from the root)
void infer::treeInference(const infer::Factor *root, bool check_consitency){
  MessagePairList msgs;
  boolA msgFlips;
  
  constructTreeMessageOrder(msgs, msgFlips, root);
  
  //cout <<order <<endl;
  
  recomputeBatchOfMessages(msgs, msgFlips, true);
  recomputeBatchOfMessages(msgs, msgFlips, false);
  
  if(check_consitency){
    for(uint i=0; i<msgs.N; i++) checkConsistency(*msgs(i));
  }
  //RAI_MSG("everything consistent :-)");
}


// =======================================================================
//
//  inference for mixture length DBNs
//

void infer::inferMixLengthUnstructured(
  arr& alpha, arr& beta, arr& PT, double& PR, double& ET,
  const arr& S, const arr& R, const arr& P, double gamma, uint Tmax,
  bool updateMode){
  DEBUG_INFER(1, cout <<RAI_HERE <<endl);
  arr a, b;
  double gt, gSum;
  uint t;
  CHECK(S.nd==1 && R.nd==1 && P.nd==2 && S.N==R.N && S.N==P.d0 && S.N==P.d1, "");
  CHECK(gamma>0. && gamma<=1., "");
  if(!updateMode){
    a = S;
    b = R;
    alpha = a;
    beta  = b * gamma;
    PT.resize(2*Tmax+1);
    PT(0) = scalarProduct(a, b);
    gt = gamma;
    for(t=1; t<=Tmax; t++){
      a = P*a;
      PT(2*t-1) = scalarProduct(a, b);
      b = b*P;
      PT(2*t) = scalarProduct(a, b);
      alpha += gt*a;
      gt *= gamma;
      beta += gt*b;
    }
    PR = 0.;
    ET = 0.;
    gt = 1.;
    gSum=0.;
    for(t=0; t<=2*Tmax; t++){
      PT(t) *= gt;
      ET += t*PT(t);
      PR += PT(t);
      gSum += gt;
      gt *= gamma;
    }
    PT /= PR;
    ET /= PR;
#ifndef RAI_NOCHECK
    //if the iteration above is only to Tmax instead of 2.*Tmax, the following check is true:
    //CHECK(fabs(scalarProduct(alpha, R)-PR)<1e-3, "");
    double dummy1 = scalarProduct(S, beta);
    double dummy2 = scalarProduct(alpha, R);
    CHECK(fabs(dummy1-gamma*dummy2)<1e-6, "");
#endif
  }else{
    CHECK(alpha.nd==1 && beta.nd==1 && alpha.N==S.N && beta.N==R.N, "");
    //alpha /= (1.-gamma); //inverse of operation on last return
    for(t=0; t<=Tmax; t++){
      alpha *= gamma;  alpha = P*alpha;  alpha += S;
      beta  = beta*P;  beta += R;        beta  *= gamma;
    }
    PT.resize(0);
    ET=-1;
    PR = scalarProduct(alpha, R);
    /*
    double dummy1 = scalarProduct(S, beta);
    double dummy2 = scalarProduct(alpha, R);
    CHECK(fabs(dummy1-gamma*dummy2)<1e-3, ""); //check holds only approx, because old alpha & beta are w.r.t old params!!
    */
  }
  //alpha *= (1.-gamma);
  PR    *= (1.-gamma);
}

void infer::inferMixLengthStructured(
  infer::Factor& alpha, infer::Factor& beta, arr& PT, double& PR, double& ET,
  const infer::VariableList& headVars, const infer::VariableList& tailVars,
  const infer::FactorList& S, const infer::FactorList& R, const infer::FactorList& P, double gamma, uint Tmax,
  bool updateMode){
  
  DEBUG_INFER(1, cout <<RAI_HERE <<endl);
  infer::Factor a, b, Shead, Rtail;
  double gt, gSum;
  uint t;
  DEBUG_INFER(2, cout <<"  headVars=" <<headVars <<"  tailVars=" <<tailVars <<endl);
  CHECK_EQ(headVars.N,tailVars.N, ""); //actually should also check that their dims are equal...
  CHECK(gamma>0. && gamma<=1., "");
  if(!updateMode){
    //get initial a and b:
    eliminationAlgorithm(a, S, headVars);
    //infer::FactorList tmp = P;  tmp.append(R);
    eliminationAlgorithm(b, R, headVars);  b.relinkTo(tailVars); //b.variables=tailVars;
    //above, we reassociate the factor to the tail variables although it
    //was initially defined over the head variables
    //we do this same ``trick'' (or hack) also below
    alpha = a;
    beta  = b;  beta.P *= gamma;
    PT.resize(2*Tmax+1);
    PT(0) = scalarProduct(a.P, b.P)*::exp(a.logP+b.logP);
    gt = gamma;
    infer::FactorList fwdList = P;  fwdList.append(&a);
    infer::FactorList bwdList;      bwdList.append(&b);  bwdList.append(P);
    for(t=1; t<=Tmax; t++){
      eliminationAlgorithm(a, fwdList, tailVars);  a.relinkTo(headVars); //a.variables=headVars;
      PT(2*t-1) = scalarProduct(a.P, b.P)*::exp(a.logP+b.logP);
      eliminationAlgorithm(b, bwdList, headVars);  b.relinkTo(tailVars); //b.variables=tailVars;
      PT(2*t)   = scalarProduct(a.P, b.P)*::exp(a.logP+b.logP);
      tensorWeightedAdd(alpha, gt, a);
      gt *= gamma;
      tensorWeightedAdd(beta , gt, b);
    }
    PR = 0.;
    ET = 0.;
    gt = 1.;
    gSum=0.;
    for(t=0; t<=2*Tmax; t++){
      PT(t) *= gt;
      ET += t*PT(t);
      PR += PT(t);
      gSum += gt;
      gt *= gamma;
    }
    PT /= PR;
    ET /= PR;
  }else{
    CHECK(alpha.variables==headVars && beta.variables==tailVars, "");
    eliminationAlgorithm(Shead, S, headVars);
    eliminationAlgorithm(Rtail, R, headVars);  Rtail.relinkTo(tailVars);  //Rtail.variables=tailVars;
    infer::FactorList fwdList = P;  fwdList.append(&alpha);
    infer::FactorList bwdList;      bwdList.append(&beta);   bwdList.append(P);
    for(t=0; t<=Tmax; t++){
      alpha.P *= gamma;
      //RAI_MSG("does that work??");
      eliminationAlgorithm(alpha, fwdList, tailVars);  alpha.relinkTo(headVars); //alpha.variables=headVars;
      tensorWeightedAdd(alpha, 1., Shead);
      
      eliminationAlgorithm(beta, bwdList, headVars);   beta.relinkTo(tailVars);  //beta.variables=tailVars;
      tensorWeightedAdd(beta, 1., Rtail);
      beta.P *= gamma;
    }
    PT.resize(0);
    ET=-1;
    PR = scalarProduct(alpha.P, Rtail.P)*::exp(alpha.logP+Rtail.logP);
  }
  //alpha.logP += ::log(1.-gamma);
  PR    *= (1.-gamma);
}
