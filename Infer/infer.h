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
/* Copyright (C) 2000, 2006  Marc Toussaint (mtoussai@inf.ed.ac.uk)
under the terms of the GNU LGPL (http://www.gnu.org/copyleft/lesser.html)
see the `util.h' file for a full copyright statement  */

#ifndef RAI_inference_h
#define RAI_inference_h

#include <Core/array.h>
#include <Core/util.h>
#include <Core/graph.h>
#include <map>
//#include <set>



//===========================================================================
//
// basic data structures: variables, factors, message pairs
//

/* Note: since variables, factors and message pairs know how they are
   linked (they store lists of adjoint objects), these types are
   already sufficient to define arbitrary models. The factor graph
   type below does nothing more than storing all variables, factors
   and messages together and simplify building connected models
   automatically */

namespace infer {
struct Variable;
struct Factor;
struct MessagePair;

//simple ``list'' types
typedef rai::Array<Variable*> VariableList;
typedef rai::Array<Factor*>   FactorList;
typedef rai::Array<MessagePair*> MessagePairList;

enum SpecialFactorType { NONE=0, AND, OR, WTA };

extern uint VarCount;
}

namespace infer {
/** a discrete random variable */
struct Variable {
  //core defining properties
  uint id;          ///< unique identifyer
  uint dim;         ///< cardinality
  
  // auxilliary & connectivity
  rai::String name;  ///< up to you...
  infer::FactorList factors;    ///< each variable knows all factors it is part of
  MessagePairList messages;  ///< each variable knows all the messages it directly connects to
  rai::Graph ats; //any convenience information (e.g. for dot);
  
  Variable();
  Variable(uint _dim, const char *_name);
  Variable(uint _dim, const char *_name, uint _id);
  ~Variable();
  operator uint() const { return id; } //TODO: REMOVE
  void write(ostream& os = std::cout) const;
};
}

namespace infer {
/** a factor (=probability table) over a tuple of variable; a list of
    factors readily defines a proper factored joint distribution

    Probability tables are multi-dimensional arrays (tensors), per
    variable a separate dimension (boolean variables: first false,
    then true). */
struct Factor {
  //core defining properties
  uintA varIds; ///< f=f(x_1, x_3, x_7) => id=[1, 3, 7]; array of variables (their id's) this factor depends on
  uintA dim;    ///< f=f(x_1, x_3, x_7) => dim=[dim(x_1), dim(x_3), dim(x_7)];  array of dimensionalities of the variables this factor depends on
  arr P;        ///< the (probability) table
  double logP;  ///< the log-scaling of the table, such that true_factor = exp(logP) * P
  
  // auxilliary & connectivity
  rai::String name;
  SpecialFactorType specialType;
  VariableList variables;
  FactorList factors;
  MessagePairList messages;          ///< each factor knows all the msg_pairs it connects to
  rai::Graph ats; //any convenience information (e.g. for dot);
  
  Factor();
  ~Factor();
  Factor(const VariableList& variables, const char *_name=NULL);
  Factor(const VariableList& variables, const arr& q, const char *_name=NULL);
  void init(const VariableList& variables);
  void relinkTo(const VariableList& variables);
  void operator=(const Factor& q);
  void setP(const arr& q);           ///< f(x) = q(x)
  void setText(const char* text);    ///< f(x) = q(x)
  void setOne();                     ///< f(x) = 1
  void setUniform();                 ///< sum(x_i) f(x_i) = 1  only if factor.P is set by hand (as a conditional)
  void setRandom();                  ///< randomize f (e.g., for automatic testing)
  void setEvidence(uint e);          ///< f(e) = 1,  f(x!=e) = 0
  bool operator==(Factor& q);        ///< check if f==q (also checks for logP)
  void getP(arr& p) const;           ///< p = exp(logP)*P
  void normalize(){ lognormScale(P, logP); logP=0.; }
  void write(std::ostream& os = std::cout, bool brief=false) const;
  void writeNice(std::ostream& os = std::cout) const;
  void writeExtremelyNice(std::ostream& os = std::cout) const;
  
  void makeLogZero(){ P *= ::exp(logP); logP=0.; }
  
  double entry(uint i);              ///< returns logP * P.elem(i)
  void checkCondNormalization(uint left=1, double tol=1e-10);
  uint numberNonZeroEntries();
};
}

namespace infer {
/** a pair of messages (fwd and bwd) which connects two
    factors. Depending on the context this can be a separator (JTA), a
    link from factor to a single-variable-factor (bi-partite factor
    graph), or a links between arbitrary factors (loopy BP) */
struct MessagePair {
  //core defining properties
  Factor m12, m21;      ///< the forward and backward message
  Factor   *f1, *f2;    ///< the first and second factor it is attached to (if at all..)
  Variable *v1, *v2;    ///< the first and second variable it is attached to (if at all..)
  Factor   *v_to_v_fac; ///< in case of variable-to-variable message: the factor it is one-to-one associated with
  VariableList variables;       ///< the variables the messages are defined over
  
  MessagePair();
  MessagePair(Factor *_f1, Factor *_f2);  //factor-to-factor message (e.g., separate in JTA)
  MessagePair(Variable *_v1, Variable *_v2, Factor *_v_to_v_fac); //var-to-var message (pair-wise net)
  MessagePair(Factor *_f1, Variable *_v2); //factor-to-var message (e.g., bi-partite graph)
  ~MessagePair();
  void init(Factor *_f1, Factor *_f2);
  void init(Variable *_v1, Variable *_v2, Factor *_v_to_v_fac);
  void init(Factor *_f1, Variable *_v2);
  void operator=(const MessagePair& s){ variables=s.variables; f1=s.f1; f2=s.f2; m12=s.m12; m21=s.m21; }
  void write(std::ostream& os) const;
  void writeIds(std::ostream& os) const;
};
}


//===========================================================================
//
// basic operations on variables, factors, & messages
//

namespace infer {
//-- operations on factors
void tensorProduct(Factor& c, const Factor& a, const Factor& b);
//[[specify remainingVars and guarantee that the output has right order!]]
void tensorProductMarginal(Factor& c, const Factor& a, const Factor& b, const VariableList& eliminateVars);
void tensorMarginal(Factor& m, const Factor& f, const VariableList& marginalVars);    //marginalVars==remaining Vars
void tensorMaxMarginal(Factor& m, const Factor& f, const VariableList& marginalVars);
void tensorMultiply(Factor& f, const Factor& m);
void tensorDivide(Factor& f, const Factor& m);
void tensorAdd(Factor& f, const Factor& m);
void tensorInvertMultiply(Factor& f, const Factor& m);
void tensorWeightedAdd(Factor& f, double w, const Factor& m);
void checkConsistent(const Factor &f);

//-- helpers for variables lists
inline uintA ids(const VariableList& vars){ uintA id(vars.N); for(uint i=0; i<id.N; i++) id(i)=vars(i)->id; return id; }

//-- operations on pure factor lists
void getJoint(Factor& joint, const FactorList& factors);
void computeEliminationOrder(VariableList& elimOrder, const FactorList& factors, const VariableList& elimVars);
void eliminateVariable(FactorList& factors, FactorList& newed_factors, Variable *var);
void eliminationAlgorithm(Factor& post, const FactorList& factors, const VariableList& remaining_vars);
void checkConsistent(const FactorList& F);

//-- operations on single messages
void collectBelief(Factor& belief, const Factor& f, const MessagePair *exclude);
void collectBelief(Factor& belief, Variable *v, const MessagePair *exclude);
void recomputeMessage_12(MessagePair& sep);
void recomputeMessage_21(MessagePair& sep);
void recomputeBatchOfMessages(MessagePairList& msgs, bool invert_order=false);
void recomputeBatchOfMessages(MessagePairList &msgs, const boolA &msgFlips, bool invert_order);
bool checkConsistency(const MessagePair& sep);
bool checkConsistencyBatch(const MessagePairList& msgs);

//-- operations on structures (after factors have been linked to messages)
void constructTreeMessageOrder(MessagePairList& msgs, boolA &msgFlips, const Factor *root);
void treeInference(const Factor *root, bool checkConsistency);

//-- LoopyBP engine
struct LoopyBP {
  ~LoopyBP();
  MessagePairList msgs;
  VariableList    vars;
  FactorList      facs;
  void clear();
  void initBipartite(const VariableList& vars, const FactorList& facs);
  void initPairwise(const VariableList& vars, const FactorList& facs);
  void getVarBeliefs(rai::Array<Factor>& beliefs, bool normalized=true);
  void getVarBelief(Factor& belief, Variable *v, bool normalized=true);
  void step();
  void step_meanfield();
  //void loopyBP_pairwise(const VariableList& vars, const FactorList& facs);
  //void loopyBP_bipartite(const VariableList& vars, const FactorList& facs);
};
void connectThemUp(VariableList& V, FactorList& F);
void getVariableBeliefs(rai::Array<arr>& post, const VariableList& vars);

//-- pipes
//inline ostream& operator<<(ostream& os, const iSpace& s)      { s.write(os); return os; }
inline ostream& operator<<(ostream& os, const Variable& v)    { v.write(os); return os; }
inline ostream& operator<<(ostream& os, const Factor& f)      { f.write(os); return os; }
inline ostream& operator<<(ostream& os, const MessagePair& s) { s.write(os); return os; }
inline ostream& operator<<(ostream& os, const VariableList& list)   { listWrite(list, os, "\n"); return os; }
inline ostream& operator<<(ostream& os, const FactorList& list)     { listWrite(list, os, "\n"); return os; }
inline ostream& operator<<(ostream& os, const MessagePairList& list){ listWrite(list, os, "\n"); return os; }
}

// =======================================================================
//
//  inference for mixture length DBNs
//

namespace infer {
void inferMixLengthUnstructured(
  arr& alpha, arr& beta, arr& PT, double& PR, double& ET,
  const arr& S, const arr& R, const arr& P, double gamma, uint Tmax,
  bool updateMode=false);
  
void inferMixLengthStructured(
  Factor& alpha, Factor& beta, arr& PT, double& PR, double& ET,
  const VariableList& headVars, const VariableList& tailVars,
  const FactorList& S, const FactorList& R, const FactorList& P, double gamma, uint Tmax,
  bool updateMode=false);
}


//=======================================================================
//
//  inference on trees
//

struct TreeNode {
  int parent;
  uint dim;
  arr P;
};
typedef rai::Array<TreeNode> Tree;

void write(Tree& tree);
void treeInference(rai::Array<arr>& posteriors, const Tree& tree);
void treeInference(rai::Array<arr>& posteriors, const Tree& forest, uintA& roots);
void randomTree(Tree& tree, uint N, uint K, uint roots=1);
std::ostream& operator<<(std::ostream& os, const TreeNode& t);


//===========================================================================
//
// Factor Graph (obsolete!!) (the new way is simply a list of variables and factors! (see LoopyBP example above))
//
namespace infer {
struct FactorGraph {
  VariableList V;
  FactorList F;   // original factors (over cliques; kept constant)
  FactorList F_v; // dummies for variable factors; only used in case of true factor graph; all 1
  MessagePairList messages;  // point to F / F_v
  FactorList B_c; // beliefs over clique factors (only for saving; optional)
  FactorList B_v; // beliefs over variable factors (only for saving; optional)
  
  ~FactorGraph(){ deleteAll(); }
  FactorGraph& operator=(const FactorGraph& M){
    RAI_MSG("das kopiert nur die pointer, erzeugt keinen neuen Factor graphen!");
    B_c=M.B_c;
    B_v=M.B_v;
    messages=M.messages;
    F=M.F;
    F_v=M.F_v;
    return *this;
  }
  
  // deletes beliefs and message pairs
  void deleteAll();
  
  void setCliqueBeliefs(const FactorList& fs_orig);
  void resetCliqueBeliefs(){setCliqueBeliefs(F);} // B_c
  void resetMessages();  // msg_pairs
  void resetVariableFactors();  // F_v
  void resetVariableBeliefs();  // B_v
  
  void checkCondNormalization_B_c(double tol=1e-10);
  
  double computeBeliefs();
  
  Factor* getBelief(Factor* f_orig);
  
  // deprecated -- nicht mehr benutzt, da auf den beliefs nicht mehr gerechnet wird
//   void checkFaithfulness(); // for whole graph; based on Marc's eq (2) / (3)

  void write(std::ostream& os = cout, bool writeBeliefs = true, bool writeMessages = true) const;
  void writeNice(std::ostream& os = cout) const;
  void writeMessagePairs(std::ostream& os) const;
  void writeVariableBeliefs(std::ostream& os = cout) const;
  
  // EFFICIENCY helpers
  
  // in case lookup for beliefs needs to be fast
  std::map<Factor*, Factor*> F2B;
  Factor* getBelief_fast(Factor* f_orig);
  void setF2Bmap();
  void addF2Bmap(Factor* f, Factor* b);
  
  // factors where variable is first (i.e., the conditioned variable)
  std::map<uint, FactorList > V2F;
  std::map<uint, Factor*> V2Fv;
  void setV2F();
  void addV2Fmap(Factor* f);
  void setV2Fv();
};
stdOutPipe(FactorGraph);
}


//===========================================================================
//
// some helper functions
//

void write(const infer::FactorList& facs, ostream& os = cout);
void writeNice(const infer::VariableList& vars);
void writeNice(const infer::FactorList& individual_factors);
void writeEdges(const infer::FactorList& individual_factors, bool withProbs = false);
void writeExtremelyNice(const infer::FactorList& facs);


// calculates marginal for given variables
void getMarginal(infer::Factor& marginal, const uintA& marginalVars, infer::FactorGraph& fg);
void getMarginal(infer::Factor& marginal, const infer::VariableList& marginalVars, infer::FactorGraph& fg);



//===========================================================================
//
// Junction Tree methods
//

namespace infer {
namespace JunctionTree {
/** triangulates graph based on factors; ensures that resulting factors are max cliques;
corresponds to UNDIRECTED_GRAPH_ELIMINATE Jordan, Chapter 3, p. 13 */
void buildTriangulatedCliques(const FactorList& factors, FactorList& triangulatedCliques);

/** Builds max spanning tree (weights of edges according to size of set of MessagePair variables */
void buildMaxSpanningTree(FactorList& factors, const VariableList& vars, FactorGraph& cliqueTree);

// main method
/** Constructs a junction tree: triangulates original graph, builds max spanning tree
and updates probabilities &*/
void constructJunctionTree(FactorGraph& junctionTree, const FactorList& factors, const VariableList& vars);

/** Update prob dist on graph by passing messages */
void collectAndDistributeInference(FactorGraph& junctionTree);

void junctionTreeInference(FactorGraph& junctionTree, const FactorList& factors, const VariableList& vars);

void checkJunctionTreeProperty(FactorGraph& junctionTree);

/** adds evidence node to graph !*/
void addEvidence(FactorGraph& junctionTree, Factor& evid);
}
}


//===========================================================================
//
//   Loopy BP
//

enum MsgCalc { WITH_DIV, NO_DIV };

namespace infer {
namespace LoopyBP_obsolete {

// with variable factors
void constructBipartiteFactorGraph(FactorGraph& fg, const FactorList& factors);

enum PassType { PARALLEL };

double passAllEdges_parallel(FactorGraph& fg);
double passAllEdges(FactorGraph& fg, PassType type);
/** computes all outgoing messages of belief */
void shoutMessages(Factor& f, MsgCalc calcMsgType = NO_DIV);

// main method
void loopy_belief_propagation(FactorGraph& fg, const FactorList& factors);
}
}










#if 1

// OLD STUFF COPIED FROM MARC'S OLD CODE
/** automatically allocate message factors for mu_fwd and mu_bwd at each edge with the
correct factor types (derived from the nodes) \ingroup infer2 */
void allocateEdgeFactors(infer::FactorGraph& G);
/** test the factor graph for equilibrium (which means converged inference) \ingroup infer2 */
bool checkEquilibrium(infer::FactorGraph& G);
void getIndices(uintA& list, const uintA& id, const uintA& dim, const uintA& mid);
void getPick(uintA& pick, const uintA& id, const uintA& mid);
void fb(infer::FactorGraph& G);
/** calls message passing for a certain selection of edges.
seq gives the indices of the edges in the first time slice.
Those indices are extrapolated over a whole DBN of time length T with Mod
edges in each time slice. A negative indes means that the esge passes backward.
dir=BACKWARD means that one starts with the last time slice of the DBN going
towards the first. \ingroup infer2 */
//void passCertainEdges(infer::FactorGraph& G, intA seq, uint T, uint Mod, FwdBwd dir);
/** \ingroup infer2 */
//void passLabeledEdges(infer::FactorGraph& G, int label, FwdBwd dir);
void clearLabeledEdges(infer::FactorGraph& G, int label);
void resetLabeledNodes(infer::FactorGraph& G, const char *name);
void resetCertainNodes(infer::FactorGraph& G, intA seq);
void resetAllNodes(infer::FactorGraph& G);
void uniformCertainNodes(infer::FactorGraph& G, intA seq, uint T, uint Mod);
/// \ingroup infer2


//////////////////////////////////////////////

// For each var there is exactly one factor with the variable as first var.
void check_exactlyOneConditional(infer::VariableList& vars, infer::FactorList& facs);
void check_atLeastOneConditional(infer::VariableList& vars, infer::FactorList& facs);

void write(const rai::Array< infer::Variable* > vars, ostream& os = cout);

//[mt] only checks the registries - not the factors
void checkMessagePairConsistency(infer::MessagePairList msg_pairs);
infer::MessagePair* getMessagePair(const infer::Factor* f1, const infer::Factor* f2);
infer::Factor* getMessage(const infer::Factor* f_from, const infer::Factor* f_to);
void writeMessage(const infer::Factor* f_from, const infer::Factor* f_to, ostream& os = cout);
void sample(infer::Factor& f, uintA& samples);///< sample (s=discrete, x=continuous)

// Accessing global objects [mt] shouldn't be visible anymore
//iSpace* get_global_space();
void get_global_vars(infer::VariableList& vars);
void print_global_vars(uintA ids);


//===========================================================================
// Algorithms on factors and message pairs and graphs

namespace infer {
// 2 types of calculating messages
/** calculate new message */
void computeMessage_withDiv(Factor& f_from, Factor& f_to); // based on Marc's eq (5) / (7)
void computeMessage_noDiv(Factor& f_from, Factor& f_to); // based on Marc's eq (4) / (6)
// using already calculated belief of f_from
void computeMessage_noDiv(Factor& f_from, Factor& b_from, Factor& f_to);

/** computes all incoming messages of belief */
void askForMessages(Factor& f, MsgCalc calcMsgType = NO_DIV); //[mt] similar to collectBelief?
// BELIEF BASED
// --> using belief factors for storage

/** Updates belief according to original factor and incoming messages */
void collectBelief(Factor& belief, const Factor& f, const MessagePair *exclude);

/** pass message from f_from to f_to and write it into b_to*/
// based on Marc's eq (2) / (3)
// if calcMsgType==with_division, the incoming msgs to f_from are not used to calc the message.
// --> important if we do a mixture of belief propagation and setting certain factors inbetween by hand
double passMessage(Factor& f_from, Factor& f_to, Factor& b_to, MsgCalc calcMsgType);

/** computes messages to all neighbors and updates these accordingly; should only
be used in case of loopy BP (for efficiency reasons) */
// NIY
// double distributeMessages(FactorGraph& fg, Factor& f, MsgCalc calcMsgType = NO_DIV);

/** Calculates posterior over the variables given in "post" using the elimination algorithm. */
//void posteriorByElimination(FactorList& factors, Factor& post);
}
#endif

//===========================================================================
//
// coda
//

#ifdef  RAI_IMPLEMENTATION
#  include "infer.cpp"
#endif

#endif



