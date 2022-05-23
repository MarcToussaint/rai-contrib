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
#ifndef RAI_mdp_help_h
#define RAI_mdp_help_h

#include <Core/array.h>

//===========================================================================
//
// fwd declarations from infer.h
//

namespace infer {
struct Variable;
struct Factor;

typedef rai::Array<Variable*> VariableList;
typedef rai::Array<Factor*>   FactorList;
}


namespace mdp {

//===========================================================================
//
// structs for problems and controllers
//

/// struct to store an MDP or POMDP (when the observation model Pyxa is given as well) \ingroup mdp
struct MDP {
  arr Px;    ///< start distribution
  arr Pxax;  ///< transition probs
  arr Pyxa;  ///< observation probs
  arr Rax;   ///< reward expectation as a function of (action, x_before)
  double gamma; ///< discounting factor
  rai::Array<uintA> neighbors; ///< state neighbor lists (for sparse worlds)
};

/// struct to store a generic structured MDP or POMDP
struct MDP_structured {
  infer::VariableList vars;
  infer::VariableList leftVars, rightVars;
  infer::VariableList obsVars, ctrlVars;
  infer::FactorList facs;
  infer::FactorList transFacs, initFacs, rewardFacs, obsFacs;
  double gamma;
  ~MDP_structured();
};

/// struct to store a standard FSC with one hidden node variable \ingroup mdp
struct FSC_lev1 {
  arr P0;    ///< initialization of node
  arr Pa0;   ///< action probs
  arr P0y0;  ///< node transition probs
};

/// struct to store a 2-level FSC with two hidden node variables \ingroup mdp
struct FSC_lev2 {
  arr P0;    ///< initialization of level 0 node
  arr P1;    ///< initialization of level 1 node
  arr Pa0;   ///< action probs
  arr P01y0; ///< level 0 transition probs P(n0'| y0, n1, n0)
  arr P1y01; ///< level 1 transition probs P(n1'| y0, n1, n0);
  
  bool hierarchical;
  arr Pe0;     //P(exit|n0) is fixed!
  arr PE_01;   //P(n0'|n1'; exit=true)
  arr P_E_0y0; //P(n0'|y', n0; exit=false)
  arr PE_1y1;  //P(n1'|y', n1; exit=true)
  arr P_E_11;  //P(n1'|n1; exit=false) is fixed!
};

struct FSC_structured {
  infer::VariableList vars, leftVars, rightVars;
  infer::FactorList   facs, initFacs, transFacs;
  ~FSC_structured();
};


//===========================================================================
//
// loading/saving problems
//

void clearMDP(MDP_structured& mdp);

void writeMDP_arr(const MDP& mdp, const char* arr_file, bool binary=false);
void writeMDP_fg(const MDP_structured& mdp, std::ostream& os, bool brief=false);
void readMDP(MDP_structured& mdp, const char* filename);
void readMDP_arr(MDP& mdp, const char* filename, bool binary=false);
void readMDP_POMDP(MDP& mdp, const char* filename);
void readMDP_fg(MDP_structured& mdp, const char* filename, bool binary=false);
void readMDP_ddgm_tabular(MDP_structured& mdp, const char *filename);

void readImageMaze(MDP& mdp, const char* ppm_file);
void reportMDP(const MDP& mdp);
void convertTonyToMDP(const char* POMDP_file);

/*void OutputAoDot(arr& Pa0, uint numObs, const char* problem,
                           std::vector<std::string>& a,
                           std::vector<std::string>& o,
                           std::map<int, std::string>& a0);  */
// Outputs a .dot file.
// To convert it to a postscript: $ dot -Tps filename.dot -o filename.ps
void OutputDotH(const char* filename, const char* problem, arr& Pa0, arr& Pe0,
                arr& PE_01, arr& P_E_0y0, arr& PE_1y1, arr& P_E_11);
                
// Output the policy to dot format
// Currently only works for 1-level policies
void OutputDot(const char* filename, const char* problem, arr& Pa0, arr& P0y0);
/*void OutputAoDot(arr& Pa0, uint numObs, const char* problem, std::vector<std::string>& a,
                      std::vector<std::string>& o,
                      std::map<int, std::string>& a0);*/

// This function simply outputs the policy in an AMPL readable format.
// It's useful to debug the value of the policy.
void SetAMPLPolicy(arr& P0y0, arr& Pa0);


//===========================================================================
//
// manipulating problems
//

extern byteA global_maze;
enum MDPProblem { tinyMaze, miniMaze, simpleMaze, heavenAndHell };
void generateStandardProblem(MDP& mdp, MDPProblem problem);
void mazeToMDP(MDP& mdp, const byteA& maze);
void addActionNoise(arr& Pxax, double eps);
void tunnel(arr& Pxax, int from=-1, int to=-1);
void tunnelToPx(arr& Pxax, int from, const arr& Px);
void tunnelRaxTo(arr& Pxax, const arr& Rax, int to);
void tunnelRaxToPx(arr& Pxax, const arr& Rax, const arr& Px);
void createNeighorList(MDP& mdp);

//checking normalization
double checkNormalization(const MDP_structured& mdp);
void checkJointNormalization(const MDP_structured& mdp, const FSC_structured& fsc);

//collapsing to a flat representation
void collapseToFlat(MDP& mdpUn, const MDP_structured& mdp);
void convert(MDP_structured& mdp, MDP& mdp_flat);


//===========================================================================
//
// manipulating FSCs
//

void clearFSC(FSC_structured& fsc);

void writeFSC_lev1(const char *arr_file, FSC_lev1& fsc, bool binary=false);
void writeFSC_lev2(const char *arr_file, FSC_lev2& fsc, bool binary=false);
void readFSC_lev1(const char *arr_file, FSC_lev1& fsc, bool binary=false);
void writeFSC_fg(const FSC_structured& fsc, std::ostream& os, bool brief=false);

//high level initialization routines
void standardInitFsc1(FSC_lev1& fsc, const MDP& mdp, uint d0);
void dirichletInitFsc1(FSC_lev1& fsc, const MDP& mdp, uint d0);
void standardInitFsc2(FSC_lev2& fsc, const MDP& mdp, uint d0, uint d1, bool hierarchical);
void standardInitFsc_structured_lev1(FSC_structured& fsc, const MDP& mdp, uint d0);
void standardInitFsc_structured_react(FSC_structured& fsc, const MDP& mdp, uint d0);
void standardInitFsc_structured_lev2(FSC_structured& fsc, const MDP& mdp, uint d0, uint d1);
void standardInitFsc_structured_levels(FSC_structured& fsc, const MDP& mdp, const uintA& levels);
void standardInitFsc_structured_levels(FSC_structured& fsc, const MDP_structured& mdp, const uintA& levels);

//lower level initialization routines
void oneNodeOneAction(arr& Pa0, uint da, uint d0, double uni=1., double noise=1., double det=100.);
void zeroNodeStart(arr& P0, uint d0, double uni=.0, double noise=.0, double zero=1.);
void generalRandTransitions(arr& P, double uni=0., double noise=.1);
void generalDirichlet(arr& P);
void generalNodeTransitions(arr& P0y0, double uni=0., double noise=.1, double stay=0.);
void generalNode0Transition(arr& P01y0,
                            uint d0, uint d1, uint dy,
                            double uni=0., double noise=.1, double stay=0.);
void generalNode1Transition(arr& P1y01,
                            uint d1, uint dy, uint d0,
                            double uni=0., double noise=.1, double stay=0.);
void generalNodeTransition(arr& P, double uni, double noise, double stay);
void initHierarchical(uint exits, FSC_lev2& fsc);
void hierarchicalToFactors(FSC_lev2& fsc);
void neutralSplit(arr& P0y0, arr& Pa0, uint i, double eps=.1);

//collapsing FSCs
void collapse2levelFSC(FSC_lev1& fsc1, const FSC_lev2& fsc2);
void collapseFSC(FSC_lev1& fsc1, const FSC_structured& fsc);
void collapseToTransitionMatrix(arr& P0x0x, arr& R0x, arr& S0x, const FSC_lev1& fsc, const MDP& mdp);


//===========================================================================
//
// Dynamic Programming in a standard MDP
//

double Qvalue(uint a, uint x, const arr&V, const MDP& mdp);
void Qfunction(arr& Q, const arr&V, const MDP& mdp);
void valueIteration(arr& V, const MDP& mdp);
void policyEvaluation(arr& V, const arr& pi, const MDP& mdp);
void maxPolicy(arr& pi, const arr& V, const MDP& mdp);
void prioritizedSweeping(arr& V, const MDP& mdp, double VerrThreshold=1e-4);


//===========================================================================
//
// Dynamic Programming in a POMDP
//

double evaluateFsc1(const FSC_lev1& fsc, const MDP& mdp, uint Tcut);


//===========================================================================
//
// inference, EM
//

enum MstepType { MstepNoisyMax, MstepExact, Mstep11Rule, MstepNone, MstepCopyExpectations };

void mdpEM(const MDP& mdp, arr& pi, arr& hatBeta, uint Tmax, float cutoffTimeFactor, MstepType mstepType, arr *Pvisited=NULL);
void calcPvisited(arr& Pvisited, const MDP& mdp);
double pomdpEM_lev1(const MDP& mdp,
                    FSC_lev1& fsc,
                    uint estepHorizon,
                    bool estepStructured,
                    bool estepIncremental,
                    MstepType mstepType,
                    double mstepRate,
                    double mstepNoise,
                    bool adaptP0,
                    arr* alpha, arr* beta,
                    std::ostream *os);
double pomdpEM_lev2(const MDP& mdp, FSC_lev2& fsm,
                    uint T, bool structuredEstep, bool maxMstep,
                    bool adaptP0, std::ostream *os=NULL);
double pomdpEM_structured(const MDP& mdp,
                          FSC_structured& fsc,
                          uint estepHorizon,
                          bool estepStructured,
                          bool estepIncremental,
                          MstepType mstepType,
                          double mstepRate,
                          double mstepNoise,
                          bool adaptP0,
                          arr* alpha, arr* beta,
                          std::ostream *os);
double pomdpEM_structured(const MDP_structured& mdp,
                          FSC_structured& fsc,
                          uint estepHorizon,
                          bool estepStructured,
                          bool estepIncremental,
                          MstepType mstepType,
                          double mstepRate,
                          double mstepNoise,
                          bool adaptP0,
                          arr* alpha, arr* beta,
                          std::ostream *os);
                          
                          
//===========================================================================
//
// MDP policy manipulation
//

void randomPolicy(arr& pi, const MDP& mdp);
void maxPolicyMap(uintA& piMap, const arr& pi);
void getPxx(arr& Pxx, const arr& Pxax, const arr& pi);
void getMaxPxxMap(uintA& PxxMap, const arr& Pxax, const arr& pi);


//===========================================================================
//
// display routines (when compiled with RAI_FREEGLUT)
//

void showMaze();
void showAB(const arr& alpha, const arr& beta);
void plotPolicyAndValue(const arr& pi, const arr& V, const MDP& mdp, bool wait);
void glDisplayGrey(const arr& x, uint d0, uint d1, bool wait, uint win);
void glDisplayRedBlue(const arr& x, uint d0, uint d1, bool wait, uint win);


//===========================================================================
//
// implementations
//

}//end of namespace mdp

#ifdef  RAI_IMPLEMENTATION
#  include "mdp.cpp"
#  include "mdp_solvers.cpp"
#  include "mdp_opengl.cpp"
#  include "mdp_tony.cpp"
#  include "mdp_pomdpEM_lev1.cpp"
#  include "mdp_pomdpEM_lev2.cpp"
#  include "mdp_pomdpEM_generic.cpp"
#  include "mdp_pomdpEM_structured.cpp"
#endif

#endif
