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
#ifdef RAI_DAI
#include "infer.h"

#define DAI_WITH_BP
#define DAI_WITH_MF
#define DAI_WITH_HAK
#define DAI_WITH_LC
#define DAI_WITH_TREEEP
#define DAI_WITH_JTREE
#define DAI_WITH_MR

#include "dai/alldai.h"



struct DaiInterface {
  dai::FactorGraph fg;
};

void createDai(DaiInterface& dai, const VariableList& vars, const FactorList& facs){
  uint i, j;
  Factor *f;
  std::vector<dai::Factor> dai_facs;
  for_list(Type,  f,  facs){
    uint nr_members = f->varIds.N;
    
    std::vector<long> labels(nr_members);
    for(j=0; j<nr_members; j++) labels[j]=f->varIds(j);
    
    std::vector<size_t> dims(nr_members);
    for(j=0; j<nr_members; j++) dims[j]=f->dim(j);
    
    dai::VarSet I_vars;
    for(j=0; j<nr_members; j++) I_vars |= dai::Var(f->varIds(j), f->dim(j));
    dai_facs.push_back(dai::Factor(I_vars, 0.0));
    
    // calculate permutation sigma (internally, members are sorted)
    std::vector<size_t> sigma(nr_members);
    dai::VarSet::iterator v = I_vars.begin();
    for(j=0; j<nr_members; j++, v++){
      long search_for = v->label();
      std::vector<long>::iterator j_loc = find(labels.begin(), labels.end(), search_for);
      sigma[j] = j_loc - labels.begin();
    }
    dai::Permute permindex(dims, sigma);
    
    //copy permuted array
    for(j=0; j<f->P.N; j++) dai_facs.back()[permindex.convert_linear_index(j)] = f->P.elem(j);
  }
  
  dai.fg = dai::FactorGraph(dai_facs);
}

#define test(METH, x, OPT) \
  dai::METH x( dai.fg, opts OPT ); \
  x.init(); \
  x.run(); \
  cout <<#METH <<" node marginals:" <<endl; \
  for( size_t i = 0; i < dai.fg.nrVars(); i++ ) \
    cout <<x.belief(dai.fg.var(i)) <<endl; \
  cout <<"Exact log partition sum: " <<x.logZ() <<endl;

/*cout <<#METH <<" factor marginals:" <<endl; \
  for( size_t I = 0; I < dai.fg.nrFactors(); I++ ) \
    cout <<x.belief(dai.fg.factor(I).vars()) <<endl; \
    */

void inference(DaiInterface& dai){
  size_t  maxiter = 10000;
  double  tol = 1e-9;
  size_t  verb = 1;
  
  dai::PropertySet opts;
  opts.Set("maxiter", maxiter);
  opts.Set("tol", tol);
  opts.Set("verbose", verb);
  
  test(JTree, jt, ("updates", std::string("HUGIN")));
  test(BP, bp, ("updates", std::string("SEQFIX"))("logdomain", false));
  test(MF, mf, ("tol", std::string("0.01"))("maxiter", std::string("100")));
  //test(HAK, hak, ("tol", std::string("0.01"))("maxiter", std::string("100"))("maxiter", std::string("100"))("maxiter", std::string("100")));
  //test(LC, lc, ("tol", std::string("0.01"))("maxiter", std::string("100")));
  //test(TreeEP, treeep, ("tol", std::string("0.01"))("maxiter", std::string("100")));
  //test(MR, mr, ("tol", std::string("0.01"))("maxiter", std::string("100")));
  
}
#endif