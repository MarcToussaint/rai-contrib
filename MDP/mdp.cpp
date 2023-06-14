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
#include "mdp.h"
#include <Infer/infer.h>
#include <vector>
#include <map>
#include <iomanip>

uintA TUP(uint a){ return uintA{a}; }
uintA TUP(uint a, uint b){ return uintA{a,b}; }
uintA TUP(uint a, uint b, uint c){ return uintA{a,b,c}; }
uintA TUP(uint a, uint b, uint c, uint d){ return uintA{a,b,c,d}; }
uintA TUP(uint  a, uint b, uint c, uint d, uint e){ return uintA{a,b,c,d,e}; }
uintA TUP(uint  a, uint b, uint c, uint d, uint e, uint f){ return uintA{a,b,c,d,e,f}; }

/** \defgroup mdp MDP and POMDP Module

 */

// @{

byteA mdp::global_maze;

//===========================================================================
//
/// \name loading/saving problems
//@{

/** @brief write an MDP into a file. For binary=false,
  the MDP arrays are dumped into the file in a fixed order in readable
  ASCII format. For binary=true the array tasgs are ASCII and
  readable, but the array contents are stored in binary double format
  (preserving full precision). */
void mdp::writeMDP_arr(const MDP& mdp, const char *filename, bool binary){
  ofstream os;
  if(binary) os.open(filename, std::ios::binary);
  else os.open(filename);
  if(!os.good()) RAI_MSG("could not open file `" <<filename <<"' for output");
  mdp.Pxax.writeTagged(os, "Pxax", binary); os <<endl;
  mdp.Pyxa.writeTagged(os, "Pyxa", binary); os <<endl;
  mdp.Px  .writeTagged(os, "Px", binary);   os <<endl;
  mdp.Rax .writeTagged(os, "Rax", binary);  os <<endl;
  arr(&mdp.gamma, 1, true).writeTagged(os, "gamma", binary);  os <<endl;
}

void mdp::writeMDP_fg(const MDP_structured& mdp, std::ostream& os, bool brief){
  os <<mdp.vars <<'\n' <<endl;
  if(!brief){
    os <<mdp.facs <<endl;
  }else{
    for(infer::Factor *f: mdp.facs){ f->write(os, true); os <<"\n"; }
  }
  os <<"\nmdp . {";
  os <<"\n  leftVars   "; listWrite(mdp.leftVars, os);
  os <<"\n  rightVars  "; listWrite(mdp.rightVars, os);
  os <<"\n  obsVars    "; listWrite(mdp.obsVars, os);
  os <<"\n  ctrlVars   "; listWrite(mdp.ctrlVars, os);
  os <<"\n  initFacs   "; listWrite(mdp.initFacs, os);
  os <<"\n  transFacs  "; listWrite(mdp.transFacs, os);
  os <<"\n  obsFacs    "; listWrite(mdp.obsFacs, os);
  os <<"\n  rewardFacs "; listWrite(mdp.rewardFacs, os);
  os <<"\n  gamma      [" <<mdp.gamma <<"]";
  os <<"\n}" <<endl;
}

void mdp::writeFSC_fg(const FSC_structured& fsc, std::ostream& os, bool brief){
  os <<fsc.vars <<'\n' <<endl;
  if(!brief){
    os <<fsc.facs <<endl;
  }else{
    for(infer::Factor *f: fsc.facs){ f->write(os, true); os <<"\n"; }
  }
  os <<"\nfsc . {";
  os <<"\n  leftVars  "; listWrite(fsc.leftVars, os);
  os <<"\n  rightVars "; listWrite(fsc.rightVars, os);
  os <<"\n  initFacs  "; listWrite(fsc.initFacs, os);
  os <<"\n  transFacs "; listWrite(fsc.transFacs, os);
  os <<"\n}" <<endl;
}

void mdp::readMDP(MDP_structured& mdp, const char *filename){
  const char *ext=filename+(strlen(filename)-6);
  if(!strcmp(ext, "dp_arr")){ MDP mdp_flat; readMDP_arr(mdp_flat, filename); convert(mdp, mdp_flat); } else if(!strcmp(ext, ".POMDP")){ MDP mdp_flat; readMDP_POMDP(mdp_flat, filename); convert(mdp, mdp_flat); } else if(!strcmp(ext, "mdp_fg")) readMDP_fg(mdp, filename);
  else if(!strcmp(ext, "abular")) readMDP_ddgm_tabular(mdp, filename);
  else HALT("don't know extension of file '" <<filename <<"'");
}

/// read and MDP from a file (see writeMDP)
void mdp::readMDP_arr(MDP& mdp, const char *filename, bool binary){
  ifstream is;
  if(binary) is.open(filename, std::ios::binary);
  else is.open(filename);
  if(!is.good()) RAI_MSG("could not open file `" <<filename <<"' for input");
  mdp.Pxax.readTagged(is, "Pxax");
  mdp.Pyxa.readTagged(is, "Pyxa");
  mdp.Px  .readTagged(is, "Px");
  mdp.Rax .readTagged(is, "Rax");
  arr(&mdp.gamma, 1, true).readTagged(is, "gamma");
}

void mdp::clearMDP(MDP_structured& mdp){
  mdp.vars.memMove=true;
  mdp.facs.memMove=true;
  listDelete(mdp.facs);
  listDelete(mdp.vars);
  mdp.leftVars.clear(); mdp.rightVars.clear(); mdp.obsVars.clear(); mdp.ctrlVars.clear();
  mdp.transFacs.clear(); mdp.initFacs.clear(); mdp.rewardFacs.clear(); mdp.obsFacs.clear();
}

mdp::FSC_structured::~FSC_structured(){ clearFSC(*this); }
mdp::MDP_structured::~MDP_structured(){ clearMDP(*this); }

void mdp::clearFSC(FSC_structured& fsc){
  fsc.vars.memMove=true;
  fsc.facs.memMove=true;
  listDelete(fsc.facs);
  listDelete(fsc.vars);
  fsc.leftVars.clear(); fsc.rightVars.clear();
  fsc.transFacs.clear(); fsc.initFacs.clear();
}

void readStringList(rai::Array<rai::String>& strings, istream& is){
  char c;
  rai::String str;
  is >>PARSE("(");
  strings.clear();
  for(;;){
    rai::skip(is);
    is.get(c);
    if(c==')') break; else is.putback(c);
    str.read(is, " ", "), \t\r\n", false);
    strings.append(str);
  }
}

template<class T> void namesToSublist(rai::Array<T*>& sub, const rai::Array<rai::String>& strings, const rai::Array<T*>& list){
  uint i;
  T *v;
  sub.clear();
  for(i=0; i<strings.N; i++){
    NIY//v=listFindByName(list, strings(i));
    CHECK(v, "variable '" <<strings(i) <<"' is unkown");
    sub.append(v);
  }
}

void mdp::convert(MDP_structured& mdp, MDP& mdp_flat){
  clearMDP(mdp);
  
  uint dx=mdp_flat.Pxax.d0, da=mdp_flat.Pxax.d1, dy=mdp_flat.Pyxa.d0;
  infer::Variable *x  = new infer::Variable(dx , "state");
  infer::Variable *a  = new infer::Variable(da , "action");
  infer::Variable *x_ = new infer::Variable(dx , "state'");
  infer::Variable *y_ = new infer::Variable(dy , "obs'");
  mdp.vars      = {y_, x_, a, x};
  mdp.leftVars  = {x};
  mdp.rightVars = {x_};
  mdp.obsVars   = {y_};
  mdp.ctrlVars  = {a};
  
  infer::Factor *Fxax = new infer::Factor({x_, a, x}, "Ftrans");  Fxax->setP(mdp_flat.Pxax);
  infer::Factor *Fyxa = new infer::Factor({y_, x_, a}, "Fobs");   Fyxa->setP(mdp_flat.Pyxa);
  infer::Factor *Fx   = new infer::Factor({x}, "Finit");        Fx->setP(mdp_flat.Px);
  infer::Factor *FRax = new infer::Factor({a, x}, "Freward");    FRax->setP(mdp_flat.Rax);
  mdp.facs       = {Fyxa, Fxax, Fx, FRax};
  mdp.transFacs  = {Fxax};
  mdp.initFacs   = {Fx};
  mdp.rewardFacs = {FRax};
  mdp.obsFacs    = {Fyxa};
  
  mdp.gamma=mdp_flat.gamma;
}

void mdp::readMDP_fg(MDP_structured& mdp, const char *filename, bool binary){
#if 1 //new code
  clearMDP(mdp);
  ifstream is;
  rai::open(is, filename);
  rai::String tag, name;
  rai::Array<rai::String> strings;
  infer::VariableList vars;
  arr P;
  uint d;
  bool insidePomdp=false;
  for(;;){
    tag.read(is, " \t\n\r", " \t\n\r({", false);
    if(!tag.N) break; //end of file
    if(tag=="variable"){
      name.read(is, " \t\n\r", " \t\n\r<({", false);
      is >>PARSE("<") >>d >>PARSE(">");
      mdp.vars.append(new infer::Variable(d, name));
    }
    if(tag=="factor"){
      name.read(is, " \t\n\r", " \t\n\r({", false);
      readStringList(strings, is);
      namesToSublist(vars, strings, mdp.vars);
      is >>P;
      infer::Factor *f=new infer::Factor(vars);
      mdp.facs.append(f);
      f->setP(P);
      f->name = name;
    }
    if(tag=="mdp"){
      name.read(is, " \t\n\r", " \t\n\r({", false);
      is >>PARSE("{");
      insidePomdp=true;
      continue;
    }
    if(tag=="}"){
      insidePomdp=false;
      break;
    }
    if(!insidePomdp) continue;
    if(tag=="leftVars"){  CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.leftVars, strings, mdp.vars); continue; }
    if(tag=="rightVars"){ CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.rightVars, strings, mdp.vars); continue; }
    if(tag=="obsVars"){   CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.obsVars, strings, mdp.vars); continue; }
    if(tag=="ctrlVars"){  CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.ctrlVars, strings, mdp.vars); continue; }
    if(tag=="initFacs"){  CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.initFacs, strings, mdp.facs); continue; }
    if(tag=="transFacs"){ CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.transFacs, strings, mdp.facs); continue; }
    if(tag=="obsFacs"){   CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.obsFacs, strings, mdp.facs); continue; }
    if(tag=="rewardFacs"){CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.rewardFacs, strings, mdp.facs); continue; }
    if(tag=="gamma"){     CHECK(insidePomdp, "");  is >>PARSE("[") >>mdp.gamma >>PARSE("]");  continue; }
    HALT("something's wrong: don't know tag '" <<tag <<"' inside mdp declaration");
  }
  
#else //old code depending on hypergraph
  clearMDP(mdp);
  Graph H;
  H <<FILE(filename);
  //cout <<"read hypergraph: " <<H <<endl;
  uint d;
  //variables
  for_list(Element, e,  H.T) if(e->type=="variable"){
    d=get<rai::String>(e->ats, "values").N;
    mdp.vars.append(new infer::Variable(d, e->name));
  }
  //cout <<"read variables:" <<endl;  listWrite(mdp.vars, cout, "\n  ");
  for_list(Type,  e,  H.T) if(e->type=="factor"){
    mdp.facs.append(new infer::Factor(e->containsIds, get<double>(e->ats, "P")));
    mdp.facs.last()->name = e->name;
  }
  //cout <<"\nread factors:" <<endl;  listWrite(mdp.facs, cout, "\n  ");
  e=listFindByType(H.T, "mdp");
  mdp.gamma = get<double>(e->ats, "discount")(0);
  //cout <<"\ngamma = " <<mdp.gamma <<endl;
  rai::Array<rai::String> S;
  S = get<rai::String>(e->ats, "initializationFacs");
  for(uint i=0; i<S.N; i++) mdp.initFacs.append(listFindByName(mdp.facs, S(i)));
  S = get<rai::String>(e->ats, "transitionFacs");
  for(uint i=0; i<S.N; i++) mdp.transFacs.append(listFindByName(mdp.facs, S(i)));
  S = get<rai::String>(e->ats, "rewardFacs");
  for(uint i=0; i<S.N; i++) mdp.rewardFacs.append(listFindByName(mdp.facs, S(i)));
  S = get<rai::String>(e->ats, "observationFacs");
  for(uint i=0; i<S.N; i++) mdp.obsFacs.append(listFindByName(mdp.facs, S(i)));
  S = get<rai::String>(e->ats, "leftVars");
  for(uint i=0; i<S.N; i++) mdp.leftVars.append(listFindByName(mdp.vars, S(i)));
  S = get<rai::String>(e->ats, "rightVars");
  for(uint i=0; i<S.N; i++) mdp.rightVars.append(listFindByName(mdp.vars, S(i)));
  S = get<rai::String>(e->ats, "observationVars");
  for(uint i=0; i<S.N; i++) mdp.obsVars.append(listFindByName(mdp.vars, S(i)));
  S = get<rai::String>(e->ats, "controlVars");
  for(uint i=0; i<S.N; i++) mdp.ctrlVars.append(listFindByName(mdp.vars, S(i)));
  
  //infer::Factor *f;
  //for_list(infer::Factor, f,  mdp.transFacs) tensorCheckCondNormalization(f->P, 1);
  //for_list(infer::Factor, f,  mdp.initFacs) tensorCheckCondNormalization(f->P, 1);
#endif
}

void mdp::readMDP_ddgm_tabular(MDP_structured& mdp, const char *filename){
  clearMDP(mdp);
  ifstream is;
  rai::open(is, filename);
  rai::String tag, name;
  rai::Array<rai::String> strings;
  rai::Array<rai::Array<rai::String> > values;
  infer::VariableList vars, rewardVars;
  infer::Factor *f;
  arr P;
  uint i;
  bool insidePomdp=false;
  for(;;){
    tag.read(is, " \t\n\r", " \t\n\r({", false);
    if(!tag.N) break; //end of file
    if(tag=="variable"){
      name.read(is, " \t\n\r", " \t\n\r({", false);
      readStringList(strings, is);
      infer::Variable *v=new infer::Variable(strings.N, name);
      mdp.vars.append(v);
      if(v->id>=values.N) values.resizeCopy(v->id+1);      values(v->id)=strings;
      name <<'\'';
      v=new infer::Variable(strings.N, name);
      mdp.vars.append(v);
      if(v->id>=values.N) values.resizeCopy(v->id+1);      values(v->id)=strings;
    }
    if(tag=="mdt"){
      name.read(is, " \t\n\r", " \t\n\r({", false);
      readStringList(strings, is);
      namesToSublist(vars, strings, mdp.vars);
      is >>P;
      mdp.facs.append(f=new infer::Factor(vars));
      f->setP(P);
      f->name = name;
    }
    if(tag=="pomdp"){
      name.read(is, " \t\n\r", " \t\n\r({", false);
      insidePomdp=true;
      continue;
    }
    if(tag=="endpomdp"){
      insidePomdp=false;
      break;
    }
    if(!insidePomdp) continue;
    if(tag=="observable"){   CHECK(insidePomdp, "");  readStringList(strings, is);  for(i=0; i<strings.N; i++) strings(i) <<'\''; /*APPEND A PRIME!*/ namesToSublist(mdp.obsVars, strings, mdp.vars); continue; }
    if(tag=="hidden"){       CHECK(insidePomdp, "");  readStringList(strings, is);  continue; }
    if(tag=="controllable"){ CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.ctrlVars, strings, mdp.vars); continue; }
    if(tag=="utility"){      CHECK(insidePomdp, "");  readStringList(strings, is);  for(i=0; i<strings.N; i++) strings(i) <<'\''; /*APPEND A PRIME!*/ namesToSublist(rewardVars, strings, mdp.vars); continue; }
    //if(tag=="separators"){   CHECK(insidePomdp, "");  readStringList(strings, is);  namesToSublist(mdp.leftVars, strings, mdp.vars); continue; }
    if(tag=="discount"){     CHECK(insidePomdp, "");  is >>PARSE("(") >>mdp.gamma >>PARSE(")");  continue; }
    if(tag=="horizon"){      CHECK(insidePomdp, "");  double z; is >>PARSE("(") >>z >>PARSE(")");  continue; }
    //a CPT declaration:
    name.read(is, " \t\n\r(", ")", true);
    infer::Variable *v=0;
    NIY//v=listFindByName(mdp.vars, tag);  CHECK(v, "");
    NIY//f=listFindByName(mdp.facs, name); CHECK(f, "");
    if(tag(tag.N-1)!='\''){
      mdp.initFacs.append(f);
      NIY//mdp.leftVars.append(listFindByName(mdp.vars, tag));
    } else if(mdp.obsVars.findValue(v)!=-1) mdp.obsFacs.append(f);
    else if(rewardVars.findValue(v)!=-1) mdp.rewardFacs.append(f);
    else if(mdp.ctrlVars.findValue(v)!=-1){ HALT("there should be no CPT for controlled vars!"); } else mdp.transFacs.append(f);
  }
  
  //get rightVars from leftVars
  for(infer::Variable *v: mdp.leftVars){
    int vid=mdp.vars.findValue(v);
    CHECK(vid!=-1, "");
    mdp.rightVars.append(mdp.vars(vid+1));
    //mdp.rightVars.append(::getvar(v->id+1));
  }
  
  //treat the reward facs special!!
  uint j;
  infer::Factor tmp;
  infer::VariableList tmpVars;
  for(infer::Variable *v: rewardVars){
    arr val(v->dim);
    strings=values(v->id);
    CHECK_EQ(strings.N,v->dim, "");
    for(j=0; j<strings.N; j++) strings(j).resetIstream() >>val(j);
    //cout <<"reward values = " <<strings <<' ' <<val <<endl;
    infer::Factor R({v});
    R.setP(val);
    for(infer::Factor *f: mdp.rewardFacs) if(f->variables.findValue(v)!=-1){
      //cout <<*f <<endl;
      tensorMultiply(*f, R);
      tmpVars=f->variables;
      tmpVars.removeValue(v);
      tensorMarginal(tmp, *f, tmpVars);
      *f = tmp;
      //cout <<*f <<endl;
    }
  }
  if(mdp.rewardFacs.N>1){ // multiple reward functions -> add them to a single factor
    tmpVars.clear();
    for(infer::Factor *f: mdp.rewardFacs) tmpVars = setUnion(tmpVars, f->variables);
    infer::Factor *ff=new infer::Factor(tmpVars);
    ff->name="REWARD_AUTO";
    ff->P.setZero();
    //cout <<*ff <<endl;
    for(infer::Factor *f: mdp.rewardFacs){
      //cout <<*f <<endl;
      tensorAdd(*ff, *f);
      //cout <<*ff <<endl;
    }
    for(infer::Factor *f: mdp.rewardFacs) mdp.facs.removeValue(f);
    listDelete(mdp.rewardFacs);
    mdp.rewardFacs.append(ff);
    mdp.facs.append(ff);
  }
}

/// .
void mdp::writeFSC_lev1(const char *filename, FSC_lev1& fsc, bool binary){
  ofstream os;
  if(binary) os.open(filename, std::ios::binary);
  else os.open(filename);
  fsc.P0  .writeTagged(os, "P0", binary);    os <<endl;
  fsc.Pa0 .writeTagged(os, "Pa0", binary);   os <<endl;
  fsc.P0y0.writeTagged(os, "P0y0", binary);  os <<endl;
}

/// .
void mdp::writeFSC_lev2(const char *filename, FSC_lev2& fsc, bool binary){
  ofstream os;
  if(binary) os.open(filename, std::ios::binary);
  else os.open(filename);
  fsc.P0   .writeTagged(os, "P0", binary);    os <<endl;
  fsc.P1   .writeTagged(os, "P1", binary);    os <<endl;
  fsc.Pa0  .writeTagged(os, "Pa0", binary);   os <<endl;
  fsc.P01y0.writeTagged(os, "P01y0", binary); os <<endl;
  fsc.P1y01.writeTagged(os, "P1y01", binary); os <<endl;
}

/// .
void mdp::readFSC_lev1(const char *filename, FSC_lev1& fsc, bool binary){
  ifstream is;
  if(binary) is.open(filename, std::ios::binary);
  else is.open(filename);
  fsc.P0  .readTagged(is, "P0");
  fsc.Pa0 .readTagged(is, "Pa0");
  fsc.P0y0.readTagged(is, "P0y0");
}

/// .
void mdp::convertTonyToMDP(const char* POMDP_file){
  MDP mdp;
  readMDP_POMDP(mdp, POMDP_file);
  writeMDP_arr(mdp, STRING(POMDP_file <<".arr"));
}

/// .
void mdp::reportMDP(const MDP& mdp){
  uint xx, y, x, a;
  //
  uintA xm(mdp.Pxax.d0);
  for(a=0; a<mdp.Pxax.d1; a++){
    xm.reshape(mdp.Pxax.d0); xm.setZero();
    for(x=0; x<mdp.Pxax.d0; x++) for(xx=0; xx<mdp.Pxax.d2; xx++)
        if(mdp.Pxax(xx, a, x)>mdp.Pxax(xm(x), a, x)) xm(x)=xx;
    xm.reshape(global_maze.d0, global_maze.d1);
    cout <<"\nmax transitions for a=" <<a <<"\n" <<xm;
  }
  
  //
  uintA ym(mdp.Pyxa.d1); ym.setZero();
  for(y=0; y<mdp.Pyxa.d0; y++)
    for(x=0; x<mdp.Pyxa.d1; x++)
      for(a=0; a<mdp.Pyxa.d2; a++) if(mdp.Pyxa(y, x, a)>mdp.Pyxa(ym(x), x, a)) ym(x)=y;
  ym.reshape(global_maze.d0, global_maze.d1);
  cout <<"\nmax observations =\n" <<ym;
}

/// .
void mdp::readImageMaze(MDP& mdp, const char* filename){
  byteA img;
  read_ppm(img, filename);
  flip_image(img);
  byteA maze(img.d0, img.d1);
  uint x, y, dx=maze.d1, dy=maze.d0;
  for(x=0; x<dx; x++) for(y=0; y<dy; y++){
      if(img(y, x, 0)==255 && img(y, x, 1)==255 && img(y, x, 2)==255) maze(y, x)=0; //white -> free
      else if(img(y, x, 0)==0   && img(y, x, 1)==0   && img(y, x, 2)==0) maze(y, x)=1; //black -> wall
      else if(img(y, x, 0)==255 && img(y, x, 1)==0   && img(y, x, 2)==0) maze(y, x)=2; //red   -> start
      else if(img(y, x, 0)==0   && img(y, x, 1)==255 && img(y, x, 2)==0) maze(y, x)=3; //green -> goal
      else if(img(y, x, 0)==0   && img(y, x, 1)==0   && img(y, x, 2)==255) maze(y, x)=4; //blue  -> avoid
      else HALT("maze image file (ppm) needs to be black/white/red/green/blue");
    }
  mazeToMDP(mdp, maze);
}

// @}

/// .
void mazeToP(const byteA& maze, arr& Px, arr& Pxax, arr& Pyxa, arr& Rax, bool trapWalls){
  CHECK_EQ(maze.nd,2, "");
  mdp::global_maze = maze;
  uint nx=maze.N, na=5, ny=16, dx=maze.d1, dy=maze.d0;
  uint x, a, i, j, k;
  
  Pxax.resize(nx, na, nx);  Pxax.setZero();
  Pyxa.resize(ny, nx, na);  Pyxa.setZero();
  Px  .resize(nx);        Px  .setZero();
  Rax .resize(na, nx);     Rax .setZero();
  
  /*
  0: free space
  1: wall
  2: start state (can be multiple)
  3: positive reward +1
  4: negative reward -1
  */
  
  for(x=0; x<nx; x++){
    i=x/dx; j=x%dx;
    if(maze(i, j)==1){
      Pxax(x, 0, x)=Pxax(x, 1, x)=Pxax(x, 2, x)=Pxax(x, 3, x)=Pxax(x, 4, x)=1.;
      for(a=0; a<na; a++) Pyxa(0, x, a)=1.;
    }else{
      k=0; //observation indicator
      if(i>0    && (trapWalls || maze(i-1, j)!=1)){ Pxax(x-dx, 0, x)=1.; }else{ k|=1; Pxax(x, 0, x)=1.; }//up
      if(j<dx-1 && (trapWalls || maze(i, j+1)!=1)){ Pxax(x+1 , 1, x)=1.; }else{ k|=2; Pxax(x, 1, x)=1.; }//right
      if(i<dy-1 && (trapWalls || maze(i+1, j)!=1)){ Pxax(x+dx, 2, x)=1.; }else{ k|=4; Pxax(x, 2, x)=1.; }//down
      if(j>0    && (trapWalls || maze(i, j-1)!=1)){ Pxax(x-1 , 3, x)=1.; }else{ k|=8; Pxax(x, 3, x)=1.; }//left
      Pxax(x, 4, x) = 1.; //stay
      for(a=0; a<na; a++) Pyxa(k, x, a)=1.; //observe `k'
    }
    if(maze(i, j)==2){
      Px(x)=1.;
    }
    if(maze(i, j)==3){
      for(a=0; a<na; a++) Rax(a, x)=1.;
    }
    if(maze(i, j)==4){
      for(a=0; a<na; a++) Rax(a, x)=-1.;
    }
  }
  
  normalizeDist(Px);
  
  //::checkNormalization(Pxax);
  //::checkNormalization(Pyxa);
  //::checkNormalization(Px);
}

/// .
void mdp::mazeToMDP(MDP& mdp, const byteA& maze){
  ::mazeToP(maze, mdp.Px, mdp.Pxax, mdp.Pyxa, mdp.Rax, true);
}

void mdp::addActionNoise(arr& Pxax, double eps){
  uint nx=Pxax.d0, na=Pxax.d1;
  uint i, j, a;
  arr Pxx;
  tensorMarginal(Pxx, Pxax, uintA{0, 2});
  Pxx /= (double)na;
  //::checkNormalization(Pxx);
  for(i=0; i<nx; i++) for(j=0; j<nx; j++) for(a=0; a<na; a++){
        Pxax(i, a, j) = (1.-eps) * Pxax(i, a, j) + eps * Pxx(i, j);
      }
  //::checkNormalization(Pxax);
}

/// .
void mdp::tunnel(arr& Pxax, int from, int to){
  uint nx=Pxax.d0, na=Pxax.d1;
  if(from<0) from+=nx;
  if(to  <0) to  +=nx;
  uint x, a;
  for(a=0; a<na; a++){
    for(x=0; x<nx; x++) Pxax(x, a, from)=0;
    Pxax(to, a, from)=1;
  }
  //::checkNormalization(Pxax);
}

/// .
void mdp::tunnelToPx(arr& Pxax, int from, const arr& Px){
  uint nx=Pxax.d0, na=Pxax.d1;
  if(from<0) from+=nx;
  uint x, a;
  for(a=0; a<na; a++){
    for(x=0; x<nx; x++) Pxax(x, a, from)=Px(x);
  }
  //::checkNormalization(Pxax);
}

/// .
void mdp::tunnelRaxToPx(arr& Pxax, const arr& Rax, const arr& Px){
  uint nx=Pxax.d0, na=Pxax.d1;
  uint x, y, a;
  for(x=0; x<nx; x++) for(a=0; a<na; a++){
      if(Rax(a, x)){
        for(y=0; y<nx; y++) Pxax(y, a, x)=Px(y);
      }
    }
  //::checkNormalization(Pxax);
}

/// .
void mdp::tunnelRaxTo(arr& Pxax, const arr& Rax, int to){
  uint nx=Pxax.d0, na=Pxax.d1;
  if(to  <0) to  +=nx;
  uint x, y, a;
  for(x=0; x<nx; x++) for(a=0; a<na; a++){
      if(Rax(a, x)){
        for(y=0; y<nx; y++) Pxax(y, a, x)=0.;
        Pxax(to, a, x)=1;
      }
    }
  //::checkNormalization(Pxax);
}

/// .
void mdp::neutralSplit(arr& P0y0, arr& Pa0, uint i, double eps){
  uint n0=P0y0.d0, ny=P0y0.d1, na=Pa0.d0;
  uint a, y, j;
  
  //make them bigger by inserting zeros
  P0y0.reshape(n0, ny*n0);
  P0y0.insRows(n0, 1);
  P0y0.reshape((n0+1)*ny, n0);
  P0y0.insColumns(n0, 1);
  Pa0.insColumns(n0, 1);
  P0y0.reshape(n0+1, ny, n0+1);
  Pa0.reshape(na, n0+1);
  
  //copy the output of the ith node indentically
  for(a=0; a<na; a++) Pa0(a, n0) = Pa0(a, i);
  for(j=0; j<n0; j++) for(y=0; y<ny; y++) P0y0(j, y, n0) = P0y0(j, y, i);
  ::checkNormalization(Pa0);
  ::checkNormalization(P0y0);
  
  //split inputs to new node between i and new node
  double p, s;
  for(j=0; j<n0; j++) for(y=0; y<ny; y++){
      p=P0y0(i, y, j);
      s=.5+rnd.uni(-eps, eps);
      P0y0(i, y, j)  = s*p;
      P0y0(n0, y, j) = (1.-s)*p;
    }
  ::checkNormalization(P0y0);
}

/// .
void generateStandardProblem(mdp::MDPProblem problem, arr& Px, arr& Pxax, arr& Pyxa, arr& Rax){
  byteA maze;
  switch(problem){
    case mdp::tinyMaze:
      STRING("[\
          2 1;\
          0 3]") >>maze;
      maze-=(byte)'0';
      
      mazeToP(maze, Px, Pxax, Pyxa, Rax, false);
      mdp::tunnel(Pxax, -1, 1);
      
      //addActionNoise(Pxax, .1);
      break;
    case mdp::miniMaze:
      STRING("[\
          0 0 0;\
          0 1 0;\
          2 1 3]") >>maze;
      maze-=(byte)'0';
      
      mazeToP(maze, Px, Pxax, Pyxa, Rax, false);
      mdp::tunnelRaxToPx(Pxax, Rax, Px);
      mdp::addActionNoise(Pxax, rai::getParameter<double>("mazeNoise"));
      mdp::showMaze();
      
      //mdp::addActionNoise(Pxax, .1);
      break;
    case mdp::simpleMaze:
      STRING("[\
          2 1 0 0 0 1 0 0 0;\
          0 1 0 1 0 1 0 1 0;\
          0 1 0 1 0 1 0 1 0;\
          0 1 0 1 0 1 0 1 0;\
          0 0 0 1 0 0 0 1 3]") >>maze;

      /*STRING("[\
      1 1 1 1 1 1 1;\
      1 3 0 0 0 0 1;\
      1 1 1 1 1 0 1;\
      1 2 0 0 0 0 1;\
      1 1 1 1 1 1 1]") >>maze;*/
      maze-=(byte)'0';
      
      mazeToP(maze, Px, Pxax, Pyxa, Rax, false);
      mdp::tunnel(Pxax, -1, -2);
      
      //addActionNoise(Pxax, .1);
      break;
    case mdp::heavenAndHell:
      STRING("[\
          3 0 4 1 4 0 3 ;\
          1 2 1 1 1 2 1 ;\
          0 0 1 1 1 0 0 ]") >>maze;
      maze-=(byte)'0';
      
      mazeToP(maze, Px, Pxax, Pyxa, Rax, false);
      mdp::tunnel(Pxax, 0, 7);
      mdp::tunnel(Pxax, 6, 13);
      
      mdp::tunnel(Pxax, 2, 9);
      mdp::tunnel(Pxax, 4, 11);
      
      mdp::tunnelToPx(Pxax, 0, Px);
      mdp::tunnelToPx(Pxax, 6, Px);
      mdp::tunnelToPx(Pxax, 2, Px);
      mdp::tunnelToPx(Pxax, 4, Px);
      
      //addActionNoise(Pxax, .1);
      break;
    default:
      HALT("don't know problem " <<problem);
  }
}

void mdp::generateStandardProblem(MDP& mdp, MDPProblem problem){
  ::generateStandardProblem(problem, mdp.Px, mdp.Pxax, mdp.Pyxa, mdp.Rax);
}

void mdp::createNeighorList(MDP& mdp){
  uint x, y, a, A=mdp.Pxax.d1, X=mdp.Px.d0;
  mdp.neighbors.resize(X);
  for(y=0; y<X; y++) for(x=0; x<=y; x++) for(a=0; a<A; a++){
        if(mdp.Pxax(y, a, x)!=0. || mdp.Pxax(x, a, y)!=0){
          mdp.neighbors(x).append(y);
          if(x!=y) mdp.neighbors(y).append(x);
          a=A;
        }
      }
}

double mdp::checkNormalization(const MDP_structured& mdp){
  infer::Factor post;
  eliminationAlgorithm(post, (mdp.transFacs, mdp.initFacs), infer::VariableList());
  arr P;
  post.getP(P);
  return P.scalar();
}

void mdp::checkJointNormalization(const MDP_structured& mdp, const FSC_structured& fsc){
  infer::Factor post;
  eliminationAlgorithm(post, (fsc.transFacs, mdp.obsFacs, mdp.transFacs, fsc.initFacs, mdp.initFacs), infer::VariableList());
  arr P;
  post.getP(P);
  //cout <<post <<endl;
  CHECK(fabs(P.scalar()-1.)<1e-10, "MDP & FSC are not jointly normalized: Z=" <<P.scalar());
}

void mdp::collapseToFlat(MDP& mdpUn, const MDP_structured& mdp){
  uint i;
  uint dx = 1;  for(i=0; i<mdp.leftVars.N; i++) dx *= mdp.leftVars(i)->dim;
  uint dy = 1;  for(i=0; i<mdp.obsVars.N; i++)  dy *= mdp.obsVars(i)->dim;
  uint da = 1;  for(i=0; i<mdp.ctrlVars.N; i++) da *= mdp.ctrlVars(i)->dim;
  
  infer::Factor post;
  eliminationAlgorithm(post, mdp.transFacs, (mdp.rightVars, mdp.ctrlVars, mdp.leftVars));
  post.getP(mdpUn.Pxax);
  mdpUn.Pxax.reshape(dx, da, dx);
  
  infer::VariableList obsvars = (mdp.obsVars, mdp.rightVars, mdp.ctrlVars);
  infer::Factor dummy(obsvars);
  infer::VariableList tmp; NIY// = {&dummy};
  eliminationAlgorithm(post, (tmp, mdp.obsFacs), obsvars);
  post.getP(mdpUn.Pyxa);
  mdpUn.Pyxa.reshape(dy, dx, da);
  
  eliminationAlgorithm(post, mdp.initFacs, mdp.leftVars);
  post.getP(mdpUn.Px);
  mdpUn.Px.reshape(dx);
  
  infer::VariableList rewardvars=(mdp.ctrlVars, mdp.leftVars);
  infer::Factor dummy2(rewardvars);
  eliminationAlgorithm(post, (/*{&dummy2}*/tmp, mdp.rewardFacs), rewardvars);
  //eliminationAlgorithm(post, mdp.rewardFacs, ids(cat(mdp.ctrlVars, mdp.leftVars)));
  //listWrite(mdp.rewardFacs, cout);
  post.getP(mdpUn.Rax);
  mdpUn.Rax.reshape(da, dx);
  
  mdpUn.gamma=mdp.gamma;
}

//===========================================================================
//
// collapsing controllers into a 1-level controller
//

void mdp::collapse2levelFSC(FSC_lev1& fsc1, const FSC_lev2& fsc2){

  uint d0=fsc2.P0.d0, d1=fsc2.P1.d0, da=fsc2.Pa0.d0, dy=fsc2.P1y01.d1;
  uint i, j, k;
  
  fsc1.P0.resize(d0, d1);
  tensorEquation(fsc1.P0, fsc2.P0, uintA{0}, fsc2.P1, uintA{1}, 0);
  fsc1.P0.reshape(d0*d1);
  
  fsc1.Pa0.resize(da, d0, d1);
  for(i=0; i<da; i++) for(j=0; j<d0; j++) for(k=0; k<d1; k++)
        fsc1.Pa0(i, j, k) = fsc2.Pa0(i, j);
  fsc1.Pa0.reshape(da, d0*d1);
  
  fsc1.P0y0.resize(TUP(d0, d1, dy, d0, d1));
  tensorEquation(fsc1.P0y0, fsc2.P01y0, uintA{0, 1, 2, 3}, fsc2.P1y01, TUP(1, 2, 3, 4), 0);
  fsc1.P0y0.reshape(d0*d1, dy, d0*d1);
  
  tensorCheckCondNormalization(fsc1.P0  , 1, 1e-10);
  tensorCheckCondNormalization(fsc1.Pa0 , 1, 1e-10);
  tensorCheckCondNormalization(fsc1.P0y0, 1, 1e-10);
};

void mdp::collapseFSC(FSC_lev1& fsc1, const FSC_structured& fsc){
  NIY;
};

void mdp::collapseToTransitionMatrix(arr& P0x0x, arr& R0x, arr& S0x, const FSC_lev1& fsc, const MDP& mdp){
  //transition matrix
  uint d0=fsc.P0.d0, dx=mdp.Pxax.d0, da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  arr P0yxa0x;
  P0yxa0x.resize(TUP(d0, dy, dx, da, d0, dx));
  P0yxa0x = 1.;
  tensorMultiply(P0yxa0x, fsc.Pa0 , TUP(3, 4));
  tensorMultiply(P0yxa0x, mdp.Pxax, TUP(2, 3, 5));
  tensorMultiply(P0yxa0x, mdp.Pyxa, TUP(1, 2, 3));
  tensorMultiply(P0yxa0x, fsc.P0y0, TUP(0, 1, 4));
  tensorMarginal(P0x0x, P0yxa0x, TUP(0, 2, 4, 5));
  P0x0x.reshape(d0*dx, d0*dx);
  
  //reward
  R0x.resize(d0, dx);
  tensorEquation(R0x, mdp.Rax, TUP(2, 1), fsc.Pa0, TUP(2, 0), 1);
  R0x.reshape(d0*dx);
  
  //start
  S0x.resize(d0, dx);
  tensorEquation(S0x, fsc.P0, TUP(0), mdp.Px, TUP(1), 0);
  S0x.reshape(d0*dx);
}

double mdp::evaluateFsc1(const FSC_lev1& fsc, const MDP& mdp, uint Tcut){
  arr P, R, S;
  collapseToTransitionMatrix(P, R, S, fsc, mdp);
  arr V(R.N);
  V=R;
  for(uint t=0; t<Tcut; t++){
    V = V*P;
    V *=mdp.gamma;
    V += R;
  }
  return scalarProduct(V, S);
}

//===========================================================================
//
// standard initializations
//

void mdp::standardInitFsc1(FSC_lev1& fsc, const MDP& mdp, uint d0){
  uint da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  fsc.P0y0.resize(d0, dy, d0);
  oneNodeOneAction(fsc.Pa0, da, d0, 1., 1., 100.);
  generalNodeTransitions(fsc.P0y0, .1, .1, 0.);
  zeroNodeStart(fsc.P0, d0);
}

void mdp::dirichletInitFsc1(FSC_lev1& fsc, const MDP& mdp, uint d0){
  uint da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  //fsc.P0.resize(d0);          generalDirichlet(fsc.P0);
  fsc.Pa0.resize(da, d0);      generalDirichlet(fsc.Pa0);
  fsc.P0y0.resize(d0, dy, d0);  generalDirichlet(fsc.P0y0);
  zeroNodeStart(fsc.P0, d0);
}

void mdp::standardInitFsc2(FSC_lev2& fsc, const MDP& mdp, uint d0, uint d1, bool hierarchical){
  //init policy
  uint da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  fsc.P01y0.resize(TUP(d0, d1, dy, d0));
  fsc.P1y01.resize(TUP(d1, dy, d0, d1));
  oneNodeOneAction(fsc.Pa0, da, d0, 1., 1., 100.);
  generalNode0Transition(fsc.P01y0, d0, d1, dy, .1, .1, 0.);
  zeroNodeStart(fsc.P0, d0);
  generalNode1Transition(fsc.P1y01, d1, dy, d0, .1, .1, 1.);
  zeroNodeStart(fsc.P1, d1);
  if(hierarchical){
    NIY//initHierarchical(min(TUP(da, d1, d0/2)), fsc);
  }else{
    fsc.hierarchical=false;
  }
}

void mdp::standardInitFsc_structured_lev1(FSC_structured& fsc, const MDP& mdp, uint d0){
  clearFSC(fsc);
  
  uint dx=mdp.Pxax.d0, da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  if(da>d0) RAI_MSG("#actions > #node0-states: that's not going to work well!");
  infer::Variable *x  = new infer::Variable(dx , "state(t)");
  infer::Variable *y  = new infer::Variable(dy , "observation(t)");
  infer::Variable *n0 = new infer::Variable(d0 , "node0(t)");
  infer::Variable *a  = new infer::Variable(da , "action(t)");
  infer::Variable *x_ = new infer::Variable(dx , "state(t+1)");
  infer::Variable *y_ = new infer::Variable(dy , "observation(t+1)");
  infer::Variable *n0_= new infer::Variable(d0 , "node0(t+1)");
  fsc.vars     = {n0_, y_, x_, a, n0, y, x};
  fsc.leftVars = {n0};
  fsc.rightVars= {n0_};
  
  infer::Factor *F0   = new infer::Factor({n0});
  infer::Factor *Fa0  = new infer::Factor({a, n0});
  infer::Factor *F0y0 = new infer::Factor({n0_, y_, n0});
  fsc.facs     = {F0y0, Fa0, F0};
  fsc.initFacs = {F0};
  fsc.transFacs= {F0y0, Fa0};
  
  oneNodeOneAction(Fa0->P, da, d0, 1., 1., 100.);
  generalNodeTransitions(F0y0->P, .1, .1, 0.);
  zeroNodeStart(F0->P, d0);
}

void mdp::standardInitFsc_structured_react(FSC_structured& fsc, const MDP& mdp, uint d0){
  clearFSC(fsc);
  
  uint dx=mdp.Pxax.d0, da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  infer::Variable *x  = new infer::Variable(dx , "state(t)");
  infer::Variable *y  = new infer::Variable(dy , "observation(t)");
  infer::Variable *n0 = new infer::Variable(d0 , "node0(t)");
  infer::Variable *a  = new infer::Variable(da , "action(t)");
  infer::Variable *x_ = new infer::Variable(dx , "state(t+1)");
  infer::Variable *y_ = new infer::Variable(dy , "observation(t+1)");
  infer::Variable *n0_= new infer::Variable(d0 , "node0(t+1)");
  fsc.vars     = {n0_, y_, x_, a, n0, y, x};
  fsc.leftVars = {n0};
  fsc.rightVars= {n0_};
  
  infer::Factor *F0    = new infer::Factor({n0});
  infer::Factor *Fa0y  = new infer::Factor({a, n0, y});
  infer::Factor *F0ya0 = new infer::Factor({n0_, y_, n0});
  fsc.facs      = {F0ya0, Fa0y, F0};
  fsc.initFacs  = {F0};
  fsc.transFacs = {F0ya0, Fa0y};
  
  generalRandTransitions(Fa0y->P, 1., 1.);
  generalNodeTransition(F0ya0->P, .1, .1, 0.);
  zeroNodeStart(F0->P, d0);
}

void mdp::standardInitFsc_structured_lev2(FSC_structured& fsc, const MDP& mdp, uint d0, uint d1){
  clearFSC(fsc);
  
  uint dx=mdp.Pxax.d0, da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  if(da>d0) RAI_MSG("#action > #node0-states: that's not going to work well!");
  infer::Variable *x  = new infer::Variable(dx , "state(t)");
  infer::Variable *y  = new infer::Variable(dy , "observation(t)");
  infer::Variable *n1 = new infer::Variable(d1 , "node1(t)");
  infer::Variable *n0 = new infer::Variable(d0 , "node0(t)");
  infer::Variable *a  = new infer::Variable(da , "action(t)");
  infer::Variable *x_ = new infer::Variable(dx , "state(t+1)");
  infer::Variable *y_ = new infer::Variable(dy , "observation(t+1)");
  infer::Variable *n1_= new infer::Variable(d1 , "node1(t+1)");
  infer::Variable *n0_= new infer::Variable(d0 , "node0(t+1)");
  fsc.vars     = {n0_, n1_, y_, x_, a, n0, n1, y, x};
  fsc.leftVars = {n0 , n1};
  fsc.rightVars= {n0_, n1_};
  
  infer::Factor *F0   = new infer::Factor({n0});
  infer::Factor *F1   = new infer::Factor({n1});
  infer::Factor *Fa0  = new infer::Factor({a, n0});
  infer::Factor *F01y0 = new infer::Factor({n0_, n1_, y_, n0});
  infer::Factor *F1y01 = new infer::Factor({n1_, y_ , n0, n1});
  fsc.facs      = {F01y0, F1y01, Fa0, F0, F1};
  fsc.initFacs  = {F0, F1};
  fsc.transFacs = {F01y0, F1y01, Fa0};
  
  oneNodeOneAction(Fa0->P, da, d0, 1., 1., 100.);
  generalNodeTransitions(F01y0->P, .1, .1, 0.);
  zeroNodeStart(F0->P, d0);
  generalNodeTransitions(F1y01->P, .1, .1, 1.);
  zeroNodeStart(F1->P, d1);
}

void mdp::standardInitFsc_structured_levels(FSC_structured& fsc, const MDP& mdp, const uintA& levels){
  clearFSC(fsc);
  
  uint dx=mdp.Pxax.d0, da=mdp.Pxax.d1, dy=mdp.Pyxa.d0;
  uint i, m=levels.N;
  if(da>levels(0)) RAI_MSG("#actions " <<da <<" > #node0-states " <<levels(0) <<" -- that's not going to work well!");
  infer::VariableList nodes(m), nodes_(m);
  infer::Variable *x  = new infer::Variable(dx , "state(t)");
  infer::Variable *y  = new infer::Variable(dy , "observation(t)");
  for(i=m; i--;) nodes(i) = new infer::Variable(levels(i) , STRING("node" <<i <<"(t)"));
  infer::Variable *a  = new infer::Variable(da , "action(t)");
  infer::Variable *x_ = new infer::Variable(dx , "state(t+1)");
  infer::Variable *y_ = new infer::Variable(dy , "observation(t+1)");
  for(i=m; i--;) nodes_(i) = new infer::Variable(levels(i) , STRING("node" <<i <<"(t+1)"));
  fsc.vars     = (nodes_, infer::VariableList{y_, x_, a}, nodes, infer::VariableList{y, x});
  fsc.leftVars = nodes ;
  fsc.rightVars= nodes_;
  
  infer::FactorList Finit(m), Ftran(m);
  for(i=m; i--;) Finit(i) = new infer::Factor({nodes(i)});
  infer::Factor *Fa0  = new infer::Factor({a, nodes(0)});
  if(m==1){
    i=0;               Ftran(i) = new infer::Factor({nodes_(i), y_, nodes(i)});
  }else{
    i=m-1;             Ftran(i) = new infer::Factor({nodes_(i), y_ , nodes(i-1), nodes(i)}); //top node
    for(i=m-2; i>0; i--) Ftran(i) = new infer::Factor({nodes_(i), nodes(i+1), y_ , nodes(i-1), nodes(i)}); //middle nodes
    i=0;               Ftran(i) = new infer::Factor({nodes_(i), nodes(i+1), y_ , nodes(i)}); //bottom node
  }
  fsc.facs      = (Ftran, infer::VariableList{/*Fa0*/}, Finit); NIY
  fsc.initFacs  = Finit;
  NIY//fsc.transFacs = (Ftran, infer::VariableList{/*Fa0*/}); NIY
  
  oneNodeOneAction(Fa0->P, da, levels(0), 1., 1., 100.);
  for(i=0; i<m; i++){
    if(i) generalNodeTransitions(Ftran(i)->P, .1, .1, 0.);
    else  generalNodeTransitions(Ftran(i)->P, .1, .1, 1.);
    zeroNodeStart(Finit(i)->P, levels(i));
  }
}

void mdp::standardInitFsc_structured_levels(FSC_structured& fsc, const MDP_structured& mdp, const uintA& levels){
#if 0
  clearFSC(fsc);
  
  uint i, m=levels.N;
  
  //----- find the variable ids for the mdp world:
  //infer::Variable *x  = listFindByName(mdp.vars, "state0");
  //infer::Variable *y  = listFindByName(mdp.vars, "observation0");
  //infer::Variable *a  = listFindByName(mdp.vars, "action");
  //infer::Variable *x_ = listFindByName(mdp.vars, "state1");
  //infer::Variable *y_ = listFindByName(mdp.vars, "observation1");
  
  uint adim=1;
  for(infer::Variable *v: mdp.ctrlVars) adim*=v->dim;
  
  if(adim>levels(0)) RAI_MSG("#actions " <<adim <<" > #node0-states " <<levels(0) <<" -- that's not going to work well!");
  infer::VariableList nodes(m), nodes_(m);
  for(i=m; i--;) nodes(i) = new infer::Variable(levels(i) , STRING("node" <<i));
  for(i=m; i--;) nodes_(i) = new infer::Variable(levels(i) , STRING("node" <<i <<"'"));
  fsc.vars     = (nodes_, nodes);
  fsc.leftVars = nodes ;
  fsc.rightVars= nodes_;
  
  infer::FactorList Finit(m), Ftran(m);
  for(i=m; i--;) Finit(i) = new infer::Factor({nodes(i)}, STRING("Finit" <<i));
  infer::Factor *Fa0  = new infer::Factor((mdp.ctrlVars, infer::FactorList{nodes(0)}), STRING("Faction"));
  if(m==1){
    i=0;               Ftran(i) = new infer::Factor(cat({nodes_(i)}, mdp.obsVars, {nodes(i)}), STRING("Ftrans" <<i));
  }else{
    i=m-1;             Ftran(i) = new infer::Factor(cat({nodes_(i)}, mdp.obsVars, {nodes(i-1), nodes(i)}), STRING("Ftrans" <<i)); //top node
    for(i=m-2; i>0; i--) Ftran(i) = new infer::Factor(cat({nodes_(i), nodes(i+1)}, mdp.obsVars, {nodes(i-1), nodes(i)}), STRING("Ftrans" <<i)); //middle nodes
    i=0;               Ftran(i) = new infer::Factor(cat({nodes_(i), nodes(i+1)}, mdp.obsVars, {nodes(i)}), STRING("Ftrans" <<i)); //bottom node
  }
  fsc.facs      = cat(Ftran, {Fa0}, Finit);
  fsc.initFacs  = Finit;
  fsc.transFacs = cat(Ftran, {Fa0});

  oneNodeOneAction(Fa0->P, adim, levels(0), 1., 1., 100.);
  for(i=0; i<m; i++){
    if(i) generalNodeTransitions(Ftran(i)->P, .1, .1, 0.);
    else  generalNodeTransitions(Ftran(i)->P, .1, .1, 1.);
    zeroNodeStart(Finit(i)->P, levels(i));
  }
#endif
NIY
  //cout <<"\nfsc initFacs:" <<endl;  listWrite(fsc.initFacs, cout, "\n  ");
  //cout <<"\nfsc transFacs:" <<endl;  listWrite(fsc.transFacs, cout, "\n  ");
  //cout <<"\nfsc facs:" <<endl;  listWrite(fsc.facs, cout, "\n  ");
}

void mdp::oneNodeOneAction(arr& Pa0, uint da, uint d0, double uni, double noise, double det){
  Pa0.resize(da, d0);
  Pa0 = uni;
  rndUniform(Pa0, 0., noise, true);
  for(uint n=1; n<d0; n++) Pa0((n-1)%da, n)+=det;
  tensorCondNormalize(Pa0, 1);
  tensorCheckCondNormalization(Pa0, 1, 1e-10);
}

/// .
void mdp::zeroNodeStart(arr& P0, uint d0, double uni, double noise, double zero){
  P0.resize(d0);
  P0 = uni;
  rndUniform(P0, .0, noise, true);
  P0(0)+=zero;
  tensorCondNormalize(P0, 1);
  tensorCheckCondNormalization(P0, 1, 1e-10);
}

void mdp::generalRandTransitions(arr& P, double uni, double noise){
  P = uni;
  rndUniform(P, .0, noise, true);
  tensorCondNormalize(P, 1);
  tensorCheckCondNormalization(P, 1, 1e-10);
}

void mdp::generalDirichlet(arr& P){
  rndNegLogUniform(P, .0, 1., false);
  tensorCondNormalize(P, 1);
  tensorCheckCondNormalization(P, 1, 1e-10);
}

void mdp::generalNodeTransition(arr& P, double uni, double noise, double stay){
  CHECK_EQ(P.dim(0),P.dim(P.nd-1), "");
  uint n=P.d0;
  uint m=P.N/(n*n);
  uint i, j;
  uintA org_dim;
  org_dim = P.dim();
  P.reshape(n, m, n);
  P = uni;
  rndUniform(P, .0, noise, true);
  for(i=0; i<n; i++) for(j=0; j<m; j++) P(i, j, i)+=stay;
  tensorCondNormalize(P, 1);
  tensorCheckCondNormalization(P, 1, 1e-10);
  P.reshape(org_dim);
}

/// .
void mdp::generalNodeTransitions(arr& P0y0, double uni, double noise, double stay){
  CHECK(P0y0.nd>=3 && P0y0.d0==P0y0.d[P0y0.nd-1], "array must be sized before this...");
  uintA oldDim;
  oldDim=P0y0.dim();
  uint d0=P0y0.d0;
  uint dy=P0y0.N/(d0*d0);
  P0y0.resize(d0, dy, d0);
  uint i, j;
  P0y0 = uni;
  rndUniform(P0y0, .0, noise, true);
  for(i=0; i<d0; i++) for(j=0; j<dy; j++) P0y0(i, j, i)+=stay;
  tensorCondNormalize(P0y0, 1);
  tensorCheckCondNormalization(P0y0, 1, 1e-10);
  P0y0.reshape(oldDim);
}

/// .
void mdp::generalNode0Transition(arr& P01y0,
                                 uint d0, uint d1, uint dy,
                                 double uni, double noise, double stay){
  P01y0.resize(TUP(d0, d1, dy, d0));
  uint i, j, k;
  P01y0 = uni;
  rndUniform(P01y0, .0, noise, true);
  for(i=0; i<d0; i++) for(j=0; j<d1; j++) for(k=0; k<dy; k++) P01y0.elem(TUP(i, j, k, i))+=stay;
  tensorCondNormalize(P01y0, 1);
  tensorCheckCondNormalization(P01y0, 1, 1e-10);
}

/// .
void mdp::generalNode1Transition(arr& P1y01,
                                 uint d1, uint dy, uint d0,
                                 double uni, double noise, double stay){
  P1y01.resize(TUP(d1, dy, d0, d1));
  uint i, j, k;
  P1y01 = uni;
  rndUniform(P1y01, .0, noise, true);
  for(i=0; i<d1; i++) for(j=0; j<dy; j++) for(k=0; k<d0; k++) P1y01.elem(TUP(i, j, k, i))+=stay;
  tensorCondNormalize(P1y01, 1);
  tensorCheckCondNormalization(P1y01, 1, 1e-10);
}


//===========================================================================
//
// strict hierarchies
//

//initialize the deterministic parts of an hierarchical controller (concerning the exit states)
void mdp::initHierarchical(uint exits, FSC_lev2& fsc){
  uint i;
  uint d0=fsc.P01y0.d0, d1=fsc.P1y01.d0, dy=fsc.P1y01.d1;
  cout <<"building hierarchical controller: d0=" <<d0 <<" d1=" <<d1 <<" dy=" <<dy
       <<" -- assuming the last " <<exits<<" level-0-nodes to be exits" <<endl;
       
  fsc.hierarchical=true;
  
  fsc.Pe0.resize(2, d0);         //P(e|n0)  is fixed!
  fsc.PE_01.resize(d0, d1);      //P(n0'|n1'; exit=true)
  fsc.P_E_0y0.resize(d0, dy, d0); //P(n0'|y', n0; exit=false)
  fsc.PE_1y1.resize(d1, dy, d1);  //P(n1'|y', n1; exit=true)
  fsc.P_E_11.resize(d1, d1);     //P(n1'|n1; exit=false) is fixed!
  
  fsc.Pe0=1.;    //no constraints on which n0 trigger an exit state?
  for(i=0; i<d0; i++){
    if(i>=d0-exits){ fsc.Pe0(1, i)=1.; fsc.Pe0(0, i)=0.; } else           { fsc.Pe0(1, i)=0.; fsc.Pe0(0, i)=1.; }
  }
  tensorCheckCondNormalization(fsc.Pe0, 1, 1e-10);
  
  fsc.P_E_11.setId();
  
  generalRandTransitions(fsc.PE_01, 1., 1.);
  generalRandTransitions(fsc.P_E_0y0, 1., 1.);
  generalRandTransitions(fsc.PE_1y1, 1., 1.);
  
  hierarchicalToFactors(fsc);
}

//compute the two-level factors (P01y0 and P1y01) from the hierarchical controller
void mdp::hierarchicalToFactors(FSC_lev2& fsc){
  CHECK(fsc.hierarchical, "");
  uint i, j, k, l;
  uint d0=fsc.P01y0.d0, d1=fsc.P1y01.d0, dy=fsc.P1y01.d1;
  
  arr P01ye0;
  P01ye0.resize(TUP(d0, d1, dy, 2, d0));
  for(i=0; i<d0; i++) for(j=0; j<d1; j++) for(k=0; k<dy; k++) for(l=0; l<d0; l++){
          P01ye0.elem(TUP(i, j, k, 0, l)) = fsc.P_E_0y0(i, k, l);
          P01ye0.elem(TUP(i, j, k, 1, l)) = fsc.PE_01(i, j);
        }
  tensorCheckCondNormalization(P01ye0, 1, 1e-10);
  
  arr P1ye1;
  P1ye1.resize(TUP(d1, dy, 2, d1));
  for(i=0; i<d1; i++) for(j=0; j<dy; j++) for(k=0; k<d1; k++){
        P1ye1.elem(TUP(i, j, 0, k)) = fsc.P_E_11(i, k);
        P1ye1.elem(TUP(i, j, 1, k)) = fsc.PE_1y1(i, j, k);
      }
  tensorCheckCondNormalization(P1ye1, 1, 1e-10);
  
  tensorEquation(fsc.P1y01, P1ye1, TUP(0, 1, 4, 3), fsc.Pe0, TUP(4, 2), 1);
  tensorEquation(fsc.P01y0, P01ye0, TUP(0, 1, 2, 4, 3), fsc.Pe0, TUP(4, 3), 1);
  tensorCheckCondNormalization(fsc.P1y01, 1, 1e-10);
  tensorCheckCondNormalization(fsc.P01y0, 1, 1e-10);
}


//===========================================================================
//===========================================================================


void OutputAoDot(arr& Pa0, uint numObs, const char* problem, std::vector<std::string>& a,
                 std::vector<std::string>& o,
                 std::map<int, std::string>& a0){
  uint i, j;
  cout <<problem <<endl;
  if(!strcmp(problem, "chainOfChains3.POMDP.arr")){
    a.push_back("doNothing"); a.push_back("wetHands"); a.push_back("dryHands");
    a.push_back("soapyHands"); a.push_back("rinseHands"); a.push_back("waterOn");
    a.push_back("waterOff");
    //  ta: taps, w: water, so: soap, to: towel, si: sink, a: away
    o.push_back("ota"); o.push_back("fta"); o.push_back("ow"); o.push_back("fw");
    o.push_back("oso"); o.push_back("fso"); o.push_back("oto"); o.push_back("fto");
    o.push_back("osi"); o.push_back("fsi"); o.push_back("oa"); o.push_back("fa");
  } else if(!strcmp(problem, "chainOfChains3.POMDP.arr")){
    a.push_back("A"); a.push_back("B"); a.push_back("C"); a.push_back("D");
    o.push_back("1");
  } else if(!strcmp(problem, "cheese-taxi-fully-observable.POMDP.arr")){
    a.push_back("N"); a.push_back("S"); a.push_back("E"); a.push_back("W");
    a.push_back("query"); a.push_back("pickup"); a.push_back("putdown");
    o.push_back("s0"); o.push_back("s1"); o.push_back("s2"); o.push_back("s3");
    o.push_back("s4"); o.push_back("s5"); o.push_back("s6"); o.push_back("s7");
    o.push_back("s8"); o.push_back("s9"); o.push_back("s10"); o.push_back("s0d1");
    o.push_back("s1d1"); o.push_back("s2d1"); o.push_back("s3d1"); o.push_back("s4d1");
    o.push_back("s5d1"); o.push_back("s6d1"); o.push_back("s7d1"); o.push_back("s8d1");
    o.push_back("s9d1"); o.push_back("s10d1"); o.push_back("s0d2"); o.push_back("s1d2");
    o.push_back("s2d2"); o.push_back("s3d2"); o.push_back("s4d2"); o.push_back("s5d2");
    o.push_back("s6d2"); o.push_back("s7d2"); o.push_back("s8d2"); o.push_back("s9d2");
    o.push_back("s10d2"); o.push_back("d0"); o.push_back("d1"); o.push_back("d2");
  }else{ // Simply put placeholders
    for(i=0; i<Pa0.d0; ++i){
      std::stringstream s;
      s <<i;
      a.push_back(s.str());
    }
    for(i=0; i<numObs; ++i){
      std::stringstream s;
      s <<i;
      o.push_back(s.str());
    }
  }
  
  // Map Nodes to Actions
  for(j=0; j<Pa0.d1; ++j){ // j is the node
    for(i=0; i<Pa0.d0; ++i){ // i is the action
      double value= Pa0(i, j);
      if(value > 0.01){
        if(value > 0.99) value = 1.0;
        //std::stringstream ss; ss<<value; std::string str; ss>>str;
        std::stringstream ss; ss<<j; std::string str; ss>>str;
        //a0[j] = a[i]+"_"+str+"_"+a0[j];
        a0[j] = str+a[i]+"_"+a0[j];
      }
    }
  }
}

// Outputs a .dot file.
// To convert it to a postscript: $ dot -Tps filename.dot -o filename.ps
void mdp::OutputDotH(const char* filename, const char* problem, arr& Pa0, arr& Pe0, arr& PE_01, arr& P_E_0y0, arr& PE_1y1,
                     arr& P_E_11){
  uint i, j;
  std::vector<std::string> a, o;
  std::map<int, std::string> a0;
  OutputAoDot(Pa0, P_E_0y0.d1, problem, a, o, a0);
  
  cout <<"Outputing graph: " <<filename <<endl;
  // Output graph
  ofstream out(filename);
  out <<"digraph " <<"ML_Controller_" <<Pa0.d1 <<"_" <<PE_1y1.d0 <<" {\n"
  <<"\trankdir=TB;\n";
  
  // Level-0
  out <<"\tsubgraph cluster_0 {\n"
  <<"\t\tstyle=filled;\n"
  <<"\t\tcolor=lightgrey;\n"
  <<"\t\tnode [style=filled, color=white];\n\t\t{rank = same;";
  for(i=0; i<Pa0.d1; ++i){
    // check if it's an exit state
    // if so, output it as a diamond....
    if(Pe0(1, i)<0.01)
      out <<"N0_" <<a0[i] <<"; ";
  }
  out <<"};\n"
  <<"\t\tlabel=\"Level 0\"\n";
  // Level-Exits
  out <<"\tsubgraph cluster_2 {\n"
  <<"\t\tstyle=filled;\n"
  <<"\t\tcolor=lightgrey;\n"
  <<"\t\tnode [style=filled, color=white, shape=diamond];";
  for(i=0; i<Pa0.d1; ++i){ if(Pe0(1, i)>0.99) out <<"E_" <<a0[i] <<" ";}
  out <<";\n"
  <<"\t\tlabel=\"Exit\"\n"
  <<"\t}\n"
  <<"\n\t}\n";
  
  // Level-1
  out <<"\tsubgraph cluster_1 {\n"
  <<"\t\tstyle=filled;\n"
  <<"\t\tcolor=lightgrey;\n"
  // This is correct but the rank=same makes dot crash once in a while, this
  // is a dot bug.
  <<"\t\tnode [style=filled, color=white];\n\t\t{rank = same;"; for(i=0; i<PE_01.d1; ++i){out <<"N1_" <<i <<"; ";}
  // <<"\t\tnode [style=filled, color=white]; "; for(i=0;i<PE_01.d1;++i){out <<"N1_" <<i <<"; ";}
  out <<"};\n"
  <<"\t\tlabel=\"Level 1\"\n"
  <<"\t}\n";
  
  
  // P_E_0y0; E is false
  for(i=0; i<P_E_0y0.d0; ++i) // incoming node is i
    for(uint k=0; k<P_E_0y0.d2; ++k){ // outgoing node is k
      bool obs = false;
      if(!(Pe0(1, k)>0.99))
        for(uint j=0; j<P_E_0y0.d1; ++j){ // outgoing obs. is j
          // Merge observations
          if(P_E_0y0(i, j, k) > 0.01){
            if(!obs){
              out <<"\tN0_" <<a0[k] <<((Pe0(1, i)>0.99)? "-> E_":"-> N0_") <<a0[i] <<" [label=\"";
              obs=true;
            }
            if(P_E_0y0(i, j, k) < 0.99){
              out <<"" <<o[j] <<":" <<std::setprecision(2)
              <<P_E_0y0(i, j, k) <<", ";
            }
          }
        }
      if(obs) out <<"\", arrowsize=0.3, labelfontsize=9.0, penwidth=0.5];\n";
    }
    
  // PE_1y1; E is true
  for(i=0; i<PE_1y1.d0; ++i) // incoming node is i
    for(uint k=0; k<PE_1y1.d2; ++k){ // outgoing node is k
      bool obs = false;
      for(uint j=0; j<PE_1y1.d1; ++j){ // outgoing obs. is j
        // Merge observations
        if(PE_1y1(i, j, k) > 0.01){
          if(!obs){
            out <<"\tN1_" <<k <<"-> N1_" <<i <<" [label = \"";
            obs=true;
          }
          if(PE_1y1(i, j, k) < 0.99){
            out <<"" <<o[j] <<":" <<std::setprecision(2)
            <<PE_1y1(i, j, k) <<", ";
          }
        }
      }
      if(obs) out <<"\", weight=50, arrowsize=0.3, labelfontsize=5.0, penwidth=0.5];\n";
    }
  // PE_01; E is true
  for(i=0; i<PE_01.d0; ++i) // incoming node is i
    for(j=0; j<PE_01.d1; ++j){ // outgoing node is j
      double value=PE_01(i, j);
      if(value > 0.01){
        if(value > 0.99)
          out <<"\tN1_" <<j <<((Pe0(1, i)>0.99)? "-> E_":"-> N0_") <<a0[i] <<" [style=dashed, arrowsize=0.3, labelfontsize=5.0, penwidth=0.5];\n";
        else
          out <<"\tN1_" <<j <<((Pe0(1, i)>0.99)? "-> E_":"-> N0_") <<a0[i] <<" [label=\"" <<std::setprecision(2) <<value <<"\", style=dashed, arrowsize=0.3, labelfontsize=5.0, penwidth=0.5];\n";
      }
    }
    
  out <<"}\n";
  out.close();
  cout <<"Done" <<endl;
}


// Output the policy to dot format
// Currently only works for 1-level policies
void mdp::OutputDot(const char* filename, const char* problem, arr& Pa0, arr& P0y0){
  std::vector<std::string> a, o;
  std::map<int, std::string> a0;
  OutputAoDot(Pa0, P0y0.d1, problem, a, o, a0);
  
  uint i, j;
  // Output graph
  ofstream out(filename);
  out <<"digraph " <<"ML_Controller_" <<Pa0.d1 <<"_" <<0 <<" {\n"
  <<"\trankdir=LR;\n"
  <<"\tsize=\"8, 5\"\n"
  <<"\tnode [shape = circle];"; for(i=0; i<P0y0.d0; ++i){out <<"N_" <<a0[i] <<" ";}
  out <<";\n"
  <<"\tnode [shape = circle];\n";
  for(i=0; i<P0y0.d0; ++i) // incoming node is i
    for(uint k=0; k<P0y0.d2; ++k){ // outgoing node is k
      bool obs = false;
      for(j=0; j<P0y0.d1; ++j){ // outgoing obs. is j
        // Merge observations
        if(P0y0(i, j, k) > 0.01){
          if(!obs){
            out <<"\tN_" <<a0[k] <<"-> N_" <<a0[i] <<"[ label = \"";
            obs=true;
          }
          if(P0y0(i, j, k) > 0.99) out <<"" <<o[j] <<", ";
          else {out <<"" <<o[j] <<":" <<std::setprecision(2) <<P0y0(i, j, k) <<", ";}
        }
      }
      if(obs) out <<"\"];\n";
    }
    
  out <<"}\n";
  out.close();
}


// This function simply outputs the policy in an AMPL readable format.
// It's useful to debug the value of the policy.
void mdp::SetAMPLPolicy(arr& P0y0, arr& Pa0){
  uint i;
  
  // re-set the number of variables (to be independent of the .data file)
  cout <<"set Nodes :=";
  for(i=0; i < P0y0.d0; i++){ cout <<" " <<i; }
  cout <<";\n";
  
  arr P0ay0;
  P0ay0.resize(TUP(P0y0.d0, Pa0.d0, P0y0.d1, P0y0.d2));
  tensorEquation(P0ay0, P0y0, TUP(0, 2, 3), Pa0, TUP(1, 3), 0);
  tensorCheckCondNormalization(P0ay0, 1, 1e-10);
  
  for(i=0; i < P0ay0.d0; i++){
    for(uint j=0; j < P0ay0.d1; j++){
      for(uint l=0; l < P0y0.d2; l++){ // switch the order of l and k since
        for(uint k=0; k < P0ay0.d2; k++){ // the AMPL input is actually P0a0y
          cout <<"\nlet x[" <<i <<", " <<j <<", " <<l <<", " <<k
               <<"] := "    <<P0ay0.elem(TUP(i, j, k, l))/*return4DNumber(P0ay0, i, j, k, l)*/ <<";";
        }
      }
    }
  }
  cout <<endl;
}


// @}
