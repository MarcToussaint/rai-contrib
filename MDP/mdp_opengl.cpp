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

#include <Gui/opengl.h>
#include <Plot/plot.h>

//#ifdef RAI_GL

OpenGL *globalGL=NULL;

void mdp::showMaze(){
  byteA &maze=global_maze;
  byteA img(maze.d0, maze.d1, 3);
  uint x, y, dx=maze.d1, dy=maze.d0;
  for(x=0; x<dx; x++) for(y=0; y<dy; y++){
      if(maze(y, x)==0){ img(y, x, 0)=255; img(y, x, 1)=255; img(y, x, 2)=255; } else if(maze(y, x)==1){ img(y, x, 0)=0;   img(y, x, 1)=0;   img(y, x, 2)=0;   } else if(maze(y, x)==2){ img(y, x, 0)=255; img(y, x, 1)=0;   img(y, x, 2)=0;   } else if(maze(y, x)==3){ img(y, x, 0)=0;   img(y, x, 1)=255; img(y, x, 2)=0;   } else if(maze(y, x)==4){ img(y, x, 0)=0;   img(y, x, 1)=0;   img(y, x, 2)=255; } else HALT("strange global maze");
    }
  flip_image(img);
  static OpenGL *gl=NULL;
  if(!gl) gl=new OpenGL;
  gl->watchImage(img, true, 10);
}

void mix(byteA& A, const byteA& B, float f=.5){
  if(f>1.) f=1.; else if(f<0.) f=0.;
  A(0)=(1.-f)*A(0)+f*B(0);
  A(1)=(1.-f)*A(1)+f*B(1);
  A(2)=(1.-f)*A(2)+f*B(2);
}

void mdp::showAB(const arr& alpha, const arr& beta){
  byteA &maze=global_maze;
  byteA img(maze.d0, maze.d1, 3);
  uint x, y, dx=maze.d1, dy=maze.d0;
  for(x=0; x<dx; x++) for(y=0; y<dy; y++){
      if(maze(y, x)==0){ img(y, x, 0)=255; img(y, x, 1)=255; img(y, x, 2)=255; } else if(maze(y, x)==1){ img(y, x, 0)=0;   img(y, x, 1)=0;   img(y, x, 2)=0;   } else if(maze(y, x)==2){ img(y, x, 0)=0;   img(y, x, 1)=0;   img(y, x, 2)=255; } else if(maze(y, x)==3){ img(y, x, 0)=255; img(y, x, 1)=0;   img(y, x, 2)=0;   } else if(maze(y, x)==4){ img(y, x, 0)=0;   img(y, x, 1)=255; img(y, x, 2)=0;   } else HALT("strange global maze");
    }
  img.reshape(alpha.N, 3);
  double aM=alpha.max(), bM=beta.max();
  for(x=0; x<alpha.N; x++) if(alpha(x)) mix(img[x](), {(byte)0, 0, 255}, alpha(x)/aM);
  for(x=0; x<alpha.N; x++) if(beta(x)) mix(img[x](), {(byte)255, 0, 0}, beta(x)/bM);
  img.reshape(maze.d0, maze.d1, 3);
  flip_image(img);
  static OpenGL *gl=NULL;
  if(!gl) gl=new OpenGL;
  gl->watchImage(img, false, 10);
}

void mdp::plotPolicyAndValue(const arr& pi, const arr& V, const MDP& mdp, bool wait){
  plot()->Opengl(false);
  plot()->colors=false;
  plot()->drawBox=true;
  
  //plotModule.grid=true;
  
  //glDisplayRedBlue(V, global_maze.d0, global_maze.d1, false);
  arr tmp;
  uintA tmpI;
  plot()->Clear();
  tmp=~V;
  tmp *= .8/tmp.max();
  tmp.reshape(global_maze.d0, global_maze.d1);
  plot()->Surface(tmp);
  
  getMaxPxxMap(tmpI, mdp.Pxax, pi);
  tmpI.reshape(global_maze.d0, global_maze.d1);
  plot()->MatrixFlow(tmpI, .25);
  plot()->update(wait);
}

// #else //RAI_GL

// #include <Core/util.h>
// void mdp::showMaze(){ RAI_MSG("display only implemented when compiling with some GL"); }
// void mdp::showAB(const arr& alpha, const arr& beta){}
// void mdp::plotPolicyAndValue(const arr& pi, const arr& V, const MDP& mdp, bool wait){}
// void mdp::glDisplayGrey(const arr &x, uint d0, uint d1, bool wait, uint win){}
// void mdp::glDisplayRedBlue(const arr &x, uint d0, uint d1, bool wait, uint win){}

// #endif
