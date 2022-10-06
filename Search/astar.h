/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include <Logic/treeSearchDomain.h>
#include <Algo/priorityQueue.h>

//===========================================================================

namespace rai {

struct AStar {
  typedef std::shared_ptr<TreeSearchNode> NodeP;
  rai::Array<NodeP> mem;
  NodeP root;
  PriorityQueue<TreeSearchNode*> queue;
  rai::Array<TreeSearchNode*> solutions;
  uint iters=0;
  int verbose=1;

  AStar(const std::shared_ptr<TreeSearchNode>& _root);

  bool step();
  void run();
  void report();
};

} //namespace

//===========================================================================

