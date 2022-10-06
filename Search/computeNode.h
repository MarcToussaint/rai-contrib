#pragma once

#include <Core/util.h>

namespace rai {

  struct ComputeNode {
    rai::String name;
    bool isComplete = false;
    bool isTerminal = false;

    virtual double compute(){
        double time = -rai::cpuTime();
        computeStep();
        time += rai::cpuTime();
        return time;
    }
    virtual void computeStep(){ HALT("need to overload"); }

    virtual int getNumDecisions() = 0;
    virtual std::shared_ptr<ComputeNode> getNewChild(uint i) = 0;

    virtual double costHeuristic(){ return 0.; }
    virtual double effortHeuristic(){ return 0.; }        //expected effort-to-go (FULL DOWN-STREAM TO LEAF NODE)
    virtual double sample(){ HALT("need to overload"); }  //get a value (at a leaf)

    virtual double groundTruthMean(){ HALT("need to overload"); }
    virtual void write(ostream& os) const{ os <<"ComputeNode"; }
  };
  stdOutPipe(ComputeNode)

}
