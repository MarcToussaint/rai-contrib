#pragma once

#include <Kin/kin.h>
#include <Kin/viewer.h>

namespace rai {

struct Animation{
  struct AnimationPart{
    StringA frameNames;
    uintA frameIDs;
    arr frameCols;
    arr X;
      
    double start = 0.;

    void write(ostream& os) const;
    void read(istream& is);
  };

  rai::Array<AnimationPart> A;

  void write(ostream& os) const;
  void read(istream& is);

  uint getT();

  void setToTime(rai::Configuration& C, const double t, const double tIgnoreBefore=0) const;
  uintA getActiveFramesAtTime(const double t);

  void play(rai::Configuration& C, bool pause);
};
stdPipes(Animation::AnimationPart)
stdPipes(Animation)

}

using AnimationPtr = std::shared_ptr<rai::Animation>;
