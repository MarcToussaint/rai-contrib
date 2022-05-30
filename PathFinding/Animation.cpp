#include <unordered_map>
#include "Animation.h"

void rai::Animation::AnimationPart::write(std::ostream& os) const{
  os <<'\n';
  arr tmp(1);
  tmp = start;
  tmp.writeTagged(os, "start");
  os <<'\n';
  frameIDs.writeTagged(os, "frameIDs");
  os <<'\n';
  frameNames.writeTagged(os, "frameNames");
  os <<'\n';
  frameCols.writeTagged(os, "frameColors");
  os <<'\n';
  X.writeTagged(os, "poses");
}

void rai::Animation::AnimationPart::read(std::istream& is) {
  arr tmp;
  tmp.readTagged(is, "start");
  start = tmp(0);
  frameIDs.readTagged(is, "frameIDs");
  frameNames.readTagged(is, "frameNames");
  frameCols.readTagged(is, "frameColors");
  X.readTagged(is, "poses");
}

void rai::Animation::read(istream& is) {
  A.readTagged(is, "animation");
}

uint rai::Animation::getT(){
  uint T=0;
  for(auto& a: A) if(a.X.d0 + a.start>T) T=a.X.d0 + a.start;
  return T;
}

void rai::Animation::setToTime(rai::Configuration& C, const double t, const double tIgnoreBefore) const{
  // TODO: make sure that the ordering is used to speed this up
  // - we know for a lot of things that they have not started yet

  auto start_time = std::chrono::high_resolution_clock::now();

  // datastructure to collect the positions
  std::unordered_map<uint, arr> poseMap;

  for(int i=A.d0-1; i>=0; --i){
    // this assumes that the newest animation part always has priority
    const auto &a = A(i);
    const double start = a.start;

    if (start < 0 || a.X.d0 == 0) {
      //LOG(-1) <<"Warning: start < 0";
      continue;
    }
    /*if (std::ceil(t) > start + a.X.d0 + 30){
      //LOG(-1) <<  "Warning: time > anim_part";
      continue;
    }*/

    // enables us to reduce the things we have to look at if we are certain that
    // animation parts before a certain time can be disregarded
    if (tIgnoreBefore > start + a.X.d0){
      continue;
    }

    arr poses;
    if(t < start){
      poses = a.X[0];
    }
    else if (t >= start + a.X.d0-1){
      poses = a.X[a.X.d0-1];
    }
    else{
      // we might want to interpolate here - quaternion interpolation could be a problem
      poses = a.X[uint(std::floor(t - start))];
    }

    uint cnt = 0;
    for (auto id: a.frameIDs){
      if (poseMap.count(id) == 0){
        poseMap[id] = poses[cnt];
      }
      cnt++;
    }
  }

  arr poses(0, 7);
  uintA ids;
  for (auto e: poseMap){
    //if (length(C.getFrameState(uintA({e.first})) - e.second) < 1e-3) {continue;}

    ids.append(e.first);
    poses.append(e.second);
  }

  C.setFrameState(poses, ids);

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end_time - start_time ).count();

  /*{
    std::ofstream f;
    f.open("./out/Anim_time.txt", std::ios_base::app);
    f << duration << std::endl;
  }*/
}

void rai::Animation::play(rai::Configuration& C, bool pause){
  rai::ConfigurationViewer V;
  V.setConfiguration(C);
  uint T=getT();

  for(uint t=0;t<T;t++){
    setToTime(C, t);
    V.setConfiguration(C, STRING("Animation t:" <<t), pause);
    rai::wait(.1);
  }
}

void rai::Animation::write(ostream& os) const {
  A.writeTagged(os, "animation");
}

