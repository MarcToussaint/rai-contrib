#pragma once

#include <Kin/kin.h>
#include <Core/thread.h>
#include <openvr/openvr.h>

namespace rai{

struct ViveController : Thread {
  
  RAI_PARAM("vivecontroller/", double, filter, .9)

  ViveController();
  ~ViveController();

  void pull(rai::Configuration& C);

  void step();

  bool gripping = false;
  rai::Transformation pose;

private:
  double origin[3] = {0, 0, 0};

  std::mutex mux;
  vr::IVRSystem* vrSystem;
  vr::VRControllerState_t controllerState;
};

} //namespace
