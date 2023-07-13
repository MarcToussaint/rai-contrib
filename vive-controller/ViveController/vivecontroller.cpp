#include "vivecontroller.h"

#ifdef RAI_VIVECONTROLLER

#include <Kin/frame.h>
double toRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

// Function to convert Euler angles to rotation matrix
void eulerToRotationMatrix(double roll, double pitch, double yaw, double (&rotationMatrix)[9]) {
    // Convert Euler angles to radians
    roll = toRadians(roll);
    pitch = toRadians(pitch);
    yaw = toRadians(yaw);

    // Calculate the elements of the rotation matrix
    double cosRoll = cos(roll);
    double sinRoll = sin(roll);
    double cosPitch = cos(pitch);
    double sinPitch = sin(pitch);
    double cosYaw = cos(yaw);
    double sinYaw = sin(yaw);

    rotationMatrix[0] = cosYaw * cosPitch;
    rotationMatrix[1] = cosYaw * sinPitch * sinRoll - sinYaw * cosRoll;
    rotationMatrix[2] = cosYaw * sinPitch * cosRoll + sinYaw * sinRoll;

    rotationMatrix[3] = sinYaw * cosPitch;
    rotationMatrix[4] = sinYaw * sinPitch * sinRoll + cosYaw * cosRoll;
    rotationMatrix[5] = sinYaw * sinPitch * cosRoll - cosYaw * sinRoll;

    rotationMatrix[6] = -sinPitch;
    rotationMatrix[7] = cosPitch * sinRoll;
    rotationMatrix[8] = cosPitch * cosRoll;
}

namespace rai{

    ViveController::ViveController() : Thread("ViveControllerThread", 0.){
        // Initialize OpenVR
        vrSystem = vr::VR_Init(nullptr, vr::VRApplication_Scene);
        if (!vrSystem) {
            LOG(-2) <<"Failed to initialize OpenVR.";
        }
        threadLoop();
    }

    ViveController::~ViveController(){
        threadClose();
        vr::VR_Shutdown();
    }

    void ViveController::pull(rai::Configuration& C) {
        rai::Frame * f = C.getFrame("viveLeft");
        if(!f){
            f = C.addFrame("viveLeft");
            f->setShape(rai::ST_marker, {.5});
        }
        {
            std::lock_guard<std::mutex> lock(mux);
            f->setPose(pose);
            if (gripping) f->setColor(arr{0, 1, 1, 1});
            else f->setColor(arr{0.5, 0.5, 0.5, 1});
        }
    }
    

  void ViveController::step(){
      std::lock_guard<std::mutex> lock(mux);

      // Update controller poses
      vr::TrackedDevicePose_t controllerPoses[vr::k_unMaxTrackedDeviceCount];
      vrSystem->GetDeviceToAbsoluteTrackingPose(vr::TrackingUniverseStanding, 0.0f, controllerPoses, vr::k_unMaxTrackedDeviceCount);

      // Iterate through tracked devices (controllers)
      for(vr::TrackedDeviceIndex_t deviceIndex = 0; deviceIndex < vr::k_unMaxTrackedDeviceCount; ++deviceIndex) {
          if (controllerPoses[deviceIndex].bPoseIsValid && vrSystem->GetTrackedDeviceClass(deviceIndex) == vr::TrackedDeviceClass_Controller) {
              // Get controller pose
              const vr::TrackedDevicePose_t& controllerPose = controllerPoses[deviceIndex];
              const vr::HmdMatrix34_t& poseMatrix = controllerPose.mDeviceToAbsoluteTracking;
              pose.pos.set(-poseMatrix.m[0][3]-origin[0], poseMatrix.m[2][3]-origin[1], poseMatrix.m[1][3]-origin[2]);

              if (deviceIndex != vr::k_unTrackedDeviceIndexInvalid) {
                  if (vrSystem->GetControllerState(deviceIndex, &controllerState, sizeof(controllerState))) {
                      // Check the button state
                      if (controllerState.ulButtonPressed >> 32 & 1) {
                          origin[0] = -poseMatrix.m[0][3];
                          origin[1] = poseMatrix.m[2][3];
                          origin[2] = poseMatrix.m[1][3];
                      }
                      gripping = (controllerState.ulButtonPressed & (1ULL << vr::k_EButton_Grip));
                  }
              }

              double alpha = atan2(poseMatrix.m[1][2], poseMatrix.m[2][2]);
              double gamma = asin(poseMatrix.m[0][2]);
              double beta = -atan2(poseMatrix.m[0][1], poseMatrix.m[0][0]);

              pose.rot.setRpy(alpha, beta, gamma);



              // Convert Euler angles to rotation matrix
              //double rotationMatrix[9];
              //eulerToRotationMatrix(alpha, beta, gamma, rotationMatrix);

              std::cout << "OpenVR Rotation Matrix:\n";
              for (int i = 0; i < 3; i++) {
                  for (int j = 0; j < 3; j++) {
                      std::cout << poseMatrix.m[i][j] << " ";
                  }
                  std::cout << std::endl;
              }

              // Print the rotation matrix



               /*double TargetMat[9] ={0,0,0,0,0,0,0,0,0};
              for (int i = 0; i < 3; i++) { // Iterate over rows of p
                  for (int j = 0; j < 3; j++) { // Iterate over columns of p
                      for (int k = 0; k < 3; k++) { // Iterate over columns of x
                          TargetMat[i * 3 + j] += poseMatrix.m[i][k] * X[k * 3 + j]; // Calculate dot product
                      }
                  }
              }*/
              // Calculate the elements of the rotation matrix
              double cosRoll = cos(alpha);
              double sinRoll = sin(alpha);
              double cosPitch = cos(beta);
              double sinPitch = sin(beta);
              double cosYaw = cos(gamma);
              double sinYaw = sin(gamma);

              rai::ArrayDouble rotationMatrix = arr{cosYaw * cosPitch, cosYaw * sinPitch * sinRoll - sinYaw * cosRoll, cosYaw * sinPitch * cosRoll + sinYaw * sinRoll,
                                                   sinYaw * cosPitch, sinYaw * sinPitch * sinRoll + cosYaw * cosRoll, sinYaw * sinPitch * cosRoll - cosYaw * sinRoll,
                                                   -sinPitch, cosPitch * sinRoll, cosPitch * cosRoll};

              /*std::cout << "Target Rotation Matrix:\n";
              for (int i = 0; i < 3; i++) {
                  for (int j = 0; j < 3; j++) {
                      std::cout << rotationMatrix.elem(i*3+j) << " ";
                  }
                  std::cout << std::endl;
              }*/

              //pose.rot.setRpy(alpha, beta, gamma);
              //pose.rot.setMatrix(arr{TargetMat[0],TargetMat[1],TargetMat[2],TargetMat[3],TargetMat[4],TargetMat[5],TargetMat[6],TargetMat[7],TargetMat[8]});
              //pose.rot.setMatrix(rotationMatrix);
          }
      }
   }
} //namespace

#else

rai::ViveController::ViveController() : Thread("ViveControllerThread") { NICO }
rai::ViveController::~ViveController(){ NICO }
void rai::ViveController::pull(rai::Configuration& C){ NICO }
void rai::ViveController::step(){ NICO }

#endif
