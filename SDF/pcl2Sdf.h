#pragma once

#include <Geo/signedDistanceFunctions.h>
#include <Kin/cameraview.h>

struct PclSdf{
  SDF_GridData pixSdf;
  SDF_GridData carSdf;

  void getPixSdfBoundingBox(const uintA& seg, const floatA& depth, uint ID, uint padding=0);

  void rayFillPixSdf(const floatA& depth, const uintA& seg, uint ID, uint depthResolution);


  void getCartesianBoundingBox(const arr& fxypxy);

  void resampleCartesianSdf(const arr& fxypxy, uint resolution, double pixScale);

  void fillCartSdf(const arr& fxypxy, const floatA& depth, const uintA& seg, uint ID, uint resolution);
};


void fillSDF(SDF_GridData& sdf, bool fwdbwd, bool depthEuclid);
