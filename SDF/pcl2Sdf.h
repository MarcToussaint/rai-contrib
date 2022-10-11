#pragma once

#include <Geo/signedDistanceFunctions.h>
#include <Kin/cameraview.h>

struct PclSdf{
  SDF_GridData pixSdf;
  SDF_GridData carSdf;

  void getPixSdfBoundingBox(const uintA& seg, const floatA& depth, uint ID);

  void rayFillPixSdf(const floatA& depth, const uintA& seg, uint ID, uint depthResolution);


  void getCartesianBoundingBox(const arr& fxypxy);

  void resampleCartesianSdf(const arr& fxypxy, uint resolution, double pixScale);
};


void fillSDF(SDF_GridData& sdf, bool fwdbwd, bool depthEuclid);
void sdf_smooth(SDF_GridData& sdf, uint width=3, uint iters=2);
