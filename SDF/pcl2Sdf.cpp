#include "pcl2Sdf.h"

#include <Geo/depth2PointCloud.h>

#include <queue>

void PclSdf::getPixSdfBoundingBox(const uintA& seg, const floatA& depth, uint ID){
    pixSdf.lo.resize(3) = 1000.;
    pixSdf.up.resize(3) = -1000.;
    for(uint py=0;py<seg.d0;py++) for(uint px=0;px<seg.d1;px++){
        if(seg(py, px)==ID){
            double dep = depth(py, px);
            if(px < pixSdf.lo(0)) pixSdf.lo(0) = px;
            if(py < pixSdf.lo(1)) pixSdf.lo(1) = py;
            if(dep < pixSdf.lo(2)) pixSdf.lo(2) = dep;
            if(px > pixSdf.up(0)) pixSdf.up(0) = px;
            if(py > pixSdf.up(1)) pixSdf.up(1) = py;
            if(dep > pixSdf.up(2)) pixSdf.up(2) = dep;
        }else{
            //        segcol(py, px, {}) = 0;
        }
    }
    pixSdf.lo(2) -= .01; pixSdf.up(2) += .01;
//    pixSdf.up(0) += 10; pixSdf.up(1) += 10;
//    pixSdf.lo(0) -= 10; pixSdf.lo(1) -= 10;
//    if(pixSdf.up(0)>seg.d1) pixSdf.up(0)=seg.d1;
//    if(pixSdf.up(1)>seg.d0) pixSdf.up(1)=seg.d0;
//    if(pixSdf.lo(0)<0) pixSdf.lo(0)=0;
//    if(pixSdf.lo(1)<0) pixSdf.lo(1)=0;
}

void PclSdf::rayFillPixSdf(const floatA& depth, const uintA& seg, uint ID, uint depthResolution){
    //-- create ray grid
    pixSdf.gridData.resize(pixSdf.up(0)-pixSdf.lo(0)+1, pixSdf.up(1)-pixSdf.lo(1)+1 , depthResolution);
    float rad = pixSdf.up(2)-pixSdf.lo(2);
    pixSdf.gridData = rad; //default initialization -> ``maximal distance to everything''

    float scale = (pixSdf.up(2)-pixSdf.lo(2))/double(pixSdf.gridData.d2-1);

    for(uint i=0;i<pixSdf.gridData.d0;i++) for(uint j=0;j<pixSdf.gridData.d1;j++){
        uint py = pixSdf.lo(1)+j;
        uint px = pixSdf.lo(0)+i;
        if(seg(py, px)==ID){
          double dep = depth(py, px);
          rai::clip(dep, pixSdf.lo(2), pixSdf.up(2));

          double thick = .4*(pixSdf.up(2)-dep);
          dep -= pixSdf.lo(2);
          for(uint k=0;k<pixSdf.gridData.d2;k++){
              double d = -( double(k)*scale - dep );
              if(d<-thick) d = -2.*thick - d;
              pixSdf.gridData(i, j, k) = d;
          }
        }
    }
    //  sdfd.animateSlices(sdfd.lo-.1, sdfd.up+.1, -1.);
}

void PclSdf::getCartesianBoundingBox(const arr& fxypxy){
    arr corners(8,3);
    for(uint i=0;i<8;i++) for(uint l=0;l<3;l++){
        corners(i, l) = ((i & (1<<l)) ? pixSdf.up(l):pixSdf.lo(l));
    }
    for(uint i=0;i<8;i++) depthData2point(&corners(i,0), fxypxy.p);
#if 1 //outer
    carSdf.lo = min(corners, 0);
    carSdf.up = max(corners, 0);
#else //inner
    corners = ~corners;
//    cout <<pixSdf.lo <<endl <<pixSdf.up <<endl <<corners <<endl;
    corners[0].sort();
    corners[1].sort();
    corners[2].sort();
    carSdf.lo = corners.col(3);
    carSdf.up = corners.col(4);
//    cout <<carSdf.lo <<endl <<carSdf.up <<endl <<corners <<endl;
#endif
}

void PclSdf::resampleCartesianSdf(const arr& fxypxy, uint resolution, double pixscale){
    uint N = resolution;
    arr X = grid(carSdf.lo, carSdf.up, {N,N,N}); //grid in cam coord
    for(uint i=0;i<X.d0;i++) point2depthData(&X(i,0), fxypxy.p);
    if(pixscale>0.){
        for(uint i=0;i<X.d0;i++){ X(i,0)*= pixscale; X(i,1)*= pixscale; }
        pixSdf.lo(0) *= pixscale; pixSdf.lo(1) *= pixscale;
        pixSdf.up(0) *= pixscale; pixSdf.up(1) *= pixscale;
    }

    carSdf.gridData = rai::convert<float>(pixSdf.eval(X)).reshape(N+1,N+1,N+1);

    if(pixscale>0.){
        pixSdf.lo(0) /= pixscale; pixSdf.lo(1) /= pixscale;
        pixSdf.up(0) /= pixscale; pixSdf.up(1) /= pixscale;
   }
}

bool updateVertex(SDF_GridData& sdf, const std::array<int,3>& p, const rai::Array<std::array<int,3>>& neighbors, const arr& dists){
    bool touched=false;
    float& dp = sdf.gridData(p[0], p[1], p[2]);
    for(uint i=0;i<neighbors.d0;i++){
        std::array<int,3> n = neighbors.p[i];
        for(uint j=0;j<3;j++) n[j] += p[j];
        if(n[0]<0 || n[1]<0 || n[2]<0 ||
                n[0]>=(int)sdf.gridData.d0 || n[1]>=(int)sdf.gridData.d1 || n[2]>=(int)sdf.gridData.d2) continue;

        float& dn = sdf.gridData(n[0], n[1], n[2]);
        float s = dists.p[i];

        if(dn>0. && dp>dn+s){ dp = dn+s; touched=true; }
        if(dn<0. && dp<dn-s){ dp = dn-s; touched=true; }
        //      if(dp==1.f || dp==-1.f){ touched=true; }
    }

    return touched;
}

void fillSDF(SDF_GridData& sdf, bool fwdbwd, bool depthEuclid){

    //grid stepsizes in 3 directions
    arr step = sdf.up-sdf.lo;
    step /= arr{(double)sdf.gridData.d0-1, (double)sdf.gridData.d1-1, (double)sdf.gridData.d2-1};
    if(depthEuclid){
        step(0) = step(2);
        step(1) = step(2);
    }

    //list of 8 neighbors
    arr dists;
    rai::Array<std::array<int,3>> neighbors;
    for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++){
        if(!i &&!j && !k) continue;
        neighbors.append( std::array<int,3>{i,j,k} );
        dists.append( length( arr{step(0)*i, step(1)*j, step(2)*k } ) );
    }

    if(fwdbwd){

        //fwd
        for(int k=sdf.gridData.d2;k--;){
            for(int i=0;i<(int)sdf.gridData.d0;i++) for(int j=0;j<(int)sdf.gridData.d1;j++){
                updateVertex(sdf, {i,j,k}, neighbors, dists);
            }
        }

        //bwd
        for(int k=0;k<(int)sdf.gridData.d2;k++){
            for(int i=0;i<(int)sdf.gridData.d0;i++) for(int j=0;j<(int)sdf.gridData.d1;j++){
                updateVertex(sdf, {i,j,k}, neighbors, dists);
            }
        }


    }else{

        //initialize queue
        std::queue<std::array<int,3>> queue;
        uintA dim=sdf.gridData.dim();
        uintA idx;
        for(uint i=0;i<sdf.gridData.N;i++){
            idx = getIndexTuple(i, dim);
            queue.push({(int)idx.p[0], (int)idx.p[1], (int)idx.p[2]});
        }

        while(queue.size()){
            //pop
            std::array<int,3> p = queue.front();
            queue.pop();

            //go over neighbors
            bool touched = updateVertex(sdf, p, neighbors, dists);

            if(touched==true){ //append all neighbors
                for(uint i=0;i<neighbors.d0;i++){
                    std::array<int,3> n = neighbors.p[i];
                    for(uint j=0;j<3;j++) n[j] += p[j];
                    if(n[0]<0 || n[1]<0 || n[2]<0 ||
                            n[0]>=(int)sdf.gridData.d0 || n[1]>=(int)sdf.gridData.d1 || n[2]>=(int)sdf.gridData.d2) continue;
                    queue.push(n);
                }
            }

            if(!(queue.size()%100000)){ cout <<p[0] <<' ' <<p[1] <<' ' <<p[2] <<" queue: " <<queue.size() <<' ' <<sdf.gridData(p[0], p[1], p[2]) <<endl; }
        }
    }
}

void sdf_smooth(SDF_GridData& sdf, uint width, uint iters){
    arr tmp;
    rai::copy(tmp, sdf.gridData);
    for(uint i=0;i<iters;i++){
        tmp = integral(tmp);
        tmp = differencing(tmp, width);
    }
    rai::copy(sdf.gridData, tmp);
}
