#ifndef MESH
#define MESH

#include "Vector.hpp"
#include "CRS.hpp"

class Mesh{
public:
    double hCoarse, hRefined, hReason;
    int matrixSize;
    int MatrixRefinedSize;
    Vector3i unityCoarse, unityRefined;

    int nzNum;
    std::vector<double> uO;
    // std::vector<double> uEGF;

    std::vector<Vector3> refined;

    double sigma;
    // double sigmaEGF;
    double deltaT;

    double * B;
    std::vector<int> rowPtr;
    int * colInd;
    double * val;
    // int * rowPtrEGF;
    // int * colIndEGF;
    // double * valEGF;
    std::vector<Vector3> pos;

};

class Mesh3D : public Mesh{
public:
    Mesh3D(Vector3 domain, double oxgD, double egfD, double oBorder, double hCoarse=10.0, double hRefined=1.0, double deltaT=1.0);
};

class Mesh2D : public Mesh{
public:
    Mesh2D(Vector3 domain, double oxgD, double egfD, double oBorder, double hCoarse=10.0, double hRefined=2.5, double deltaT=1.0);
};

#endif /* end of include guard: MESH */
