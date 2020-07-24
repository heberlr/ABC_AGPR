#ifndef MACRO
#define MACRO

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "Cell.hpp"
#include "Mesh.hpp"
#include "Frame.hpp"
#include "ConfigHandler.hpp"
#include "CRS.hpp"
#include "mgmres.hpp"
#include "FileFactory.hpp"

class Macro{
protected:
    double * oUptake; // remove
    // double * egfSource; // remove
    Mesh* mesh;
    Frame* frame;
    ConfigHandler* config;
public:
  double DeadConfluence;
	double LiveConfluence;
    Macro(Mesh* mesh, Frame* frame, ConfigHandler* config)
    {
        this->mesh = mesh;
        this->frame = frame;
        this->config = config;
        this->oUptake = new double[this->mesh->matrixSize];
        // this->egfSource = new double[this->mesh->matrixSize];
        this->DeadConfluence = 0.0;
        this->LiveConfluence = 0.0;
    }
    virtual void reaction() = 0;
    virtual void diference() = 0;
    // virtual void diferenceEGF() = 0;
};

class Macro3D : public Macro{
public:
    Macro3D(Mesh* mesh, Frame* frame, ConfigHandler* config) : Macro(mesh, frame, config){}
    void reaction();
    void diference();
    // void diferenceEGF();
};

class Macro2D : public Macro{
public:
    Macro2D(Mesh* mesh, Frame* frame, ConfigHandler* config) : Macro(mesh, frame, config){}
    void reaction();
    void diference();
    // void diferenceEGF();
};

#endif /* end of include guard: MACRO */
