#ifndef FRAME_FACTORY
#define FRAME_FACTORY

#include <fstream>
#include <cmath>
#include <string>
#include <iostream>
#include "Frame.hpp"
#include "ConfigHandler.hpp"

class FrameFactory{
public:
    static Frame* makeFrame2D(std::string fileName,double oConsumption);

    static Frame* makeFrame3D(std::string fileName,double oConsumption);
};

#endif /* end of include guard: FRAME_FACTORY */
