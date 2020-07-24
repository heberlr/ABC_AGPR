#ifndef CELL_MANIPULATION
#define CELL_MANIPULATION

#include "NR3.hpp"
#include "Ran.hpp"

#include <vector>
#include <string>
#include <cmath>

#include "Cell.hpp"
#include "Vector.hpp"
#include "Frame.hpp"
#include "ConfigHandler.hpp"
#include "Mesh.hpp"

#include <iostream>
#include <algorithm>

class CellManipulation{
public:
    static Cell divide(Cell* cell, double rand1, double rand2);
    static void force(Frame *frame, ConfigHandler *config);
    static double norma(Vector3 pos);
    static Vector3 normal(Vector3 coordinates, double domainRadius);
    static Vector3 func_var_phi(Vector3 normal, double actionRadius, int n);
    static Vector3 func_var_psi(Vector3 normal, double nucleusRadius, double radius, double M, int m);
    static void calculateCellSigmas(Cell *cell, Mesh *mesh);
    static void updateFrame(Frame *frame, ConfigHandler *config, Mesh *mesh, Ran *ran);
};

class CellManipulation2D : public CellManipulation{
public:
    static Vector3 normal(Vector3 coordinates, double domainRadius);
    static void force(Frame *frame, ConfigHandler *config);
    static Cell divide(Cell* cell, double rand1);
    static void calculateCellSigmas(Cell *cell, Mesh *mesh);
    static void updateFrame(Frame *frame, ConfigHandler *config, Mesh *mesh, Ran *ran);

    // static double norma(Vector3 pos);
    // static Vector3 func_var_phi(Vector3 normal, double actionRadius, int n);
    // static Vector3 func_var_psi(Vector3 normal, double nucleusRadius, double radius, double M, int m);

};

#endif /* end of include guard: CELL_MANIPULATION */
