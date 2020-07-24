#ifndef FRAME
#define FRAME

#include <vector>
#include "Cell.hpp"

class Frame
{
public:
    int time, tumorCells, outCells;
    Vector3 domain;
    std::vector<Cell> cells;

    Frame(Vector3 domain = Vector3(), int time = 0, int outCells = 0, int tumorCells = 0, std::vector<Cell> cells = std::vector<Cell>())
    :
        domain(domain),
        time(time),
        outCells(outCells),
        tumorCells(tumorCells),
        cells(cells)
    {}

    std::string to_string()
    {
        std::string out = "============================================\n";
        out += "Time = " + std::to_string(this->time) + " h\n";
        out += "Tumor cells = " + std::to_string(this->tumorCells) + "\n";
        out += "Normal cells = " + std::to_string(this->cells.size()-this->tumorCells) + "\n";
        out += "Total agents = " + std::to_string(this->cells.size()) + "\n";
        out += "Total cells left = " + std::to_string(this->outCells) + "\n";
        out += "============================================\n";
        return out;
    }

};

#endif /* end of include guard: FRAME */
