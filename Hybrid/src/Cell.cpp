#include "Cell.hpp"

std::string Cell::to_string() {
    std::string out = "Pos = " + this->coordinates.to_string();
    out += " Radius | Nuclear = " + std::to_string(this->nucleusRadius) + ", Cell = " + std::to_string(this->radius) + ", Max = " + std::to_string(this->actionRadius);
    out += " Time = " + std::to_string(this->lifetime);
    out += " V = " + this->speed.to_string();

    switch (this->state) {
        case 1:
            out += " State = Quiescent";
            break;
        case 2:
            out += " State = Proliferative";
            break;
        case 4:
            out += " State = Apoptotic";
            break;
        case 5:
            out += " State = G1";
            break;
        default:
            out += " State = " + this->state;
    }
    return out;
}
