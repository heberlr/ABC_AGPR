#ifndef CELL
#define CELL

#include "Vector.hpp"

class Cell
{
public:
	int state;

	Vector3 coordinates,
			speed,
			force;

	double	nucleusRadius,
			radius,
			actionRadius,
			oConsumption,
			radiusRate;
			// egfSource,
			// calcification;

	int previousState,
		lifetime;

	double 	sigmaO;
			// sigmaEGF;

	Cell(int state = -1, Vector3 coordinates = Vector3(), double nucleusRadius = 0.0,
		double radius = 0.0, double actionRadius = 0.0, int lifetime = 0,
		int previousState = 0, double oConsumption = 0.0, double egfSource = 0.0,
		double calcification = 0.0, double sigmaO = 0.0, double sigmaEGF = 0.0, Vector3 speed = Vector3(), Vector3 force = Vector3())
		:
			state(state),
			coordinates(coordinates),
			nucleusRadius(nucleusRadius),
			radius(radius),
			actionRadius(actionRadius),
			lifetime(lifetime),
			previousState(previousState),
			oConsumption(oConsumption),
			// egfSource(egfSource),
			// calcification(calcification),
			radiusRate(0.0),
			speed(speed),
			sigmaO(sigmaO),
			// sigmaEGF(sigmaEGF),
			force(force)
	{}

	std::string to_string();

};

#endif /* end of include guard: CELL */
