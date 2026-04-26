#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include "vector.h"

constexpr int NX = 128;
constexpr int NY = 128;
constexpr float DX = 1.0f;

using VectorField = std::vector<Vector2D>;
using ScalarField = std::vector<float>;
using WallField = std::vector<unsigned char>;

struct ScalarSample {
    float value;
    float minValue;
    float maxValue;
};

inline int IX(int x, int y) {
    return y * NX + x;
}

void initializeWalls(WallField& walls);
void clearWallCells(VectorField& field, const WallField& walls);
void clearWallCells(ScalarField& field, const WallField& walls);
void addForce(VectorField& velocity, const VectorField& force, float dt, const WallField& walls);
void advectVector(VectorField& dst, const VectorField& src, const VectorField& velocity, float dt, const WallField& walls);
void advectScalar(ScalarField& dst, const ScalarField& src, const VectorField& velocity, float dt, const WallField& walls);
void advectScalarSL(ScalarField& dst, const ScalarField& src, const VectorField& velocity, float dt, const WallField& walls, ScalarField* clampMin = nullptr, ScalarField* clampMax = nullptr);
Vector2D sampleVector(const VectorField& field, float x, float y, const WallField& walls);
ScalarSample sampleScalarWithRange(const ScalarField& field, float x, float y, const WallField& walls);
float sampleScalar(const ScalarField& field, float x, float y, const WallField& walls);
void computeDivergence(const VectorField& velocity, ScalarField& divergence, const WallField& walls);
void solvePressure(ScalarField& pressure, const ScalarField& divergence, const WallField& walls);
void subtractPressureGradient(VectorField& velocity, const ScalarField& pressure, const WallField& walls);
void addSource(ScalarField& smoke, ScalarField& temperature, const WallField& walls);
void project(VectorField& velocity, ScalarField& pressure, ScalarField& divergence, const WallField& walls);
void writeSmokePGM(const ScalarField& smoke, int frame);
void writeWallMaskPGM(const WallField& walls, int frame);
void calculateBuoyancy(VectorField &buoyancy, const ScalarField &temperature, const ScalarField &smoke, const WallField& walls);
void scalarDecay(ScalarField& target, const float decayRate, const float dt);
void applyDomainBoundary(VectorField& velocity, const WallField& walls);
void computeCurl(ScalarField& curl, const VectorField& velocity, const WallField& walls);
void computeVorticityForce(VectorField& vorticityForce, const ScalarField& curl, const WallField& walls);
#endif
