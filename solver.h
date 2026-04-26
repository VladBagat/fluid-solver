#ifndef SOLVER_H
#define SOLVER_H

#include <iosfwd>
#include <vector>

#include "vector.h"

constexpr int NX = 128;
constexpr int NY = 128;
constexpr float DX = 1.0f;

using VectorField = std::vector<Vector2D>;
using ScalarField = std::vector<float>;
using WallField = std::vector<unsigned char>;

struct SimParams {
    float dt = 0.15f;
    float smokeDecay = 0.01f;
    float temperatureDecay = 0.01f;
    float buoyancyStrength = 0.2f;
    float smokeWeight = 0.05f;
    float vorticityStrength = 0.1f;
    int pressureIterations = 35;
};

struct SimDiagnostics {
    float maxVelocity = 0.0f;
    float maxSmoke = 0.0f;
    float maxDivergenceBeforeProjection = 0.0f;
    float maxDivergenceAfterProjection = 0.0f;
};

struct ScalarSample {
    float value;
    float minValue;
    float maxValue;
};

inline int IX(int x, int y) {
    return y * NX + x;
}

class FluidSim {
public:
    explicit FluidSim(const SimParams& params = SimParams{});

    void step();

    const ScalarField& smokeField() const;
    const WallField& wallField() const;
    const SimDiagnostics& diagnostics() const;

private:
    SimParams params_;
    SimDiagnostics diagnostics_;

    VectorField velocity_;
    VectorField velocityPrev_;
    VectorField buoyancy_;
    VectorField vorticityForce_;

    ScalarField smoke_;
    ScalarField smokePrev_;
    ScalarField temperature_;
    ScalarField temperaturePrev_;
    ScalarField pressure_;
    ScalarField divergence_;
    ScalarField curl_;

    ScalarField scalarForward_;
    ScalarField scalarBackward_;
    ScalarField scalarClampMin_;
    ScalarField scalarClampMax_;

    WallField walls_;

    void initializeWalls();
    void clearWallCells(VectorField& field) const;
    void clearWallCells(ScalarField& field) const;

    void addSource();
    void advectVelocity();
    void applyForces();
    void project();
    void advectScalars();

    void addForce(const VectorField& force);
    void advectVector(VectorField& dst, const VectorField& src, const VectorField& velocity) const;
    void advectScalar(ScalarField& dst, const ScalarField& src);
    void advectScalarSL(
        ScalarField& dst,
        const ScalarField& src,
        const VectorField& velocity,
        float dt,
        ScalarField* clampMin = nullptr,
        ScalarField* clampMax = nullptr
    ) const;

    Vector2D sampleVector(const VectorField& field, float x, float y) const;
    ScalarSample sampleScalarWithRange(const ScalarField& field, float x, float y) const;

    void computeDivergence();
    void solvePressure();
    void subtractPressureGradient();
    void applyDomainBoundary();
    void scalarDecay(ScalarField& target, float decayRate) const;
    void calculateBuoyancy();
    void computeCurl();
    void computeVorticityForce();

    float maxVelocity() const;
    float maxScalar(const ScalarField& field) const;
    float maxAbsScalar(const ScalarField& field) const;
};

void writeSmokePGM(const ScalarField& smoke, int frame);
void writeWallMaskPGM(const WallField& walls, int frame);
void writeDiagnosticsHeader(std::ostream& out);
void writeDiagnosticsRow(std::ostream& out, int frame, const SimDiagnostics& diagnostics);

#endif
