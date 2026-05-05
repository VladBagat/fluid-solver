#include "solver.h"

#include <algorithm>
#include <cmath>

FluidSim::FluidSim(const SimParams& params)
    : params_(params),
      buoyancy_(NX * NY),
      vorticityForce_(NX * NY),
      u((NX + 1) * NY),
      uPrev_((NX + 1) * NY),
      v(NX * (NY + 1)),
      vPrev_(NX * (NY + 1)),
      smoke_(NX * NY),
      smokePrev_(NX * NY),
      temperature_(NX * NY),
      temperaturePrev_(NX * NY),
      pressure_(NX * NY),
      divergence_(NX * NY),
      curl_(NX * NY),
      scalarForward_(NX * NY),
      scalarBackward_(NX * NY),
      scalarClampMin_(NX * NY),
      scalarClampMax_(NX * NY) {}

void FluidSim::step() {
    addSource();
    advectVelocity();
    applyForces();
    project();
    advectScalars();

    diagnostics_.maxVelocity = maxVelocity();
    diagnostics_.maxSmoke = maxScalar(smoke_);
}

const ScalarField& FluidSim::smokeField() const {
    return smoke_;
}

const SimDiagnostics& FluidSim::diagnostics() const {
    return diagnostics_;
}

void FluidSim::addSource() {
    const int cx = NX / 2;
    const int cy = NY / 4;

    for (int y = cy - 3; y <= cy + 3; y++) {
        for (int x = cx - 3; x <= cx + 3; x++) {
            int id = IX(x, y);
            smoke_[id] = 1.0f;
            temperature_[id] = 5.0f;
        }
    }
}

void FluidSim::advectVelocity() {
    uPrev_ = u;
    vPrev_ = v;

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x <= NX; x++) {
            float faceX = float(x) - 0.5f;
            float faceY = float(y);
            Vector2D velocity = sampleMacVelocity(uPrev_, vPrev_, faceX, faceY);
            float prevX = faceX - params_.dt * velocity.x;
            float prevY = faceY - params_.dt * velocity.y;

            u[U_INDEX(x, y)] = sampleU(uPrev_, prevX, prevY);
        }
    }

    for (int y = 0; y <= NY; y++) {
        for (int x = 0; x < NX; x++) {
            float faceX = float(x);
            float faceY = float(y) - 0.5f;
            Vector2D velocity = sampleMacVelocity(uPrev_, vPrev_, faceX, faceY);
            float prevX = faceX - params_.dt * velocity.x;
            float prevY = faceY - params_.dt * velocity.y;

            v[V_INDEX(x, y)] = sampleV(vPrev_, prevX, prevY);
        }
    }

    applyDomainBoundary();
}

void FluidSim::applyForces() {
    calculateBuoyancy();
    addForce(buoyancy_);
    applyDomainBoundary();

    computeCurl();
    computeVorticityForce();
    addForce(vorticityForce_);
    applyDomainBoundary();
}

void FluidSim::project() {
    computeDivergence();
    diagnostics_.maxDivergenceBeforeProjection = maxAbsScalar(divergence_);

    solvePressure();
    subtractPressureGradient();
    applyDomainBoundary();

    computeDivergence();
    diagnostics_.maxDivergenceAfterProjection = maxAbsScalar(divergence_);
}

void FluidSim::advectScalars() {
    smokePrev_ = smoke_;
    advectScalar(smoke_, smokePrev_);

    temperaturePrev_ = temperature_;
    advectScalar(temperature_, temperaturePrev_);

    scalarDecay(smoke_, params_.smokeDecay);
    scalarDecay(temperature_, params_.temperatureDecay);
}

void FluidSim::addForce(const VectorField& force) {
    for (int y = 0; y < NY; y++) {
        for (int x = 1; x < NX; x++) {
            float forceX = 0.5f * (force[IX(x - 1, y)].x + force[IX(x, y)].x);
            u[U_INDEX(x, y)] += forceX * params_.dt;
        }
    }

    for (int y = 1; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            float forceY = 0.5f * (force[IX(x, y - 1)].y + force[IX(x, y)].y);
            v[V_INDEX(x, y)] += forceY * params_.dt;
        }
    }
}

void FluidSim::advectScalar(ScalarField& dst, const ScalarField& src) {
    advectScalarSL(scalarForward_, src, params_.dt, &scalarClampMin_, &scalarClampMax_);
    advectScalarSL(scalarBackward_, scalarForward_, -params_.dt);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);

            float corrected = scalarForward_[id] + 0.5f * (src[id] - scalarBackward_[id]);
            dst[id] = std::clamp(corrected, scalarClampMin_[id], scalarClampMax_[id]);
        }
    }
}

void FluidSim::advectScalarSL(
    ScalarField& dst,
    const ScalarField& src,
    float dt,
    ScalarField* clampMin,
    ScalarField* clampMax
) const {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            Vector2D velocity = sampleMacVelocity(float(x), float(y));
            float prevX = float(x) - dt * velocity.x;
            float prevY = float(y) - dt * velocity.y;

            ScalarSample sample = sampleScalarWithRange(src, prevX, prevY);
            dst[id] = sample.value;
            if (clampMin) {
                (*clampMin)[id] = sample.minValue;
            }
            if (clampMax) {
                (*clampMax)[id] = sample.maxValue;
            }
        }
    }
}

float FluidSim::sampleField(const ScalarField& field, int width, int height, float x, float y) const {
    x = std::clamp(x, 0.0f, float(width - 1));
    y = std::clamp(y, 0.0f, float(height - 1));

    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));

    int x1 = std::min(x0 + 1, width - 1);
    int y1 = std::min(y0 + 1, height - 1);

    float tx = x - float(x0);
    float ty = y - float(y0);

    float bottomLeft = field[y0 * width + x0];
    float bottomRight = field[y0 * width + x1];
    float topLeft = field[y1 * width + x0];
    float topRight = field[y1 * width + x1];

    float bottom = bottomLeft * (1.0f - tx) + bottomRight * tx;
    float top = topLeft * (1.0f - tx) + topRight * tx;

    return bottom * (1.0f - ty) + top * ty;
}

float FluidSim::sampleU(const ScalarField& field, float x, float y) const {
    return sampleField(field, NX + 1, NY, x + 0.5f, y);
}

float FluidSim::sampleV(const ScalarField& field, float x, float y) const {
    return sampleField(field, NX, NY + 1, x, y + 0.5f);
}

Vector2D FluidSim::sampleMacVelocity(float x, float y) const {
    return sampleMacVelocity(u, v, x, y);
}

Vector2D FluidSim::sampleMacVelocity(const ScalarField& uField, const ScalarField& vField, float x, float y) const {
    return Vector2D(sampleU(uField, x, y), sampleV(vField, x, y));
}

ScalarSample FluidSim::sampleScalarWithRange(const ScalarField& field, float x, float y) const {
    x = std::clamp(x, 0.0f, float(NX - 1));
    y = std::clamp(y, 0.0f, float(NY - 1));

    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));

    int x1 = std::min(x0 + 1, NX - 1);
    int y1 = std::min(y0 + 1, NY - 1);

    float tx = x - float(x0);
    float ty = y - float(y0);

    float result = 0.0f;
    float weightSum = 0.0f;
    float minValue = 0.0f;
    float maxValue = 0.0f;
    bool foundSample = false;

    auto addSample = [&](int sx, int sy, float weight) {
        int id = IX(sx, sy);
        if (weight > 0.0f) {
            float value = field[id];
            result += value * weight;
            weightSum += weight;
            if (!foundSample) {
                minValue = value;
                maxValue = value;
                foundSample = true;
            } else {
                minValue = std::min(minValue, value);
                maxValue = std::max(maxValue, value);
            }
        }
    };

    addSample(x0, y0, (1.0f - tx) * (1.0f - ty));
    addSample(x1, y0, tx * (1.0f - ty));
    addSample(x0, y1, (1.0f - tx) * ty);
    addSample(x1, y1, tx * ty);

    if (weightSum == 0.0f) {
        return {0.0f, 0.0f, 0.0f};
    }

    return {result / weightSum, minValue, maxValue};
}

void FluidSim::computeDivergence() {
    std::fill(divergence_.begin(), divergence_.end(), 0.0f);

    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX - 1; x++) {
            int id = IX(x, y);

            float rightX = u[U_INDEX(x + 1, y)];
            float leftX = u[U_INDEX(x, y)];
            float topY = v[V_INDEX(x, y + 1)];
            float bottomY = v[V_INDEX(x, y)];

            divergence_[id] = (rightX - leftX + topY - bottomY) / DX;
        }
    }
}

void FluidSim::solvePressure() {
    std::fill(pressure_.begin(), pressure_.end(), 0.0f);

    for (int iter = 0; iter < params_.pressureIterations; iter++) {
        for (int y = 1; y < NY - 1; y++) {
            for (int x = 1; x < NX - 1; x++) {
                int id = IX(x, y);

                float right = pressure_[IX(x + 1, y)];
                float left = pressure_[IX(x - 1, y)];
                float top = pressure_[IX(x, y + 1)];
                float bottom = pressure_[IX(x, y - 1)];

                pressure_[id] = (right + left + top + bottom - divergence_[id]) / 4.0f;
            }
        }
    }
}

void FluidSim::subtractPressureGradient() {
    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX; x++) {
            float gradient = (pressure_[IX(x, y)] - pressure_[IX(x - 1, y)]) / DX;
            u[U_INDEX(x, y)] -= gradient;
        }
    }

    for (int y = 1; y < NY; y++) {
        for (int x = 1; x < NX - 1; x++) {
            float gradient = (pressure_[IX(x, y)] - pressure_[IX(x, y - 1)]) / DX;
            v[V_INDEX(x, y)] -= gradient;
        }
    }
}

void FluidSim::applyDomainBoundary() {
    for (int y = 0; y < NY; y++) {
        u[U_INDEX(0, y)] = 0.0f;
        u[U_INDEX(NX, y)] = 0.0f;
    }

    for (int x = 0; x < NX; x++) {
        v[V_INDEX(x, 0)] = 0.0f;
        v[V_INDEX(x, NY)] = 0.0f;
    }
}

void FluidSim::scalarDecay(ScalarField& target, float decayRate) const {
    float factor = std::exp(-decayRate * params_.dt);
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            target[IX(x, y)] *= factor;
        }
    }
}

void FluidSim::calculateBuoyancy() {
    constexpr float ambientTemperature = 0.0f;

    std::fill(buoyancy_.begin(), buoyancy_.end(), Vector2D(0.0f, 0.0f));

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            buoyancy_[id].y =
                params_.buoyancyStrength * (temperature_[id] - ambientTemperature) -
                params_.smokeWeight * smoke_[id];
        }
    }
}

void FluidSim::computeCurl() {
    std::fill(curl_.begin(), curl_.end(), 0.0f);

    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX - 1; x++) {
            float dvdx = (sampleMacVelocity(float(x + 1), float(y)).y -
                          sampleMacVelocity(float(x - 1), float(y)).y) * 0.5f / DX;
            float dudy = (sampleMacVelocity(float(x), float(y + 1)).x -
                          sampleMacVelocity(float(x), float(y - 1)).x) * 0.5f / DX;

            curl_[IX(x, y)] = dvdx - dudy;
        }
    }
}

void FluidSim::computeVorticityForce() {
    std::fill(vorticityForce_.begin(), vorticityForce_.end(), Vector2D(0.0f, 0.0f));

    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX - 1; x++) {
            int id = IX(x, y);

            float dwdx = (std::abs(curl_[IX(x + 1, y)]) - std::abs(curl_[IX(x - 1, y)])) * 0.5f / DX;
            float dwdy = (std::abs(curl_[IX(x, y + 1)]) - std::abs(curl_[IX(x, y - 1)])) * 0.5f / DX;

            float length = std::sqrt(dwdx * dwdx + dwdy * dwdy) + 1e-5f;

            float nx = dwdx / length;
            float ny = dwdy / length;

            float fx = ny * curl_[id];
            float fy = -nx * curl_[id];

            vorticityForce_[id].x = params_.vorticityStrength * fx;
            vorticityForce_[id].y = params_.vorticityStrength * fy;
        }
    }
}

float FluidSim::maxVelocity() const {
    float result = 0.0f;
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            Vector2D velocity = sampleMacVelocity(float(x), float(y));
            float magnitude = std::sqrt(velocity.x * velocity.x + velocity.y * velocity.y);
            result = std::max(result, magnitude);
        }
    }
    return result;
}

float FluidSim::maxScalar(const ScalarField& field) const {
    float result = 0.0f;
    for (float value : field) {
        result = std::max(result, value);
    }
    return result;
}

float FluidSim::maxAbsScalar(const ScalarField& field) const {
    float result = 0.0f;
    for (float value : field) {
        result = std::max(result, std::abs(value));
    }
    return result;
}
