#include "solver.h"

#include <algorithm>
#include <cmath>

FluidSim::FluidSim(const SimParams& params)
    : params_(params),
      velocity_(NX * NY),
      velocityPrev_(NX * NY),
      buoyancy_(NX * NY),
      vorticityForce_(NX * NY),
      u((NX + 1) * NY),
      v(NX * (NY + 1)),
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
    velocityPrev_ = velocity_;
    advectVector(velocity_, velocityPrev_, velocityPrev_);
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
    for (int i = 0; i < NX * NY; i++) {
        velocity_[i] += force[i] * params_.dt;
    }
}

void FluidSim::advectVector(VectorField& dst, const VectorField& src, const VectorField& velocity) const {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            Vector2D v = velocity[id];
            float prevX = x - params_.dt * v.x;
            float prevY = y - params_.dt * v.y;

            dst[id] = sampleVector(src, prevX, prevY);
        }
    }
}

void FluidSim::advectScalar(ScalarField& dst, const ScalarField& src) {
    advectScalarSL(scalarForward_, src, velocity_, params_.dt, &scalarClampMin_, &scalarClampMax_);
    advectScalarSL(scalarBackward_, scalarForward_, velocity_, -params_.dt);

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
    const VectorField& velocity,
    float dt,
    ScalarField* clampMin,
    ScalarField* clampMax
) const {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            Vector2D v = velocity[id];
            float prevX = x - dt * v.x;
            float prevY = y - dt * v.y;

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

Vector2D FluidSim::sampleVector(const VectorField& field, float x, float y) const {
    x = std::clamp(x, 0.0f, float(NX - 1));
    y = std::clamp(y, 0.0f, float(NY - 1));

    int x0 = static_cast<int>(std::floor(x));
    int y0 = static_cast<int>(std::floor(y));

    int x1 = std::min(x0 + 1, NX - 1);
    int y1 = std::min(y0 + 1, NY - 1);

    float tx = x - float(x0);
    float ty = y - float(y0);

    Vector2D result(0.0f, 0.0f);
    float weightSum = 0.0f;

    auto addSample = [&](int sx, int sy, float weight) {
        int id = IX(sx, sy);
        if (weight > 0.0f) {
            result += field[id] * weight;
            weightSum += weight;
        }
    };

    addSample(x0, y0, (1.0f - tx) * (1.0f - ty));
    addSample(x1, y0, tx * (1.0f - ty));
    addSample(x0, y1, (1.0f - tx) * ty);
    addSample(x1, y1, tx * ty);

    if (weightSum == 0.0f) {
        return Vector2D(0.0f, 0.0f);
    }

    return result / weightSum;
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
        for (int x = 1; x < NX - 1; x++) {
            int id = IX(x, y);

            float right = pressure_[IX(x + 1, y)];
            float left = pressure_[IX(x - 1, y)];
            float top = pressure_[IX(x, y + 1)];
            float bottom = pressure_[IX(x, y - 1)];

            float gradX = 0.5f * (right - left) / DX;
            float gradY = 0.5f * (top - bottom) / DX;

            velocity_[id].x -= gradX;
            velocity_[id].y -= gradY;
        }
    }
}

void FluidSim::applyDomainBoundary() {
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
            float dvdx = (velocity_[IX(x + 1, y)].y - velocity_[IX(x - 1, y)].y) * 0.5f / DX;
            float dudy = (velocity_[IX(x, y + 1)].x - velocity_[IX(x, y - 1)].x) * 0.5f / DX;

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
    for (const Vector2D& velocity : velocity_) {
        float magnitude = std::sqrt(velocity.x * velocity.x + velocity.y * velocity.y);
        result = std::max(result, magnitude);
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
