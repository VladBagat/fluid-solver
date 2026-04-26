#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include "solver.h"

int main() {
    VectorField velocity(NX * NY);
    VectorField velocityPrev(NX * NY);

    VectorField buoyancy(NX * NY);

    ScalarField smoke(NX * NY);
    ScalarField smokePrev(NX * NY);

    ScalarField temperature(NX * NY);
    ScalarField temperaturePrev(NX * NY);

    ScalarField pressure(NX * NY);
    ScalarField divergence(NX * NY);

    ScalarField curl(NX * NY);
    VectorField vorticityForce(NX * NY);

    WallField walls(NX * NY);
    initializeWalls(walls);
    writeWallMaskPGM(walls, 0);

    const float dt = 0.15f;

    const float smokeDecay = 0.01f;
    const float temperatureDecay = 0.01f;

    for (int frame = 0; frame < 500; frame++) {

        addSource(smoke, temperature, walls);
        
        velocityPrev = velocity;
        advectVector(velocity, velocityPrev, velocityPrev, dt, walls);
        applyDomainBoundary(velocity, walls);
        
        calculateBuoyancy(buoyancy, temperature, smoke, walls);
        addForce(velocity, buoyancy, dt, walls);
        applyDomainBoundary(velocity, walls);
        
        computeCurl(curl, velocity, walls);
        computeVorticityForce(vorticityForce, curl, walls);
        addForce(velocity, vorticityForce, dt, walls);
        applyDomainBoundary(velocity, walls);

        project(velocity, pressure, divergence, walls);
        applyDomainBoundary(velocity, walls);

        smokePrev = smoke;
        advectScalar(smoke, smokePrev, velocity, dt, walls);
        clearWallCells(smoke, walls);

        temperaturePrev = temperature;
        advectScalar(temperature, temperaturePrev, velocity, dt, walls);
        clearWallCells(temperature, walls);

        scalarDecay(smoke, smokeDecay, dt);
        scalarDecay(temperature, temperatureDecay, dt);

        clearWallCells(smoke, walls);
        clearWallCells(temperature, walls);

        writeSmokePGM(smoke, frame);
    }

    return 0;
}

void computeCurl(ScalarField& curl, const VectorField& velocity, const WallField& walls) {
    std::fill(curl.begin(), curl.end(), 0.0f);

    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX - 1; x++) {
            if (walls[IX(x, y)]) {
                continue;
            }

            float dvdx = (velocity[IX(x + 1, y)].y - velocity[IX(x - 1, y)].y) * 0.5f / DX;
            float dudy = (velocity[IX(x, y + 1)].x - velocity[IX(x, y - 1)].x) * 0.5f / DX;

            curl[IX(x, y)] = dvdx - dudy;
        }
    }
}

void computeVorticityForce(VectorField& vorticityForce, const ScalarField& curl, const WallField& walls) {
    std::fill(vorticityForce.begin(), vorticityForce.end(), Vector2D(0.0f, 0.0f));

    const float epsilon = 0.1f;

    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX - 1; x++) {
            int id = IX(x,y);

            if (walls[id]) {
                continue;
            }

            float dwdx = (std::abs(curl[IX(x+1,y)]) - std::abs(curl[IX(x-1,y)])) * 0.5f / DX;
            float dwdy = (std::abs(curl[IX(x,y+1)]) - std::abs(curl[IX(x,y-1)])) * 0.5f / DX;

            float length = sqrt(dwdx*dwdx + dwdy*dwdy) + 1e-5f;

            float Nx = dwdx / length;
            float Ny = dwdy / length;

            auto Fx =  Ny * curl[id];
            auto Fy = -Nx * curl[id];

            vorticityForce[id].x = epsilon * Fx;
            vorticityForce[id].y = epsilon * Fy;
        }
    }
}

void initializeWalls(WallField& walls) {
    std::fill(walls.begin(), walls.end(), 0);

    for (int y = 0; y < NY; y++) {
        walls[IX(0, y)] = 1;
        walls[IX(NX - 1, y)] = 1;
    }

    for (int x = 0; x < NX; x++) {
        walls[IX(x, 0)] = 1;
        walls[IX(x, NY - 1)] = 1;
    }

    const int shelfX0 = NX / 4;
    const int shelfX1 = (3 * NX) / 4;
    const int shelfCenterY = (2 * NY) / 3;
    const int shelfY0 = shelfCenterY - 1;
    const int shelfY1 = shelfCenterY + 1;

    for (int y = shelfY0; y <= shelfY1; y++) {
        for (int x = shelfX0; x <= shelfX1; x++) {
            walls[IX(x, y)] = 1;
        }
    }
}

void clearWallCells(VectorField& field, const WallField& walls) {
    for (int i = 0; i < NX * NY; i++) {
        if (walls[i]) {
            field[i] = Vector2D(0.0f, 0.0f);
        }
    }
}

void clearWallCells(ScalarField& field, const WallField& walls) {
    for (int i = 0; i < NX * NY; i++) {
        if (walls[i]) {
            field[i] = 0.0f;
        }
    }
}

void scalarDecay(ScalarField& target, const float decayRate, const float dt) {
    float factor = exp(-decayRate * dt);
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            target[IX(x,y)] *= factor;
        }
    }    
}

void calculateBuoyancy(VectorField &buoyancy, const ScalarField &temperature, const ScalarField &smoke, const WallField& walls) {
    const float buoyancyStrength = 0.2f;
    const float ambientTemperature = 0.0f;
    const float smokeWeight = 0.05f;

    std::fill(buoyancy.begin(), buoyancy.end(), Vector2D(0.0f, 0.0f));

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            if (!walls[id]) {
                buoyancy[id].y = buoyancyStrength * (temperature[id] - ambientTemperature) - smokeWeight * smoke[id];
            }
        }
    }   
}

void addSource(ScalarField& smoke, ScalarField& temperature, const WallField& walls) {
    int cx = NX / 2;
    int cy = NY / 4;

    for (int y = cy - 3; y <= cy + 3; y++) {
        for (int x = cx - 3; x <= cx + 3; x++) {
            int id = IX(x, y);
            if (walls[id]) {
                continue;
            }

            smoke[id] = 1.0f;
            temperature[id] = 5.0f;
        }
    }
}

void addForce(VectorField& velocity, const VectorField& force, float dt, const WallField& walls) {
    for (int i = 0; i < NX * NY; i++) {
        if (walls[i]) {
            velocity[i] = Vector2D(0.0f, 0.0f);
        } else {
            velocity[i] += (force[i] * dt);
        }
    }
}

void advectVector(VectorField& dst, const VectorField& src, const VectorField& velocity, float dt, const WallField& walls) {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            if (walls[id]) {
                dst[id] = Vector2D(0.0f, 0.0f);
                continue;
            }

            Vector2D v = velocity[id];

            float prev_x = x - dt * v.x;
            float prev_y = y - dt * v.y;

            dst[id] = sampleVector(src, prev_x, prev_y, walls);
        }
    }
}

void advectScalar(ScalarField& dst, const ScalarField& src, const VectorField& velocity, float dt, const WallField& walls) {
    ScalarField forward(NX * NY);
    ScalarField backward(NX * NY);
    ScalarField clampMin(NX * NY);
    ScalarField clampMax(NX * NY);

    advectScalarSL(forward, src, velocity, dt, walls, &clampMin, &clampMax);
    advectScalarSL(backward, forward, velocity, -dt, walls);

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            auto id = IX(x,y);

            if (walls[id]) {
                dst[id] = 0.0f;
                continue;
            }
            
            float corrected = forward[id] + 0.5f * (src[id] - backward[id]);
            dst[id] = std::clamp(corrected, clampMin[id], clampMax[id]);
        }
    }
}

void advectScalarSL(ScalarField& dst, const ScalarField& src, const VectorField& velocity, float dt, const WallField& walls, ScalarField* clampMin, ScalarField* clampMax) {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);
            if (walls[id]) {
                dst[id] = 0.0f;
                if (clampMin) {
                    (*clampMin)[id] = 0.0f;
                }
                if (clampMax) {
                    (*clampMax)[id] = 0.0f;
                }
                continue;
            }

            Vector2D v = velocity[id];

            float prev_x = x - dt * v.x;
            float prev_y = y - dt * v.y;

            ScalarSample sample = sampleScalarWithRange(src, prev_x, prev_y, walls);
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

Vector2D sampleVector(const VectorField& field, float x, float y, const WallField& walls) {
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
        if (weight > 0.0f && !walls[id]) {
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

ScalarSample sampleScalarWithRange(const ScalarField& field, float x, float y, const WallField& walls) {
    x = std::clamp(x, 0.0f, float(NX - 1));
    y = std::clamp(y, 0.0f, float(NY - 1));

    int x0 = int(std::floor(x));
    int y0 = int(std::floor(y));

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
        if (weight > 0.0f && !walls[id]) {
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

float sampleScalar(const ScalarField& field, float x, float y, const WallField& walls) {
    return sampleScalarWithRange(field, x, y, walls).value;
}

void computeDivergence(const VectorField &velocity, ScalarField &divergence, const WallField& walls)
{ 
    std::fill(divergence.begin(), divergence.end(), 0.0f);
    for (int y = 1; y < NY - 1; y++){
        for (int x = 1; x < NX - 1; x++) {
            int id = IX(x, y);
            if (walls[id]) {
                continue;
            }

            float rightX = walls[IX(x + 1, y)] ? 0.0f : velocity[IX(x + 1, y)].x;
            float leftX = walls[IX(x - 1, y)] ? 0.0f : velocity[IX(x - 1, y)].x;
            float topY = walls[IX(x, y + 1)] ? 0.0f : velocity[IX(x, y + 1)].y;
            float bottomY = walls[IX(x, y - 1)] ? 0.0f : velocity[IX(x, y - 1)].y;

            float div = 0.5f * (
                rightX - leftX +
                topY - bottomY
            ) / DX;
            divergence[id] = div;
        }
    }
    
}

void solvePressure(ScalarField& pressure, const ScalarField& divergence, const WallField& walls) {
    std::fill(pressure.begin(), pressure.end(), 0.0f);

    for (int iter = 0; iter < 35; iter++) {
        for (int y = 1; y < NY - 1; y++) {
            for (int x = 1; x < NX - 1; x++) {
                int id = IX(x, y);
                if (walls[id]) {
                    pressure[id] = 0.0f;
                    continue;
                }

                float center = pressure[id];
                float right = walls[IX(x + 1, y)] ? center : pressure[IX(x + 1, y)];
                float left = walls[IX(x - 1, y)] ? center : pressure[IX(x - 1, y)];
                float top = walls[IX(x, y + 1)] ? center : pressure[IX(x, y + 1)];
                float bottom = walls[IX(x, y - 1)] ? center : pressure[IX(x, y - 1)];

                pressure[id] =
                    (right +
                     left +
                     top +
                     bottom -
                     divergence[id]) / 4.0f;
            }
        }
    }

    clearWallCells(pressure, walls);
}

void subtractPressureGradient(VectorField& velocity, const ScalarField& pressure, const WallField& walls) {
    for (int y = 1; y < NY - 1; y++) {
        for (int x = 1; x < NX - 1; x++) {
            int id = IX(x, y);
            if (walls[id]) {
                velocity[id] = Vector2D(0.0f, 0.0f);
                continue;
            }

            float center = pressure[id];
            float right = walls[IX(x + 1, y)] ? center : pressure[IX(x + 1, y)];
            float left = walls[IX(x - 1, y)] ? center : pressure[IX(x - 1, y)];
            float top = walls[IX(x, y + 1)] ? center : pressure[IX(x, y + 1)];
            float bottom = walls[IX(x, y - 1)] ? center : pressure[IX(x, y - 1)];

            float gradX =
                0.5f * (right - left) / DX;

            float gradY =
                0.5f * (top - bottom) / DX;

            velocity[id].x -= gradX;
            velocity[id].y -= gradY;

            if (walls[IX(x - 1, y)] && velocity[id].x < 0.0f) {
                velocity[id].x = 0.0f;
            }
            if (walls[IX(x + 1, y)] && velocity[id].x > 0.0f) {
                velocity[id].x = 0.0f;
            }
            if (walls[IX(x, y - 1)] && velocity[id].y < 0.0f) {
                velocity[id].y = 0.0f;
            }
            if (walls[IX(x, y + 1)] && velocity[id].y > 0.0f) {
                velocity[id].y = 0.0f;
            }
        }
    }

    clearWallCells(velocity, walls);
}

void project(VectorField& velocity, ScalarField& pressure, ScalarField& divergence, const WallField& walls) {
    computeDivergence(velocity, divergence, walls);
    solvePressure(pressure, divergence, walls);
    subtractPressureGradient(velocity, pressure, walls);
}

void applyDomainBoundary(VectorField& velocity, const WallField& walls) {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int id = IX(x, y);

            if (walls[id]) {
                velocity[id] = Vector2D(0.0f, 0.0f);
                continue;
            }

            if (x > 0 && walls[IX(x - 1, y)] && velocity[id].x < 0.0f) {
                velocity[id].x = 0.0f;
            }

            if (x + 1 < NX && walls[IX(x + 1, y)] && velocity[id].x > 0.0f) {
                velocity[id].x = 0.0f;
            }

            if (y > 0 && walls[IX(x, y - 1)] && velocity[id].y < 0.0f) {
                velocity[id].y = 0.0f;
            }

            if (y + 1 < NY && walls[IX(x, y + 1)] && velocity[id].y > 0.0f) {
                velocity[id].y = 0.0f;
            }
        }
    }
}

void writeSmokePGM(const ScalarField& smoke, int frame) {
    std::filesystem::create_directories("animation");

    std::ostringstream name;
    name << "animation/frame_"
         << std::setw(5) << std::setfill('0') << frame
         << ".pgm";

    std::ofstream out(name.str(), std::ios::binary);

    out << "P5\n" << NX << " " << NY << "\n255\n";

    for (int y = NY - 1; y >= 0; y--) {
        for (int x = 0; x < NX; x++) {
            float d = std::clamp(smoke[IX(x, y)], 0.0f, 1.0f);
            unsigned char pixel = static_cast<unsigned char>(d * 255.0f);
            out.write(reinterpret_cast<char*>(&pixel), 1);
        }
    }
}

void writeWallMaskPGM(const WallField& walls, int frame) {
    std::filesystem::create_directories("animation");

    std::ostringstream name;
    name << "animation/wall_mask_"
         << std::setw(5) << std::setfill('0') << frame
         << ".pgm";

    std::ofstream out(name.str(), std::ios::binary);

    out << "P5\n" << NX << " " << NY << "\n255\n";

    for (int y = NY - 1; y >= 0; y--) {
        for (int x = 0; x < NX; x++) {
            unsigned char pixel = walls[IX(x, y)] ? 255 : 0;
            out.write(reinterpret_cast<char*>(&pixel), 1);
        }
    }
}
