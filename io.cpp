#include "solver.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace {
void ensureAnimationDirectory() {
    std::filesystem::create_directories("animation");
}
}

void writeSmokePGM(const ScalarField& smoke, int frame) {
    ensureAnimationDirectory();

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
    ensureAnimationDirectory();

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

void writeDiagnosticsHeader(std::ostream& out) {
    out << "frame,"
        << "max_velocity,"
        << "max_smoke,"
        << "max_divergence_before_projection,"
        << "max_divergence_after_projection\n";
}

void writeDiagnosticsRow(std::ostream& out, int frame, const SimDiagnostics& diagnostics) {
    out << frame << ','
        << diagnostics.maxVelocity << ','
        << diagnostics.maxSmoke << ','
        << diagnostics.maxDivergenceBeforeProjection << ','
        << diagnostics.maxDivergenceAfterProjection << '\n';
}
