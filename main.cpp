#include "solver.h"

#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
    FluidSim sim;

    writeWallMaskPGM(sim.wallField(), 0);

    std::ofstream diagnosticsCsv("animation/diagnostics.csv");
    writeDiagnosticsHeader(diagnosticsCsv);

    for (int frame = 0; frame <= 500; frame++) {
        sim.step();

        if (frame % 10 == 0) {
            const SimDiagnostics& diagnostics = sim.diagnostics();
            writeDiagnosticsRow(diagnosticsCsv, frame, diagnostics);

            if (frame % 100 == 0) {
                std::cout << "frame " << std::setw(4) << frame
                        << " max_velocity=" << diagnostics.maxVelocity
                        << " max_smoke=" << diagnostics.maxSmoke
                        << " max_divergence_before=" << diagnostics.maxDivergenceBeforeProjection
                        << " max_divergence_after=" << diagnostics.maxDivergenceAfterProjection
                        << '\n';
            }
        }
        writeSmokePGM(sim.smokeField(), frame);
    }

    return 0;
}
