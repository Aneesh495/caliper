/**
 * Caliper Metrology Engine: Unit Tests for C Bridge & Pipeline Harness
 * File: test_c_bridge.c
 *
 * Verifies end-to-end native execution of Signed Distance Field inspection,
 * deviation metrics, and GD&T tolerance checks.
 */

#include "caliper_bridge.h"
#include "caliper_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

static bool test_bridge_inspection(void) {
    CaliperEngine *engine = caliper_engine_create(4);
    assert(engine != NULL);

    /* Construct reference CAD triangle: (0,0,0) - (10,0,0) - (0,10,0) in XY plane */
    CaliperMesh *cad = caliper_mesh_create(1);
    cad->triangles[0].normal = (CaliperVec3){0.0, 0.0, 1.0};
    cad->triangles[0].v0 = (CaliperVec3){0.0, 0.0, 0.0};
    cad->triangles[0].v1 = (CaliperVec3){10.0, 0.0, 0.0};
    cad->triangles[0].v2 = (CaliperVec3){0.0, 10.0, 0.0};

    /* Construct 4 scan points above triangle at Z = +0.020 mm */
    size_t n = 4;
    CaliperPointCloud *scan = caliper_pointcloud_create(n, true, true, true);
    for (size_t i = 0; i < n; ++i) {
        scan->x[i] = 1.0 + (double)i * 1.5;
        scan->y[i] = 1.0 + (double)i * 1.0;
        scan->z[i] = 0.020; /* 20 micron positive deviation */
    }

    CaliperInspectionSummary summary;
    bool ok = caliper_run_inspection(engine, scan, cad, 0.050, &summary);
    assert(ok);
    assert(summary.evaluated_points == n);
    assert(fabs(summary.mean_deviation_mm - 0.020) < 1e-4);
    assert(summary.in_tolerance_pct == 100.0);
    assert(summary.inspection_passed);

    caliper_engine_free(engine);
    caliper_mesh_free(cad);
    caliper_pointcloud_free(scan);
    return true;
}

static bool test_bridge_gdt(void) {
    CaliperEngine *engine = caliper_engine_create(4);
    assert(engine != NULL);

    size_t n = 12;
    CaliperPointCloud *scan = caliper_pointcloud_create(n, true, true, false);
    for (size_t i = 0; i < n; ++i) {
        double theta = (double)i * (2.0 * 3.141592653589793 / (double)n);
        scan->x[i] = 10.0 * cos(theta);
        scan->y[i] = 10.0 * sin(theta);
        scan->z[i] = 0.010; /* Flatness = 0.0, Cylindricity = 0.0 */
    }

    CaliperGdtSummary gdt;
    bool ok = caliper_evaluate_gdt_sample(engine, scan, 0.025, 0.030, &gdt);
    assert(ok);
    assert(gdt.flatness_passed);
    assert(gdt.overall_gdt_passed);

    caliper_engine_free(engine);
    caliper_pointcloud_free(scan);
    return true;
}

int main(void) {
    printf("Running C Bridge Unit Tests...\n");
    bool insp_ok = test_bridge_inspection();
    printf("  Inspection pipeline: %s\n", insp_ok ? "PASSED" : "FAILED");
    bool gdt_ok = test_bridge_gdt();
    printf("  GD&T evaluation:     %s\n", gdt_ok ? "PASSED" : "FAILED");

    if (insp_ok && gdt_ok) {
        printf("All C Bridge tests passed successfully.\n");
        return 0;
    }
    return 1;
}
