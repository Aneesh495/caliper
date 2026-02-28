/**
 * Caliper Metrology Engine: C Bridge and Execution Harness
 * File: caliper_bridge.h
 *
 * High-level native C orchestration interface connecting file I/O,
 * data buffers, and the data-parallel Futhark computational kernels.
 */

#ifndef CALIPER_BRIDGE_H
#define CALIPER_BRIDGE_H

#include "caliper_io.h"
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double mean_deviation_mm;
    double std_dev_mm;
    double rmse_mm;
    double max_positive_mm;
    double max_negative_mm;
    double in_tolerance_pct;
    size_t evaluated_points;
    bool inspection_passed;
} CaliperInspectionSummary;

typedef struct {
    double measured_flatness_mm;
    double flatness_tolerance_mm;
    bool flatness_passed;

    double measured_cylindricity_mm;
    double cylindricity_tolerance_mm;
    bool cylindricity_passed;

    double measured_perpendicularity_mm;
    double perpendicularity_tolerance_mm;
    bool perpendicularity_passed;

    double measured_position_mm;
    double position_tolerance_mm;
    double position_bonus_mm;
    bool position_passed;

    bool overall_gdt_passed;
} CaliperGdtSummary;

typedef struct CaliperEngine CaliperEngine;

/* Engine lifecycle */
CaliperEngine* caliper_engine_create(int num_threads);
void caliper_engine_free(CaliperEngine *engine);

/* Inspection pipelines */
bool caliper_run_inspection(
    CaliperEngine *engine,
    CaliperPointCloud *scan_pc,
    const CaliperMesh *cad_mesh,
    double tolerance_mm,
    CaliperInspectionSummary *summary_out
);

bool caliper_evaluate_gdt_sample(
    CaliperEngine *engine,
    const CaliperPointCloud *scan_pc,
    double flatness_tol_mm,
    double cylindricity_tol_mm,
    CaliperGdtSummary *gdt_out
);

#ifdef __cplusplus
}
#endif

#endif /* CALIPER_BRIDGE_H */
