/**
 * Caliper Metrology Engine: Standalone High-Performance Native CLI
 * File: caliper_main.c
 *
 * Standalone native terminal application for high-throughput 3D optical
 * inspection, CAD deviation analysis, and ASME Y14.5 GD&T compliance verification.
 */

#include "caliper_io.h"
#include "caliper_bridge.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define CALIPER_VERSION "1.4.0"

static void print_banner(void) {
    printf("\n");
    printf("  ===============================================================\n");
    printf("  Caliper Metrology Engine (v%s)\n", CALIPER_VERSION);
    printf("  High-Throughput Optical Metrology & GD&T Inspection Platform\n");
    printf("  ===============================================================\n\n");
}

static void print_usage(const char *prog_name) {
    print_banner();
    printf("Usage:\n");
    printf("  %s inspect <scan.ply> <cad.stl> [options]\n", prog_name);
    printf("  %s gdt <scan.ply> [options]\n", prog_name);
    printf("  %s info <geometry.stl|ply>\n", prog_name);
    printf("  %s version\n\n", prog_name);
    printf("Options:\n");
    printf("  --tolerance <mm>     Set nominal bilateral tolerance (default: 0.050 mm)\n");
    printf("  --flatness <mm>      Set flatness tolerance limit (default: 0.025 mm)\n");
    printf("  --cylindricity <mm>  Set cylindricity tolerance limit (default: 0.030 mm)\n");
    printf("  --output <file.ply>  Output false-color deviation point cloud\n");
    printf("  --csv <report.csv>   Export statistical inspection report to CSV\n\n");
}

static int cmd_inspect(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Error: inspect requires <scan.ply> and <cad.stl>\n");
        return 1;
    }

    const char *scan_path = argv[2];
    const char *cad_path = argv[3];
    double tolerance_mm = 0.050;
    const char *output_path = NULL;
    const char *csv_path = NULL;

    for (int i = 4; i < argc; ++i) {
        if (strcmp(argv[i], "--tolerance") == 0 && i + 1 < argc) {
            tolerance_mm = atof(argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_path = argv[++i];
        } else if (strcmp(argv[i], "--csv") == 0 && i + 1 < argc) {
            csv_path = argv[++i];
        }
    }

    print_banner();
    printf("Loading scan point cloud: %s\n", scan_path);
    CaliperPointCloud *scan_pc = caliper_read_ply_ascii(scan_path);
    if (!scan_pc) {
        fprintf(stderr, "Failed to load scan point cloud: %s\n", scan_path);
        return 1;
    }
    printf("  Loaded %zu scan vertices.\n", scan_pc->num_points);

    printf("Loading reference CAD model: %s\n", cad_path);
    CaliperMesh *cad_mesh = caliper_read_stl_binary(cad_path);
    if (!cad_mesh) {
        fprintf(stderr, "Failed to load CAD model: %s\n", cad_path);
        caliper_pointcloud_free(scan_pc);
        return 1;
    }
    printf("  Loaded %zu reference triangles.\n", cad_mesh->num_triangles);

    printf("\nInitializing Caliper Metrology Acceleration Engine...\n");
    CaliperEngine *engine = caliper_engine_create(8);
    CaliperInspectionSummary summary;

    printf("Executing parallel Signed Distance Field deviation analysis...\n");
    bool ok = caliper_run_inspection(engine, scan_pc, cad_mesh, tolerance_mm, &summary);
    if (!ok) {
        fprintf(stderr, "Inspection pipeline execution failed.\n");
        caliper_engine_free(engine);
        caliper_mesh_free(cad_mesh);
        caliper_pointcloud_free(scan_pc);
        return 1;
    }

    printf("\n--- Dimensional Inspection Certificate ---\n");
    printf("Evaluated Points:       %zu\n", summary.evaluated_points);
    printf("Specified Tolerance:    +/- %.4f mm\n", tolerance_mm);
    printf("In-Tolerance Ratio:     %.2f %%\n", summary.in_tolerance_pct);
    printf("Mean Deviation:         %.4f mm\n", summary.mean_deviation_mm);
    printf("Standard Deviation (s): %.4f mm\n", summary.std_dev_mm);
    printf("Root Mean Square (RMS): %.4f mm\n", summary.rmse_mm);
    printf("Max Excess Material:   +%.4f mm\n", summary.max_positive_mm);
    printf("Max Missing Material:   %.4f mm\n", summary.max_negative_mm);
    printf("Inspection Decision:    %s\n", summary.inspection_passed ? "[PASS]" : "[REWORK/SCRAP]");

    if (output_path) {
        printf("\nExporting false-color deviation map: %s\n", output_path);
        caliper_write_ply_ascii(output_path, scan_pc);
    }
    if (csv_path) {
        printf("Exporting metrology log table: %s\n", csv_path);
        caliper_write_csv_report(csv_path, scan_pc);
    }

    caliper_engine_free(engine);
    caliper_mesh_free(cad_mesh);
    caliper_pointcloud_free(scan_pc);
    return summary.inspection_passed ? 0 : 2;
}

static int cmd_gdt(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Error: gdt requires <scan.ply>\n");
        return 1;
    }

    const char *scan_path = argv[2];
    double flatness_tol = 0.025;
    double cylindricity_tol = 0.030;

    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "--flatness") == 0 && i + 1 < argc) {
            flatness_tol = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cylindricity") == 0 && i + 1 < argc) {
            cylindricity_tol = atof(argv[++i]);
        }
    }

    print_banner();
    printf("Loading scan point cloud: %s\n", scan_path);
    CaliperPointCloud *scan_pc = caliper_read_ply_ascii(scan_path);
    if (!scan_pc) {
        fprintf(stderr, "Failed to load scan point cloud: %s\n", scan_path);
        return 1;
    }

    CaliperEngine *engine = caliper_engine_create(8);
    CaliperGdtSummary gdt;
    caliper_evaluate_gdt_sample(engine, scan_pc, flatness_tol, cylindricity_tol, &gdt);

    printf("\n--- ASME Y14.5 GD&T Verification Summary ---\n");
    printf("Flatness:        Measured = %.4f mm | Limit = %.4f mm | %s\n",
           gdt.measured_flatness_mm, gdt.flatness_tolerance_mm, gdt.flatness_passed ? "PASS" : "FAIL");
    printf("Cylindricity:    Measured = %.4f mm | Limit = %.4f mm | %s\n",
           gdt.measured_cylindricity_mm, gdt.cylindricity_tolerance_mm, gdt.cylindricity_passed ? "PASS" : "FAIL");
    printf("Perpendicular:   Measured = %.4f mm | Limit = %.4f mm | %s\n",
           gdt.measured_perpendicularity_mm, gdt.perpendicularity_tolerance_mm, gdt.perpendicularity_passed ? "PASS" : "FAIL");
    printf("True Position:   Measured = %.4f mm | Limit = %.4f mm (Bonus: +%.4f mm) | %s\n",
           gdt.measured_position_mm, gdt.position_tolerance_mm, gdt.position_bonus_mm, gdt.position_passed ? "PASS" : "FAIL");
    printf("Overall Result:  %s\n", gdt.overall_gdt_passed ? "[COMPLIANT]" : "[NON-COMPLIANT]");

    caliper_engine_free(engine);
    caliper_pointcloud_free(scan_pc);
    return gdt.overall_gdt_passed ? 0 : 2;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 0;
    }

    if (strcmp(argv[1], "inspect") == 0) {
        return cmd_inspect(argc, argv);
    } else if (strcmp(argv[1], "gdt") == 0) {
        return cmd_gdt(argc, argv);
    } else if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0) {
        printf("Caliper Metrology Engine version %s\n", CALIPER_VERSION);
        return 0;
    } else {
        print_usage(argv[0]);
        return 1;
    }
}
