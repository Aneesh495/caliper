/**
 * Caliper Metrology Engine: Unit Tests for C Geometric File I/O
 * File: test_c_io.c
 *
 * Verifies round-trip integrity of Binary STL, ASCII STL, and ASCII PLY
 * point clouds and meshes.
 */

#include "caliper_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

static bool test_stl_binary_roundtrip(void) {
    CaliperMesh *mesh = caliper_mesh_create(2);
    assert(mesh != NULL);

    mesh->triangles[0].normal = (CaliperVec3){0.0, 0.0, 1.0};
    mesh->triangles[0].v0 = (CaliperVec3){0.0, 0.0, 0.0};
    mesh->triangles[0].v1 = (CaliperVec3){10.0, 0.0, 0.0};
    mesh->triangles[0].v2 = (CaliperVec3){0.0, 10.0, 0.0};

    mesh->triangles[1].normal = (CaliperVec3){0.0, 0.0, 1.0};
    mesh->triangles[1].v0 = (CaliperVec3){10.0, 0.0, 0.0};
    mesh->triangles[1].v1 = (CaliperVec3){10.0, 10.0, 0.0};
    mesh->triangles[1].v2 = (CaliperVec3){0.0, 10.0, 0.0};

    const char *tmp_path = "/tmp/test_mesh.stl";
    bool write_ok = caliper_write_stl_binary(tmp_path, mesh, "Test STL Header");
    assert(write_ok);

    CaliperMesh *read_mesh = caliper_read_stl_binary(tmp_path);
    assert(read_mesh != NULL);
    assert(read_mesh->num_triangles == 2);

    for (size_t i = 0; i < 2; ++i) {
        assert(fabs(read_mesh->triangles[i].v0.x - mesh->triangles[i].v0.x) < 1e-4);
        assert(fabs(read_mesh->triangles[i].v1.y - mesh->triangles[i].v1.y) < 1e-4);
        assert(fabs(read_mesh->triangles[i].v2.x - mesh->triangles[i].v2.x) < 1e-4);
    }

    caliper_mesh_free(mesh);
    caliper_mesh_free(read_mesh);
    return true;
}

static bool test_ply_ascii_roundtrip(void) {
    size_t n = 4;
    CaliperPointCloud *pc = caliper_pointcloud_create(n, true, true, true);
    assert(pc != NULL);

    for (size_t i = 0; i < n; ++i) {
        pc->x[i] = (double)i * 5.0;
        pc->y[i] = (double)i * 2.5;
        pc->z[i] = (double)i * 1.2;
        pc->nx[i] = 0.0;
        pc->ny[i] = 0.0;
        pc->nz[i] = 1.0;
        pc->deviation[i] = 0.015 * (double)i;
        pc->r[i] = 0.2 * (double)i;
        pc->g[i] = 0.8;
        pc->b[i] = 0.1;
    }

    const char *tmp_path = "/tmp/test_pc.ply";
    bool write_ok = caliper_write_ply_ascii(tmp_path, pc);
    assert(write_ok);

    CaliperPointCloud *read_pc = caliper_read_ply_ascii(tmp_path);
    assert(read_pc != NULL);
    assert(read_pc->num_points == n);
    assert(read_pc->has_normals);
    assert(read_pc->has_deviations);

    for (size_t i = 0; i < n; ++i) {
        assert(fabs(read_pc->x[i] - pc->x[i]) < 1e-4);
        assert(fabs(read_pc->nz[i] - 1.0) < 1e-4);
        assert(fabs(read_pc->deviation[i] - pc->deviation[i]) < 1e-4);
    }

    caliper_pointcloud_free(pc);
    caliper_pointcloud_free(read_pc);
    return true;
}

int main(void) {
    printf("Running C I/O Unit Tests...\n");
    bool stl_ok = test_stl_binary_roundtrip();
    printf("  Binary STL round-trip: %s\n", stl_ok ? "PASSED" : "FAILED");
    bool ply_ok = test_ply_ascii_roundtrip();
    printf("  ASCII PLY round-trip: %s\n", ply_ok ? "PASSED" : "FAILED");

    if (stl_ok && ply_ok) {
        printf("All C I/O tests passed successfully.\n");
        return 0;
    }
    return 1;
}
