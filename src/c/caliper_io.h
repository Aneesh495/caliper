/**
 * Caliper Metrology Engine: High-Speed Geometric File I/O
 * File: caliper_io.h
 *
 * Binary and ASCII parsers and writers for standard industrial CAD and
 * coordinate metrology interchange formats: Binary/ASCII STL, PLY, and CSV.
 */

#ifndef CALIPER_IO_H
#define CALIPER_IO_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x, y, z;
} CaliperVec3;

typedef struct {
    CaliperVec3 v0;
    CaliperVec3 v1;
    CaliperVec3 v2;
    CaliperVec3 normal;
} CaliperTriangle;

typedef struct {
    size_t num_points;
    double *x;
    double *y;
    double *z;
    double *nx;
    double *ny;
    double *nz;
    double *deviation;
    double *r;
    double *g;
    double *b;
    bool has_normals;
    bool has_deviations;
    bool has_colors;
} CaliperPointCloud;

typedef struct {
    size_t num_triangles;
    CaliperTriangle *triangles;
} CaliperMesh;

/* Point cloud memory management */
CaliperPointCloud* caliper_pointcloud_create(size_t capacity, bool with_normals, bool with_deviations, bool with_colors);
void caliper_pointcloud_free(CaliperPointCloud *pc);

/* Mesh memory management */
CaliperMesh* caliper_mesh_create(size_t num_triangles);
void caliper_mesh_free(CaliperMesh *mesh);

/* STL I/O */
bool caliper_write_stl_binary(const char *filepath, const CaliperMesh *mesh, const char *header_info);
bool caliper_write_stl_ascii(const char *filepath, const CaliperMesh *mesh, const char *solid_name);
CaliperMesh* caliper_read_stl_binary(const char *filepath);

/* PLY I/O */
bool caliper_write_ply_ascii(const char *filepath, const CaliperPointCloud *pc);
bool caliper_write_ply_binary(const char *filepath, const CaliperPointCloud *pc);
CaliperPointCloud* caliper_read_ply_ascii(const char *filepath);

/* CSV Metrology Log I/O */
bool caliper_write_csv_report(const char *filepath, const CaliperPointCloud *pc);

#ifdef __cplusplus
}
#endif

#endif /* CALIPER_IO_H */
