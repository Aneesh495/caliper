/**
 * Caliper Metrology Engine: High-Speed Geometric File I/O Implementation
 * File: caliper_io.c
 *
 * Implements high-throughput binary and ASCII parsers and serializers for
 * STL, PLY, and CSV formats.
 */

#include "caliper_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

CaliperPointCloud* caliper_pointcloud_create(size_t capacity, bool with_normals, bool with_deviations, bool with_colors) {
    CaliperPointCloud *pc = (CaliperPointCloud*)calloc(1, sizeof(CaliperPointCloud));
    if (!pc) return NULL;

    pc->num_points = capacity;
    pc->has_normals = with_normals;
    pc->has_deviations = with_deviations;
    pc->has_colors = with_colors;

    if (capacity > 0) {
        pc->x = (double*)malloc(capacity * sizeof(double));
        pc->y = (double*)malloc(capacity * sizeof(double));
        pc->z = (double*)malloc(capacity * sizeof(double));
        if (!pc->x || !pc->y || !pc->z) {
            caliper_pointcloud_free(pc);
            return NULL;
        }

        if (with_normals) {
            pc->nx = (double*)malloc(capacity * sizeof(double));
            pc->ny = (double*)malloc(capacity * sizeof(double));
            pc->nz = (double*)malloc(capacity * sizeof(double));
            if (!pc->nx || !pc->ny || !pc->nz) {
                caliper_pointcloud_free(pc);
                return NULL;
            }
        }

        if (with_deviations) {
            pc->deviation = (double*)malloc(capacity * sizeof(double));
            if (!pc->deviation) {
                caliper_pointcloud_free(pc);
                return NULL;
            }
        }

        if (with_colors) {
            pc->r = (double*)malloc(capacity * sizeof(double));
            pc->g = (double*)malloc(capacity * sizeof(double));
            pc->b = (double*)malloc(capacity * sizeof(double));
            if (!pc->r || !pc->g || !pc->b) {
                caliper_pointcloud_free(pc);
                return NULL;
            }
        }
    }
    return pc;
}

void caliper_pointcloud_free(CaliperPointCloud *pc) {
    if (!pc) return;
    free(pc->x);
    free(pc->y);
    free(pc->z);
    free(pc->nx);
    free(pc->ny);
    free(pc->nz);
    free(pc->deviation);
    free(pc->r);
    free(pc->g);
    free(pc->b);
    free(pc);
}

CaliperMesh* caliper_mesh_create(size_t num_triangles) {
    CaliperMesh *mesh = (CaliperMesh*)calloc(1, sizeof(CaliperMesh));
    if (!mesh) return NULL;

    mesh->num_triangles = num_triangles;
    if (num_triangles > 0) {
        mesh->triangles = (CaliperTriangle*)calloc(num_triangles, sizeof(CaliperTriangle));
        if (!mesh->triangles) {
            free(mesh);
            return NULL;
        }
    }
    return mesh;
}

void caliper_mesh_free(CaliperMesh *mesh) {
    if (!mesh) return;
    free(mesh->triangles);
    free(mesh);
}

/* ----------------------------------------------------------------------------
 * Binary STL I/O
 * ---------------------------------------------------------------------------- */

#pragma pack(push, 1)
typedef struct {
    float normal[3];
    float v0[3];
    float v1[3];
    float v2[3];
    uint16_t attr_byte_count;
} StlBinaryTriangle;
#pragma pack(pop)

bool caliper_write_stl_binary(const char *filepath, const CaliperMesh *mesh, const char *header_info) {
    if (!filepath || !mesh) return false;

    FILE *f = fopen(filepath, "wb");
    if (!f) return false;

    char header[80];
    memset(header, 0, sizeof(header));
    if (header_info) {
        strncpy(header, header_info, sizeof(header) - 1);
    } else {
        snprintf(header, sizeof(header), "Caliper Metrology Engine Binary STL Export");
    }
    fwrite(header, 1, 80, f);

    uint32_t count = (uint32_t)mesh->num_triangles;
    fwrite(&count, sizeof(uint32_t), 1, f);

    for (size_t i = 0; i < mesh->num_triangles; ++i) {
        const CaliperTriangle *t = &mesh->triangles[i];
        StlBinaryTriangle st;
        st.normal[0] = (float)t->normal.x;
        st.normal[1] = (float)t->normal.y;
        st.normal[2] = (float)t->normal.z;

        st.v0[0] = (float)t->v0.x;
        st.v0[1] = (float)t->v0.y;
        st.v0[2] = (float)t->v0.z;

        st.v1[0] = (float)t->v1.x;
        st.v1[1] = (float)t->v1.y;
        st.v1[2] = (float)t->v1.z;

        st.v2[0] = (float)t->v2.x;
        st.v2[1] = (float)t->v2.y;
        st.v2[2] = (float)t->v2.z;

        st.attr_byte_count = 0;
        fwrite(&st, sizeof(StlBinaryTriangle), 1, f);
    }

    fclose(f);
    return true;
}

CaliperMesh* caliper_read_stl_binary(const char *filepath) {
    if (!filepath) return NULL;

    FILE *f = fopen(filepath, "rb");
    if (!f) return NULL;

    char header[80];
    if (fread(header, 1, 80, f) != 80) {
        fclose(f);
        return NULL;
    }

    uint32_t count = 0;
    if (fread(&count, sizeof(uint32_t), 1, f) != 1) {
        fclose(f);
        return NULL;
    }

    CaliperMesh *mesh = caliper_mesh_create(count);
    if (!mesh) {
        fclose(f);
        return NULL;
    }

    for (uint32_t i = 0; i < count; ++i) {
        StlBinaryTriangle st;
        if (fread(&st, sizeof(StlBinaryTriangle), 1, f) != 1) {
            caliper_mesh_free(mesh);
            fclose(f);
            return NULL;
        }

        CaliperTriangle *t = &mesh->triangles[i];
        t->normal.x = (double)st.normal[0];
        t->normal.y = (double)st.normal[1];
        t->normal.z = (double)st.normal[2];

        t->v0.x = (double)st.v0[0];
        t->v0.y = (double)st.v0[1];
        t->v0.z = (double)st.v0[2];

        t->v1.x = (double)st.v1[0];
        t->v1.y = (double)st.v1[1];
        t->v1.z = (double)st.v1[2];

        t->v2.x = (double)st.v2[0];
        t->v2.y = (double)st.v2[1];
        t->v2.z = (double)st.v2[2];
    }

    fclose(f);
    return mesh;
}

/* ----------------------------------------------------------------------------
 * ASCII STL Writer
 * ---------------------------------------------------------------------------- */

bool caliper_write_stl_ascii(const char *filepath, const CaliperMesh *mesh, const char *solid_name) {
    if (!filepath || !mesh) return false;

    FILE *f = fopen(filepath, "w");
    if (!f) return false;

    const char *name = solid_name ? solid_name : "caliper_cad_model";
    fprintf(f, "solid %s\n", name);

    for (size_t i = 0; i < mesh->num_triangles; ++i) {
        const CaliperTriangle *t = &mesh->triangles[i];
        fprintf(f, "  facet normal %.6e %.6e %.6e\n", t->normal.x, t->normal.y, t->normal.z);
        fprintf(f, "    outer loop\n");
        fprintf(f, "      vertex %.6e %.6e %.6e\n", t->v0.x, t->v0.y, t->v0.z);
        fprintf(f, "      vertex %.6e %.6e %.6e\n", t->v1.x, t->v1.y, t->v1.z);
        fprintf(f, "      vertex %.6e %.6e %.6e\n", t->v2.x, t->v2.y, t->v2.z);
        fprintf(f, "    endloop\n");
        fprintf(f, "  endfacet\n");
    }

    fprintf(f, "endsolid %s\n", name);
    fclose(f);
    return true;
}

/* ----------------------------------------------------------------------------
 * PLY Point Cloud I/O
 * ---------------------------------------------------------------------------- */

bool caliper_write_ply_ascii(const char *filepath, const CaliperPointCloud *pc) {
    if (!filepath || !pc) return false;

    FILE *f = fopen(filepath, "w");
    if (!f) return false;

    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
    fprintf(f, "comment Caliper Metrology Engine Inspection Point Cloud\n");
    fprintf(f, "element vertex %zu\n", pc->num_points);
    fprintf(f, "property double x\n");
    fprintf(f, "property double y\n");
    fprintf(f, "property double z\n");

    if (pc->has_normals) {
        fprintf(f, "property double nx\n");
        fprintf(f, "property double ny\n");
        fprintf(f, "property double nz\n");
    }
    if (pc->has_deviations) {
        fprintf(f, "property double deviation_mm\n");
    }
    if (pc->has_colors) {
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
    }
    fprintf(f, "end_header\n");

    for (size_t i = 0; i < pc->num_points; ++i) {
        fprintf(f, "%.6f %.6f %.6f", pc->x[i], pc->y[i], pc->z[i]);
        if (pc->has_normals) {
            fprintf(f, " %.6f %.6f %.6f", pc->nx[i], pc->ny[i], pc->nz[i]);
        }
        if (pc->has_deviations) {
            fprintf(f, " %.6f", pc->deviation[i]);
        }
        if (pc->has_colors) {
            int r_int = (int)(pc->r[i] * 255.0);
            int g_int = (int)(pc->g[i] * 255.0);
            int b_int = (int)(pc->b[i] * 255.0);
            if (r_int < 0) r_int = 0; if (r_int > 255) r_int = 255;
            if (g_int < 0) g_int = 0; if (g_int > 255) g_int = 255;
            if (b_int < 0) b_int = 0; if (b_int > 255) b_int = 255;
            fprintf(f, " %d %d %d", r_int, g_int, b_int);
        }
        fprintf(f, "\n");
    }

    fclose(f);
    return true;
}

CaliperPointCloud* caliper_read_ply_ascii(const char *filepath) {
    if (!filepath) return NULL;

    FILE *f = fopen(filepath, "r");
    if (!f) return NULL;

    char line[256];
    size_t num_verts = 0;
    bool has_norm = false;
    bool has_dev = false;
    bool has_col = false;

    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "element vertex", 14) == 0) {
            sscanf(line, "element vertex %zu", &num_verts);
        } else if (strncmp(line, "property double nx", 18) == 0 || strncmp(line, "property float nx", 17) == 0) {
            has_norm = true;
        } else if (strstr(line, "deviation")) {
            has_dev = true;
        } else if (strstr(line, "red")) {
            has_col = true;
        } else if (strncmp(line, "end_header", 10) == 0) {
            break;
        }
    }

    if (num_verts == 0) {
        fclose(f);
        return NULL;
    }

    CaliperPointCloud *pc = caliper_pointcloud_create(num_verts, has_norm, has_dev, has_col);
    if (!pc) {
        fclose(f);
        return NULL;
    }

    for (size_t i = 0; i < num_verts; ++i) {
        double px, py, pz;
        if (fscanf(f, "%lf %lf %lf", &px, &py, &pz) != 3) break;
        pc->x[i] = px;
        pc->y[i] = py;
        pc->z[i] = pz;

        if (has_norm) {
            double nx, ny, nz;
            if (fscanf(f, "%lf %lf %lf", &nx, &ny, &nz) == 3) {
                pc->nx[i] = nx;
                pc->ny[i] = ny;
                pc->nz[i] = nz;
            }
        }
        if (has_dev) {
            double dev;
            if (fscanf(f, "%lf", &dev) == 1) {
                pc->deviation[i] = dev;
            }
        }
        if (has_col) {
            int r, g, b;
            if (fscanf(f, "%d %d %d", &r, &g, &b) == 3) {
                pc->r[i] = (double)r / 255.0;
                pc->g[i] = (double)g / 255.0;
                pc->b[i] = (double)b / 255.0;
            }
        }
    }

    fclose(f);
    return pc;
}

/* ----------------------------------------------------------------------------
 * CSV Quality Report Writer
 * ---------------------------------------------------------------------------- */

bool caliper_write_csv_report(const char *filepath, const CaliperPointCloud *pc) {
    if (!filepath || !pc) return false;

    FILE *f = fopen(filepath, "w");
    if (!f) return false;

    fprintf(f, "point_index,x_mm,y_mm,z_mm");
    if (pc->has_normals) {
        fprintf(f, ",nx,ny,nz");
    }
    if (pc->has_deviations) {
        fprintf(f, ",deviation_mm");
    }
    fprintf(f, "\n");

    for (size_t i = 0; i < pc->num_points; ++i) {
        fprintf(f, "%zu,%.6f,%.6f,%.6f", i, pc->x[i], pc->y[i], pc->z[i]);
        if (pc->has_normals) {
            fprintf(f, ",%.6f,%.6f,%.6f", pc->nx[i], pc->ny[i], pc->nz[i]);
        }
        if (pc->has_deviations) {
            fprintf(f, ",%.6f", pc->deviation[i]);
        }
        fprintf(f, "\n");
    }

    fclose(f);
    return true;
}
