#ifndef IMPORTER_H
#define IMPORTER_H

#include <vector>

#include "aarect.h"
#include "bvh.h"
#include "hittable_list.h"
#include "material.h"
#include "transform.h"
#include "triangle.h"
#include "utility.h"
#include "sphere.h"

// Importers. Header only libraries taken from Github
#include <iostream>
#include <fstream>
#include <string>
#include "../ImportedResource/npy.hpp"
#include "../ImportedResource/stl_reader_sebastian.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "../ImportedResource/tiny_obj_loader.h"

extern bool TrustNormal;

shared_ptr<TriangleMesh> obj_loader (std :: string fname, Transform tran) {

    float minX = 100000, minY = 100000, minZ = 100000;
    float maxX = -100000, maxY = -100000 ,maxZ = -100000;

    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(fname)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }
    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();

    // In our cases, there is usually only a single shape. Thus we can simply call shapes[0]
    const int number_of_shapes = int(shapes.size());
    const int number_of_triangles = int(shapes[0].mesh.num_face_vertices.size());
    const int number_of_vertices = int(shapes[0].mesh.indices.size());

    std :: vector <int> vertices_index (number_of_vertices);
    std :: vector <vec3> p (number_of_vertices);
    std :: vector <vec3> n (number_of_vertices); 

    size_t index_offset = 0;
    for (size_t f = 0; f < number_of_triangles; f++) {
        // Loop over vertices in the face.
        size_t fv = 3;
        for (size_t v = 0; v < fv; v++) {
            // access to vertex
            tinyobj::index_t idx = shapes[0].mesh.indices[index_offset + v];
            vertices_index[index_offset + v] = int(idx.vertex_index);

            tinyobj::real_t vx = attrib.vertices[3*size_t(idx.vertex_index)+0];
            tinyobj::real_t vy = attrib.vertices[3*size_t(idx.vertex_index)+1];
            tinyobj::real_t vz = attrib.vertices[3*size_t(idx.vertex_index)+2];
            p[idx.vertex_index] = vec3(float(vx), float(vy), float(vz));

            // Size Checker
            if (vx > maxX) maxX = vx;
            if (vy > maxY) maxY = vy;
            if (vz > maxZ) maxZ = vz;
            if (vx < minX) minX = vx;
            if (vy < minY) minY = vy;
            if (vz < minZ) minZ = vz;

            // Check if normal_index is zero or positive. negative = no normal data
            if (idx.normal_index >= 0 && TrustNormal) {
                tinyobj::real_t nx = attrib.normals[3*size_t(idx.normal_index)+0];
                tinyobj::real_t ny = attrib.normals[3*size_t(idx.normal_index)+1];
                tinyobj::real_t nz = attrib.normals[3*size_t(idx.normal_index)+2];
                n[idx.normal_index] = vec3(float(nx), float(ny), float(nz));
            }
        }
        index_offset += fv;
    }
    auto mesh = std :: shared_ptr<TriangleMesh> (new TriangleMesh(tran, number_of_triangles, &vertices_index[0], number_of_vertices, &p[0], &n[0]));

    std :: cout << "* Minimum and maximum extent of each model" << std :: endl;
    std :: cout << "min point : " << point3(minX, minY, minZ) << std :: endl;
    std :: cout << "max point : " << point3(maxX, maxY, maxZ) << std :: endl;

    return mesh;
}

hittable_list obj_world_builder (std :: vector <std :: string> fnames, std :: vector <shared_ptr<material>> mats, 
std :: vector <Transform> trans, std :: vector <shared_ptr<hittable>> lights, std :: vector <shared_ptr<hittable>> walls, unsigned int& seed) {
    
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;

    for (int i = 0; i < fnames.size(); i++) {
        const auto& fname = fnames[i];
        shared_ptr<material> mat = mats[i];

        auto mesh = obj_loader(fname, trans[i]);
        for (int j = 0; j < mesh -> nTriangles; j++) {
            shared_ptr<hittable> tri = make_shared<Triangle>(mesh, j, mat);
            world.add(tri);
        }
    }

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    for (const auto& i : lights) {
        world_bvh.add(i);
    }
    for (const auto& i : walls) {
        world_bvh.add(i);
    }

    return world_bvh;
}

hittable_list stl_world_builder (std :: vector <std :: string> fnames, std :: vector <shared_ptr<material>> mats, 
std :: vector <Transform> trans, std :: vector <shared_ptr<hittable>> lights, std :: vector <shared_ptr<hittable>> walls, unsigned int& seed) {
    
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;
    vec3 v0, v1, v2;

    // // Just in case of completely empty model data
    // bool empty_model = true;

    // On send thought, these maximum&minimum probably should be universal
    float minX = 100000, minY = 100000, minZ = 100000;
    float maxX = -100000, maxY = -100000 ,maxZ = -100000;

    for (int i = 0; i < fnames.size(); i++) {
        const auto& fname = fnames[i];
        shared_ptr<material> mat = mats[i];

        try {
            stl_reader :: StlMesh <float, unsigned int> mesh (fname);
            for(size_t itri = 0; itri < mesh.num_tris(); ++itri) {
                const float* c0 = mesh.tri_corner_coords (itri, 0);
                v0 = vec3(c0[0], c0[1], c0[2]);
                const float* c1 = mesh.tri_corner_coords (itri, 1);
                v1 = vec3(c1[0], c1[1], c1[2]);
                const float* c2 = mesh.tri_corner_coords (itri, 2);
                v2 = vec3(c2[0], c2[1], c2[2]);

                // Size Checker
                if (c0[0] > maxX) maxX = c0[0];
                if (c0[1] > maxY) maxY = c0[1];
                if (c0[2] > maxZ) maxZ = c0[2];
                if (c0[0] < minX) minX = c0[0];
                if (c0[1] < minY) minY = c0[1];
                if (c0[2] < minZ) minZ = c0[2];
                if (c1[0] > maxX) maxX = c1[0];
                if (c1[1] > maxY) maxY = c1[1];
                if (c1[2] > maxZ) maxZ = c1[2];
                if (c1[0] < minX) minX = c1[0];
                if (c1[1] < minY) minY = c1[1];
                if (c1[2] < minZ) minZ = c1[2];
                if (c2[0] > maxX) maxX = c2[0];
                if (c2[1] > maxY) maxY = c2[1];
                if (c2[2] > maxZ) maxZ = c2[2];
                if (c2[0] < minX) minX = c2[0];
                if (c2[1] < minY) minY = c2[1];
                if (c2[2] < minZ) minZ = c2[2];

                // std :: cout << "v0 : " << v0 << std :: endl;
                // std :: cout << "v1 : " << v1 << std :: endl;
                // std :: cout << "v2 : " << v2 << std :: endl;

                if (TrustNormal) {
                    const float* n = mesh.tri_normal (itri);
                    vec3 norm = vec3(n[0], n[1], n[2]);
                    world.add(make_shared<triangle>(v0, v1, v2, norm, mat));
                }
                else {
                    world.add(make_shared<triangle>(v0, v1, v2, mat));
                    // std :: cout << "A triangle got made/added" << std :: endl;
                }
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));
    // std :: cout << "Max of worldbox : " << world_bvh.WorldBound().max() << std :: endl;
    // std :: cout << "Min of worldbox : " << world_bvh.WorldBound().min() << std :: endl;

    std :: cout << "* Minimum and maximum extent of the whole model" << std :: endl;
    std :: cout << "min point : " << point3(minX, minY, minZ) << std :: endl;
    std :: cout << "max point : " << point3(maxX, maxY, maxZ) << std :: endl;

    for (const auto& i : lights) {
        world_bvh.add(i);
    }
    for (const auto& i : walls) {
        world_bvh.add(i);
    }

    return world_bvh;
}

hittable_list stl_world_builder_with_assets (std :: vector <std :: string> fnames, std :: vector <shared_ptr<material>> mats, 
std :: vector <Transform> trans, std :: vector <shared_ptr<hittable>> lights, std :: vector <shared_ptr<hittable>> walls, 
std :: vector <std :: string> asset_names, std :: vector <shared_ptr<material>> asset_mats, std :: vector <Transform> asset_trans, unsigned int& seed) {
    
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;
    vec3 v0, v1, v2;

    // On send thought, these maximum&minimum probably should be universal
    float minX = 100000, minY = 100000, minZ = 100000;
    float maxX = -100000, maxY = -100000 ,maxZ = -100000;

    // On second thought, these maximum&minimum probably should be universal
    float minAssetsX = 100000, minAssetsY = 100000, minAssetsZ = 100000;
    float maxAssetsX = -100000, maxAssetsY = -100000 ,maxAssetsZ = -100000;

    // If it is coredump here that probably is wrong model name
    for (int i = 0; i < fnames.size(); i++) {
        const auto& fname = fnames[i];
        shared_ptr<material> mat = mats[i];

        try {
            stl_reader :: StlMesh <float, unsigned int> mesh (fname);
            for(size_t itri = 0; itri < mesh.num_tris(); ++itri) {
                const float* c0 = mesh.tri_corner_coords (itri, 0);
                v0 = vec3(c0[0], c0[1], c0[2]);
                const float* c1 = mesh.tri_corner_coords (itri, 1);
                v1 = vec3(c1[0], c1[1], c1[2]);
                const float* c2 = mesh.tri_corner_coords (itri, 2);
                v2 = vec3(c2[0], c2[1], c2[2]);

                // Size Checker
                if (c0[0] > maxX) maxX = c0[0];
                if (c0[1] > maxY) maxY = c0[1];
                if (c0[2] > maxZ) maxZ = c0[2];
                if (c0[0] < minX) minX = c0[0];
                if (c0[1] < minY) minY = c0[1];
                if (c0[2] < minZ) minZ = c0[2];
                if (c1[0] > maxX) maxX = c1[0];
                if (c1[1] > maxY) maxY = c1[1];
                if (c1[2] > maxZ) maxZ = c1[2];
                if (c1[0] < minX) minX = c1[0];
                if (c1[1] < minY) minY = c1[1];
                if (c1[2] < minZ) minZ = c1[2];
                if (c2[0] > maxX) maxX = c2[0];
                if (c2[1] > maxY) maxY = c2[1];
                if (c2[2] > maxZ) maxZ = c2[2];
                if (c2[0] < minX) minX = c2[0];
                if (c2[1] < minY) minY = c2[1];
                if (c2[2] < minZ) minZ = c2[2];


                if (TrustNormal) {
                    const float* n = mesh.tri_normal (itri);
                    vec3 norm = vec3(n[0], n[1], n[2]);
                    world.add(make_shared<triangle>(v0, v1, v2, norm, mat));
                }
                else {
                    world.add(make_shared<triangle>(v0, v1, v2, mat));
                    // std :: cout << "A triangle got made/added" << std :: endl;
                }
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    std::cout << "Pass checkpoint 1" << std::endl;

    for (int i = 0; i < asset_names.size(); i++) {
        const auto& fname = asset_names[i];
        shared_ptr<material> mat = asset_mats[i];

        std::cout << "Assets got made" << std::endl;
        try {
            stl_reader :: StlMesh <float, unsigned int> mesh (fname);
            for(size_t itri = 0; itri < mesh.num_tris(); ++itri) {
                const float* c0 = mesh.tri_corner_coords (itri, 0);
                v0 = vec3(c0[0], c0[1], c0[2]);
                const float* c1 = mesh.tri_corner_coords (itri, 1);
                v1 = vec3(c1[0], c1[1], c1[2]);
                const float* c2 = mesh.tri_corner_coords (itri, 2);
                v2 = vec3(c2[0], c2[1], c2[2]);
                // std::cout << "Asset points (x, y, z) :" << v0 << ", " << v1 << ", " << v2 << std::endl;

                // Size Checker
                if (c0[0] > maxAssetsX) maxAssetsX = c0[0];
                if (c0[1] > maxAssetsY) maxAssetsY = c0[1];
                if (c0[2] > maxAssetsZ) maxAssetsZ = c0[2];
                if (c0[0] < minAssetsX) minAssetsX = c0[0];
                if (c0[1] < minAssetsY) minAssetsY = c0[1];
                if (c0[2] < minAssetsZ) minAssetsZ = c0[2];
                if (c1[0] > maxAssetsX) maxAssetsX = c1[0];
                if (c1[1] > maxAssetsY) maxAssetsY = c1[1];
                if (c1[2] > maxAssetsZ) maxAssetsZ = c1[2];
                if (c1[0] < minAssetsX) minAssetsX = c1[0];
                if (c1[1] < minAssetsY) minAssetsY = c1[1];
                if (c1[2] < minAssetsZ) minAssetsZ = c1[2];
                if (c2[0] > maxAssetsX) maxAssetsX = c2[0];
                if (c2[1] > maxAssetsY) maxAssetsY = c2[1];
                if (c2[2] > maxAssetsZ) maxAssetsZ = c2[2];
                if (c2[0] < minAssetsX) minAssetsX = c2[0];
                if (c2[1] < minAssetsY) minAssetsY = c2[1];
                if (c2[2] < minAssetsZ) minAssetsZ = c2[2];

                if (TrustNormal) {
                    const float* n = mesh.tri_normal (itri);
                    vec3 norm = vec3(n[0], n[1], n[2]);
                    world.add(make_shared<triangle>(v0, v1, v2, norm, mat));
                }
                else {
                    world.add(make_shared<triangle>(v0, v1, v2, mat));
                    // std :: cout << "A triangle got made/added" << std :: endl;
                }
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model(asset) or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    // Will this work ?
    // for (int i = 0; i < asset_names.size(); i++) {
    //     const auto& fname = asset_names[i];
    //     shared_ptr<material> mat = asset_mats[i];

    //     auto mesh = obj_loader(fname, asset_trans[i]);
    //     for (int j = 0; j < mesh -> nTriangles; j++) {
    //         shared_ptr<hittable> tri = make_shared<Triangle>(mesh, j, mat);
    //         world.add(tri);
    //     }
    // }

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    std :: cout << "* Minimum and maximum extent of the fluid model" << std :: endl;
    std :: cout << "min point : " << point3(minX, minY, minZ) << std :: endl;
    std :: cout << "max point : " << point3(maxX, maxY, maxZ) << std :: endl;

    std :: cout << "* Minimum and maximum extent of the asset models" << std :: endl;
    std :: cout << "min point : " << point3(minAssetsX, minAssetsY, minAssetsZ) << std :: endl;
    std :: cout << "max point : " << point3(maxAssetsX, maxAssetsY, maxAssetsZ) << std :: endl;

    for (const auto& i : lights) {
        world_bvh.add(i);
    }
    for (const auto& i : walls) {
        world_bvh.add(i);
    }

    return world_bvh;
}

hittable_list num_array_world_builder (std :: vector <std :: string> fnames, std :: vector <shared_ptr<material>> mats, std :: vector<unsigned long> NumOfDataPoints,
std :: vector<float> RadOfSphere, std :: vector <shared_ptr<hittable>> lights, std :: vector <shared_ptr<hittable>> walls, unsigned int& seed) {
    
    std :: cerr << std :: endl << "[The importation process]" << std :: endl;
    hittable_list world;

    float minX = 100000, minY = 100000, minZ = 100000;
    float maxX = -100000, maxY = -100000 ,maxZ = -100000;

    bool fortran_order = true;
    unsigned long data_points_dimension = 3;
    
    for (int i = 0; i < fnames.size(); i++) {
        std :: cerr << fnames[i] << std :: endl;
        float particle_rad = RadOfSphere[i];
        auto mat = mats[i];
        unsigned long data_points_number = NumOfDataPoints[i];

        std::vector<unsigned long> shape {data_points_number, data_points_dimension};
        std::vector<float> data_list;
        npy :: LoadArrayFromNumpy(fnames[i], shape, fortran_order, data_list);

        float x_axis, y_axis, z_axis;
        int j, starting_index;

        for (j = 0; j < data_points_number; j++) {
            // Since there are 3 data points for each point in 3D space
            starting_index = 3 * j;
            x_axis = data_list[starting_index];
            y_axis = data_list[starting_index + 1];
            z_axis = data_list[starting_index + 2];
            point3 center(x_axis, y_axis, z_axis);
            
            // Size Checker
            if (x_axis > maxX) maxX = x_axis;
            if (y_axis > maxY) maxY = y_axis;
            if (z_axis > maxZ) maxZ = z_axis;
            if (x_axis < minX) minX = x_axis;
            if (y_axis < minY) minY = y_axis;
            if (z_axis < minZ) minZ = z_axis;

            auto particle_instance = make_shared<sphere>(center, particle_rad, mat);
            world.add(particle_instance);
        }
    }

    std :: cout << "* Minimum and maximum extent of each model" << std :: endl;
    std :: cout << "min point : " << point3(minX, minY, minZ) << std :: endl;
    std :: cout << "max point : " << point3(maxX, maxY, maxZ) << std :: endl;

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    for (const auto& i : lights) {
        world_bvh.add(i);
    }
    for (const auto& i : walls) {
        world_bvh.add(i);
    }

    return world_bvh;
}
hittable_list mix_world_builder (std :: vector <std :: string> sphere_names, std :: vector <shared_ptr<material>> sphere_mats, 
std :: vector<unsigned long> NumOfDataPoints, std :: vector<float> RadOfSphere, 
std :: vector <std :: string> obj_names, std :: vector <shared_ptr<material>> obj_mats, std :: vector <Transform> trans, 
std :: vector <shared_ptr<hittable>> lights, std :: vector <shared_ptr<hittable>> walls, unsigned int& seed) {
    
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;

    float minX = 100000, minY = 100000, minZ = 100000;
    float maxX = -100000, maxY = -100000 ,maxZ = -100000;

    bool fortran_order = true;
    unsigned long data_points_dimension = 3;
    
    for (int i = 0; i < sphere_names.size(); i++) {
        std :: cerr << sphere_names[i] << std :: endl;
        float particle_rad = RadOfSphere[i];
        auto mat = sphere_mats[i];
        unsigned long data_points_number = NumOfDataPoints[i];

        std::vector<unsigned long> shape {data_points_number, data_points_dimension};
        std::vector<float> data_list;
        npy :: LoadArrayFromNumpy(sphere_names[i], shape, fortran_order, data_list);

        float x_axis, y_axis, z_axis;
        int j, starting_index;

        for (j = 0; j < data_points_number; j++) {
            
            // Since there are 3 data points for each point in 3D space
            starting_index = 3 * j;
            x_axis = data_list[starting_index];
            y_axis = data_list[starting_index + 1];
            z_axis = data_list[starting_index + 2];
            point3 center(x_axis, y_axis, z_axis);
            
            // Size Checker
            if (x_axis > maxX) maxX = x_axis;
            if (y_axis > maxY) maxY = y_axis;
            if (z_axis > maxZ) maxZ = z_axis;
            if (x_axis < minX) minX = x_axis;
            if (y_axis < minY) minY = y_axis;
            if (z_axis < minZ) minZ = z_axis;

            auto particle_instance = make_shared<sphere>(center, particle_rad, mat);
            world.add(particle_instance);
        }
    }
    for (int i = 0; i < obj_names.size(); i++) {
        const auto& fname = obj_names[i];
        shared_ptr<material> mat = obj_mats[i];

        auto mesh = obj_loader(fname, trans[i]);
        for (int j = 0; j < mesh -> nTriangles; j++) {
            shared_ptr<hittable> tri = make_shared<Triangle>(mesh, j, mat);
            world.add(tri);
        }
    }

    std :: cout << "* Minimum and maximum extent of each model" << std :: endl;
    std :: cout << "min point : " << point3(minX, minY, minZ) << std :: endl;
    std :: cout << "max point : " << point3(maxX, maxY, maxZ) << std :: endl;

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    for (const auto& i : lights) {
        world_bvh.add(i);
    }
    for (const auto& i : walls) {
        world_bvh.add(i);
    }

    return world_bvh;
}
#endif