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

#include "setup.h"

// Importers. Header only libraries taken from Github
#include <iostream>
#include <fstream>
#include <string>
#include "../ImportedResource/npy.hpp"
#include "../ImportedResource/stl_reader_sebastian.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "../ImportedResource/tiny_obj_loader.h"

hittable_list stl_world_builder (std :: vector <std :: string> fnames, std :: vector <shared_ptr<material>> mats, 
std :: vector <std :: string> asset_names, std :: vector <shared_ptr<material>> asset_mats, unsigned int& seed, vec3& worldMin, vec3& worldMax) {
    
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;
    vec3 v0, v1, v2;

    float minX = 100000, minY = 100000, minZ = 100000;
    float maxX = -100000, maxY = -100000 ,maxZ = -100000;

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

                world.add(make_shared<triangle>(v0, v1, v2, mat));
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    for (int i = 0; i < asset_names.size(); i++) {
        const auto& fname = asset_names[i];
        shared_ptr<material> mat = asset_mats[i];

        std::cout << "Assets being made :" << fname << std::endl;
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

                world.add(make_shared<triangle>(v0, v1, v2, mat));
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model(asset) or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    std::cout << "Just before entering world bvh building" << std::endl;

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    std :: cout << "* Minimum and maximum extent of the world" << std :: endl;
    std :: cout << "min point : " << point3(minX, minY, minZ) << std :: endl;
    std :: cout << "max point : " << point3(maxX, maxY, maxZ) << std :: endl;

    // Returning the minimum and maximum point of the world back as well
    worldMin = vec3(minX, minY, minZ);
    worldMax = vec3(maxX, maxY, maxZ);

    return world_bvh;
}

hittable_list sphere_test (unsigned int& seed, float radius, vec3& worldMin, vec3& worldMax) {
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;

    world.add(make_shared<sphere>(vec3(60.0f, 0.0f, 0.0f), radius, make_shared<lambertian>(vec3(1.0f, 0.2f, 0.2f))));
    // world.add(make_shared<sphere>(vec3(), radius, make_shared<metal>(vec3(0.5f, 0.5f, 0.5f), 0.0f)));
    world.add(make_shared<sphere>(vec3(), radius, make_shared<metal>(vec3(0.5f, 0.5f, 0.5f), 0.0f)));
    // world.add(make_shared<sphere>(vec3(0.0f, 21.0f, 0.0f), radius, make_shared<metal>(vec3(0.5f, 0.5f, 0.5f), 0.0f)));
    world.add(make_shared<sphere>(vec3(0.0f, 60.0f, 0.0f), radius, make_shared<beer_lambert_dielectric>(1.5f, 1.0, vec3(0.1f, 0.1f, 0.1f))));

    world.add(make_shared<sphere>(vec3(0.0f, 0.0f, -1000.0f), 950, make_shared<lambertian>(vec3(0.5f, 0.5f, 0.5f))));
    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    // Returning the minimum and maximum point of the world back as well
    worldMin = vec3(-radius, -radius, -radius);
    worldMax = vec3(radius, radius, radius);

    return world_bvh;
}

hittable_list sphere_world_builder (std::vector<sphereInputHandler> sphereInputList, 
std :: vector <std :: string> asset_names, std :: vector <shared_ptr<material>> asset_mats, unsigned int& seed, vec3& worldMin, vec3& worldMax) {
    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;

    for (size_t sphereDataNo = 0; sphereDataNo < sphereInputList.size(); sphereDataNo++)
    {
        auto spherePtrList = sphereInputList[sphereDataNo].getCompleteSphereList();
        for (auto& spherePtr : spherePtrList) {
            world.add(spherePtr);
        }

        // Will think about multi-sphere input later
        sphereInputList[sphereDataNo].getLimit(worldMin, worldMax);
    }

    vec3 v0, v1, v2;
    for (int i = 0; i < asset_names.size(); i++) {
        const auto& fname = asset_names[i];
        shared_ptr<material> mat = asset_mats[i];

        std::cout << "Assets being made :" << fname << std::endl;
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
                if (c0[0] > worldMax.x()) worldMax.x(c0[0]);
                if (c0[1] > worldMax.y()) worldMax.y(c0[1]);
                if (c0[2] > worldMax.z()) worldMax.z(c0[2]);
                if (c0[0] < worldMin.x()) worldMin.x(c0[0]);
                if (c0[1] < worldMin.y()) worldMin.y(c0[1]);
                if (c0[2] < worldMin.z()) worldMin.z(c0[2]);
                if (c1[0] > worldMax.x()) worldMax.x(c1[0]);
                if (c1[1] > worldMax.y()) worldMax.y(c1[1]);
                if (c1[2] > worldMax.z()) worldMax.z(c1[2]);
                if (c1[0] < worldMin.x()) worldMin.x(c1[0]);
                if (c1[1] < worldMin.y()) worldMin.y(c1[1]);
                if (c1[2] < worldMin.z()) worldMin.z(c1[2]);
                if (c2[0] > worldMax.x()) worldMax.x(c2[0]);
                if (c2[1] > worldMax.y()) worldMax.y(c2[1]);
                if (c2[2] > worldMax.z()) worldMax.z(c2[2]);
                if (c2[0] < worldMin.x()) worldMin.x(c2[0]);
                if (c2[1] < worldMin.y()) worldMin.y(c2[1]);
                if (c2[2] < worldMin.z()) worldMin.z(c2[2]);

                world.add(make_shared<triangle>(v0, v1, v2, mat));
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model(asset) or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    // Debug
    std::cout << "The minimum point is : " << worldMin << std::endl;
    std::cout << "The maximum point is : " << worldMax << std::endl;

    return world_bvh;
}

hittable_list all_world_builder (std :: vector <std :: string> fnames, std :: vector <shared_ptr<material>> mats,
std::vector<sphereInputHandler> sphereInputList, std :: vector <std :: string> asset_names, std :: vector <shared_ptr<material>> asset_mats, 
unsigned int& seed, vec3& worldMin, vec3& worldMax) {

    std :: cerr << "[The importation process]" << std :: endl;
    hittable_list world;

    for (size_t sphereDataNo = 0; sphereDataNo < sphereInputList.size(); sphereDataNo++)
    {
        auto spherePtrList = sphereInputList[sphereDataNo].getCompleteSphereList();
        for (auto& spherePtr : spherePtrList) {
            world.add(spherePtr);
        }

        // Will think about multi-sphere input later
        sphereInputList[sphereDataNo].getLimit(worldMin, worldMax);
    }

    vec3 v0, v1, v2;
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
                if (c0[0] > worldMax.x()) worldMax.x(c0[0]);
                if (c0[1] > worldMax.y()) worldMax.y(c0[1]);
                if (c0[2] > worldMax.z()) worldMax.z(c0[2]);
                if (c0[0] < worldMin.x()) worldMin.x(c0[0]);
                if (c0[1] < worldMin.y()) worldMin.y(c0[1]);
                if (c0[2] < worldMin.z()) worldMin.z(c0[2]);
                if (c1[0] > worldMax.x()) worldMax.x(c1[0]);
                if (c1[1] > worldMax.y()) worldMax.y(c1[1]);
                if (c1[2] > worldMax.z()) worldMax.z(c1[2]);
                if (c1[0] < worldMin.x()) worldMin.x(c1[0]);
                if (c1[1] < worldMin.y()) worldMin.y(c1[1]);
                if (c1[2] < worldMin.z()) worldMin.z(c1[2]);
                if (c2[0] > worldMax.x()) worldMax.x(c2[0]);
                if (c2[1] > worldMax.y()) worldMax.y(c2[1]);
                if (c2[2] > worldMax.z()) worldMax.z(c2[2]);
                if (c2[0] < worldMin.x()) worldMin.x(c2[0]);
                if (c2[1] < worldMin.y()) worldMin.y(c2[1]);
                if (c2[2] < worldMin.z()) worldMin.z(c2[2]);

                world.add(make_shared<triangle>(v0, v1, v2, mat));
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    for (int i = 0; i < asset_names.size(); i++) {
        const auto& fname = asset_names[i];
        shared_ptr<material> mat = asset_mats[i];

        std::cout << "Assets being made :" << fname << std::endl;
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
                if (c0[0] > worldMax.x()) worldMax.x(c0[0]);
                if (c0[1] > worldMax.y()) worldMax.y(c0[1]);
                if (c0[2] > worldMax.z()) worldMax.z(c0[2]);
                if (c0[0] < worldMin.x()) worldMin.x(c0[0]);
                if (c0[1] < worldMin.y()) worldMin.y(c0[1]);
                if (c0[2] < worldMin.z()) worldMin.z(c0[2]);
                if (c1[0] > worldMax.x()) worldMax.x(c1[0]);
                if (c1[1] > worldMax.y()) worldMax.y(c1[1]);
                if (c1[2] > worldMax.z()) worldMax.z(c1[2]);
                if (c1[0] < worldMin.x()) worldMin.x(c1[0]);
                if (c1[1] < worldMin.y()) worldMin.y(c1[1]);
                if (c1[2] < worldMin.z()) worldMin.z(c1[2]);
                if (c2[0] > worldMax.x()) worldMax.x(c2[0]);
                if (c2[1] > worldMax.y()) worldMax.y(c2[1]);
                if (c2[2] > worldMax.z()) worldMax.z(c2[2]);
                if (c2[0] < worldMin.x()) worldMin.x(c2[0]);
                if (c2[1] < worldMin.y()) worldMin.y(c2[1]);
                if (c2[2] < worldMin.z()) worldMin.z(c2[2]);

                world.add(make_shared<triangle>(v0, v1, v2, mat));
            }
        }
        catch (std :: exception& e) {
            std :: cout << "WARNING! : " << e.what() << std :: endl;
            std :: cout << "This may indicate either a really problematic model(asset) or just a harmless blank model." << std :: endl << std :: endl;
        }
    }

    hittable_list world_bvh;
    world_bvh.add(make_shared<bvh_node>(world, seed));

    // Debug
    std::cout << "The minimum point is : " << worldMin << std::endl;
    std::cout << "The maximum point is : " << worldMax << std::endl;

    return world_bvh;

}

#endif