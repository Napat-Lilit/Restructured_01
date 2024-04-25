#include "source/HomemadeResource/color.h"
#include "source/HomemadeResource/aabb.h"
#include "source/HomemadeResource/aarect.h"
#include "source/HomemadeResource/bvh.h"
#include "source/HomemadeResource/camera.h"
#include "source/HomemadeResource/color.h"
#include "source/HomemadeResource/efloat.h"
#include "source/HomemadeResource/hit_record.h"
#include "source/HomemadeResource/hittable_list.h"
#include "source/HomemadeResource/material.h"
#include "source/HomemadeResource/sphere.h"
#include "source/HomemadeResource/triangle.h"
#include "source/HomemadeResource/importer.h"
#include "source/HomemadeResource/material.h"
#include "source/HomemadeResource/pdf.h"
#include "source/HomemadeResource/tracer.h"
#include "source/HomemadeResource/setup.h"

// String stream
#include <sstream>
#include <fstream>
#include <iomanip>

#include <vector>
#include <string>
#include <omp.h>
#include <iostream>
#include <chrono>

bool TrustNormal;
bool naive;
bool HasLights;

std :: string out_file_name;
int starting_index;
int ending_index;

int parallel_start;
int parallel_end;

std :: vector<std :: string> model_names;
std :: vector<std :: string> sphere_names;
std :: vector<int> model_numbers;
std :: vector<std :: string> model_types;
std :: vector<shared_ptr<material>> model_mats;
std :: vector<shared_ptr<material>> sphere_mats;
std :: vector<shared_ptr<hittable>> lights;
std :: vector<shared_ptr<hittable>> walls;
std :: vector<unsigned long> sphere_num;
std :: vector<float> sphere_rad;

// To determined what type of importor/world builder should be used
std :: string CombinationType;

// Camera. Be careful about the relation between vup and our viewing direction.
vec3 vup;
point3 lookfrom;
point3 lookat;

int samples_per_pixel;
int sample_max_depth;

// Background Color
color BackgroundColor;

// Assets
bool HasAssets;
std :: vector<std :: string> asset_names;
std :: vector<shared_ptr<material>> asset_mats;
// Transformation to the models
std :: vector <Transform> asset_trans;

int main(int argc, char **argv) {
    // First argv (excluding the name of the program) should be num_threads 
    int num_threads = atoi(argv[1]);
    omp_set_num_threads(num_threads);

    // Image Setup
    const auto aspect_ratio = 1.0;
    const int image_width = 1500;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    int padding_amount = 3;

    // Timing
    auto start = std :: chrono :: high_resolution_clock::now();
    Initialization();
    auto ray_color = (naive)? naive_tracer : explicit_tracer;

    auto vfov = 40.0;
    auto aperture = 0.0;
    auto distance_to_focus = 2.0;
    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, distance_to_focus);

    // World Setter
    unsigned init_seed;

    // Transformation to the models
    std :: vector <Transform> trans;
    Transform scaling_tran_01 = Scale(1.0, 1.0, 1.0);
    trans.push_back(scaling_tran_01);

    // std :: cout << "the value of starting_index is " << starting_index << std :: endl;

    for (int i = starting_index; i <= ending_index ; i++) {
        // Additional Setups
        std :: ofstream myfile;
        std :: stringstream ss;
        ss << "./results/" << out_file_name << std :: setw(padding_amount) << std :: setfill('0') << i << ".ppm";   // Output image file name
        std :: string fname = ss.str();
        myfile.open(fname);

        std :: vector<std :: string> current_model_names;
        std :: vector<std :: string> current_sphere_names;
        std :: vector<shared_ptr<hittable>> current_lights;
        for (int j = 0; j < model_names.size(); j++) {
            std :: stringstream ss02;
            if (CombinationType == "STL") {
                ss02 << model_names[j] << std :: setw(padding_amount) << std :: setfill('0') << i << ".stl";
            }
            else {
                ss02 << model_names[j] << std :: setw(padding_amount) << std :: setfill('0') << i << ".obj";
            }
            current_model_names.push_back(ss02.str());
        }
        for (int j = 0; j < sphere_names.size(); j++) {
            std :: stringstream ss02;
            ss02 << sphere_names[j] << std :: setw(padding_amount) << std :: setfill('0') << i << ".npy";
            current_sphere_names.push_back(ss02.str());
        }
        // Only support aarect lights. In case of many light sources, make the x component slightly different between lights
        hittable_list world;
        if (CombinationType == "IsoSurface") 
            world = obj_world_builder(current_model_names, model_mats, trans, lights, walls, init_seed);
        else if (CombinationType == "Spheres") 
            world = num_array_world_builder(current_sphere_names, sphere_mats, sphere_num, sphere_rad, lights, walls, init_seed);
        else if (CombinationType == "STL") {
            world = stl_world_builder(current_model_names, model_mats, trans, lights, walls, init_seed);
            // world = stl_world_builder_with_assets(current_model_names, model_mats, trans, lights, walls, asset_names, asset_mats, asset_trans, init_seed);
        }
        else if (CombinationType == "IsoSurfaceAndSphere") {
            // How to deal with this again
            world = mix_world_builder(current_sphere_names, sphere_mats, sphere_num, sphere_rad, 
            current_model_names, model_mats, trans, lights, walls, init_seed);
        }
        else 
            std :: cerr << "Unknown combination type" << std :: endl;
        
        auto end = std :: chrono :: high_resolution_clock :: now();
        auto total_duration = std :: chrono :: duration_cast<std :: chrono :: seconds>(end - start);
        std :: cerr << "\nFinished reading up the " << i << "-th file in : " << total_duration.count() << " seconds" << std :: endl;
        std :: cerr << "There are " << world.objects.size() << " objects in world (models counted as one + number of lights and other objects such as walls in the world)" 
        << std :: endl << std :: endl;
        std :: cerr << "The importation and world buiding completed with no issues, proceed to the main rendering process" << std :: endl << std :: endl;
        std :: cerr << "---------------------------------------------------------------------------------" << std :: endl;
        std :: cerr << std :: endl << "[The rendering process log]" << std :: endl << "* Ideally, there should be nothing here" << std :: endl;

        myfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        static color picture_map[image_height][image_width];

        #pragma omp parallel
        {
            unsigned seed = omp_get_thread_num();

            #pragma omp for schedule(dynamic, 8)
            for (int j = image_height - 1; j >= 0; j--){
                color tmp[image_width];
                for (int i = 0; i < image_width; i++){
                    color pixel_color(0, 0, 0);
                    for (int s = 0; s < samples_per_pixel; s++) {

                        auto u = (i + random_float(seed)) / (image_width - 1);
                        auto v = (j + random_float(seed)) / (image_height - 1);
                        ray r = cam.get_ray(u, v, seed);

                        pixel_color += ray_color(r, BackgroundColor, world, sample_max_depth, lights, true, seed);
                    }
                    tmp[i] = pixel_color;
                }
                std :: copy(std :: begin(tmp), std :: end(tmp), std :: begin(picture_map[j]));
            }
        }
        end = std :: chrono :: high_resolution_clock :: now();
        total_duration = std :: chrono :: duration_cast<std :: chrono :: seconds>(end - start);
        std :: cerr << std :: endl << "Finish the rendering loop in : " << total_duration.count() << " seconds ----> Proceed to the image saving loop" << std :: endl;

        for (int j = image_height - 1; j>=0; j--){
            for (int i = 0; i < image_width; i++){
                write_color(myfile, picture_map[j][i], samples_per_pixel);
            }
        }

        end = std :: chrono :: high_resolution_clock :: now();
        total_duration = std :: chrono :: duration_cast<std :: chrono :: seconds>(end - start);
        std :: cerr << "Finish the whole process in : " << total_duration.count() 
        << " seconds ----> Image written to " << fname << std :: endl << std :: endl;
        std :: cerr << "--------------------------------------------------------------------------------- Finish processing file no. " << i << std :: endl;

        myfile.close();
    }
}