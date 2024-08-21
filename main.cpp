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
#include <cmath>

/*
    Undergoing refactoring + restructuring, many functionalities are unusable at the moment.
    1. Light source importance sampling
    2. OBJ models (most notable, trust_normal is removed -> Not sure if we should reinstated in the future or not...)
    2. Sphere or any other magic models
    3. Dynamic ambient light
    4. Denoising is removed, will make it into a separate sub-project
*/

std :: string out_file_name;
int starting_index;
int ending_index;

int parallel_start;
int parallel_end;

std :: vector<std :: string> model_names;
std :: vector<int> model_numbers;
std :: vector<shared_ptr<material>> model_mats;
std :: vector<shared_ptr<hittable>> lights;

// Camera. Be careful about the relation between vup and our viewing direction.
vec3 vup;
point3 lookfrom;
point3 lookat;

int samples_per_pixel;
int sample_max_depth;

// Background Color
color BackgroundColorForIllumination;
color BackgroundColorForOutput;

// Assets
std :: vector<std :: string> asset_names;
std :: vector<shared_ptr<material>> asset_mats;

int main(int argc, char **argv) {
    // First argv (excluding the name of the program) should be num_threads 
    int num_threads = atoi(argv[1]);
    omp_set_num_threads(num_threads);

    // Image Setup
    const auto aspect_ratio = 1.7777;
    const int image_width = 1920;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int denoising_sample = 10;

    int padding_amount = 3;

    // Timing
    auto start = std :: chrono :: high_resolution_clock::now();
    Initialization();
    // auto ray_color = (naive)? naive_tracer : explicit_tracer;
    auto ray_color = naive_tracer; // On the process of refactoring, temporary forcing users to use naive rendering for simplicity

    auto vfov = 40.0;
    auto aperture = 0.0;
    auto distance_to_focus = 2.0;
    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, distance_to_focus);

    // World Setter
    unsigned init_seed;

    for (int i = starting_index; i <= ending_index ; i++) {
        // Denoise Test
        std :: ofstream myfile_normal;
        std :: stringstream ss_normal;
        ss_normal << "./results/" << out_file_name << "_normalMap_" << std :: setw(padding_amount) << std :: setfill('0') << i << ".ppm";   // Output image file name
        std :: string fname_normal = ss_normal.str();
        myfile_normal.open(fname_normal);
        std :: ofstream myfile_position;
        std :: stringstream ss_position;
        ss_position << "./results/" << out_file_name << "_positionMap_" << std :: setw(padding_amount) << std :: setfill('0') << i << ".ppm";   // Output image file name
        std :: string fname_position = ss_position.str();
        myfile_position.open(fname_position);
        std :: ofstream myfile_standardized;
        std :: stringstream ss_standardized;
        ss_standardized << "./results/" << out_file_name << "_raytracedMap_" << std :: setw(padding_amount) << std :: setfill('0') << i << ".ppm";   // Output image file name
        std :: string fname_standardized = ss_standardized.str();
        myfile_standardized.open(fname_standardized);

        std :: vector<std :: string> current_model_names;
        for (int j = 0; j < model_names.size(); j++) {
            std :: stringstream ss02;
            ss02 << model_names[j] << std :: setw(padding_amount) << std :: setfill('0') << i << ".stl";
            current_model_names.push_back(ss02.str());
        }
        hittable_list world;
        // These will be used to make positional mapping, which is necessary for denoising
        vec3 worldMin;
        vec3 worldMax;
        world = stl_world_builder(current_model_names, model_mats, asset_names, asset_mats, init_seed, worldMin, worldMax);
        
        auto end = std :: chrono :: high_resolution_clock :: now();
        auto total_duration = std :: chrono :: duration_cast<std :: chrono :: seconds>(end - start);
        std :: cerr << "\nFinished reading up the " << i << "-th file in : " << total_duration.count() << " seconds" << std :: endl;
        std :: cerr << "There are " << world.objects.size() << " models in world" 
        << std :: endl << std :: endl;
        std :: cerr << "The importation and world buiding completed with no issues, proceed to the main rendering process" << std :: endl << std :: endl;
        std :: cerr << "---------------------------------------------------------------------------------" << std :: endl;
        std :: cerr << std :: endl << "[The rendering process log]" << std :: endl << "* Ideally, there should be nothing here" << std :: endl;

        // Denoise Test
        myfile_normal << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        myfile_position << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        myfile_standardized << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        std::vector<std::vector<vec3>> normal_map(image_height, std::vector<vec3>(image_width));
        std::vector<std::vector<vec3>> position_map(image_height, std::vector<vec3>(image_width));
        std::vector<std::vector<vec3>> standardized_picture_map(image_height, std::vector<vec3>(image_width));

        #pragma omp parallel
        {
            unsigned seed = omp_get_thread_num();

            #pragma omp for schedule(dynamic, 8)
            for (int j = image_height - 1; j >= 0; j--){
                color tmp_normal[image_width];
                color tmp_position[image_width];
                color tmp_standardized[image_width];

                for (int i = 0; i < image_width; i++){
                    color pixel_color(0, 0, 0);
                    for (int s = 0; s < samples_per_pixel; s++) {

                        auto u = (i + random_float(seed)) / (image_width - 1);
                        auto v = (j + random_float(seed)) / (image_height - 1);
                        ray r = cam.get_ray(u, v, seed);

                        // true here refers to the first shot of ray tracing -> used for determining which background color to use
                        pixel_color += ray_color(r, BackgroundColorForIllumination, BackgroundColorForOutput, world, sample_max_depth, true, seed);
                    }
                    tmp_standardized[i] = pixel_color/samples_per_pixel;
                }
                // Denoiser Test
                for (int i = 0; i < image_width; i++){
                    color pixel_color_normal(0, 0, 0);
                    color pixel_color_position(0, 0, 0);
                    for (int s = 0; s < denoising_sample; s++) {

                        auto u = (i + random_float(seed)) / (image_width - 1);
                        auto v = (j + random_float(seed)) / (image_height - 1);
                        ray r = cam.get_ray(u, v, seed);
                        
                        pixel_color_normal += normal_color(r, world);
                        pixel_color_position += position_color(r, world, worldMin, worldMax, vec3());
                    }
                    tmp_normal[i] = pixel_color_normal/(static_cast<float>(denoising_sample));
                    tmp_position[i] = pixel_color_position/(static_cast<float>(denoising_sample));
                }

                std :: copy(std :: begin(tmp_normal), std :: end(tmp_normal), std :: begin(normal_map[j]));
                std :: copy(std :: begin(tmp_position), std :: end(tmp_position), std :: begin(position_map[j]));
                std :: copy(std :: begin(tmp_standardized), std :: end(tmp_standardized), std :: begin(standardized_picture_map[j]));
            }
        }
        end = std :: chrono :: high_resolution_clock :: now();
        total_duration = std :: chrono :: duration_cast<std :: chrono :: seconds>(end - start);
        std :: cerr << std :: endl << "Finish the rendering loop in : " << total_duration.count() << " seconds ----> Proceed to the image saving loop" << std :: endl;

        for (int j = image_height - 1; j>=0; j--){
            for (int i = 0; i < image_width; i++){
                write_color(myfile_normal, normal_map[j][i], 1);
                write_color(myfile_position, position_map[j][i], 1);
                write_color(myfile_standardized, standardized_picture_map[j][i], 1);
            }
        }
        std :: cerr << "Normal image written to " << fname_normal << std :: endl;
        std :: cerr << "Position image written to " << fname_position << std :: endl;
        std :: cerr << "Standard image written to " << fname_standardized << std :: endl;
        myfile_normal.close();
        myfile_position.close();
        myfile_standardized.close();

        end = std :: chrono :: high_resolution_clock :: now();
        total_duration = std :: chrono :: duration_cast<std :: chrono :: seconds>(end - start);
        std :: cerr << std :: endl << "Finish saving results in : " << total_duration.count() << " seconds" << std :: endl << std :: endl;
        std :: cerr << "--------------------------------------------------------------------------------- Finish processing file no. " << i << std :: endl;
    }
}