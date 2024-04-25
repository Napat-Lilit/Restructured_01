#ifndef SETUP_H
#define SETUP_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "material.h"
#include "aarect.h"
#include "transform.h"

extern std :: string CombinationType;

extern std :: string out_file_name;
extern int starting_index;
extern int ending_index;

extern int parallel_start;
extern int parallel_end;

extern std :: vector<std :: string> model_names;
extern std :: vector<std :: string> sphere_names;
extern std :: vector<shared_ptr<material>> model_mats;
extern std :: vector<shared_ptr<material>> sphere_mats;
extern std :: vector <shared_ptr<hittable>> lights;
extern std :: vector <shared_ptr<hittable>> walls;
extern std :: vector <unsigned long> sphere_num;
extern std :: vector <float> sphere_rad;

extern bool TrustNormal;
extern bool naive;
extern bool HasLights;
extern bool DynamicLight;

extern vec3 vup;
extern point3 lookfrom;
extern point3 lookat;

extern int samples_per_pixel;
extern int sample_max_depth;

extern color BackgroundColor;

// Try to put assets in
extern bool HasAssets;
extern std :: vector<std :: string> asset_names;
extern std :: vector<shared_ptr<material>> asset_mats;
extern std :: vector <Transform> asset_trans;

bool InitializeModels();
bool InitializeYesOrNo();
bool InitializeLookFrom();
bool InitializeLights();
bool InitializeIO();
bool InitializeAssets();

bool Initialization() {
    std :: cerr << std :: endl << "[The initialization process]" << std :: endl;
    if (!(InitializeModels())) {
        std :: cout << "Something went wrong with models initialization" << std :: endl;
        return false;
    } 
    if (!(InitializeYesOrNo())) {
        std :: cout << "Something went wrong with Yes/No initialization" << std :: endl;
        return false;
    } 
    if (!(InitializeLookFrom())) {
        std :: cout << "Something went wrong with look initialization" << std :: endl;
        return false;
    } 
    if (!(InitializeLights())) {
        std :: cout << "Something went wrong with lights initialization" << std :: endl;
        return false;
    }
    if (!(InitializeIO())) {
        std :: cout << "Something went wrong with IO initialization" << std :: endl;
        return false;
    }
    if (!(InitializeAssets())) {
        std :: cout << "Something went wrong with Assets initialization" << std :: endl;
        return false;
    }

    std :: cout << std :: endl << "The initialization completed with no issues, proceed to the importation of models" << std :: endl << std :: endl;
    std :: cout << "---------------------------------------------------------------------------------" << std :: endl << std :: endl;
    return true;
}

bool InitializeIO() {
    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/IOSetup.txt");

    std :: string myline;
    std :: string myword;
    if (SetUp.is_open()) {
        while (SetUp) {
            std :: getline (SetUp, myline);
            std :: istringstream iss(myline);

            std :: vector<std::string> words;
            while (std :: getline(iss, myword, ' '))
                words.push_back(myword);

            if (words.size() == 0) continue;
            else if (words[0] == "OutputName:") {
                out_file_name = words[1];
            }
            else if (words[0] == "StartingIndex:") {
                starting_index = std :: stoi(words[1]);
                // std :: cout << "Staring index in setup = " << starting_index << std :: endl;
            }
            else if (words[0] == "EndingIndex:") {
                ending_index = std :: stoi(words[1]);
            }
            else if (words[0] == "ParallelStartingIndex:") {
                parallel_start = std :: stoi(words[1]);
            }
            else if (words[0] == "ParallelEndingIndex:") {
                parallel_end = std :: stoi(words[1]);
            }
            else if (words[0] == "#") {
                continue;
            }
            else {
                return false;
            }
        }
    }
    else {
        std :: cout << "Couldn't open IOSetup.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}

bool InitializeModels() {
    std :: vector<shared_ptr<material>> WallMats;

    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/ModelsSetup.txt");

    // Used to keep track of where the wall geo should be mapped to wall color
    int WallCounters = 0;

    std :: string myline;
    std :: string myword;
    if (SetUp.is_open()) {
        while (SetUp) {
            std :: getline (SetUp, myline);
            std :: istringstream iss(myline);

            std :: vector<std::string> words;
            while (std :: getline(iss, myword, ' '))
                words.push_back(myword);

            // Empty line, just continue to the next line
            if (words.size() == 0) continue;
            else if (words[0] == "CombinationType:") {
                CombinationType = words[1];
            }
            else if (words[0] == "OBJModelName:") {
                for (int i = 1; i < words.size(); i++) {
                    model_names.push_back(words[i]);
                    std :: cerr << "OBJModelNames : " << words[i] << std :: endl;
                }
            }
            else if (words[0] == "OBJModelMats:") {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "lam") {
                        float r = std :: stof(parameters[1]);
                        float g = std :: stof(parameters[2]);
                        float b = std :: stof(parameters[3]);
                        color temp(r, g, b);
                        auto mat = make_shared<lambertian>(temp);
                        model_mats.push_back(mat);
                    }
                    else if (parameters[0] == "dielectric") {
                        float index_refrac = std :: stof(parameters[1]);
                        auto mat = make_shared<dielectric>(index_refrac);
                        model_mats.push_back(mat);
                    }
                    else if (parameters[0] == "beer") {
                        float index_refrac = std :: stof(parameters[1]);
                        float density = std :: stof(parameters[2]);
                        float r = std :: stof(parameters[3]);
                        float g = std :: stof(parameters[4]);
                        float b = std :: stof(parameters[5]);
                        color temp(r, g, b);
                        auto mat = make_shared<beer_lambert_dielectric>(index_refrac, density, temp);
                        model_mats.push_back(mat);
                    }
                    else if (parameters[0] == "metal") {
                        float fuzz = std :: stof(parameters[1]);
                        float r = std :: stof(parameters[2]);
                        float g = std :: stof(parameters[3]);
                        float b = std :: stof(parameters[4]);
                        color temp(r, g, b);
                        auto mat = make_shared<metal>(temp, fuzz);
                        model_mats.push_back(mat);
                    }
                    else {
                        std :: cerr << "Find an unrecognized materials in OBJModelmats, please check for spelling." << std :: endl;
                        return false;
                    }
                }
            }
            else if (words[0] == "SphereModelName:") {
                for (int i = 1; i < words.size(); i++) {
                    sphere_names.push_back(words[i]);
                    std :: cerr << "SphereModelNames : " << words[i] << std :: endl;
                }
            }
            else if (words[0] == "SphereModelMats:") {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "lam") {
                        float r = std :: stof(parameters[1]);
                        float g = std :: stof(parameters[2]);
                        float b = std :: stof(parameters[3]);
                        color temp(r, g, b);
                        auto mat = make_shared<lambertian>(temp);
                        sphere_mats.push_back(mat);
                    }
                    else if (parameters[0] == "beer") {
                        float index_refrac = std :: stof(parameters[1]);
                        float density = std :: stof(parameters[2]);
                        float r = std :: stof(parameters[3]);
                        float g = std :: stof(parameters[4]);
                        float b = std :: stof(parameters[5]);
                        color temp(r, g, b);
                        auto mat = make_shared<beer_lambert_dielectric>(index_refrac, density, temp);
                        sphere_mats.push_back(mat);
                    }
                    else if (parameters[0] == "metal") {
                        float fuzz = std :: stof(parameters[1]);
                        float r = std :: stof(parameters[2]);
                        float g = std :: stof(parameters[3]);
                        float b = std :: stof(parameters[4]);
                        color temp(r, g, b);
                        auto mat = make_shared<metal>(temp, fuzz);
                        sphere_mats.push_back(mat);
                    }
                    else {
                        std :: cerr << "Find an unrecognized materials in SphereModelMats, please check for spelling." << std :: endl;
                        return false;
                    }
                }
            }
            else if (words[0] == "SphereModelNum:") {
                for (int i = 1; i < words.size(); i++) {
                    sphere_num.push_back(std :: stoi(words[i]));
                }
            }
            else if (words[0] == "SphereModelRad:") {
                for (int i = 1; i < words.size(); i++) {
                    sphere_rad.push_back(std :: stof(words[i]));
                }
            }
            else if (words[0] == "WallMats:") {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "lam") {
                        float r = std :: stof(parameters[1]);
                        float g = std :: stof(parameters[2]);
                        float b = std :: stof(parameters[3]);
                        color temp(r, g, b);
                        auto mat = make_shared<lambertian>(temp);
                        WallMats.push_back(mat);
                    }
                    else if (parameters[0] == "beer") {
                        float index_refrac = std :: stof(parameters[1]);
                        float density = std :: stof(parameters[2]);
                        float r = std :: stof(parameters[3]);
                        float g = std :: stof(parameters[4]);
                        float b = std :: stof(parameters[5]);
                        color temp(r, g, b);
                        auto mat = make_shared<beer_lambert_dielectric>(index_refrac, density, temp);
                        WallMats.push_back(mat);
                    }
                    else if (parameters[0] == "metal") {
                        float fuzz = std :: stof(parameters[1]);
                        float r = std :: stof(parameters[2]);
                        float g = std :: stof(parameters[3]);
                        float b = std :: stof(parameters[4]);
                        color temp(r, g, b);
                        auto mat = make_shared<metal>(temp, fuzz);
                        WallMats.push_back(mat);
                    }
                    else {
                        std :: cerr << "Find an unrecognized materials, please check for spelling." << std :: endl;
                        return false;
                    }
                }
            }
            else if (words[0] == "WallGeo:") {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    // std :: cerr << words[i] << std :: endl;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "xy" || parameters[0] == "yx") {
                        std :: cerr << "xy wall got made" << std :: endl;
                        walls.push_back(make_shared<xy_rect>(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]), std :: stof(parameters[4]), std :: stof(parameters[5]), WallMats[WallCounters]));
                    }
                    else if (parameters[0] == "xz" || parameters[0] == "zx") {
                        std :: cerr << "xz wall got made" << std :: endl;
                        walls.push_back(make_shared<xz_rect>(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]), std :: stof(parameters[4]), std :: stof(parameters[5]), WallMats[WallCounters]));
                    }
                    else if (parameters[0] == "yz" || parameters[0] == "zy") {
                        std :: cerr << "yz wall got made" << std :: endl;
                        walls.push_back(make_shared<yz_rect>(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]), std :: stof(parameters[4]), std :: stof(parameters[5]), WallMats[WallCounters]));
                    }
                    else {
                        return false;
                    }
                    WallCounters++;
                }
            }
            else if (words[0] == "#") {
                continue;
            }
            else {
                std :: cerr << "Find an unrecognized leading word, please check for spelling." << std :: endl;
                return false;
            }
        }
    }
    else {
        std :: cout << "Couldn't open ModelsSetup.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}

bool InitializeLookFrom() {
    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/LookFrom.txt");

    std :: string myline;
    std :: string myword;
    if (SetUp.is_open()) {
        while (SetUp) {
            std :: getline (SetUp, myline);
            std :: istringstream iss(myline);

            std :: vector<std::string> words;
            while (std :: getline(iss, myword, ' '))
                words.push_back(myword);

            if (words.size() == 0) continue;
            else if (words[0] == "vup:") {
                vup = vec3(std :: stof(words[1]), std :: stof(words[2]), std :: stof(words[3]));
            }
            else if (words[0] == "lookfrom:") {
                lookfrom = vec3(std :: stof(words[1]), std :: stof(words[2]), std :: stof(words[3]));
            }
            else if (words[0] == "lookat:") {
                lookat = vec3(std :: stof(words[1]), std :: stof(words[2]), std :: stof(words[3]));
            }
            else if (words[0] == "samples_per_pixel:") {
                samples_per_pixel = std :: stoi(words[1]);
            }
            else if (words[0] == "sample_max_depth:") {
                sample_max_depth = std :: stoi(words[1]);
            }
            else if (words[0] == "#") {
                continue;
            }
            else {
                return false;
            }
        }
    }
    else {
        std :: cout << "Couldn't open LookFrom.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}
bool InitializeLights() {
    std :: vector<shared_ptr<material>> LightsMat;

    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/LightsSetup.txt");

    std :: string myline;
    std :: string myword;
    if (SetUp.is_open()) {
        while (SetUp) {
            std :: getline (SetUp, myline);
            std :: istringstream iss(myline);

            std :: vector<std::string> words;
            while (std :: getline(iss, myword, ' '))
                words.push_back(myword);
            // std :: cerr << myline << std :: endl;

            if (words.size() == 0) continue;
            else if (words[0] == "LightGeo:" && HasLights) {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    // std :: cerr << words[i] << std :: endl;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "xy") {
                        std :: cerr << "xy light got made" << std :: endl;
                        lights.push_back(make_shared<xy_rect>(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]), std :: stof(parameters[4]), std :: stof(parameters[5]), LightsMat[i-1]));
                    }
                    else if (parameters[0] == "xz") {
                        std :: cerr << "xz light got made" << std :: endl;
                        lights.push_back(make_shared<xz_rect>(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]), std :: stof(parameters[4]), std :: stof(parameters[5]), LightsMat[i-1]));
                    }
                    else if (parameters[0] == "yz") {
                        std :: cerr << "yz light got made" << std :: endl;
                        lights.push_back(make_shared<yz_rect>(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]), std :: stof(parameters[4]), std :: stof(parameters[5]), LightsMat[i-1]));
                    }
                    else {
                        return false;
                    }
                }
            }
            else if (words[0] == "LightMats:" && HasLights) {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "diff") {
                        LightsMat.push_back(make_shared<diffuse_light>(color(std :: stof(parameters[1]), std :: stof(parameters[2]), 
                        std :: stof(parameters[3]))));
                    }
                    else {
                        return false;
                    }
                }
            }
            else if (words[0] == "BackgroundColor:") {
                std :: istringstream iss(words[1]);
                std :: vector<std::string> parameters;
                while (std :: getline(iss, myword, ','))
                    parameters.push_back(myword);
                BackgroundColor = color(std :: stof(parameters[0]), std :: stof(parameters[1]), std :: stof(parameters[2]));
                std :: cerr << "Ambient light set to : " << BackgroundColor << std :: endl;
            }
            else if (words[0] == "#") {
                continue;
            }
            else {
                if (words[0] == "LightGeo:" || words[0] == "LightMats:"){
                    continue;
                }
                std :: cerr << "Found something else : " << words[0] << std :: endl;
                return false;
            }
        }
    }
    else {
        std :: cout << "Couldn't open LightsSetup.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}
bool InitializeYesOrNo() {
    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/YesOrNo.txt");

    std :: string myline;
    std :: string myword;
    if (SetUp.is_open()) {
        while (SetUp) {
            std :: getline (SetUp, myline);
            std :: istringstream iss(myline);

            std :: vector<std::string> words;
            while (std :: getline(iss, myword, ' '))
                words.push_back(myword);

            if (words.size() == 0) continue;
            else if (words[0] == "TrustNormal:") {
                if (words[1] == "true" || words[1] == "True") TrustNormal = true;
                else if (words[1] == "false" || words[1] == "False") TrustNormal = false;
                else return false;
            }
            else if (words[0] == "Naive:") {
                if (words[1] == "true" || words[1] == "True") naive = true;
                else if (words[1] == "false" || words[1] == "False") naive = false;
                else return false;
            }
            else if (words[0] == "HasLights:") {
                if (words[1] == "true" || words[1] == "True") HasLights= true;
                else if (words[1] == "false" || words[1] == "False") HasLights = false;
                else return false;
            }
            else if (words[0] == "DynamicLight:") {
                if (words[1] == "true" || words[1] == "True") DynamicLight = true;
                else if (words[1] == "false" || words[1] == "False") DynamicLight = false;
                else return false;
            }
            else if (words[0] == "#") {
                continue;
            }
            else {
                return false;
            }
        }
    }
    else {
        std :: cout << "Couldn't open YesOrNo.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}
bool InitializeAssets() {
    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/Assets.txt");

    std :: string myline;
    std :: string myword;
    if (SetUp.is_open()) {
        while (SetUp) {
            std :: getline (SetUp, myline);
            std :: istringstream iss(myline);

            std :: vector<std::string> words;
            while (std :: getline(iss, myword, ' '))
                words.push_back(myword);

            // Empty line, just continue to the next line
            if (words.size() == 0) continue;
            else if (words[0] == "HasAssets:") {
                if (words[1] == "false") HasAssets=false;
                else if (words[1] == "true") HasAssets=true;
                else {
                    std :: cerr << "Find an unrecognized argument in HasAssets, please check for spelling." << std :: endl;
                    return false;
                }
            }
            else if (words[0] == "OBJModelName:") {
                for (int i = 1; i < words.size(); i++) {
                    asset_names.push_back(words[i]);
                    std :: cerr << "OBJModelNames(Assets) : " << words[i] << std :: endl;
                }
            }
            else if (words[0] == "OBJModelMats:") {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);
                    if (parameters[0] == "lam") {
                        float r = std :: stof(parameters[1]);
                        float g = std :: stof(parameters[2]);
                        float b = std :: stof(parameters[3]);
                        color temp(r, g, b);
                        auto mat = make_shared<lambertian>(temp);
                        asset_mats.push_back(mat);
                    }
                    else if (parameters[0] == "dielectric") {
                        float index_refrac = std :: stof(parameters[1]);
                        auto mat = make_shared<dielectric>(index_refrac);
                        asset_mats.push_back(mat);
                    }
                    else if (parameters[0] == "beer") {
                        float index_refrac = std :: stof(parameters[1]);
                        float density = std :: stof(parameters[2]);
                        float r = std :: stof(parameters[3]);
                        float g = std :: stof(parameters[4]);
                        float b = std :: stof(parameters[5]);
                        color temp(r, g, b);
                        auto mat = make_shared<beer_lambert_dielectric>(index_refrac, density, temp);
                        asset_mats.push_back(mat);
                    }
                    else if (parameters[0] == "metal") {
                        float fuzz = std :: stof(parameters[1]);
                        float r = std :: stof(parameters[2]);
                        float g = std :: stof(parameters[3]);
                        float b = std :: stof(parameters[4]);
                        color temp(r, g, b);
                        auto mat = make_shared<metal>(temp, fuzz);
                        asset_mats.push_back(mat);
                    }
                    else {
                        std :: cerr << "Find an unrecognized materials in OBJModelmats(Assets), please check for spelling." << std :: endl;
                        return false;
                    }
                }
            }
            else if (words[0] == "OBJModelTran:") {
                for (int i = 1; i < words.size(); i++) {
                    std :: istringstream iss(words[i]);
                    std :: vector<std::string> parameters;
                    while (std :: getline(iss, myword, ','))
                        parameters.push_back(myword);

                    if (parameters[0] == "1") {
                        if (asset_trans.size() < 1) {
                            Transform Temp = Scale(1.0, 1.0, 1.0);
                            asset_trans.push_back(Temp);
                        }
                        Transform Temp = asset_trans[0];
                        Transform NewTran;

                        if (parameters[1] == "scale") {
                            float x = std :: stof(parameters[2]);
                            float y = std :: stof(parameters[3]);
                            float z = std :: stof(parameters[4]);
                            NewTran = Scale(x, y, z);
                        }
                        else if (parameters[1] == "translate") {
                            float x = std :: stof(parameters[2]);
                            float y = std :: stof(parameters[3]);
                            float z = std :: stof(parameters[4]);
                            NewTran = Translate(vec3(x, y, z));
                        }

                        Transform TotalTran = NewTran * Temp;
                        asset_trans[0] = TotalTran;
                    }
                    else if (parameters[0] == "2") {
                        if (asset_trans.size() < 2) {
                            Transform Temp = Scale(1.0, 1.0, 1.0);
                            asset_trans.push_back(Temp);
                        }
                        Transform Temp = asset_trans[1];
                        Transform NewTran;

                        if (parameters[1] == "scale") {
                            float x = std :: stof(parameters[2]);
                            float y = std :: stof(parameters[3]);
                            float z = std :: stof(parameters[4]);
                            NewTran = Scale(x, y, z);
                        }
                        else if (parameters[1] == "translate") {
                            float x = std :: stof(parameters[2]);
                            float y = std :: stof(parameters[3]);
                            float z = std :: stof(parameters[4]);
                            NewTran = Translate(vec3(x, y, z));
                        }

                        Transform TotalTran = NewTran * Temp;
                        asset_trans[1] = TotalTran;
                    }
                }
            }
            else if (words[0] == "#") {
                continue;
            }
            else {
                return false;
            }
        }
    }
    else {
        std :: cout << "Couldn't open Assets.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}

#endif