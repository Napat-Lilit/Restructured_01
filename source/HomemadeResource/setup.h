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

extern std :: string out_file_name;
extern int starting_index;
extern int ending_index;

extern std :: vector<std :: string> model_names;
extern std :: vector<shared_ptr<material>> model_mats;

extern vec3 vup;
extern point3 lookfrom;
extern point3 lookat;

extern int samples_per_pixel;
extern int sample_max_depth;
extern color BackgroundColorForIllumination;
extern color BackgroundColorForOutput;

// Assets
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
    // if (!(InitializeYesOrNo())) {
    //     std :: cout << "Something went wrong with Yes/No initialization" << std :: endl;
    //     return false;
    // } 
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
    SetUp.open("./SetupConfig/DynamicModelsSetup.txt");

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
                    else if (parameters[0] == "diff") {
                        float r = std :: stof(parameters[1]);
                        float g = std :: stof(parameters[2]);
                        float b = std :: stof(parameters[3]);
                        color temp(r, g, b);
                        auto mat = make_shared<diffuse_light>(temp);
                        model_mats.push_back(mat);
                    }
                    else {
                        std :: cerr << "Find an unrecognized materials in OBJModelmats, please check for spelling." << std :: endl;
                        return false;
                    }
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
            else if (words[0] == "BackgroundColorForIllumination:") {
                std :: istringstream iss(words[1]);
                std :: vector<std::string> parameters;
                while (std :: getline(iss, myword, ','))
                    parameters.push_back(myword);
                BackgroundColorForIllumination = color(std :: stof(parameters[0]), std :: stof(parameters[1]), std :: stof(parameters[2]));
                std :: cerr << "Ambient light set to : " << BackgroundColorForIllumination << std :: endl;
            }
            else if (words[0] == "BackgroundColorForOutput:") {
                std :: istringstream iss(words[1]);
                std :: vector<std::string> parameters;
                while (std :: getline(iss, myword, ','))
                    parameters.push_back(myword);
                BackgroundColorForOutput = color(std :: stof(parameters[0]), std :: stof(parameters[1]), std :: stof(parameters[2]));
                std :: cerr << "Output background set to : " << BackgroundColorForOutput << std :: endl;
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
        std :: cout << "Couldn't open LightsSetup.txt\n";
        return false;
    }
    SetUp.close();
    return true;
}
bool InitializeAssets() {
    std :: ifstream SetUp;
    SetUp.open("./SetupConfig/StaticModelsSetup.txt");

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
                    else if (parameters[0] == "diff") {
                        float r = std :: stof(parameters[1]);
                        float g = std :: stof(parameters[2]);
                        float b = std :: stof(parameters[3]);
                        color temp(r, g, b);
                        auto mat = make_shared<diffuse_light>(temp);
                        asset_mats.push_back(mat);
                    }
                    else {
                        std :: cerr << "Find an unrecognized materials in OBJModelmats(Assets), please check for spelling." << std :: endl;
                        return false;
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
