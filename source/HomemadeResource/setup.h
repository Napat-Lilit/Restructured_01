#ifndef SETUP_H
#define SETUP_H

#include <sstream>
// #include <stdexcept>
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

// Camera Batchs
extern std::vector<vec3> vupList;
extern std::vector<point3> lookfromList;
extern std::vector<point3> lookatList;

std::shared_ptr<material> getMatFromTxt(std::string matString) {

    std :: istringstream iss(matString);
    std :: string myword;
    std :: vector<std::string> parameters;

    while (std :: getline(iss, myword, ','))
        parameters.push_back(myword);
    if (parameters[0] == "lam") {
        float r = std :: stof(parameters[1]);
        float g = std :: stof(parameters[2]);
        float b = std :: stof(parameters[3]);
        color temp(r, g, b);
        auto mat = make_shared<lambertian>(temp);
        return mat;
    }
    else if (parameters[0] == "dielectric") {
        float index_refrac = std :: stof(parameters[1]);
        auto mat = make_shared<dielectric>(index_refrac);
        return mat;
    }
    else if (parameters[0] == "beer") {
        float index_refrac = std :: stof(parameters[1]);
        float density = std :: stof(parameters[2]);
        float r = std :: stof(parameters[3]);
        float g = std :: stof(parameters[4]);
        float b = std :: stof(parameters[5]);
        color temp(r, g, b);
        auto mat = make_shared<beer_lambert_dielectric>(index_refrac, density, temp);
        return mat;
    }
    else if (parameters[0] == "metal") {
        float fuzz = std :: stof(parameters[1]);
        float r = std :: stof(parameters[2]);
        float g = std :: stof(parameters[3]);
        float b = std :: stof(parameters[4]);
        color temp(r, g, b);
        auto mat = make_shared<metal>(temp, fuzz);
        return mat;
    }
    else if (parameters[0] == "diff") {
        float r = std :: stof(parameters[1]);
        float g = std :: stof(parameters[2]);
        float b = std :: stof(parameters[3]);
        color temp(r, g, b);
        auto mat = make_shared<diffuse_light>(temp);
        return mat;
    }

    // std :: cerr << "Find an unrecognized materials in SphereMats configuration, please check for spelling." << std :: endl;
    throw std::runtime_error("Find an unrecognized materials in SphereMats configuration, please check for spelling.");
}
// Sphere handling
class sphereInputHandler {
    public:
        sphereInputHandler(int _type) : type(_type) {};
        ~sphereInputHandler() {};

        void processData() {
            if (type == 1)
            {
                auto mat = getMatFromTxt(ShphereMats[0]);
                std :: ifstream SetUp;
                SetUp.open(SphereNameCurrent);

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

                        vec3 pos(std::stof(words[0]), std::stof(words[1]), std::stof(words[2]));
                        CompleteSphereList.push_back(make_shared<sphere>(pos, SphereRadius[0], mat));

                        // Limit
                        updateLimit(std::stof(words[0])+SphereRadius[0], std::stof(words[1])+SphereRadius[0], std::stof(words[2])+SphereRadius[0]);
                        updateLimit(std::stof(words[0])-SphereRadius[0], std::stof(words[1])-SphereRadius[0], std::stof(words[2])-SphereRadius[0]);

                        // Debug
                        std::cout << "The following sphere got made :" << pos 
                        << " + " << SphereRadius[0] << " + " << ShphereMats[0] << "\n";
                    }
                }
                else {
                    std :: cout << "Couldn't open sphere file :" << SphereNameCurrent << "\n";
                }
                SetUp.close();
            }
            else if (type == 2)
            {
                auto mat = getMatFromTxt(ShphereMats[0]);
                std :: ifstream SetUp;
                SetUp.open(SphereNameCurrent);

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

                        vec3 pos(std::stof(words[0]), std::stof(words[1]), std::stof(words[2]));
                        CompleteSphereList.push_back(make_shared<sphere>(pos, std::stof(words[3]), mat));

                        // Limit
                        updateLimit(std::stof(words[0])+std::stof(words[3]), std::stof(words[1])+std::stof(words[3]), std::stof(words[2])+std::stof(words[3]));
                        updateLimit(std::stof(words[0])-std::stof(words[3]), std::stof(words[1])-std::stof(words[3]), std::stof(words[2])-std::stof(words[3]));
                    
                        // Debug
                        std::cout << "The following sphere got made :" << pos 
                        << " + " << std::stof(words[3]) << " + " << ShphereMats[0] << "\n";
                    }
                }
                else {
                    std :: cout << "Couldn't open sphere file :" << SphereNameCurrent << "\n";
                }
                SetUp.close();
            }
            else if (type == 3)
            {
                std::vector<std::shared_ptr<material>> matList;
                for (size_t typeNo = 0; typeNo < getTypesNumber(); typeNo++)
                {
                    auto mat = getMatFromTxt(ShphereMats[typeNo]);
                    matList.push_back(mat);
                }
                
                std :: ifstream SetUp;
                SetUp.open(SphereNameCurrent);

                std :: string myline;
                std :: string myword;
                if (SetUp.is_open()) {
                    while (SetUp) {
                        std :: getline (SetUp, myline);
                        std :: istringstream iss(myline);

                        std :: vector<std::string> words;
                        while (iss >> myword) {
                            words.push_back(myword);
                        }
                        if (words.size() == 0) continue;

                        vec3 pos(std::stof(words[0]), std::stof(words[1]), std::stof(words[2]));
                        CompleteSphereList.push_back
                        (make_shared<sphere>(pos, SphereRadius[static_cast<int>(std::stof(words[3]))], matList[static_cast<int>(std::stof(words[3]))]));
                        // (make_shared<sphere>(pos, SphereRadius[static_cast<int>(std::stof(words[3])) - 1], matList[static_cast<int>(std::stof(words[3]) - 1.f)]));
                    
                        // Limit
                        updateLimit(std::stof(words[0])+SphereRadius[std::stoi(words[3])], std::stof(words[1])+SphereRadius[std::stoi(words[3])], std::stof(words[2])+SphereRadius[std::stoi(words[3])]);
                        updateLimit(std::stof(words[0])-SphereRadius[std::stoi(words[3])], std::stof(words[1])-SphereRadius[std::stoi(words[3])], std::stof(words[2])-SphereRadius[std::stoi(words[3])]);
                    }
                }
                else {
                    std :: cout << "Couldn't open sphere file :" << SphereNameCurrent << "\n";
                }
                SetUp.close();
            }
            else if (type == 4)
            {
                // Debug
                std::cout << "Just before sphere processing occurs" << std::endl;

                std::vector<std::shared_ptr<material>> matList;
                for (size_t typeNo = 0; typeNo < getTypesNumber(); typeNo++)
                {
                    auto mat = getMatFromTxt(ShphereMats[typeNo]);
                    matList.push_back(mat);
                }
                
                std :: ifstream SetUp;
                SetUp.open(SphereNameCurrent);

                std :: string myline;
                std :: string myword;
                if (SetUp.is_open()) {
                    while (SetUp) {
                        std :: getline (SetUp, myline);
                        std :: istringstream iss(myline);

                        std :: vector<std::string> words;
                        while (iss >> myword) {
                            words.push_back(myword);
                        }
                        if (words.size() == 0) continue;

                        vec3 pos(std::stof(words[0]), std::stof(words[1]), std::stof(words[2]));
                        CompleteSphereList.push_back
                        (make_shared<sphere>(pos, std::stof(words[3]), matList[static_cast<int>(std::stof(words[4]))]));
                        // (make_shared<sphere>(pos, std::stof(words[3]), matList[static_cast<int>(std::stof(words[4])) - 1]));

                        // Limit
                        updateLimit(std::stof(words[0])+std::stof(words[3]), std::stof(words[1])+std::stof(words[3]), std::stof(words[2])+std::stof(words[3]));
                        updateLimit(std::stof(words[0])-std::stof(words[3]), std::stof(words[1])-std::stof(words[3]), std::stof(words[2])-std::stof(words[3]));
                    }
                }
                else {
                    std :: cout << "Couldn't open sphere file :" << SphereNameCurrent << "\n";
                }
                SetUp.close();

                // Debug
                std::cout << "Just after sphere processing occurs" << std::endl;
            }
            else {
                std :: cout << "Something went wrong with processing sphere data" << std :: endl;
            }
        }
        // Debug
        void printInfo() {
            std::cout << "Type is " << type << std::endl;
            std::cout << "Sphere's name is " << SphereName << std::endl;
            for (size_t i = 0; i < SphereRadius.size(); i++)
            {
                std::cout << "Radius no." << i+1 << " is " << SphereRadius[i] << std::endl;
            }
            for (size_t i = 0; i < ShphereMats.size(); i++)
            {
                std::cout << "Material no." << i+1 << " is " << ShphereMats[i] << std::endl;
            }
        }
        int getTypesNumber() {
            return ShphereMats.size();
        }
        void resetLimit() {
            minX = 10000.f;
            minY = 10000.f;
            minZ = 10000.f;
            maxX = -10000.f;
            maxY = -10000.f;
            maxZ = -10000.f;
        }
        void updateLimit(float xPos, float yPos, float zPos) {
            if (xPos < minX) minX = xPos;
            if (yPos < minY) minY = yPos;
            if (zPos < minZ) minZ = zPos;

            if (xPos > maxX) maxX = xPos;
            if (yPos > maxY) maxY = yPos;
            if (zPos > maxZ) maxZ = zPos;
        }
        void getLimit(vec3& minLimit, vec3&maxLimit) {
            minLimit = vec3(minX, minY, minZ);
            maxLimit = vec3(maxX, maxY, maxZ);
        }

        void updateFileNumber(int num, int padding_amount) {
            // Clear all old data
            CompleteSphereList.clear();
            resetLimit();

            std :: string current_model_names;
            std :: stringstream ss02;
            ss02 << SphereName << std :: setw(padding_amount) << std :: setfill('0') << num << ".gt";
            SphereNameCurrent = ss02.str();
        }
        std::vector<std::shared_ptr<sphere>> getCompleteSphereList() {
            return CompleteSphereList;
        }

        std::string SphereNameCurrent;
        std::string SphereName;
        std::vector<float> SphereRadius;
        std::vector<std::string> ShphereMats;
        std::vector<std::shared_ptr<sphere>> CompleteSphereList;

        int fileNo;

    private:
        // Type of the input
        const int type = 0;
        // Limit of models
        float minX = 10000.f;
        float minY = 10000.f;
        float minZ = 10000.f;
        float maxX = -10000.f;
        float maxY = -10000.f;
        float maxZ = -10000.f;
};
extern std::vector<sphereInputHandler> sphereInputList;

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

// EXR
extern bool hasExr;
extern bool hasOverlayBackgroundColor;
extern std :: string exrFile;

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
            else if (words[0] == "SphereInputType:") {
                if (words[1] == "1")
                {
                    sphereInputHandler sphereConfigReader(1);
                    std :: getline (SetUp, myline);
                    std :: istringstream issName(myline);
                    std :: vector<std::string> sphereWords;
                    while (std :: getline(issName, myword, ' '))
                        sphereWords.push_back(myword);
                    sphereConfigReader.SphereName = sphereWords[1];

                    std :: getline (SetUp, myline);
                    std :: istringstream issRad(myline);
                    sphereWords.clear();
                    while (std :: getline(issRad, myword, ' '))
                        sphereWords.push_back(myword);
                    sphereConfigReader.SphereRadius.push_back(std::stof(sphereWords[1]));

                    std :: getline (SetUp, myline);
                    std :: istringstream issMat(myline);
                    sphereWords.clear();
                    while (std :: getline(issMat, myword, ' '))
                        sphereWords.push_back(myword);
                    for (size_t i = 1; i < sphereWords.size(); i++)
                    {
                        sphereConfigReader.ShphereMats.push_back(sphereWords[i]);   
                    }

                    sphereInputList.push_back(sphereConfigReader);
                }
                else if (words[1] == "2")
                {
                    sphereInputHandler sphereConfigReader(2);
                    std :: getline (SetUp, myline);
                    std :: istringstream issName(myline);
                    std :: vector<std::string> sphereWords;
                    while (std :: getline(issName, myword, ' '))
                        sphereWords.push_back(myword);
                    sphereConfigReader.SphereName = sphereWords[1];

                    std :: getline (SetUp, myline);
                    std :: istringstream issMat(myline);
                    sphereWords.clear();
                    while (std :: getline(issMat, myword, ' '))
                        sphereWords.push_back(myword);
                    for (size_t i = 1; i < sphereWords.size(); i++)
                    {
                        sphereConfigReader.ShphereMats.push_back(sphereWords[i]);   
                    }

                    sphereInputList.push_back(sphereConfigReader);
                }
                else if (words[1] == "3")
                {
                    sphereInputHandler sphereConfigReader(3);
                    std :: getline (SetUp, myline);
                    std :: istringstream issName(myline);
                    std :: vector<std::string> sphereWords;
                    while (std :: getline(issName, myword, ' '))
                        sphereWords.push_back(myword);
                    sphereConfigReader.SphereName = sphereWords[1];

                    std :: getline (SetUp, myline);
                    std :: istringstream issRad(myline);
                    sphereWords.clear();
                    while (std :: getline(issRad, myword, ' '))
                        sphereWords.push_back(myword);
                    for (size_t i = 1; i < sphereWords.size(); i++)
                    {
                        sphereConfigReader.SphereRadius.push_back(std::stof(sphereWords[i])); 
                    }

                    std :: getline (SetUp, myline);
                    std :: istringstream issMat(myline);
                    sphereWords.clear();
                    while (std :: getline(issMat, myword, ' '))
                        sphereWords.push_back(myword);
                    for (size_t i = 1; i < sphereWords.size(); i++)
                    {
                        sphereConfigReader.ShphereMats.push_back(sphereWords[i]);   
                    }

                    sphereInputList.push_back(sphereConfigReader);
                }
                else if (words[1] == "4")
                {
                    sphereInputHandler sphereConfigReader(4);
                    std :: getline (SetUp, myline);
                    std :: istringstream issName(myline);
                    std :: vector<std::string> sphereWords;
                    while (std :: getline(issName, myword, ' '))
                        sphereWords.push_back(myword);
                    sphereConfigReader.SphereName = sphereWords[1];

                    std :: getline (SetUp, myline);
                    std :: istringstream issMat(myline);
                    sphereWords.clear();
                    while (std :: getline(issMat, myword, ' '))
                        sphereWords.push_back(myword);
                    for (size_t i = 1; i < sphereWords.size(); i++)
                    {
                        sphereConfigReader.ShphereMats.push_back(sphereWords[i]);   
                    }

                    sphereInputList.push_back(sphereConfigReader);
                }
                else {
                    std :: cerr << "Find an unrecognized sphere input type" << std :: endl;
                    return false;
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
            else if (words[0] == "cameras_batch:") {
                std::string movementFile = words[1];
                std::ifstream cameraSetup;
                cameraSetup.open(movementFile);

                std :: string cameraline;
                std :: string cameraword;

                if(cameraSetup.is_open()) {
                    while (cameraSetup)
                    {
                        std :: getline (cameraSetup, cameraline);
                        std :: istringstream iss(cameraline);

                        std :: vector<std::string> camerawords;
                        while (std :: getline(iss, cameraword, ' '))
                            camerawords.push_back(cameraword);

                        if (camerawords.size() == 0) continue;

                        vec3 tempVup(std::stof(camerawords[0]), std::stof(camerawords[1]), std::stof(camerawords[2]));
                        vec3 tempLookfrom(std::stof(camerawords[3]), std::stof(camerawords[4]), std::stof(camerawords[5]));
                        vec3 tempLookat(std::stof(camerawords[6]), std::stof(camerawords[7]), std::stof(camerawords[8]));

                        vupList.push_back(tempVup);
                        lookfromList.push_back(tempLookfrom);
                        lookatList.push_back(tempLookat);
                    }
                }

                // Debug
                // for (size_t iVp = 0; iVp < vupList.size(); iVp++)
                // {
                //     std::cout << "Debugging -> Vp:" << vupList[iVp] << std::endl;
                // }
                // for (size_t iLookfrom = 0; iLookfrom < lookfromList.size(); iLookfrom++)
                // {
                //     std::cout << "Debugging -> Lookfrom:" << lookfromList[iLookfrom] << std::endl;
                // }
                // for (size_t iLookat = 0; iLookat < lookatList.size(); iLookat++)
                // {
                //     std::cout << "Debugging -> Lookat:" << lookatList[iLookat] << std::endl;
                // }
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
            else if (words[0] == "BackgroundExrForIllumination:") {
                exrFile = words[1];
                hasExr = true;
                std :: cerr << "Ambient light set with file : " << exrFile << std :: endl;
            }
            else if (words[0] == "BackgroundColorForOutput:") {
                hasOverlayBackgroundColor = true;
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