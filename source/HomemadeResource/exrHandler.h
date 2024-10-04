#ifndef EXR_H
#define EXR_H

#define TINYEXR_IMPLEMENTATION
#include "../ImportedResource/tinyexr.h"

#include "vec3.h"
#include "ray.h"

class exrHandler {
    public:
        exrHandler() {}
        void setFile (std::string _exrFileName){
            exrFileName = _exrFileName;
            const char* err = NULL;

            int ret = LoadEXR(&outMap, &width, &height, exrFileName.data(), &err);

            if (ret != TINYEXR_SUCCESS) {
                if (err) {
                    fprintf(stderr, "ERR : %s\n", err);
                    FreeEXRErrorMessage(err); // release memory of error message.
                }
            }
        }
        ~exrHandler() {free(outMap);}

        color getSampleFromRay(const ray& lightRay) const {
            vec3 rayDirection = lightRay.direction();
            float theta = atan2(rayDirection.y(), rayDirection.x());
            float phi = acos(rayDirection.z() / rayDirection.length());

            float u = (theta + Pi)/(2.f * Pi);
            float v = phi/Pi;

            int xPosition = static_cast<int>(u*width);
            int yPosition = static_cast<int>(v*height);
            int indexNum = 4 * (yPosition * width + xPosition);

            // Debug
            // std::cout << color(outMap[indexNum], outMap[indexNum + 1], outMap[indexNum + 2]) << std::endl;

            return color(outMap[indexNum], outMap[indexNum + 1], outMap[indexNum + 2]);
        }

    private:
        std::string exrFileName;
        int width;
        int height;
        float* outMap = nullptr;
};

#endif