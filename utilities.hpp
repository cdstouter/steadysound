#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <stdint.h>

float calculateRMS(float *data, long long samples);

float VtoDB(float v);

float DBtoV(float db);

struct gainPoint {
    float gain;
    int64_t position;
};

#endif
