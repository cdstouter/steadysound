#include <cmath>
#include <sndfile.h>

float calculateRMS(float *data, long long samples) {
    float mean = 0.0;
    for (long long i = 0; i<samples; i++) {
        float square = pow(data[i], 2.0);
        mean += square;
    }
    mean = pow(mean / (float)samples, 0.5);
    return mean;
}

float VtoDB(float v) {
    return 20.0f * log10(v);
}

float DBtoV(float db) {
    return pow(10.0, db / 20.0);
}

struct gainPoint {
    float gain;
    sf_count_t position;
};

