#include <stdint.h>
#include "utilities.hpp"

class Limiter {
    public:
        Limiter(int attack, float release);
        ~Limiter();
        void process(float *data, int64_t frames, int channels, float *lookaheadData, int64_t lookaheadFrames);
        int64_t getTotalFrames() {return totalFrames;}
        int64_t getLimitedFrames() {return limitedFrames;}
    private:
        float attack, release;
        float gain;
        float *attackCurve;
        int64_t limitedFrames;
        int64_t totalFrames;
};
