#include "limiter.hpp"

Limiter::Limiter(int attack, float release) {
    this->attack = attack;
    this->release = release;

    // set up the attack curve
    attackCurve = new float[attack];
    for (int i=0; i<attack; i++) {
        float fade = 1.0 - ((float)i / (float)attack);
        // apply a power curve to the fade?
        //fade = std::pow(fade, 0.5);
        attackCurve[i] = fade;
    }

    gain = 0.0;
    limitedFrames = 0;
    totalFrames = 0;
}

Limiter::~Limiter() {
    delete[] attackCurve;
}

void Limiter::process(float *data, int64_t frames, int channels, float *lookaheadData, int64_t lookaheadFrames) {
    for (int64_t frame = 0; frame < frames; frame++) {
        totalFrames++;
        // calculate attack
        float max = 0.0;
        for (int64_t ahead = 0; ahead < attack; ahead++) {
            float sample;
            int64_t aheadFrame = frame + ahead;
            if (aheadFrame < frames) {
                for (int channel = 0; channel < channels; channel++) {
                    sample = data[aheadFrame * channels + channel];
                    sample = std::abs(sample * attackCurve[ahead]);
                    if (sample > max) max = sample;
                }
            } else if (lookaheadData && (aheadFrame - frames) < lookaheadFrames) {
                aheadFrame = aheadFrame - frames;
                for (int channel = 0; channel < channels; channel++) {
                    sample = lookaheadData[aheadFrame * channels + channel];
                    sample = std::abs(sample * attackCurve[ahead]);
                    if (sample > max) max = sample;
                }
            }
        }
        float neededGain = -(VtoDB(max) + 0.001);
        if (neededGain < gain) {
            gain = neededGain;
        }
        // apply the gain
        float mult = DBtoV(gain);
        for (int channel = 0; channel < channels; channel++) {
            data[frame * channels + channel] *= mult;
        }
        if (mult < 1.0) limitedFrames++;
        // release
        gain += release;
        if (gain > 0) gain = 0.0;
    }
}
