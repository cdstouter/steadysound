#include <sndfile.h>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

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

int main(void) {
    std::string filename = "input.wav";
    int blockSize = 1024;
    float minDb = -45;
    float targetDb = -20;
    int windowSize = 20;
    bool limiter = true;
    float limiterRelease = 0.0006; // dB per sample
    float limiterAttack = 1024; // number of samples

    SF_INFO info;
    info.format = 0;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &info);

    if (!infile) {
        std::cout << "Error opening file: " << sf_strerror(NULL) << std::endl;
        return 1;
    }


    // the file was opened, let's calculate some RMS!
    std::cout << "Calculating RMS values..." << std::endl;
    std::vector<float> rmsBlocks;
    sf_count_t frames;
    float *data = new float[blockSize * info.channels];
    do {
        frames = sf_readf_float(infile, data, blockSize);
        //std::cout << samples << std::endl;
        float rms = calculateRMS(data, frames * info.channels);
        float rmsDb = VtoDB(rms);
        rmsBlocks.push_back(rmsDb);
    } while (frames == blockSize);
    std::cout << rmsBlocks.size() << " RMS blocks calculated." << std::endl;
    delete[] data;

    // use a ring buffer to calculate a moving average over the RMS blocks
    int ringBufferSize = windowSize * 2 + 1;
    float ringBuffer[ringBufferSize];
    for (int i=0; i<ringBufferSize; i++) {
        ringBuffer[i] = -1000.0;
    }
    int ringBufferPos = 0;
    std::vector<float>::size_type currentBlock = 0;
    // start filling the ring buffer
    for (int i=0; i<windowSize; i++) {
        if (currentBlock < rmsBlocks.size()) {
            ringBuffer[ringBufferPos] = rmsBlocks[currentBlock];
        } else {
            ringBuffer[ringBufferPos] = -1000.0;
        }
        currentBlock++;
        ringBufferPos++;
    }
    // go through all of our blocks
    std::cout << "Calculating gain points..." << std::endl;
    std::vector<gainPoint> gainPoints;
    int minAverageBlocks = windowSize * 1.5;
    for (std::vector<float>::size_type i=0; i<rmsBlocks.size(); i++) {
        if (currentBlock < rmsBlocks.size()) {
            ringBuffer[ringBufferPos] = rmsBlocks[currentBlock];
        } else {
            ringBuffer[ringBufferPos] = -1000.0;
        }
        currentBlock++;
        ringBufferPos++;
        if (ringBufferPos >= ringBufferSize) ringBufferPos = 0;
        // calculate the average RMS if we have enough blocks that aren't silent
        int numBlocks = 0;
        float rms = 0;
        for (int j=0; j<ringBufferSize; j++) {
            if (ringBuffer[j] > minDb) {
                numBlocks++;
                rms += ringBuffer[j];
            }
        }
        if (numBlocks > minAverageBlocks) {
            rms = rms / (float)numBlocks;
            float gain = targetDb - rms;
            gainPoint gp;
            gp.gain = gain;
            gp.position = (i * blockSize) + (blockSize / 2);
            gainPoints.push_back(gp);
        }
    }
    std::cout << gainPoints.size() << " gain points calculated." << std::endl;
    if (gainPoints.size() == 0) {
        std::cout << "Error: no gain points found, nothing to do." << std::endl;
        sf_close(infile);
        return 1;
    }

    SF_INFO outinfo;
    outinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    outinfo.samplerate = info.samplerate;
    outinfo.channels = info.channels;
    SNDFILE *outfile = sf_open("output.wav", SFM_WRITE, &outinfo);

    if (!outfile) {
        std::cout << "Error opening file: " << sf_strerror(NULL) << std::endl;
        sf_close(infile);
        return 1;
    }

    // now we go through the file and apply the gain!
    float limiterGain = 0.0;
    std::cout << "Applying gain..." << std::endl;
    int processSize = 1024;
    if (limiter && limiterAttack > processSize) processSize = limiterAttack;
    data = new float[processSize * info.channels];
    float *oldData = new float[processSize * info.channels];
    sf_count_t oldFrames;
    sf_seek(infile, 0, SEEK_SET);
    sf_count_t currentFrame = 0;
    sf_count_t clipped = 0;
    int currentPerc = 0;
    int blockNum = 0;
    int totalBlocks = info.frames / processSize;
    gainPoint firstGainPoint = gainPoints[0];
    gainPoint lastGainPoint = gainPoints[gainPoints.size() - 1];
    int currentGainPointNumber = 0;
    gainPoint currentGainPoint = firstGainPoint;
    gainPoint nextGainPoint;
    if (gainPoints.size() > 1) nextGainPoint = gainPoints[1];
    do {
        // flip the buffers and read new data
        float *temp = oldData;
        data = oldData;
        oldData = temp;
        oldFrames = frames;
        frames = sf_readf_float(infile, data, processSize);
        for (sf_count_t i=0; i<frames; i++) {
            // calculate the gain
            float gain = 0.0;
            if (currentFrame <= firstGainPoint.position) {
                gain = firstGainPoint.gain;
            } else if (currentFrame >= lastGainPoint.position) {
                gain = lastGainPoint.gain;
            } else {
                // interpolate between gain points
                if (currentFrame >= nextGainPoint.position) {
                    currentGainPointNumber++;
                    currentGainPoint = gainPoints[currentGainPointNumber];
                    nextGainPoint = gainPoints[currentGainPointNumber + 1];
                }
                float pos = (float)(currentFrame - currentGainPoint.position) / (float)(nextGainPoint.position - currentGainPoint.position);
                gain = ((1.0f - pos) * currentGainPoint.gain) + (pos * nextGainPoint.gain);
            }
            //float gain = interpolateGain(gainPoints, currentFrame);
            if (limiter) gain += limiterGain;
            float gainMult = DBtoV(gain);
            // apply the gain
            float highestValue = 0.0;
            for (int j=0; j<info.channels; j++) {
                float sample = data[i * info.channels + j];
                sample = sample * gainMult;
                if (!limiter && sample > 1.0) clipped++;
                if (sample > highestValue) highestValue = sample;
                data[i * info.channels + j] = sample;
            }
            if (highestValue > 1.0 && limiter) {
                // kick the limiter in
                highestValue += .001;
                float additionalLimiterGain = -VtoDB(highestValue);
                limiterGain += additionalLimiterGain;
                // and drop the attack frames down
                for (int attack=0; attack<limiterAttack; attack++) {
                    int frame = i - attack;
                    float gain = additionalLimiterGain * ((float)(limiterAttack - attack) / (float)limiterAttack);
                    float mult = DBtoV(gain);
                    if (frame >= 0) {
                        for (int j=0; j<info.channels; j++) {
                            data[frame * info.channels + j] *= mult;
                        }
                    } else {
                        int oldFrame = frame + oldFrames;
                        for (int j=0; j<info.channels; j++) {
                            oldData[oldFrame * info.channels + j] *= mult;
                        }
                    }
                }
                // and drop our current frame down
                for (int j=0; j<info.channels; j++) {
                    data[i * info.channels + j] /= highestValue;
                }
                
            }
            currentFrame++;
            limiterGain += limiterRelease;
            if (limiterGain > 0) limiterGain = 0;
        }
        // write the frames to the new file
        if (blockNum > 0) sf_writef_float(outfile, oldData, oldFrames);
        blockNum++;
        if (totalBlocks > 0) {
            int thisPerc = (blockNum * 100) / totalBlocks;
            if (thisPerc > currentPerc) {
                currentPerc = thisPerc;
                std::cout << currentPerc << "% done" << std::endl;
            }
        }
    } while (frames == processSize);
    sf_writef_float(outfile, data, frames);
    if (clipped) std::cout << "WARNING: " << clipped << " samples clipped and limiter disabled" << std::endl;

    sf_close(infile);
    sf_close(outfile);

    return 0;
}
