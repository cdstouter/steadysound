#include <sndfile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "utilities.hpp"
#include "limiter.hpp"

namespace po = boost::program_options;

int main(int argc, char **argv) {
    std::string inputFilename = "";
    std::string outputFilename = "";
    int blockSize = 1024;
    float minDb = -45;
    float targetDb = -20;
    int windowSize = 20;
    bool limiterUsed = true;
    float limiterRelease = 0.0006;
    int limiterAttack = 1024;
    bool checkSilence = false;
    float predictive = 0.5;

    po::options_description prettyDesc("Options"); // this version doesn't include the positional arguments
    prettyDesc.add_options()
        ("help", "Show this help information")
        ("target,t", po::value<float>(&targetDb)->required(), "Set target dB level")
        ("silence,s", po::value<float>(&minDb)->default_value(-45.0), "Minimum dB value for audio content, anything below this is considered silence.")
        ("check-silence", "Skips dynamics processing and cuts out silence from the output file. Use to test the silence level.")
        ("window,w", po::value<int>(&windowSize)->default_value(20), "Window size for smoothing volume changes.")
        ("predictive", po::value<float>(&predictive)->default_value(0.5), "Predictive factor. 1.0 is fully predictive, 0.0 is fully reactive.")
        ("block-size", po::value<int>(&blockSize)->default_value(1024), "Block size for calculating RMS values.")
        ("limiter-attack", po::value<int>(&limiterAttack)->default_value(4), "Limiter attack time in samples. Increasing this value will directly increase processing time.")
        ("limiter-release", po::value<float>(&limiterRelease)->default_value(0.0006), "Limiter release time in dB/sample.")
        ("limiter-disable", "Disable the limiter completely - may cause clipping.")
    ;
    po::options_description desc("Options"); // this one is actually used to parse, don't use required()
    desc.add_options()
        ("help", "Show this help information")
        ("input-file", po::value<std::string>(&inputFilename))
        ("output-file", po::value<std::string>(&outputFilename))
        ("target,t", po::value<float>(&targetDb), "Set target dB level")
        ("silence,s", po::value<float>(&minDb)->default_value(-45.0), "Minimum dB value for audio content, anything below this is considered silence.")
        ("check-silence", "Skips dynamics processing and cuts out silence from the output file. Use to test the silence level.")
        ("window,w", po::value<int>(&windowSize)->default_value(20), "Window size for smoothing volume changes.")
        ("predictive", po::value<float>(&predictive)->default_value(0.5), "Predictive factor. 1.0 is fully predictive, 0.0 is fully reactive.")
        ("block-size", po::value<int>(&blockSize)->default_value(1024), "Block size for calculating RMS values.")
        ("limiter-attack", po::value<int>(&limiterAttack)->default_value(4), "Limiter attack time in samples. Increasing this value will directly increase processing time.")
        ("limiter-release", po::value<float>(&limiterRelease)->default_value(0.0006), "Limiter release time in dB/sample.")
        ("disable-limiter", "Disable the limiter completely - may cause clipping.")
    ;
    po::positional_options_description p;
    p.add("input-file", 1).add("output-file", 1);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);
    } catch(std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    if (vm.count("help")) {
        std::cout << "Usage: " << argv[0] << " input.wav output.wav OPTIONS" << std::endl << std::endl;
        std::cout << prettyDesc << std::endl;
        return 1;
    }

    // do we have the required values?
    if (!vm.count("input-file") || !vm.count("output-file")) {
        std::cout << "The input and output filenames are required." << std::endl;
        return 1;
    }
    if (vm.count("check-silence")) {
        checkSilence = true;
    }
    if (!vm.count("target") && !checkSilence) {
        std::cout << "The target dB level (-t, --target)is required." << std::endl;
        return 1;
    }
    if (vm.count("disable-limiter")) {
        limiterUsed = false;
    }
    // check that everything's within a valid range
    if (blockSize < 32) {
        std::cout << "Block size must be >= 32." << std::endl;
        return 1;
    }
    if (minDb >= 0) {
        std::cout << "Silence level must be < 0 dB." << std::endl;
        return 1;
    }
    if (targetDb >= 0) {
        std::cout << "Target level must be < 0 dB." << std::endl;
        return 1;
    }
    if (windowSize < 1) {
        std::cout << "Window size must be >= 1." << std::endl;
        return 1;
    }
    if (predictive < 0.0 && predictive > 1.0) {
        std::cout << "Predictive factor must be between 0 and 1." << std::endl;
        return 1;
    }
    if (limiterRelease <= 0) {
        std::cout << "Limiter release must be > 0." << std::endl;
        return 1;
    }
    if (limiterAttack < 1) {
        std::cout << "Limiter attack must be > 0." << std::endl;
        return 1;
    }

    SF_INFO info;
    info.format = 0;
    SNDFILE* infile = sf_open(inputFilename.c_str(), SFM_READ, &info);

    if (!infile) {
        std::cout << "Error opening input file: " << sf_strerror(NULL) << std::endl;
        return 1;
    }

    SNDFILE *outfile = NULL;
    if (checkSilence) {
        SF_INFO outinfo;
        outinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
        outinfo.samplerate = info.samplerate;
        outinfo.channels = info.channels;
        outfile = sf_open(outputFilename.c_str(), SFM_WRITE, &outinfo);

        if (!outfile) {
            std::cout << "Error opening output file: " << sf_strerror(NULL) << std::endl;
            sf_close(infile);
            return 1;
        }
    }

    // the file was opened, let's calculate some RMS!
    std::cout << "Calculating RMS values..." << std::endl;
    std::vector<float> rmsBlocks;
    sf_count_t frames;
    float *data = new float[blockSize * info.channels];
    do {
        frames = sf_readf_float(infile, data, blockSize);
        float rms = calculateRMS(data, frames * info.channels);
        float rmsDb = VtoDB(rms);
        rmsBlocks.push_back(rmsDb);
        if (checkSilence && rmsDb > minDb) {
            sf_writef_float(outfile, data, frames);
        }
    } while (frames == blockSize);
    std::cout << rmsBlocks.size() << " RMS blocks calculated." << std::endl;
    delete[] data;

    if (checkSilence) {
        std::cout << "Done." << std::endl;

        sf_close(infile);
        sf_close(outfile);
        return 0;
    }

    // use a ring buffer to calculate a moving average over the RMS blocks
    int ringBufferSize = windowSize * 2 + 1;
    float ringBuffer[ringBufferSize];
    for (int i=0; i<ringBufferSize; i++) {
        ringBuffer[i] = -1000.0;
    }
    int ringBufferPos = 0;
    std::vector<float>::size_type currentBlock = 0;
    // start filling the ring buffer
    int ringBufferPrefill = 1 + (int)(predictive * (float)windowSize * 2.0);
    for (int i=0; i<ringBufferPrefill; i++) {
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
    float maximumGain = 0;
    float minimumGain = 0;
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
            if (gainPoints.size() == 1 || gain > maximumGain) maximumGain = gain;
            if (gainPoints.size() == 1 || gain < minimumGain) minimumGain = gain;
        }
    }
    std::cout << gainPoints.size() << " gain points calculated." << std::endl;
    if (gainPoints.size() == 0) {
        std::cout << "Error: no gain points found, nothing to do." << std::endl;
        sf_close(infile);
        return 1;
    }
    std::cout << "Minimum gain: " << minimumGain << " dB." << std::endl;
    std::cout << "Maximum gain: " << maximumGain << " dB." << std::endl;

    SF_INFO outinfo;
    outinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    outinfo.samplerate = info.samplerate;
    outinfo.channels = info.channels;
    outfile = sf_open(outputFilename.c_str(), SFM_WRITE, &outinfo);

    if (!outfile) {
        std::cout << "Error opening output file: " << sf_strerror(NULL) << std::endl;
        sf_close(infile);
        return 1;
    }

    // set up the limiter
    Limiter *limiter = NULL;
    if (limiterUsed) limiter = new Limiter(limiterAttack, limiterRelease);

    // now we go through the file and apply the gain!
    std::cout << "Applying gain..." << std::endl;
    int processSize = 1024;
    if (limiterUsed && limiterAttack > processSize) processSize = limiterAttack;
    float *backingData1 = new float[processSize * info.channels];
    float *backingData2 = new float[processSize * info.channels];
    data = backingData1;
    float *oldData = backingData2;
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
            float gainMult = DBtoV(gain);
            // apply the gain
            float highestValue = 0.0;
            for (int j=0; j<info.channels; j++) {
                float sample = data[i * info.channels + j];
                sample = sample * gainMult;
                if (!limiterUsed && sample > 1.0) clipped++;
                if (sample > highestValue) highestValue = sample;
                data[i * info.channels + j] = sample;
            }
            currentFrame++;
        }
        // apply the limiter
        limiter->process(oldData, oldFrames, info.channels, data, frames);
        // write the frames to the new file
        if (blockNum > 0) sf_writef_float(outfile, oldData, oldFrames);
        blockNum++;
        if (totalBlocks > 0) {
            int thisPerc = (blockNum * 100) / totalBlocks;
            if (thisPerc > currentPerc) {
                currentPerc = thisPerc;
                std::cout << ">" << std::setw(3) << currentPerc << "% done" << "\r" << std::flush;
            }
        }
    } while (frames == processSize);
    std::cout << "Done.     " << std::endl;

    limiter->process(data, frames, info.channels, NULL, 0);
    sf_writef_float(outfile, data, frames);
    if (clipped) std::cout << "WARNING: " << clipped << " samples clipped and limiter disabled" << std::endl;

    if (limiterUsed) {
        int limitedPercent = 100.0 * ((float)limiter->getLimitedFrames() / (float)limiter->getTotalFrames());
        if (limiter->getLimitedFrames() > 0 && limitedPercent == 0) {
            std::cout << "Limiter applied to <1% of audio." << std::endl;
        } else {
            std::cout << "Limiter applied to " << limitedPercent << "% of audio." << std::endl;
        }
        if (limitedPercent > 10) {
            std::cout << "WARNING: You may want to lower the target dB level so that less limiting is applied." << std::endl;
        }
        delete limiter;
    }
    delete[] backingData1;
    delete[] backingData2;

    sf_close(infile);
    sf_close(outfile);

    return 0;
}
