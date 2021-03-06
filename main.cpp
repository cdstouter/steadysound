#include <sndfile.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "version.hpp"
#include "utilities.hpp"
#include "limiter.hpp"
#include "ringbuffer.hpp"

namespace po = boost::program_options;

int main(int argc, char **argv) {
    std::string inputFilename = "";
    std::string outputFilename = "";
    int blockSize = 1024;
    float minDb = -80;
    bool targetDbSet = false;
    float targetDb = -20.0;
    bool limiterUsed = true;
    float limiterRelease = 25;
    int limiterAttack = 1024;
    bool checkSilence = false;
    bool analyze = false;
    float correction = 1.0;
    int lookAhead, lookBehind;

    po::options_description optionsMain("Main options"); // this version doesn't include the positional arguments
    optionsMain.add_options()
        ("help", "Show this help information")
        ("target,t", po::value<float>(&targetDb), "Set target dB level. Set to the median RMS dB level if not specified.")
        ("silence,s", po::value<float>(&minDb)->default_value(-80.0), "Minimum dB value for audio content, anything below this is considered silence.")
        ("correction,c", po::value<float>(&correction)->default_value(1.0), "Correction factor. 0 means no volume levelling effect, 1 means full effect.")
        ("block-size,b", po::value<int>(&blockSize)->default_value(500), "Block size in msec for calculating RMS values.")
        ("lookahead,A", po::value<int>(&lookAhead)->default_value(10), "Number of RMS blocks to look ahead at.")
        ("lookbehind,B", po::value<int>(&lookBehind)->default_value(10), "Number of RMS blocks to look behind at.")
    ;
    po::options_description optionsAnalysis("Analysis options");
    optionsAnalysis.add_options()
        ("analyze", "Skips dynamics processing and outputs a GnuPlot file with an RMS graph of the input file.")
        ("check-silence", "Skips dynamics processing and cuts out silence from the output file. Use to test the silence level.")
    ;
    po::options_description optionsPeakLimiter("Peak limiter options");
    optionsPeakLimiter.add_options()
        ("limiter-attack", po::value<int>(&limiterAttack)->default_value(4), "Limiter attack time in samples. Increasing this value will directly increase processing time.")
        ("limiter-release", po::value<float>(&limiterRelease)->default_value(25.0), "Limiter release time in dB/second.")
        ("limiter-disable,L", "Disable the limiter completely - may cause clipping.")
    ;
    po::options_description optionsPositionals("Options"); // this one is actually used to parse, don't use required()
    optionsPositionals.add_options()
        ("input-file", po::value<std::string>(&inputFilename))
        ("output-file", po::value<std::string>(&outputFilename))
    ;
    po::positional_options_description optionsPositionalsDescription;
    optionsPositionalsDescription.add("input-file", 1).add("output-file", 1);

    //prettyDesc.add(positionals);
    po::options_description optionsAll;
    optionsAll.add(optionsMain);
    optionsAll.add(optionsAnalysis);
    optionsAll.add(optionsPeakLimiter);
    optionsAll.add(optionsPositionals);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(optionsAll).positional(optionsPositionalsDescription).run(), vm);
        po::notify(vm);
    } catch(std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    if (argc == 1 || vm.count("help")) {
        std::cout << "SteadySound v" << STEADYSOUND_MAJOR_VERSION << "." << STEADYSOUND_MINOR_VERSION << std::endl << std::endl;
        std::cout << "Usage: " << argv[0] << " input output OPTIONS" << std::endl << std::endl;
        std::cout << optionsMain << std::endl;
        std::cout << optionsAnalysis << std::endl;
        std::cout << optionsPeakLimiter << std::endl;
        std::cout << "The input file can be in any format that libsndfile supports, which includes ";
        std::cout << "WAV, AIFF, FLAC, and OGG. Any samplerate and number of channels is supported. ";
        std::cout << "The output file will be saved in the same format as the input file, no matter what extension is given. ";
        std::cout << "Note: MP3 is not supported." << std::endl;
        return 1;
    }

    // do we have the required values?
    if (vm.count("analyze")) {
        analyze = true;
    }
    if (!vm.count("input-file")) {
        std::cout << "The input filename is required." << std::endl;
        return 1;
    }
    if (!vm.count("output-file") && !analyze) {
        std::cout << "The output filename is required." << std::endl;
        return 1;
    }
    if (vm.count("check-silence")) {
        checkSilence = true;
    }
    if (vm.count("target")) {
        targetDbSet = true;
    }
    if (vm.count("disable-limiter")) {
        limiterUsed = false;
    }
    // check that everything's within a valid range
    if (blockSize < 10) {
        std::cout << "Block size must be >= 10." << std::endl;
        return 1;
    }
    if (minDb >= 0) {
        std::cout << "Silence level must be < 0 dB." << std::endl;
        return 1;
    }
    if (correction < 0.0 || correction > 1.0) {
        std::cout << "Correction factor must be between 0 and 1." << std::endl;
        return 1;
    }
    if (targetDb >= 0) {
        std::cout << "Target level must be < 0 dB." << std::endl;
        return 1;
    }
    if (lookAhead < 0) {
        std::cout << "Lookahead must be >= 0." << std::endl;
        return 1;
    }
    if (lookBehind < 0) {
        std::cout << "Lookbehind must be >= 0." << std::endl;
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
    std::ofstream gnuplot;
    if (checkSilence) {
        SF_INFO outinfo;
        outinfo.format = info.format;
        outinfo.samplerate = info.samplerate;
        outinfo.channels = info.channels;
        outfile = sf_open(outputFilename.c_str(), SFM_WRITE, &outinfo);

        if (!outfile) {
            std::cout << "Error opening output file: " << sf_strerror(NULL) << std::endl;
            sf_close(infile);
            return 1;
        }
    }

    if (analyze) {
        gnuplot.open((inputFilename + ".gp").c_str());
        if (!gnuplot.is_open()) {
            std::cout << "Error opening output GnuPlot file." << std::endl;
            sf_close(infile);
            if (checkSilence) sf_close(outfile);
            return 1;
        }
        gnuplot << "# Time RMS" << std::endl << "$data << EOD" << std::endl;
    }

    // convert blockSize from msec to samples
    blockSize = (blockSize * info.samplerate) / 1000;

    // the file was opened, let's calculate some RMS!
    std::cout << "Calculating RMS values..." << std::endl;
    std::vector<float> rmsBlocks;
    sf_count_t frames;
    float *data = new float[blockSize * info.channels];
    float medianRMS = 0.0;
    int64_t medianRMSCount = 0;
    sf_count_t currentFrame = 0;
    float currentSeconds = 0;
    do {
        frames = sf_readf_float(infile, data, blockSize);
        float rms = calculateRMS(data, frames * info.channels);
        float rmsDb = VtoDB(rms);
        rmsBlocks.push_back(rmsDb);
        if (rmsDb > minDb) {
            medianRMS += rmsDb;
            medianRMSCount++;
            if (checkSilence) {
                sf_writef_float(outfile, data, frames);
            }
        }
        if (analyze) {
            currentSeconds = (float)currentFrame / (float)info.samplerate;
            gnuplot << currentSeconds << " " << rmsDb << std::endl;
            currentFrame += frames;
        }
    } while (frames == blockSize);
    std::cout << rmsBlocks.size() << " RMS blocks calculated." << std::endl;
    medianRMS = medianRMS / medianRMSCount;
    std::cout << "RMS median: " << medianRMS << " dB." << std::endl;
    // if the target dB isn't set, use the median as the target
    if (!targetDbSet) {
        targetDb = medianRMS;
        std::cout << "Target dB set to " << targetDb << " dB." << std::endl;
    }
    delete[] data;

    if (checkSilence || analyze) {
        std::cout << "Done." << std::endl;
    }
    if (checkSilence) {
        std::cout << "The input audio has been written to " << outputFilename << " with all detected silence removed." << std::endl;
        sf_close(infile);
        sf_close(outfile);
    }
    if (analyze) {
        gnuplot << "EOD" << std::endl << std::endl;
        gnuplot << "set title \"" << inputFilename << "\"" << std::endl;
        gnuplot << "set xlabel \"Time (s)\"" << std::endl;
        gnuplot << "set xrange [0:" << currentSeconds << "]" << std::endl;
        gnuplot << "set ylabel \"RMS (dB)\"" << std::endl;
        gnuplot << "set style line 1 lw 0.5 lc variable" << std::endl;
        gnuplot << "plot \"$data\" using 1:2:2 with lines lw 0.5 palette title \"\"" << std::endl;
        gnuplot.close();

        std::cout << std::endl << "An RMS graph of the input file has been written to the GnuPlot file " << inputFilename << ".gp. ";
        std::cout << "To view the RMS graph, run GnuPlot like this:" << std::endl << std::endl;
        std::cout << "gnuplot " << inputFilename << ".gp -" << std::endl;
    }
    if (analyze || checkSilence) return 0;

    // go through all of our blocks
    std::cout << "Calculating gain points..." << std::endl;
    std::vector<gainPoint> gainPoints;
    float maximumGain = 0;
    float minimumGain = 0;
    int minAverageBlocks = (int)((float)(lookAhead + lookBehind) * 0.75);
    if (minAverageBlocks == 0) minAverageBlocks = 1;
    for (std::vector<float>::size_type i=0; i<rmsBlocks.size(); i++) {
        float rms = -1000.0;
        int numBlocks = 0;
        // lookbehind
        for (int j=1; j<=lookBehind; j++) {
            int64_t pos = (int64_t)i - (int64_t)j;
            if (pos < 0) break;
            float dB = rmsBlocks[pos];
            if (dB > minDb) {
                float gain = 1.0 - ((float)j / (float)(lookBehind + 1));
                gain = std::pow(gain, 0.75);
                gain = VtoDB(gain);
                dB = dB + gain;
                if (dB > rms) rms = dB;
                numBlocks++;
            }
        }
        // lookahead
        for (int j=0; j<=lookAhead; j++) {
            int64_t pos = (int64_t)i + (int64_t)j;
            if (pos >= (int64_t)rmsBlocks.size()) break;
            float dB = rmsBlocks[pos];
            if (dB > minDb) {
                float gain = 1.0 - ((float)j / (float)(lookAhead + 1));
                gain = std::pow(gain, 0.75);
                gain = VtoDB(gain);
                dB = dB + gain;
                if (dB > rms) rms = dB;
                numBlocks++;
            }
        }
        if (numBlocks > minAverageBlocks) {
            float correctedGain = targetDb - rms;
            float uncorrectedGain = targetDb - medianRMS;
            float gain = (correction * correctedGain) + ((1.0 - correction) * uncorrectedGain);
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
    outinfo.format = info.format;
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
    if (limiterUsed) limiter = new Limiter(limiterAttack, limiterRelease / (float)info.samplerate);

    // now we go through the file and apply the gain!
    std::cout << "Applying gain..." << std::endl;
    int processSize = 1024;
    if (limiterUsed && limiterAttack > processSize) processSize = limiterAttack;
    float *backingData1 = new float[processSize * info.channels];
    float *backingData2 = new float[processSize * info.channels];
    sf_count_t frames1, frames2;
    sf_seek(infile, 0, SEEK_SET);
    currentFrame = 0;
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
        // read new data
        frames1 = sf_readf_float(infile, backingData1, processSize);
        for (sf_count_t i=0; i<frames1; i++) {
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
                float sample = backingData1[i * info.channels + j];
                sample = sample * gainMult;
                if (!limiterUsed && sample > 1.0) clipped++;
                if (sample > highestValue) highestValue = sample;
                backingData1[i * info.channels + j] = sample;
            }
            currentFrame++;
        }
        // flip the buffers
        std::swap(backingData1, backingData2);
        std::swap(frames1, frames2);
        // apply the limiter
        if (blockNum > 0) limiter->process(backingData1, frames1, info.channels, backingData2, frames2);
        // write the frames to the new file
        if (blockNum > 0) sf_writef_float(outfile, backingData1, frames1);
        blockNum++;
        if (totalBlocks > 0) {
            int thisPerc = (blockNum * 100) / totalBlocks;
            if (thisPerc > currentPerc) {
                currentPerc = thisPerc;
                std::cout << ">" << std::setw(3) << currentPerc << "% done" << "\r" << std::flush;
            }
        }
    } while (frames2 == processSize);
    std::cout << "Done.     " << std::endl;

    limiter->process(backingData2, frames2, info.channels, NULL, 0);
    sf_writef_float(outfile, backingData2, frames2);
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
