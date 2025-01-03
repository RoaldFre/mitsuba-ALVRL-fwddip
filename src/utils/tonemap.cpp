/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/util.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>
#if defined(WIN32)
# include <mitsuba/core/getopt.h>
#endif
#ifdef MTS_OPENMP
# include <omp.h>
#endif

MTS_NAMESPACE_BEGIN

class Tonemap : public Utility {
public:
    void help() {
        cout << endl;
        cout << "Synopsis: Loads one or more EXR/RGBE images and writes tonemapped 8-bit PNG/JPGs";
        cout << endl;
        cout << "Usage: mtsutil tonemap [options] <EXR/RGBE file (s)>" << endl;
        cout << "Options/Arguments:" << endl;
        cout << "   -h             Display this help text" << endl << endl;
        cout << "   -g gamma       Specify the gamma value (The default is -1 => sRGB)" << endl << endl;
        cout << "   -m multiplier  Multiply the pixel values by 'multiplier' (Default = 1)" << endl << endl;
        cout << "   -b r,g,b       Color balance: apply the specified per-channel multipliers" << endl << endl;
        cout << "   -c x,y,w,h     Crop: tonemap a given rectangle instead of the entire image" << endl << endl;
        cout << "   -s w,h         Resize the output image to the specified resolution" << endl << endl;
        cout << "   -F name        Name of the resampling filter used by -s (Default = lanczos)" << endl << endl;
        cout << "   -r x,y,w,h,i   Add a rectangle at the specified position and intensity, e.g." << endl;
        cout << "                  to make paper figures. The intensity should be in [0, 255]." << endl << endl;
        cout << "   -f fmt         Request a certain output format (png/jpg, default:png)" << endl << endl;
        cout << "   -a             Require the output image to have an alpha channel" << endl << endl;
        cout << "   -p key,burn    Run Reinhard et al.'s photographic tonemapping operator. 'key'" << endl;
        cout << "                  between [0, 1] chooses between low and high-key images and" << endl;
        cout << "                  'burn' (also [0, 1]) controls how much highlights may burn out" << endl << endl;
        cout << "   -B fov         Apply a bloom filter that simulates scattering in the human" << endl;
        cout << "                  eye. Requires the approx. field of view of the images to be" << endl;
        cout << "                  processed in order to compute a point spread function." << endl << endl;
        cout << "   -x             Temporal coherence mode: activate this flag when tonemapping " << endl;
        cout << "                  frames of an animation using the '-p' option to avoid flicker" << endl << endl;
        cout << "   -o file        Save the output with a given filename" << endl << endl;
        cout << "   -t             Multithreaded: process several files in parallel" << endl << endl;
        cout << "   -M             Merge all images into one final, averaged image" << endl << endl;
        cout << "   -C num         Cutoff: if using -M, only merge <num> images and ignore the rest" << endl << endl;
        cout << "   -R val         When merging images with -M, first merge them in to this" << endl;
        cout << "                  many bins and drop the extremal bin values." << endl;
        cout << "                  (legacy: val = robustness fraction: discard this fraction" << endl;
        cout << "                  of extremal values (e.g. -R 0.05 discards the 5\% lowest and" << endl;
        cout << "                  5\% highest samples of the given files to merge)" << endl << endl;
        cout << "   -n             Don't produce a LDR tonemapped image, useful for combined use with -M" << endl << endl;
        cout << " The operations are ordered as follows: 1. crop, 2. bloom, 3. resize, 4. color" << endl;
        cout << " balance, 5. tonemap, 6. annotate. To simply process a directory full of EXRs" << endl;
        cout << " in parallel, run the following: 'mtsutil tonemap -t path-to-directory/*.exr'" << endl;
    }

    typedef struct {
        int r[5];
    } Rect;

    /**
     * Computes a bloom filter based on
     *
     * "Physically-Based Glare Effects for Digital Images" by
     * Greg Spencer, Peter Shirley, Kurt Zimmerman and Donald P. Greenberg
     * SIGGRAPH 1995
     */
    ref<Bitmap> computeBloomFilter(int size, Float fov) {
        ref<Bitmap> bitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat, Vector2i(size));

        Float scale       = 2.f / (size - 1),
              halfLength  = std::tan(.5f * degToRad(fov));

        Float *ptr = bitmap->getFloatData();
        double sum = 0;

        for (int y=0; y<size; ++y) {
            for (int x=0; x<size; ++x) {
                Float xf = x*scale - 1,
                      yf = y*scale - 1,
                      r = std::sqrt(xf*xf+yf*yf),
                      angle = radToDeg(std::atan(r * halfLength)),
                      tmp   = angle + 0.02f,
                      f0 = 2.61e6f * math::fastexp(-2500*angle*angle),
                      f1 = 20.91f / (tmp*tmp*tmp),
                      f2 = 72.37f / (tmp*tmp),
                      f  = 0.384f*f0 + 0.478f*f1 + 0.138f*f2;

                *ptr++ = f;
                sum += f;
            }
        }
        ptr = bitmap->getFloatData();
        Float normalization = (Float) (1/sum);
        for (int i=0; i<size*size; ++i)
            *ptr++ *= normalization;

        return bitmap;
    }


    void mergeBitmaps(const std::vector<ref<Bitmap>> &bitmaps, ref<Bitmap> &result, size_t numToDrop=0,
            size_t beginIdx = 0, size_t endIdx = (size_t)-1) {
        size_t numBitmaps = bitmaps.size();

        if (endIdx == (size_t)-1)
            endIdx = numBitmaps;

        Log(EInfo, "Got %d bitmaps, merge %d - %d, drop %d",
                numBitmaps, beginIdx, endIdx, numToDrop);

        if (beginIdx >= endIdx)
            Log(EError, "Cannot merge <= 0 bitmaps!");
        size_t numToMerge = endIdx - beginIdx;

        if (numToDrop >= numToMerge)
            Log(EError, "Cannot drop %d samples when I only have to merge %d bitmaps!",
                    numToDrop, numToMerge);

        // make sure we don't lose precision in accumulation -> convert to doubles
        result = result->convert(result->getPixelFormat(), Bitmap::EFloat64);

        size_t nEntries =
            (size_t) bitmaps[0]->getSize().x *
            (size_t) bitmaps[0]->getSize().y *
            bitmaps[0]->getChannelCount();

        std::vector<double> samples(numToMerge);
        for (size_t i=0; i<nEntries; ++i) {
            for (size_t b=0; b<numToMerge; b++)
                samples[b] = bitmaps[beginIdx + b]->getFloat64Data()[i];
            std::sort(samples.begin(), samples.end()); 
            /* note: full sort is somewhat overkill for cutting the 
             * top and bottom part, but it *does* make the following 
             * addition numerically robust! */
            double avg = 0;
            for (size_t j = numToDrop; j < numToMerge-numToDrop; j++) {
                avg += samples[j];
            }
            avg = avg / (numToMerge - 2*numToDrop);

            result->getFloat64Data()[i] = avg;
        }
    }

    int run(int argc, char **argv) {
        ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
        int optchar;
        char *end_ptr = NULL;
        optind = 1;
        Float gamma = -1, multiplier = 1;
        Bitmap::EFileFormat format = Bitmap::EPNG;
        Float cbal[] = {1, 1, 1};
        int crop[] = {0, 0, -1, -1};
        int resize[] = {-1, -1};
        Float tonemapper[] = {-1, -1};
        bool temporalCoherence = false;
        std::vector<Rect> rects;
        std::string outputFilename;
        Bitmap::EPixelFormat pixelFormat = Bitmap::ERGB;
        Float logAvgLuminance = 0, maxLuminance = 0;
        bool runParallel = false;
        ReconstructionFilter *rfilter = NULL;
        Float bloomFov = 0;
        std::string rfilterName = "lanczos";
        bool merge = false;
        bool noLDR = false;
        Float robustFraction = 0;
        size_t maxNumBitmapsToMerge = 0;

        /* Parse command-line arguments */
        while ((optchar = getopt(argc, argv, "htMC:R:xag:m:f:r:b:c:o:p:s:B:F:n")) != -1) {
            switch (optchar) {
                case 'h': {
                        help();
                        return 0;
                    }
                    break;

                case 'g':
                    gamma = (Float) strtod(optarg, &end_ptr);
                    if (*end_ptr != '\0')
                        SLog(EError, "Could not parse the gamma value!");
                    break;

                case 'x':
                    temporalCoherence = true;
                    break;

                case 'f': {
                      std::string fmt = boost::to_lower_copy(std::string(optarg));
                      if (fmt == "png")
                          format = Bitmap::EPNG;
                      else if (fmt == "jpg" || fmt == "jpeg")
                          format = Bitmap::EJPEG;
                      else
                          SLog(EError, "Unknown format! (must be png/jpg)");
                    }
                    break;

                case 'F':
                    rfilterName = std::string(optarg);
                    break;

                case 'B':
                    bloomFov = (Float) strtod(optarg, &end_ptr);
                    #if !defined(MTS_HAS_FFTW)
                        Log(EWarn, "Applying a bloom filter without FFTW support compiled into "
                            "Mitsuba is likely going to be very, very slow!");
                    #endif
                    if (*end_ptr != '\0')
                        SLog(EError, "Could not parse the bloom field of view!");
                    break;


                case 'm':
                    multiplier = (Float) strtod(optarg, &end_ptr);
                    if (*end_ptr != '\0')
                        SLog(EError, "Could not parse the multiplier!");
                    break;

                case 'a':
                    pixelFormat = Bitmap::ERGBA;
                    break;

                case 'b': {
                        std::vector<std::string> tokens = tokenize(optarg, ", ");
                        if (tokens.size() != 3)
                            Log(EError, "Invalid color balancing parameter!");
                        for (int i=0; i<3; ++i) {
                            cbal[i] = (Float) std::strtod(tokens[i].c_str(), &end_ptr);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse floating point number "
                                    "in color balancing parameter!");
                        }
                    }
                    break;

                case 'c': {
                        std::vector<std::string> tokens = tokenize(optarg, ", ");
                        if (tokens.size() != 4)
                            Log(EError, "Invalid crop parameter!");
                        for (int i=0; i<4; ++i) {
                            crop[i] = (int) std::strtol(tokens[i].c_str(), &end_ptr, 10);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse integer in crop parameter!");
                        }
                    }
                    break;

                case 'p': {
                        std::vector<std::string> tokens = tokenize(optarg, ", ");
                        if (tokens.size() != 2)
                            Log(EError, "Invalid tone mapper parameter!");
                        for (int i=0; i<2; ++i) {
                            tonemapper[i] = (Float) std::strtod(tokens[i].c_str(), &end_ptr);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse tone mapper parameters!");
                        }
                    }
                    break;

                case 's': {
                        std::vector<std::string> tokens = tokenize(optarg, ", ");
                        if (tokens.size() != 2)
                            Log(EError, "Invalid resize parameter!");
                        for (int i=0; i<2; ++i) {
                            resize[i] = (int) std::strtol(tokens[i].c_str(), &end_ptr, 10);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse integer in resize parameter!");
                        }
                    }
                    break;


                case 'r': {
                        std::vector<std::string> tokens = tokenize(optarg, ", ");
                        if (tokens.size() != 5)
                            Log(EError, "Invalid rectangle parameter!");
                        Rect r;
                        for (int i=0; i<5; ++i) {
                            r.r[i] = (int) std::strtol(tokens[i].c_str(), &end_ptr, 10);
                            if (*end_ptr != '\0')
                                Log(EError, "Cannot parse integer in rectangle parameter!");
                        }
                        rects.push_back(r);
                    }
                    break;

                case 'o':
                    outputFilename = optarg;
                    break;

                case 't':
                    runParallel = true;
                    break;

                case 'M':
                    merge = true;
                    break;

                case 'C':
                    maxNumBitmapsToMerge = strtol(optarg, &end_ptr, 10);
                    if (*end_ptr != '\0')
                        SLog(EError, "Could not parse the merge Cutoff number!");
                    break;

                case 'R':
                    robustFraction = (Float) strtod(optarg, &end_ptr);
                    if (*end_ptr != '\0')
                        SLog(EError, "Could not parse the Robust fraction value!");
                    break;

                case 'n':
                    noLDR = true;
                    break;
            }
        }

        if (bloomFov != 0 && (bloomFov <= 0 || bloomFov >= 180))
            Log(EError, "Bloom field of view value must be between 0 and 180!");

        if (runParallel) {
            if (outputFilename != "" || temporalCoherence || merge) {
                Log(EWarn, "Requested multithreaded tonemapping along with incompatible options, disabling threading..");
                runParallel = false;
            } else {
                Log(EInfo, "Performing multithreaded tonemapping ..");
            }
        }

        if (optind == argc) {
            help();
            return 0;
        }

        if (pixelFormat == Bitmap::ERGBA && format == Bitmap::EJPEG)
            Log(EError, "JPEG images do not support an alpha channel!");

        if (resize[0] != -1) {
            /* A resampling operation was requested; use a Lanczos Sinc reconstruction filter by default */
            rfilter = static_cast<ReconstructionFilter *> (PluginManager::getInstance()->
                    createObject(MTS_CLASS(ReconstructionFilter), Properties(rfilterName)));
            rfilter->configure();
        }

        if (runParallel) {
            ref<Logger> logger = Thread::getThread()->getLogger();
            std::vector<std::string> messages;

            #if defined(MTS_OPENMP)
                #pragma omp parallel for schedule(static)
            #endif
            for (int i=optind; i<argc; ++i) {
                Thread *thread = Thread::getThread();
                if (!thread) {
                    thread = Thread::registerUnmanagedThread("omp");
                    thread->setLogger(logger);
                }
                try {
                    fs::path inputFile = fileResolver->resolve(argv[i]);
                    Log(EInfo, "Loading image \"%s\" ..", inputFile.string().c_str());
                    ref<FileStream> is = new FileStream(inputFile, FileStream::EReadOnly);
                    ref<Bitmap> input = new Bitmap(Bitmap::EAuto, is);

                    if (crop[2] != -1 && crop[3] != -1)
                        input = input->crop(Point2i(crop[0], crop[1]), Vector2i(crop[2], crop[3]));

                    if (bloomFov != 0) {
                        int maxDim = std::max(input->getWidth(), input->getHeight());
                        if (maxDim % 2 == 0)
                            ++maxDim;

                        ref<Bitmap> bloomFilter = computeBloomFilter(maxDim, bloomFov);

                        if (input->getComponentFormat() != Bitmap::EFloat)
                            input = input->convert(input->getPixelFormat(), Bitmap::EFloat);

                        Log(EInfo, "Convolving image with bloom filter ..");
                        input->convolve(bloomFilter);
                    }

                    if (resize[0] != -1)
                        input = input->resample(rfilter, ReconstructionFilter::EClamp,
                                ReconstructionFilter::EClamp, Vector2i(resize[0], resize[1]));

                    if (cbal[0] != 1 || cbal[1] != 1 || cbal[2] != 1)
                        input->colorBalance(cbal[0], cbal[1], cbal[2]);

                    if (tonemapper[0] != -1) {
                        Float logAvgLuminance = 0, maxLuminance = 0;
                        input->tonemapReinhard(logAvgLuminance, maxLuminance, tonemapper[0], tonemapper[1]);
                        Log(EInfo, "Tonemapper reports: log-average luminance = %f, max. luminance = %f",
                            logAvgLuminance, maxLuminance);
                    }

                    ref<Bitmap> output = input->convert(pixelFormat, Bitmap::EUInt8, gamma, multiplier);

                    for (size_t i=0; i<rects.size(); ++i) {
                        int *r = rects[i].r;
                        output->drawRect(Point2i(r[0], r[1]), Vector2i(r[2], r[3]), Spectrum(r[4]/255.0f));
                    }

                    fs::path outputFile = inputFile;
                    if (format == Bitmap::EPNG)
                        outputFile.replace_extension(".png");
                    else if (format == Bitmap::EJPEG)
                        outputFile.replace_extension(".jpg");
                    else
                        Log(EError, "Unknown target format!");

                    Log(EInfo, "Writing tonemapped image to \"%s\" ..", outputFile.string().c_str());

                    ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
                    output->write(format, os);
                } catch (const std::exception &e) {
                    #pragma omp critical
                    messages.push_back(e.what());
                }
            }
            if (!messages.empty()) {
                Log(EWarn, "The tonemapping worker threads encountered several issues:");
                for (size_t i=0; i<messages.size(); ++i)
                    Log(EWarn, "Exception %i: %s", (int) i, messages[i].c_str());
            }
        } else if (merge) {
            if (outputFilename == "") {
                Log(EError, "When merging images, an explicit output file is required!");
            }
            ref<Bitmap> bloomFilter;

            size_t n = 0;
            fs::path inputFile;
            ref<FileStream> is;
            ref<Bitmap> input;
            Bitmap::EComponentFormat originalComponentFormat;
            // Load first image that we can as base (ignore errors)
            do {
                try {
                    inputFile = fileResolver->resolve(argv[optind]);
                    Log(EInfo, "Loading image \"%s\" ..", inputFile.string().c_str());
                    is = new FileStream(inputFile, FileStream::EReadOnly);
                    input = new Bitmap(Bitmap::EAuto, is);
                } catch (const std::exception &e) {
                    Log(EWarn, "Problem loading file \"%s\".", inputFile.string().c_str());
                    Log(EWarn, "Error was: %s.", e.what());
                    optind++;
                    if (optind >= argc) {
                        Log(EError, "Could not load any file!");
                        exit(1);
                    }
                    continue;
                }
                n++;
                originalComponentFormat = input->getComponentFormat();
            } while (n == 0);
            // make sure we don't lose precision in accumulation -> convert to doubles
            input = input->convert(input->getPixelFormat(), Bitmap::EFloat64);


            
            if (robustFraction == 0) {
                /* Nothing special requested, can just accumulate */
                for (int i=optind+1; i<argc; ++i) {
                    inputFile = fileResolver->resolve(argv[i]);
                    Log(EInfo, "Loading image \"%s\" ..", inputFile.string().c_str());
                    ref<Bitmap> thisInput;
                    try {
                        is = new FileStream(inputFile, FileStream::EReadOnly);
                        thisInput = new Bitmap(Bitmap::EAuto, is);
                        thisInput = thisInput->convert(thisInput->getPixelFormat(), Bitmap::EFloat64);
                    } catch (const std::exception &e) {
                        Log(EWarn, "Problem loading file \"%s\".", inputFile.string().c_str());
                        Log(EWarn, "Error was: %s.", e.what());
                        continue;
                    }
                    input->accumulate(thisInput.get());
                    n++;
                    if (maxNumBitmapsToMerge > 0 && n >= maxNumBitmapsToMerge) {
                        Log(EInfo, "Stopping because cutoff number of bitmaps (%ld) was reached", maxNumBitmapsToMerge);
                        break;
                    }
                }
                input->scale(1.0/n);
            } else {
                /* Robust merging requested */

                if (maxNumBitmapsToMerge > 0 && (size_t)argc < maxNumBitmapsToMerge) {
                    Log(EError, "Requested merge cutoff of %zu bitmaps, but only %d input files provided. Bailing out!",
                            maxNumBitmapsToMerge, argc);
                }

                /* Memory is cheap: get everything in-mem */
                std::vector<ref<Bitmap> > bitmaps;
                for (int i=optind; i<argc; ++i) {
                    inputFile = fileResolver->resolve(argv[i]);
                    Log(EInfo, "Loading image \"%s\" ..", inputFile.string().c_str());
                    try {
                        is = new FileStream(inputFile, FileStream::EReadOnly);
                        ref<Bitmap> thisInput = new Bitmap(Bitmap::EAuto, is);
                        thisInput = thisInput->convert(thisInput->getPixelFormat(), Bitmap::EFloat64);
                        bitmaps.push_back(thisInput);
                    } catch (const std::exception &e) {
                        Log(EWarn, "Problem loading file \"%s\".", inputFile.string().c_str());
                        Log(EWarn, "Error was: %s.", e.what());
                        continue;
                    }
                    if (maxNumBitmapsToMerge > 0 && bitmaps.size() >= maxNumBitmapsToMerge) {
                        Log(EInfo, "Stopped loading input files because cutoff "
                                "num files reached (%zu)", maxNumBitmapsToMerge);
                        break;
                    }
                }
                size_t numBitmaps = bitmaps.size();
                size_t numToDrop;

                /* We have the least amount of bias when our pixel 
                 * distributions are as symmetrical as possible (as close 
                 * to a gaussian as possible), according to CLT, this 
                 * happens in the limit of N->infty. If we get 
                 * robustFraction == -1, then we automatically merge all 
                 * bitmaps into three sets from which we only use the 
                 * median value. We can hack that in here by pre-merging 
                 * the bitmaps into 3 merged ones. */
                if (robustFraction == -1 || robustFraction > 1) {
                    if (robustFraction == -1)
                        robustFraction = 3; // legacy
                    size_t numBins = robustFraction + 0.5;
                    numToDrop = 1;
                    if (numBins > numBitmaps) {
                        Log(EWarn, "Requested %d bins, but only got %d bitmaps -> using %d bins",
                                numBins, numBitmaps, numBitmaps);
                    } else {
                        Log(EInfo, "Merging all %d bitmaps into %d bins and droping extremal bins.",
                                numBitmaps, numBins);
                        std::vector<ref<Bitmap> > mergedBitmaps(numBins);
                        for (size_t i = 0; i < numBins; i++) {
                            mergedBitmaps[i] = new Bitmap(
                                    input->getPixelFormat(), Bitmap::EFloat64,
                                    input->getSize(), input->getChannelCount());
                            size_t startIdx = std::round(((double) i  ) * numBitmaps / numBins);
                            size_t endIdx   = std::round(((double) i+1) * numBitmaps / numBins);
                            mergeBitmaps(bitmaps, mergedBitmaps[i], 0, startIdx, endIdx);
                        }
                        /* Just pretend like we were given our merged bitmaps */
                        bitmaps = mergedBitmaps;
                        numBitmaps = numBins;
                    }
                } else {
                    /* 'Regular' dropping of samples out of all input bitmaps
                     * [OLD, LEGACY, More bias than 'numBins' approach] */
                    numToDrop = std::max((size_t)1,(size_t)(0.5 + numBitmaps * robustFraction));
                }

                if (numBitmaps < 3) {
                    Log(EWarn, "Requested robust merging, but could not "
                            "load at least 3 bitmaps! (loaded %d bitmaps)! "
                            "WILL NOT DROP ANY OUTLIERS!",
                            numBitmaps);
                    numToDrop = 0;
                }

                Assert(numBitmaps - 2*numToDrop > 0);

                Log(EInfo, "Requested robust merging: dropping %d lowest and %d highest samples out of %d",
                        numToDrop, numToDrop, numBitmaps);

                mergeBitmaps(bitmaps, input, numToDrop); // original input gets overwritten 
                size_t nEntries =
                    (size_t) input->getSize().x *
                    (size_t) input->getSize().y *
                    input->getChannelCount();

                std::vector<double> samples(numBitmaps);
                for (size_t i=0; i<nEntries; ++i) {
                    for (size_t b=0; b<numBitmaps; b++)
                        samples[b] = bitmaps[b]->getFloat64Data()[i];
                    std::sort(samples.begin(), samples.end()); 
                    /* note: full sort is somewhat overkill for cutting the 
                     * top and bottom part, but it *does* make the following 
                     * addition numerically robust! */
                    double avg = 0;
                    for (size_t j = numToDrop; j < numBitmaps-numToDrop; j++) {
                        avg += samples[j];
                    }
                    avg = avg / (numBitmaps - 2*numToDrop);
                }
            }


            ref<Bitmap> outputExr = input->convert(pixelFormat, originalComponentFormat);
            fs::path exrOutputFile = outputFilename;
            exrOutputFile.replace_extension(".exr");
            Log(EInfo, "Writing merged image to \"%s\" ..", exrOutputFile.string().c_str());
            ref<FileStream> exros = new FileStream(exrOutputFile, FileStream::ETruncReadWrite);
            outputExr->write(Bitmap::EOpenEXR, exros);


            if (!noLDR) {
                if (crop[2] != -1 && crop[3] != -1)
                    input = input->crop(Point2i(crop[0], crop[1]), Vector2i(crop[2], crop[3]));

                if (bloomFov != 0) {
                    int maxDim = std::max(input->getWidth(), input->getHeight());
                    if (maxDim % 2 == 0)
                        ++maxDim;

                    if (bloomFilter == NULL || bloomFilter->getWidth() != maxDim)
                        bloomFilter = computeBloomFilter(maxDim, bloomFov);

                    if (input->getComponentFormat() != Bitmap::EFloat)
                        input = input->convert(input->getPixelFormat(), Bitmap::EFloat);

                    Log(EInfo, "Convolving image with bloom filter ..");
                    input->convolve(bloomFilter);
                }

                if (resize[0] != -1)
                    input = input->resample(rfilter, ReconstructionFilter::EClamp,
                        ReconstructionFilter::EClamp, Vector2i(resize[0], resize[1]));

                if (cbal[0] != 1 || cbal[1] != 1 || cbal[2] != 1)
                    input->colorBalance(cbal[0], cbal[1], cbal[2]);

                if (tonemapper[0] != -1) {
                    input->tonemapReinhard(logAvgLuminance, maxLuminance, tonemapper[0], tonemapper[1]);
                    Log(EInfo, "Tonemapper reports: log-average luminance = %f, max. luminance = %f",
                        logAvgLuminance, maxLuminance);
                    if (!temporalCoherence) {
                        logAvgLuminance = 0;
                        maxLuminance = 0;
                    }
                }

                ref<Bitmap> output = input->convert(pixelFormat, Bitmap::EUInt8, gamma, multiplier);

                for (size_t i=0; i<rects.size(); ++i) {
                    int *r = rects[i].r;
                    output->drawRect(Point2i(r[0], r[1]), Vector2i(r[2], r[3]), Spectrum(r[4]/255.0f));
                }

                fs::path outputFile = outputFilename;
                Log(EInfo, "Writing tonemapped image to \"%s\" ..", outputFile.string().c_str());

                ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
                output->write(format, os);
            }
        } else {
            ref<Bitmap> bloomFilter;
            for (int i=optind; i<argc; ++i) {
                fs::path inputFile = fileResolver->resolve(argv[i]);
                Log(EInfo, "Loading image \"%s\" ..", inputFile.string().c_str());
                ref<FileStream> is = new FileStream(inputFile, FileStream::EReadOnly);
                ref<Bitmap> input = new Bitmap(Bitmap::EAuto, is);

                if (crop[2] != -1 && crop[3] != -1)
                    input = input->crop(Point2i(crop[0], crop[1]), Vector2i(crop[2], crop[3]));

                if (bloomFov != 0) {
                    int maxDim = std::max(input->getWidth(), input->getHeight());
                    if (maxDim % 2 == 0)
                        ++maxDim;

                    if (bloomFilter == NULL || bloomFilter->getWidth() != maxDim)
                        bloomFilter = computeBloomFilter(maxDim, bloomFov);

                    if (input->getComponentFormat() != Bitmap::EFloat)
                        input = input->convert(input->getPixelFormat(), Bitmap::EFloat);

                    Log(EInfo, "Convolving image with bloom filter ..");
                    input->convolve(bloomFilter);
                }

                if (resize[0] != -1)
                    input = input->resample(rfilter, ReconstructionFilter::EClamp,
                        ReconstructionFilter::EClamp, Vector2i(resize[0], resize[1]));

                if (cbal[0] != 1 || cbal[1] != 1 || cbal[2] != 1)
                    input->colorBalance(cbal[0], cbal[1], cbal[2]);

                if (tonemapper[0] != -1) {
                    input->tonemapReinhard(logAvgLuminance, maxLuminance, tonemapper[0], tonemapper[1]);
                    Log(EInfo, "Tonemapper reports: log-average luminance = %f, max. luminance = %f",
                        logAvgLuminance, maxLuminance);
                    if (!temporalCoherence) {
                        logAvgLuminance = 0;
                        maxLuminance = 0;
                    }
                }

                ref<Bitmap> output = input->convert(pixelFormat, Bitmap::EUInt8, gamma, multiplier);

                for (size_t i=0; i<rects.size(); ++i) {
                    int *r = rects[i].r;
                    output->drawRect(Point2i(r[0], r[1]), Vector2i(r[2], r[3]), Spectrum(r[4]/255.0f));
                }

                fs::path outputFile = inputFile;
                if (outputFilename == "") {
                    if (format == Bitmap::EPNG)
                        outputFile.replace_extension(".png");
                    else if (format == Bitmap::EJPEG)
                        outputFile.replace_extension(".jpg");
                    else
                        Log(EError, "Unknown target format!");
                } else {
                    outputFile = outputFilename;
                }

                Log(EInfo, "Writing tonemapped image to \"%s\" ..", outputFile.string().c_str());

                ref<FileStream> os = new FileStream(outputFile, FileStream::ETruncReadWrite);
                output->write(format, os);
            }
        }
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Tonemap, "Command line batch tonemapper")
MTS_NAMESPACE_END
