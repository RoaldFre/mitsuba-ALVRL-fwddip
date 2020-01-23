#include <mitsuba/core/plugin.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/truncnorm.h>

MTS_NAMESPACE_BEGIN

class TruncNorm : public Utility {
public:

    int run(int argc, char **argv) {
        if (argc != 6)
            Log(EError, "Expected 5 arguments, but got %d! Need:\n"
                    "<mean>\n"
                    "<stddev>\n"
                    "<lowerBound>\n"
                    "<upperBound>\n"
                    "<numSamples>",
                    argc - 1
                    );

        char *endPtr = NULL;

        Float mean = strtod(argv[1], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read mean");

        Float stddev = strtod(argv[2], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read stddev");

        Float lowerBound = strtod(argv[3], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read lowerBound");

        Float upperBound = strtod(argv[4], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read upperBound");

        if (lowerBound > upperBound)
            Log(EError, "Incompatible lower and upper bounds!");

        long long numSamples = llround(strtod(argv[5], &endPtr));
        if (*endPtr != '\0')
            Log(EError, "could not read num samples");

        Properties samplerProps("independent");
        Sampler *sampler = static_cast<Sampler *>(
                PluginManager::getInstance()->createObject(
                    MTS_CLASS(Sampler), samplerProps));

        for (long long i = 0; i < numSamples; i++) {
            double z = truncnorm(mean, stddev, lowerBound, upperBound, sampler);
            double pdf = truncnormPdf(mean, stddev, lowerBound, upperBound, z);
            cerr << z << " " << pdf << endl;
        }

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(TruncNorm, "Utility for generating truncated normal values")
MTS_NAMESPACE_END
