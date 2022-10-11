#include <mitsuba/core/plugin.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/dss.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

class DssDirectionAnalyser : public Utility {
public:

    int run(int argc, char **argv) {
        if (argc != 12)
            Log(EError, "Expected 11 arguments, but got %d! Need:\n"
                    "<scene.xml> <numEval> <|R|> <RcosTheta> <nOutCosTheta> <nOutPhi> <dOutCosTheta> <dOutPhi> <dInCosTheta> <dInPhi> <outputFilePrefix>\n"
                    "Geometrical setup:\n"
                    "  n_in = (0,0,1)\n"
                    "  R = |R|*('RsinTheta',0,RcosTheta)\n"
                    "  d_out, d_in and n_out are given as spherical coordinates with theta & phi\n"
                    "Set <dInCosTheta> or <dInPhi> to -1 to sample that dimension\n",
                    argc - 1
                    );

        char *endPtr = NULL;

        char *sceneFile = argv[1];

        long long numEval = llround(strtod(argv[2], &endPtr));
        if (*endPtr != '\0')
            Log(EError, "could not read num samples");

        Float lRl = strtod(argv[3], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read |R|");

        Float RcosTheta = strtod(argv[4], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read RcosTheta");

        Float nOutCosTheta = strtod(argv[5], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read nOutCosTheta");

        Float nOutPhi = strtod(argv[6], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read nOutPhi");

        Float dOutCosTheta = strtod(argv[7], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read dOutCosTheta");

        Float dOutPhi = strtod(argv[8], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read dOutPhi");

        Float dInCosTheta = strtod(argv[9], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read dInCosTheta");

        Float dInPhi = strtod(argv[10], &endPtr);
        if (*endPtr != '\0')
            Log(EError, "could not read dInPhi");

        std::string outputFilePrefix(argv[11]);
        std::string outputFileActual(outputFilePrefix + "_actual");
        std::string outputFileIdeal(outputFilePrefix + "_ideal");


        if (dInCosTheta != -1 && dInPhi != -1)
            Log(EError, "Need at least one of dInCosTheta or dInPhi equal to -1");
        bool fullGridSampling = (dInCosTheta == -1 && dInPhi == -1);

        Properties samplerProps("independent");
        Sampler *sampler = static_cast<Sampler *>(
                PluginManager::getInstance()->createObject(
                    MTS_CLASS(Sampler), samplerProps));

        ref<Scene> scene(loadScene(sceneFile));
        scene->initialize();

        if (scene->getSubsurfaceIntegrators().size() != 1)
            Log(EError, "Expected to find exactly one SubsurfaceIntegrator in the scene file!, but got %d",
                    scene->getSubsurfaceIntegrators().size());
        if (!scene->getSubsurfaceIntegrators()[0]->getClass()->derivesFrom(MTS_CLASS(DirectSamplingSubsurface)))
            Log(EError, "Expected to find a DirectSamplingSubsurface SubsurfaceIntegrator!");
        DirectSamplingSubsurface *dss = static_cast<DirectSamplingSubsurface*>(
                scene->getSubsurfaceIntegrators()[0].get());

        if (scene->getShapes().size() != 1)
            Log(EError, "Expected to find exactly one Shape in the scene file!, but got %d",
                    scene->getShapes().size());
        Shape *shape = scene->getShapes()[0];

#if 0
        if (!shape->hasBSDF())
            Log(EError, "Expected to find a BSDF associated with the shape!");
        BSDF *bsdf = shape->getBSDF();
#endif

        long long numIntSamples = std::max((long long)100000, numEval);

        Intersection its_out, its_in;
        Vector n_in(0,0,1);
        Frame inFrame(n_in);

        Point p_in(0,0,0);
        Point p_out(lRl*math::safe_sqrt(1-math::square(RcosTheta)), 0, lRl*RcosTheta);

        Vector n_out(cylindricalDirection(nOutCosTheta, nOutPhi));
        Vector d_out(cylindricalDirection(dOutCosTheta, dOutPhi));

        its_out.shFrame = Frame(n_out);
        its_out.geoFrame = its_out.shFrame;
        its_out.p = p_out;
        its_out.shape = shape;

        its_in.shFrame = inFrame;
        its_in.geoFrame = inFrame;
        its_in.p = p_in;
        its_in.shape = shape;

        EMeasure bsdfMeasure = EDiscrete;

        void *extraParams = malloc(dss->extraParamsSize());
        if (dss->sampleExtraParams(scene, its_out, d_out, its_in, NULL, Spectrum(1.0f), extraParams, sampler).isZero())
            Log(EError, "Could not sample extra parameters!");

        //auto [integral, integralErr] = dss->computeDirectionsIntegral(
        auto thePair = dss->computeDirectionsIntegral(
                scene, Spectrum(1.0f),
                its_out, d_out, false,
                its_in,
                extraParams,
                bsdfMeasure, sampler, true /* absify! */,
                numIntSamples, 1.0);
        Float integral = thePair.first;
        Float integralErr = thePair.second;


        cout << integral << "+-"<<integralErr<<endl;
        if (integral == 0)
            return 1;

        std::fstream fActual;
        std::fstream fIdeal;
        fActual.open(outputFileActual, std::fstream::out);
        fIdeal.open(outputFileIdeal, std::fstream::out);
        int numEvalPerDim = (fullGridSampling ? int(sqrt(numEval)+0.5) : numEval);
        int numEvalCos = dInCosTheta == -1 ? numEvalPerDim : 1;
        int numEvalPhi =      dInPhi == -1 ? numEvalPerDim : 1;
        for (int i = 0; i < numEvalCos; i++) {
            Float x = (i+0.5) / numEvalPerDim;
            Float d_in_cosTheta = (dInCosTheta == -1 ? x : dInCosTheta);

            for (int j = 0; j < numEvalPhi; j++) {
                Float y = (j+0.5) / numEvalPerDim;
                Float d_in_phi      = (dInPhi == -1 ? TWO_PI * y : dInPhi);

                Float d_in_sinTheta = math::safe_sqrt(1 - math::square(d_in_cosTheta));
                Float d_in_cosPhi, d_in_sinPhi;
                math::sincos(d_in_phi, &d_in_sinPhi, &d_in_cosPhi);
                Vector d_in = d_in_sinTheta*d_in_cosPhi * inFrame.s
                            + d_in_sinTheta*d_in_sinPhi * inFrame.t
                            - d_in_cosTheta * inFrame.n;

                /* Evaluate the sampling pdf */
                Vector rec_wi = -d_in; // hardcoded for null brdf
                its_in.wi = inFrame.toLocal(rec_wi);
                //EMeasure bsdfMeasure = EDiscrete;
                Float thePdf = dss->pdfBssrdfDirection(
                        scene, its_out, d_out, its_in, d_in,
                        extraParams, Spectrum(1.0f));

                /* Evaluate BSSRDF */
                Spectrum bssrdfVal = dss->bssrdf(scene, p_in, d_in, n_in, p_out,
                        d_out, n_out, extraParams);
                Float theBssrdf = bssrdfVal.averageNonNan();
                Assert(dot(d_in, n_in) <= 0);
                Float cosTheta = -dot(d_in, n_in);

                fActual << thePdf;
                fIdeal << theBssrdf*cosTheta/integral;
                if (j != numEvalPhi-1) {
                    fActual<< " ";
                    fIdeal << " ";
                }
            }
            fActual << endl;
            fIdeal << endl;
        }
        fActual.close();
        fIdeal.close();

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(DssDirectionAnalyser, "Utility for analysing the direction sampling of a Direct Sampling Subsurface")
MTS_NAMESPACE_END
