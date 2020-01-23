#pragma once
#if !defined(__MITSUBA_MEDIUM_FWDSCAT_H_)
#define __MITSUBA_MEDIUM_FWDSCAT_H_

#include <mitsuba/render/dipoleModel.h>

MTS_NAMESPACE_BEGIN

/* This ensures that the pdf calculation from a set of ('numerically
 * rounded') directions doesn't become badly conditioned when compared to
 * the calculation of the pdf during sampling. */
// WARNING: Single precision has not been tested! Value is just a guess...
#ifdef SINGLE_PRECISION
# define MTS_FWDSCAT_DIRECTION_MIN_MU 1e-3
#else
# define MTS_FWDSCAT_DIRECTION_MIN_MU 1e-4
#endif




class MTS_EXPORT_RENDER FwdScat final : public DipoleModel {
public:
    FwdScat(Float sigS, Float sigA, Float g, Float eta,
            int channel, const Properties &props) :
                DipoleModel(sigS, sigA, g, eta, channel, props),
                p(0.5 * sigS * (1 - g)) { 
        if (g < 0 || g >= 1) {
            Log(EError, "Valid values for g are in [0,1). "
                    "Sensible values are close to 1.");
        }
    }

    FwdScat(Stream *stream, InstanceManager *manager) :
            DipoleModel(stream, manager),
            p(stream->readFloat()) { }

    void serialize(Stream *stream, InstanceManager *manager) const {
        DipoleModel::serialize(stream, manager);
        stream->writeFloat(p);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "FwdScat[sigma_s="<<m_sigS
                <<", sigma_a="<<m_sigA
                <<", p="<<p
                <<", eta="<<m_eta
                <<"]";
        return oss.str();
    }

    size_t extraParamsSize() const {
        return sizeof(Float); // only param is length
    }

    virtual Float getRequestedDirectionalCosineHemisphereWeight() const {
        return 0.1;
    }

    virtual Float evalMonopole(const Monopole &m) const;

    Float evalMonopole(Vector u0, Vector uL, Vector R, Float length) const;

    Float evalPlaneSource(Vector u0, Vector uL,
            Vector n, Float Rz, Float length) const;


    virtual Float sampleExtraParamsMonopole(
            const Monopole &m, void *extraParams, Sampler *sampler) const;

    virtual Float pdfExtraParamsMonopole(const Monopole &m) const;

    virtual Float sampleDirectionMonopole(Monopole &m, Sampler *sampler) const;

    virtual Float pdfDirectionMonopole(const Monopole &m) const;

    virtual Float realSourceWeight_margOverParamsAndDirections(
            const Monopole &real, const Monopole &virt) const {
        return 0.5; // TODO
    }

    virtual Float realSourceWeight_margOverParams(
            const Monopole &real, const Monopole &virt) const {
        return 0.5; // TODO
    }

    virtual Float realSourceWeight_margOverDirections(
            const Monopole &real, const Monopole &virt) const {
        return 0.5; // TODO
    }



protected:
    void calcValues(double length, double &C, double &D, double &E, double &F,
            double *Z=NULL) const;
    double absorptionAndNormalizationConstant(Float theLength) const;

    /// Returns the pdf
    Float sampleLengthShortLimit(
            Vector R, const Vector *u0, Vector uL, Float &s, Sampler *sampler) const;
    Float pdfLengthShortLimit(
            Vector R, const Vector *u0, Vector uL, Float s) const;
    void implLengthShortLimit(
            Vector R, const Vector *u0, Vector uL, Float &s, Sampler *sampler, Float *pdf) const;
    void implLengthShortLimitKnownU0(
            Vector R, Vector u0, Vector uL, Float &s, Sampler *sampler, Float *pdf) const;
    void implLengthShortLimitMargOverU0(
            Vector R, Vector uL, Float &s, Sampler *sampler, Float *pdf) const;
    void implLengthShortLimitMargOverU0_internal(
            Vector R, Vector uL, Float &s, Sampler *sampler, Float *pdf, Float safetyFac) const;

    /// Returns the pdf
    Float sampleLengthLongLimit(
            Vector R, Vector uL, Float &s, Sampler *sampler) const;
    Float pdfLengthLongLimit(
            Vector R, Vector uL, Float s) const;

    /// Returns the pdf
    Float sampleLengthAbsorption(
            Float &s, Sampler *sampler) const;
    Float pdfLengthAbsorption(
            Float s) const;


    /// Returns the pdf
    Float sampleDirectionBoundaryAwareMonopole_BRDF(
            Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, Sampler *sampler) const;
    Float pdfDirectionBoundaryAwareMonopole_BRDF(
            const Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s) const;
    void implDirectionBoundaryAwareMonopole_BRDF(
            Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, Sampler *sampler, Float *pdf) const;

    /// Returns the pdf
    Float sampleDirectionBoundaryAwareMonopole(
            Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, bool useEffectiveBRDF, Sampler *sampler) const;
    Float pdfDirectionBoundaryAwareMonopole(
            const Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, bool useEffectiveBRDF) const;

    /// Returns the pdf
    Float sampleDirectionBoundaryAwareMonopole_orig(
            Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, Sampler *sampler) const;
    Float pdfDirectionBoundaryAwareMonopole_orig(
            const Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s) const;

    /// Returns the pdf
    Float sampleDirectionBoundaryAwareMonopole_bis(
            Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, Sampler *sampler) const;
    Float pdfDirectionBoundaryAwareMonopole_bis(
            const Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s) const;
    void implDirectionBoundaryAwareMonopole_bis(
            Vector &u0, const Vector &n0, const Vector &uL, const Vector &R,
            Float s, Sampler *sampler, Float *pdf) const;

    const Float p;  /// Inverse length scale of forward scattering model

    MTS_DECLARE_CLASS();
};

MTS_NAMESPACE_END


#include "fwdscat_impl.h"


#endif /* __MITSUBA_MEDIUM_FWDSCAT_H_ */
