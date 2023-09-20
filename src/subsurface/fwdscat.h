#pragma once
#if !defined(__MITSUBA_MEDIUM_FWDSCAT_H_)
#define __MITSUBA_MEDIUM_FWDSCAT_H_

#include <mitsuba/render/dipoleModel.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT FwdScat final : public DipoleModel {
public:
    FwdScat(Float sigS, Float sigA, Float g, Float eta, const Properties &props) :
                DipoleModel(sigS, sigA, g, eta, props),
                p(0.5 * sigS * (1 - g)) { 
        if (g < 0 || g >= 1) {
            Log(EError, "Valid values for g are in [0,1). "
                    "Sensible values are close to 1.");
        }

        m_direction_min_mu = props.getFloat("direction_min_mu", -1);
        if (m_direction_min_mu == -1) {
#ifdef SINGLE_PRECISION
            // WARNING: Single precision has not been tested! Value is just a guess...
            m_direction_min_mu = 1e-3;
#else
            m_direction_min_mu = 1e-4;
#endif
        }

        m_exactMargOverDirectionWeight = props.getBoolean("exactMargOverDirectionWeight", true);

        m_debug_override_ps = props.getFloat("debug_override_ps", -1);
        m_debug_uniform_ps_min = props.getFloat("debug_uniform_ps_min", -1);
        m_debug_uniform_ps_max = props.getFloat("debug_uniform_ps_max", -1);
        m_debug_requested_hemi_weight = props.getFloat("debug_requested_hemi_weight", -1);

        m_useRayDirSurfSampler = props.getBoolean("rayDirSurfSampler", true); /* TODO Disable by default because very high ray cost for not much benefit unless very sparse media (in which case subsurf model isn't the best choice anyway) */
        m_useBidirectionalRayDirSurfSampler = props.getBoolean("bidirectionalRayDirSurfSampler", true);
        m_rayDirSurfSamplerStrategy = props.getString("rayDirSurfSamplerStrategy", "pt"); // TODO change default to ss or om?

        Log(EInfo, "FWDSCAT DEBUG: ps[%f|%f..%f], hemi %f", 
                m_debug_override_ps,
                m_debug_uniform_ps_min,
                m_debug_uniform_ps_max,
                m_debug_requested_hemi_weight);

        Log(EInfo, "Loaded FwdScat medium with p = %f, g %f; %s)",
                p, g, toString().c_str());
    }

    FwdScat(Stream *stream, InstanceManager *manager) :
            DipoleModel(stream, manager),
            p(stream->readFloat()) {
        m_direction_min_mu = stream->readFloat();
        m_exactMargOverDirectionWeight = stream->readBool();
        m_debug_override_ps = stream->readFloat();
        m_debug_uniform_ps_min = stream->readFloat();
        m_debug_uniform_ps_max = stream->readFloat();
        m_debug_requested_hemi_weight = stream->readFloat();
        m_useRayDirSurfSampler = stream->readBool();
        m_useBidirectionalRayDirSurfSampler = stream->readBool();
        m_rayDirSurfSamplerStrategy = stream->readString();
        Log(EInfo, "unser FwdScat: %s", toString().c_str());
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        DipoleModel::serialize(stream, manager);
        stream->writeFloat(p);
        stream->writeFloat(m_direction_min_mu);
        stream->writeBool(m_exactMargOverDirectionWeight);
        stream->writeFloat(m_debug_override_ps);
        stream->writeFloat(m_debug_uniform_ps_min);
        stream->writeFloat(m_debug_uniform_ps_max);
        stream->writeFloat(m_debug_requested_hemi_weight);
        stream->writeBool(m_useRayDirSurfSampler);
        stream->writeBool(m_useBidirectionalRayDirSurfSampler);
        stream->writeString(m_rayDirSurfSamplerStrategy);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "FwdScat[sigma_s="<<m_sigS
                <<", sigma_a="<<m_sigA
                <<", p="<<p
                <<", minMu="<<m_direction_min_mu
                <<", eta="<<m_eta
                <<", exactMargOverDirectionWeight="<<m_exactMargOverDirectionWeight
                <<", rayDirSurfSampler="<<m_useRayDirSurfSampler
                <<", bidirRayDirSurfSamp="<<m_useBidirectionalRayDirSurfSampler
                <<", rayDirSurfSamplerStrategy="<<m_rayDirSurfSamplerStrategy
                <<"]";
        return oss.str();
    }

    size_t extraParamsSize() const {
        return sizeof(Float); // only param is length
    }

    virtual Float getRequestedDirectionalCosineHemisphereWeight() const {
        if (m_debug_requested_hemi_weight < 0)
            return 0.1;
        Log(EInfo, "RETURNING DEBUG HEMI WEIGHT %f", m_debug_requested_hemi_weight);
        return m_debug_requested_hemi_weight;
    }

    virtual Float evalMonopole(const Monopole &m) const;

    Float evalMonopole(Vector u0, Vector uL, Vector R, Float length) const;

    Float evalIsotropicPointSource(Vector uL, Vector R, Float length) const;

    Float evalPlaneSource(Vector u0, Vector uL,
            Vector n, Float Rz, Float length) const;


    virtual Float sampleExtraParamsMonopole(
            const Monopole &m, void *extraParams, Sampler *sampler) const;

    virtual Float pdfExtraParamsMonopole(const Monopole &m) const;

    virtual Float sampleDirectionMonopole(Monopole &m, Sampler *sampler) const;

    virtual Float pdfDirectionMonopole(const Monopole &m) const;

    virtual Float monopoleWeight_margOverParamsAndDirections(
            const Monopole &m) const {
        return 1; // TODO
    }

    virtual Float monopoleWeight_margOverParams(
            const Monopole &m) const {
        return 1; // TODO
    }

    virtual Float monopoleWeight_margOverDirections(
            const Monopole &m) const {
        if (!m_exactMargOverDirectionWeight)
            return 1;

        Float s = *(static_cast<const Float*>(m.extraParams));
        /* Warning: this does not take into account any restriction to proper
         * *incoming* directions wrt the incoming normal! */
        return evalIsotropicPointSource(m.d_out, m.R, s);
    }



protected:
    void calcValues(double length, double &C, double &D, double &E, double &F,
            double *Z=NULL) const;
    double absorptionAndNormalizationConstant(Float theLength) const;

    /// Returns the pdf
    Float sampleLengthShortLimit(
            Vector R, const Vector *u0, Vector uL, Float &s, Sampler *sampler, bool isPlaneSource) const;
    Float pdfLengthShortLimit(
            Vector R, const Vector *u0, Vector uL, Float s, bool isPlaneSource) const;
    void implLengthShortLimit(
            Vector R, const Vector *u0, Vector uL, Float &s, Sampler *sampler, Float *pdf, bool isPlaneSource) const;
    void implLengthShortLimitKnownU0(
            Vector R, Vector u0, Vector uL, Float &s, Sampler *sampler, Float *pdf, bool isPlaneSource) const;
    void implLengthShortLimitMargOverU0(
            Vector R, Vector uL, Float &s, Sampler *sampler, Float *pdf, bool isPlaneSource) const;
    void implLengthShortLimitMargOverU0_internal(
            Vector R, Vector uL, Float &s, Sampler *sampler, Float *pdf, Float safetyFac, bool isPlaneSource) const;
    Float sampleLengthLongLimit_BRDF(Vector R, Vector uL, Float &s, Sampler *sampler) const;
    Float pdfLengthLongLimit_BRDF(Vector R, Vector uL, Float s) const;
    Float implLengthLongLimit_BRDF(Vector R, Vector uL, Float &s, Sampler *sampler) const;

    /// Returns the pdf
    Float sampleLengthLongLimit(
            Vector R, Vector uL, Float &s, Sampler *sampler, bool isPlaneSource) const;
    Float pdfLengthLongLimit(
            Vector R, Vector uL, Float s, bool isPlaneSource) const;

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

    /**
     * Lower limit on the effective 'mu' of the directional part.
     * This ensures that the pdf calculation from a set of ('numerically
     * rounded') directions doesn't become badly conditioned when compared to
     * the calculation of the pdf during sampling.
     */
    float m_direction_min_mu;

    bool m_exactMargOverDirectionWeight;

    Float m_debug_override_ps;
    Float m_debug_uniform_ps_min;
    Float m_debug_uniform_ps_max;
    Float m_debug_requested_hemi_weight;

    MTS_DECLARE_CLASS();

    /* Not really needed at this level, but added here so we can query it
     * easily when wrapping this in a DipoleModel to set up the surface
     * samplers. Quick and dirty, not very clean, needs restructure (TODO). */
    bool m_useRayDirSurfSampler;
    bool m_useBidirectionalRayDirSurfSampler;
    std::string m_rayDirSurfSamplerStrategy;
public:
    bool useRayDirSurfSampler() const {
        return m_useRayDirSurfSampler;
    }
    bool useBidirectionalrayDirSurfSampler() const {
        return m_useBidirectionalRayDirSurfSampler;
    }
    std::string getRayDirSurfSamplerStrategy() const {
        return m_rayDirSurfSamplerStrategy;
    }
};

MTS_NAMESPACE_END


#include "fwdscat_impl.h"


#endif /* __MITSUBA_MEDIUM_FWDSCAT_H_ */
