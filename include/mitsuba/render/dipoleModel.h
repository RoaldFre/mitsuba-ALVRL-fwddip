#if !defined(__DIPOLE_MODEL_H)
#define __DIPOLE_MODEL_H

#include <mitsuba/render/dss.h>
#include "../../../src/medium/materials.h" // XXX ehhh, not the cleanest dependency...

/* Reject incoming directions that come from within the actual geometry
 * (i.e. w.r.t. the actual local normal at the incoming point instead of,
 * for instance, the modified tangent plane normal)? */
#define MTS_DIPOLE_MODEL_REJECT_INCOMING_WRT_TRUE_SURFACE_NORMAL true

MTS_NAMESPACE_BEGIN

enum TangentPlaneMode {
    EUnmodifiedIncoming,
    EUnmodifiedOutgoing,
    EFrisvadEtAl,
    EFrisvadEtAlWithMeanNormal,
};

enum DipoleMode {
    EReal = 1,
    EVirt = 2,
    ERealAndVirt = EReal | EVirt,
};

enum ZvMode {
    EClassicDiffusion, /// As in the original Jensen et al. dipole
    EBetterDipoleZv,   /// As in the better dipole model of d'Eon
    EFrisvadEtAlZv,    /// As in the directional dipole model of Frisvad et al.
};


struct MTS_EXPORT_RENDER Monopole {
    Normal n_in;
    Vector d_in;
    Normal n_out;
    Vector d_out;
    Vector R; // p_out - p_in
    const void *extraParams;
    bool planeSource;

    Monopole() :
            n_in(Vector(0./0.)),  d_in(Vector(0./0.)),
            n_out(Vector(0./0.)), d_out(Vector(0./0.)),
            R(Vector(0./0.)), extraParams(NULL) { }

    Monopole(Normal n_in, Vector d_in, Normal n_out, Vector d_out, Vector R,
                const void *extraParams, bool isPlaneSource) :
            n_in(n_in),   d_in(d_in),
            n_out(n_out), d_out(d_out),
            R(R), extraParams(extraParams),
            planeSource(isPlaneSource) { }

    Monopole(Normal n_in, const Vector *d_in, Normal n_out, Vector d_out, Vector R,
                const void *extraParams, bool isPlaneSource) :
            n_in(n_in),   d_in(d_in ? *d_in : Vector(0./0.)),
            n_out(n_out), d_out(d_out),
            R(R), extraParams(extraParams),
            planeSource(isPlaneSource) { }

    Monopole(Vector d_out, Normal n_out, Vector R, const void *extraParams, bool isPlaneSource) :
            n_in(Vector(0./0.)), d_in(Vector(0./0.)),
            n_out(n_out),        d_out(d_out),
            R(R), extraParams(extraParams),
            planeSource(isPlaneSource) { }


    bool hasDin() const { return d_in.isFinite(); }
    bool hasNin() const { return n_in.isFinite(); }
    bool isPlaneSource() const { return planeSource; }
};


struct MTS_EXPORT_RENDER DipoleConfig {
    bool rejectInternalIncoming; /// reject 'internal' incoming direction wrt the effective normal?
    bool reciprocal;
    TangentPlaneMode tangentMode;
    ZvMode zvMode;
    bool useEffectiveBRDF;
    DipoleMode dipoleMode;
    /**
     * Bit of a hack for index-MISmatched dipole configurations. This makes
     * the dipole refract its directions and change the virtual source
     * displacement (as determined by the Zvmode). This is an 'implicit'
     * boundary condition, as opposed to an explicit 'index matched'
     * (eta = 1) coupling to a proper BSDF as boundary.  */
    Float eta;

    DipoleConfig(const Properties &props);

    DipoleConfig(Stream *stream) {
        rejectInternalIncoming = stream->readBool();
        reciprocal = stream->readBool();
        tangentMode = static_cast<TangentPlaneMode>(stream->readInt());
        zvMode = static_cast<ZvMode>(stream->readInt());
        dipoleMode = static_cast<DipoleMode>(stream->readInt());
        useEffectiveBRDF = stream->readBool();
        eta = stream->readFloat();
    }

    void save(Stream *stream) const {
        stream->writeBool(rejectInternalIncoming);
        stream->writeBool(reciprocal);
        stream->writeInt(tangentMode);
        stream->writeInt(zvMode);
        stream->writeInt(dipoleMode);
        stream->writeBool(useEffectiveBRDF);
        stream->writeFloat(eta);
    }

    /// Should we reject the d_in of m wrt the given effective n_in?
    bool shouldRejectDin(const Monopole &m, const Vector &n_in_effective) const {
        if (rejectInternalIncoming && m.hasDin() && dot(n_in_effective, m.d_in) > 0)
            return true;
        return false;
    }
};




/**
 * Abstract class to construct (monochromatic) dipole models from monopoles.
 */
class MTS_EXPORT_RENDER DipoleModel : public SerializableObject {
public:
    DipoleModel(Float sigS, Float sigA, Float g, Float eta, const Properties &props) :
                m_sigS(sigS), m_sigA(sigA), m_g(g), m_eta(eta) {
        if (g < -1 || g > 1)
            Log(EError, "Invalid value for g: %f, should be in (-1,1)", m_g);
    }

    DipoleModel(Stream *stream, InstanceManager *manager) :
            SerializableObject(stream, manager) {
        Log(EInfo, "unser DipoleModel");
        m_sigS = stream->readFloat();
        m_sigA = stream->readFloat();
        m_g = stream->readFloat();
        m_eta = stream->readFloat();
    }

    void serialize(Stream *stream, InstanceManager *manger) const {
        stream->writeFloat(m_sigS);
        stream->writeFloat(m_sigA);
        stream->writeFloat(m_g);
        stream->writeFloat(m_eta);
    }

    virtual Float evalMonopole(const Monopole &m) const = 0;

    virtual size_t extraParamsSize() const = 0;

    /** Within a DipoleDSS, supplement the our direction sampling with a 
     *  simple directional cosine lobe with this weight */
    virtual Float getRequestedDirectionalCosineHemisphereWeight() const = 0;


    /// Returns the pdf
    virtual Float sampleExtraParamsMonopole(
            const Monopole &m, void *extraParams, Sampler *sampler) const = 0;

    virtual Float pdfExtraParamsMonopole(const Monopole &m) const = 0;

    /// Returns the pdf
    virtual Float sampleDirectionMonopole(Monopole &m, Sampler *sampler) const = 0;

    virtual Float pdfDirectionMonopole(const Monopole &m) const = 0;

    /** Relative weight of real vs virtual source, marginalized over 
     * extraParams and incoming directions */
    virtual Float realSourceWeight_margOverParamsAndDirections(
            const Monopole &real, const Monopole &virt) const = 0;

    /** Relative weight of real vs virtual source, marginalized over 
     * extraParams */
    virtual Float realSourceWeight_margOverParams(
            const Monopole &real, const Monopole &virt) const = 0;

    /** Relative weight of real vs virtual source, marginalized over 
     * incoming directions */
    virtual Float realSourceWeight_margOverDirections(
            const Monopole &real, const Monopole &virt) const = 0;




    /// Returns the pdf
    Float sampleExtraParamsDipole(
            const Normal &n_in,  const Vector *d_in,
            const Normal &n_out, const Vector &d_out,
            const Vector &R, void *extraParams,
            const DipoleConfig &dipConf, Sampler *sampler) const;

    Float pdfExtraParamsDipole(
            const Normal &n_in,  const Vector *d_in,
            const Normal &n_out, const Vector &d_out,
            const Vector &R, const void *extraParams,
            const DipoleConfig &dipConf) const;

    /// Returns the pdf
    Float sampleDirectionDipole(
            const Normal &n_in,  Vector       &d_in,
            const Normal &n_out, const Vector &d_out,
            const Vector &R, const void *extraParams,
            const DipoleConfig &dipConf, Sampler *sampler) const;

    Float pdfDirectionDipole(
            const Normal &n_in,  const Vector &d_in,
            const Normal &n_out, const Vector &d_out,
            const Vector &R, const void *extraParams,
            const DipoleConfig &dipConf) const;

    bool getDipoleParameters(
            const Monopole &m, const DipoleConfig &dipConf,
            Vector &n_in_effective, Float &zv) const;

    Float evalDipole(
            const Normal &n0, const Vector &u0_external,
            const Normal &nL, const Vector &uL_external,
            const Vector &R, const void *extraParams,
            const DipoleConfig &dipConf) const;

    /// Compute the virtual monopole corresponding to the given real monopole.
    bool realToVirt(
            const DipoleConfig &dipConf, const Monopole &r, Monopole &v) const;

    /// Compute the real monopole corresponding to the given virtual monopole.
    bool virtToReal(
            const Vector &R_real, const Vector &n_in_real,
            const DipoleConfig &dipConf, const Monopole &v, Monopole &r) const;

    MTS_DECLARE_CLASS();

protected:
    Float evalDipole_internal(
            const Normal &n0, const Vector &u0_external,
            const Normal &nL, const Vector &uL_external,
            const Vector &R, const void *extraParams,
            const DipoleConfig &dipConf) const;

    bool shouldRejectRealDin(
            const DipoleConfig &dipConf, const Monopole &r) const;

    /* Params to calculate dipole parameters (Note: can also demand a fixed 
     * dipoleConfig and cache e.g. zv...) */
    Float m_sigS; /// sigma_s
    Float m_sigA; /// sigma_a
    Float m_g;
    Float m_eta;
};



/**
 * Construct Dipole Model based DirectSamplingSubsurface materials. 
 * Internally typically has a set of DipoleModels, one for for each 
 * spectral channel.
 *
 * Usage: instantiate me and implement config() to setup your surface 
 * samplers as required.
 */
template <typename DipMod> class MTS_EXPORT_RENDER DipoleDSS final : public DirectSamplingSubsurface {
public:
    /// Setup and register your surface samplers here with registerSampler().
    virtual void configure();

    DipoleDSS(const Properties &props) :
            DirectSamplingSubsurface(props),
            m_dipConf(props) {
        Float eta;
        lookupMaterial(props, m_sigmaS, m_sigmaA, m_g, &eta);
        if (m_sigmaS.min() == m_sigmaS.max()
                && m_sigmaA.min() == m_sigmaA.max()
                && m_g.min() == m_g.max()) {
            // Monochrome!
            m_dipoles.resize(1);
        } else {
            m_dipoles.resize(SPECTRUM_SAMPLES);
        }
        for (size_t i = 0; i < m_dipoles.size(); i++) {
            m_dipoles[i] = new DipMod(m_sigmaS[i], m_sigmaA[i], m_g[i], eta, props);
        }

        m_dirHemiWeight = m_dipoles[0]->getRequestedDirectionalCosineHemisphereWeight();
    }

    DipoleDSS(Stream *stream, InstanceManager *manager) : 
            DirectSamplingSubsurface(stream, manager),
            m_dipConf(stream),
            m_sigmaS(stream),
            m_sigmaA(stream),
            m_g(stream) {
        m_dipoles.resize(stream->readSize());
        for (size_t i = 0; i < m_dipoles.size(); i++) {
            Log(EInfo, "Unserializing DipMod %d",i);
            m_dipoles[i] = new DipMod(stream, manager);
        }
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        DirectSamplingSubsurface::serialize(stream, manager);
        m_dipConf.save(stream);
        m_sigmaS.serialize(stream);
        m_sigmaA.serialize(stream);
        m_g.serialize(stream);
        stream->writeSize(m_dipoles.size());
        for (size_t i = 0; i < m_dipoles.size(); i++) {
            m_dipoles[i]->serialize(stream, manager);
        }
    }

    /** 
     * The extra parameters are stored in memory in an array of packets of the form:
     *    struct ExtraParamPacket {
     *        char individualExtraParamData[individualExtraParamSize()];
     *        bool isValid; // Was the sampling of these parameters successful?
     *    };
     * There is one such packet for each spectral channel.
     */
    inline size_t individualExtraParamSize() const {
        return m_dipoles[0]->extraParamsSize();
    }
    inline size_t extraParamsPacketSize() const {
        return individualExtraParamSize() + sizeof(bool);
    }
    inline size_t extraParamsSize() const {
        return SPECTRUM_SAMPLES * extraParamsPacketSize();
    }
    /// const version
    inline const void *getIndividualExtraParams(const void *extraParams, int i) const {
        return static_cast<const char*>(extraParams) + i*extraParamsPacketSize();
    }
    /// non-const version
    inline void *getIndividualExtraParams(void *extraParams, int i) const {
        return static_cast<char*>(extraParams) + i*extraParamsPacketSize();
    }
    inline bool isValidExtraParam(const void *extraParams, int i) const {
        const char *packet = static_cast<const char*>(extraParams)
                    + i*extraParamsPacketSize();
        return *reinterpret_cast<const bool*>(packet + individualExtraParamSize());
    }
    inline void setValidExtraParam(void *extraParams, int i, bool isValid) const {
        char *packet = static_cast<char*>(extraParams)
                    + i*extraParamsPacketSize();
        *reinterpret_cast<bool*>(packet + individualExtraParamSize()) = isValid;
    }


    inline virtual Spectrum bssrdf(const Scene *scene,
            const Point &p_in,  const Vector &d_in,  const Normal &n_in,
            const Point &p_out, const Vector &d_out, const Normal &n_out,
            const void *extraParams) const {
        Assert(MTS_DSS_ALLOW_INTERNAL_INCOMING_DIR || dot(d_in, n_in) <= 0);
        Assert(m_allowIncomingOutgoingDirections || dot(d_out, n_out) >= 0);
        Spectrum result;
        for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
            // Shortcut for when the given spectra are effectively 1D:
            if (m_dipoles.size() == 1 && i > 0) {
                result[i] = result[0];
                continue;
            }
            if (!isValidExtraParam(extraParams, i)) {
                result[i] = 0;
                continue;
            }

            const void *myExtraParams = getIndividualExtraParams(extraParams, i);
            result[i] = m_dipoles[i]->evalDipole(
                    n_in, d_in, n_out, d_out, p_out - p_in,
                    myExtraParams, m_dipConf);
        }
        return result;
    }

    inline virtual Spectrum sampleExtraParams(const Scene *scene,
            const Intersection &its_out, const Vector &d_out,
            const Intersection &its_in,  const Vector *d_in,
            const Spectrum &throughput, void *extraParams,
            Sampler *sampler) const {
        Vector n_in = its_in.shFrame.n;
        Vector n_out = its_out.shFrame.n;
        Point p_in = its_in.p;
        Point p_out = its_out.p;
        Vector R = p_out - p_in;
        Assert(dot(d_out, n_out) >= -Epsilon);
        Assert(!d_in || dot(*d_in, n_in) <= Epsilon);
        Assert(!m_dipConf.useEffectiveBRDF || n_in == n_out);
        Assert(!m_dipConf.useEffectiveBRDF || p_in == p_out);
        Assert(!throughput.isZero());

        Spectrum pdf;
        if (m_dipoles.size() == 1) {
            Assert(throughput[0] != 0);
            pdf = Spectrum(m_dipoles[0]->sampleExtraParamsDipole(
                    n_in, d_in, n_out, d_out, R, extraParams, m_dipConf, sampler));
            setValidExtraParam(extraParams, 0, pdf[0] != 0);
        } else {
            for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
                if (throughput[i] == 0) {
                    pdf[i] = 0;
                } else {
                    void *myExtraParams = getIndividualExtraParams(extraParams, i);
                    pdf[i] = m_dipoles[i]->sampleExtraParamsDipole(
                            n_in, d_in, n_out, d_out, R, myExtraParams, m_dipConf, sampler);
                }
                setValidExtraParam(extraParams, i, pdf[i] != 0);
            }
        }

#ifdef MTS_DIPOLE_DEBUG
        if (pdf.isZero())
            return pdf;
        Spectrum pdf2 = pdfExtraParams(scene, its_out, d_out, its_in,
                d_in, throughput, extraParams);
        if (!(fabs((pdf - pdf2).average()/(pdf + pdf2).average()) < 1e-3))
            SLog(EWarn, "Inconsistent pdfs: %s vs %s, rel diff %s",
                    pdf.toString().c_str(), pdf2.toString().c_str(),
                    ((pdf-pdf2)/(pdf+pdf2)).toString().c_str());
#endif

        return pdf;
    }

    virtual Spectrum pdfExtraParams(const Scene *scene,
            const Intersection &its_out, const Vector &d_out,
            const Intersection &its_in,  const Vector *d_in,
            const Spectrum &throughput, const void *extraParams) const {
        Vector n_in = its_in.shFrame.n;
        Vector n_out = its_out.shFrame.n;
        Point p_in = its_in.p;
        Point p_out = its_out.p;
        Vector R = p_out - p_in;
        Assert(dot(d_out, n_out) >= -Epsilon);
        Assert(!d_in || dot(*d_in, n_in) <= Epsilon);
        Assert(!m_dipConf.useEffectiveBRDF || n_in == n_out);
        Assert(!m_dipConf.useEffectiveBRDF || p_in == p_out);

        Spectrum pdf;
        if (m_dipoles.size() == 1) {
            Assert(throughput[0] != 0);
            pdf = Spectrum(m_dipoles[0]->pdfExtraParamsDipole(
                    n_in, d_in, n_out, d_out, R, extraParams, m_dipConf));
        } else {
            for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
                if (throughput[i] == 0) {
                    pdf[i] = 0;
                } else {
                    const void *myExtraParams = getIndividualExtraParams(extraParams, i);
                    pdf[i] = m_dipoles[i]->pdfExtraParamsDipole(
                            n_in, d_in, n_out, d_out, R, myExtraParams, m_dipConf);
                }
            }
        }
        return pdf;
    }




    /**
     * MIS weighting of importance sampling the transport and sampling the
     * (cosine) hemisphere */
    virtual Float sampleBssrdfDirection(const Scene *scene,
            const Intersection &its_out, const Vector &d_out,
            Intersection &its_in,        Vector       &d_in,
            const void *extraParams, const Spectrum &throughput,
            Sampler *sampler) const {
        Float pdfHemi, pdfImp;
        if (m_dirHemiWeight == 1 || sampler->next1D() < m_dirHemiWeight) {
            // Default implementation does cosine hemisphere sampling
            pdfHemi = DirectSamplingSubsurface::sampleBssrdfDirection(
                    scene, its_out, d_out, its_in, d_in, extraParams, throughput, sampler);
            if (pdfHemi == 0)
                return 0;
            if (m_dirHemiWeight == 1)
                return pdfHemi;
            pdfImp = pdfDirectionImportance(
                    scene, its_out, d_out, its_in, d_in, extraParams, throughput);
        } else {
            pdfImp = sampleDirectionImportance(
                    scene, its_out, d_out, its_in, d_in, extraParams, throughput, sampler);
            if (pdfImp == 0)
                return 0;
            pdfHemi = DirectSamplingSubsurface::pdfBssrdfDirection(
                    scene, its_out, d_out, its_in, d_in, extraParams, throughput);
        }
        Float pdf = m_dirHemiWeight * pdfHemi + (1 - m_dirHemiWeight) * pdfImp;

#ifdef MTS_DIPOLE_DEBUG
        if (pdf == 0)
            return 0;
        Float pdf2 = pdfBssrdfDirection(scene, its_out, d_out, its_in,
                d_in, extraParams, throughput);
        if (!(math::abs((pdf - pdf2)/(pdf + pdf2)) < 1e-3))
            SLog(EWarn, "Inconsistent pdfs: %e vs %e, rel diff %e",
                    pdf, pdf2, (pdf-pdf2)/(pdf+pdf2));
#endif
        return pdf;

    }

    virtual Float pdfBssrdfDirection(const Scene *scene,
            const Intersection &its_out, const Vector &d_out,
            const Intersection &its_in,  const Vector &d_in,
            const void *extraParams, const Spectrum &throughput) const {
#if !MTS_DSS_ALLOW_INTERNAL_INCOMING_DIR
        Assert(-dot(d_in, its_in.shFrame.n) >= 0);
#endif
        Float pdf = 0;
        pdf += m_dirHemiWeight
                * DirectSamplingSubsurface::pdfBssrdfDirection(
                    scene, its_out, d_out, its_in, d_in, extraParams, throughput);
        if (m_dirHemiWeight == 1)
            return pdf;
        pdf += (1 - m_dirHemiWeight)
                * pdfDirectionImportance(
                    scene, its_out, d_out, its_in, d_in, extraParams, throughput);
        return pdf;
    }


    /* Can only sample one direction, because otherwise we would be
     * branching into degenerate single-spectral-channel transport.
     *
     * We sample a non-zero throughput channel uniformly and the effective
     * pdf becomes averaged pdf over all (non-zero-throughput and
     * non-(-1)-length) channels (essentially like the MIS balance
     * heuristic). Returns the pdf */
    virtual Float sampleDirectionImportance(const Scene *scene,
            const Intersection &its_out, const Vector &d_out,
            Intersection &its_in,        Vector       &d_in,
            const void *extraParams, const Spectrum &throughput,
            Sampler *sampler) const {
        Vector n_in = its_in.shFrame.n;
        Vector n_out = its_out.shFrame.n;
        Point p_in = its_in.p;
        Point p_out = its_out.p;
        Vector R = p_out - p_in;
        Assert(dot(d_out, n_out) >= -Epsilon);
        Assert(!m_dipConf.useEffectiveBRDF || n_in == n_out);
        Assert(!m_dipConf.useEffectiveBRDF || p_in == p_out);
        Assert(!throughput.isZero());

        Float pdf;
        if (m_dipoles.size() == 1) {
            pdf = m_dipoles[0]->sampleDirectionDipole(
                    n_in, d_in, n_out, d_out, R,
                    extraParams, m_dipConf, sampler);
            if (pdf == 0)
                return 0;
        } else {
            /* Only consider nonzero throughput channels */
            int i = throughput.sampleNonZeroChannelUniform(sampler);
            const void *myExtraParams = getIndividualExtraParams(extraParams, i);
            pdf = m_dipoles[i]->sampleDirectionDipole(
                    n_in, d_in, n_out, d_out, R,
                    myExtraParams, m_dipConf, sampler);
            if (pdf == 0)
                return 0;
            int N = 1; // number of nonzero throughput channels
            for (int j = 0; j < SPECTRUM_SAMPLES; j++) {
                if (j == i || throughput[j] == 0)
                    continue;
                const void *myExtraParams = getIndividualExtraParams(extraParams, j);
                pdf += m_dipoles[j]->pdfDirectionDipole(
                        n_in, d_in, n_out, d_out, R,
                        myExtraParams, m_dipConf);
                N++;
            }
            pdf /= N;
        }
        Assert(dot(d_in, n_in) <= Epsilon);
        its_in.wi = its_in.toLocal(d_in);
        return pdf;
    }

    virtual Float pdfDirectionImportance(const Scene *scene,
            const Intersection &its_out, const Vector &d_out,
            const Intersection &its_in,  const Vector &d_in,
            const void *extraParams, const Spectrum &throughput) const {
        Vector n_in = its_in.shFrame.n;
        Vector n_out = its_out.shFrame.n;
        Point p_in = its_in.p;
        Point p_out = its_out.p;
        Vector R = p_out - p_in;
        Assert(dot(d_out, n_out) >= -Epsilon);
        Assert(dot(d_in, n_in) <= Epsilon);
        Assert(!m_dipConf.useEffectiveBRDF || n_in == n_out);
        Assert(!m_dipConf.useEffectiveBRDF || p_in == p_out);
        Assert(!throughput.isZero());

        if (m_dipoles.size() == 1) {
            return m_dipoles[0]->pdfDirectionDipole(
                    n_in, d_in, n_out, d_out, R, extraParams, m_dipConf);
        } else {
            Float pdf = 0;
            int N = 0; // number of nonzero-throughput channels
            for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
                if (throughput[i] == 0)
                    continue;
                const void *myExtraParams = getIndividualExtraParams(extraParams, i);
                pdf += m_dipoles[i]->pdfDirectionDipole(
                        n_in, d_in, n_out, d_out, R,
                        myExtraParams, m_dipConf);
                N++;
            }
            if (N > 0)
                pdf /= N;
            return pdf;
        }
    }

    MTS_DECLARE_CLASS()

protected:
    DipoleConfig m_dipConf;
    ref_vector<DipMod> m_dipoles;
    Float m_dirHemiWeight;

    // Added for convenience, e.g. for use in configure():
    Spectrum m_sigmaS;
    Spectrum m_sigmaA;
    Spectrum m_g;
};








inline DipoleConfig::DipoleConfig(const Properties &props) {
    rejectInternalIncoming = props.getBoolean(
            "rejectInternalIncoming", false);
    reciprocal = props.getBoolean("reciprocal", false);

    std::string zvModeStr = props.getString("zvMode", "diff");
    if (zvModeStr == "diff") {
        zvMode = EClassicDiffusion;
    } else if (zvModeStr == "better") {
        zvMode = EBetterDipoleZv;
    } else if (zvModeStr == "Frisvad") {
        zvMode = EFrisvadEtAlZv;
    } else {
        SLog(EError, "Unknown zvMode: %s", zvModeStr.c_str());
    }

    std::string tangentModeStr = props.getString("tangentMode", "Frisvad");
    if (tangentModeStr == "incoming") {
        tangentMode = EUnmodifiedIncoming;
    } else if (tangentModeStr == "outgoing") {
        tangentMode = EUnmodifiedOutgoing;
    } else if (tangentModeStr == "Frisvad") {
        tangentMode = EFrisvadEtAl;
    } else if (tangentModeStr == "FrisvadMean") {
        tangentMode = EFrisvadEtAlWithMeanNormal;
    } else {
        SLog(EError, "Unknown tangentMode: %s", tangentModeStr.c_str());
    }

    bool onlyReal = props.getBoolean("onlyReal", false);
    bool onlyVirt = props.getBoolean("onlyVirt", false);
    if (onlyReal && onlyVirt)
        SLog(EError, "Requested both *only* real and *only* virtual "
                "contributions!");
    if (onlyReal) {
        dipoleMode = EReal;
        SLog(EInfo, "Requested only contributions from real source");
    } else if (onlyVirt) {
        dipoleMode = EVirt;
        SLog(EInfo, "Requested only contributions from (positive) "
                "virtual source");
    } else {
        dipoleMode = ERealAndVirt;
    }

    useEffectiveBRDF = props.getBoolean("useEffectiveBRDF", false);

    Spectrum _sigmaS, _sigmaA, _g;
    lookupMaterial(props, _sigmaS, _sigmaA, _g, &eta);

    if (eta != 1) {
        SLog(EWarn, "This Dipole Config is configured for "
                "non-index matched, 'implicit' boundaries "
                "internally (eta = %f, presumably coupled to a NULL or "
                "index-matched BSDF). Although this is possible, it is "
                "much advised to keep eta = 1 here and then couple this "
                "'index matched' subsurface scattering model to a "
                "refractive BSDF for more correct, 'explicit' boundary "
                "conditions!", eta);
    }
}




inline Float dEon_C1(const Float n) {
    Float r;
    if (n > 1.0) {
        r = -9.23372 + n * (22.2272 + n * (-20.9292 + n * (10.2291 + n * (-2.54396 + 0.254913 * n))));
    } else {
        r = 0.919317 + n * (-3.4793 + n * (6.75335 + n *  (-7.80989 + n *(4.98554 - 1.36881 * n))));
    }
    return r / 2.0;
}
inline Float dEon_C2(const Float n) {
    Float r = -1641.1 + n * (1213.67 + n * (-568.556 + n * (164.798 + n * (-27.0181 + 1.91826 * n))));
    r += (((135.926 / n) - 656.175) / n + 1376.53) / n;
    return r / 3.0;
}

inline Float dEon_A(const Float eta) {
    return (1 + 3*dEon_C2(eta)) / (1 - 2*dEon_C1(eta));
}

FINLINE bool DipoleModel::getDipoleParameters(
        const Monopole &m, const DipoleConfig &dipConf,
        Vector &n_in_effective, Float &zv) const {
    /* n_in_effective */
    switch (dipConf.tangentMode) {
    case EFrisvadEtAl:
        /* Use the modified tangent plane of the directional dipole model
         * of Frisvad et al */
        if (m.R.length() == 0) {
            n_in_effective = m.n_in;
        } else {
            if (cross(m.n_in, m.R).length() == 0)
                return false;
            n_in_effective = cross(normalize(m.R), normalize(cross(m.n_in, m.R)));
            Assert(dot(n_in_effective, m.n_in) > -Epsilon);
        }
        break;
    case EFrisvadEtAlWithMeanNormal: {
        /* Like the tangent plane of Frisvad et al, but based on an
         * 'average' normal at incoming and outgoing point instead of on
         * the incoming normal. This should immediately give reciprocity as
         * a bonus. */
        Vector sumNormal = m.n_in + m.n_out;
        if (m.R.length() == 0) {
            n_in_effective = m.n_in;
        } else {
            if (cross(sumNormal,m.R).length() == 0)
                return false;
            n_in_effective = cross(normalize(m.R), normalize(cross(sumNormal, m.R)));
        }
        break; }
    case EUnmodifiedIncoming:
        n_in_effective = m.n_in; break;
    case EUnmodifiedOutgoing:
        n_in_effective = m.n_out; break;
    default:
        SLog(EError, "Unknown tangentMode: %d", dipConf.tangentMode);
        return false; // keep compiler happy
    }
    if (!n_in_effective.isFinite()) {
        SLog(EWarn, "Non-finite n_in_effective: %s", n_in_effective.toString().c_str());
        return false;
    }

    Assert(math::abs(n_in_effective.length() - 1) < Epsilon);

    /* Zv */
    Float sigma_sp = m_sigS * (1 - m_g);
    Float sigma_tp = sigma_sp + m_sigA;
    switch (dipConf.zvMode) {
    case EFrisvadEtAlZv: {
        if (sigma_tp == 0 || sigma_sp == 0)
            return false;
        Float D = 1./(3.*sigma_tp);
        Float alpha_p = sigma_sp / sigma_tp;
        Float d_e = 2.131 * D / sqrt(alpha_p);
        Float A = dEon_A(m_eta);
        zv = 2*A*d_e;
        break; }
    case EBetterDipoleZv: {
        if (sigma_tp == 0)
            return false;
        Float D = (2*m_sigA + sigma_sp)/(3*math::square(sigma_tp));
        Float A = dEon_A(m_eta);
        zv = 4*A*D;
        break; }
    case EClassicDiffusion: {
        if (sigma_tp == 0)
            return false;
        Float Fdr = fresnelDiffuseReflectance(1 / m_eta);
        Float A = (1 + Fdr) / (1 - Fdr);
        Float D = 1./(3*sigma_tp);
        zv = 4*A*D;
        break; }
    default:
        SLog(EError, "Unknown VirtSourceHeight mode %d", dipConf.zvMode);
        return false;
    }

    return true;
}



FINLINE bool DipoleModel::realToVirt(
        const DipoleConfig &dipConf, const Monopole &r, Monopole &v) const {
    Vector n_in_effective;
    Float zv;
    if (!getDipoleParameters(r, dipConf, n_in_effective, zv))
        return false;
#ifdef MTS_DIPOLE_DEBUG
    if (dipConf.shouldRejectDin(r, n_in_effective)) {
        Log(EWarn, "Transforming real to virtual, but real should have been rejected!");
    }
#endif

    v = r;
    v.R = r.R - zv * n_in_effective;
    /* The 'virtual normal' n_in_virt also gets transformed so that upon 
     * transforming a virtual direction d_in_virt to its corresponding 
     * d_in, that d_in is on the correct side of the actual boundary as 
     * determined by n_in. */
    v.d_in = r.d_in  -  2*dot(n_in_effective, r.d_in) * n_in_effective;
    v.n_in = r.n_in  -  2*dot(n_in_effective, r.n_in) * n_in_effective;

    return true;
}

FINLINE bool DipoleModel::shouldRejectRealDin(
        const DipoleConfig &dipConf, const Monopole &r) const {
    Vector n_in_effective;
    Float _zv;
    if (!getDipoleParameters(r, dipConf, n_in_effective, _zv))
        return true; // can't even get the parameters, so reject

    return dipConf.shouldRejectDin(r, n_in_effective);
}


FINLINE bool DipoleModel::virtToReal(
        const Vector &R_real, const Vector &n_in_real,
        const DipoleConfig &dipConf, const Monopole &v, Monopole &r) const {
    r = v;
    r.R = R_real;
    r.n_in = n_in_real;

    Vector n_in_effective;
    Float zv;
    if (!getDipoleParameters(r, dipConf, n_in_effective, zv))
        return false;

    r.d_in = v.d_in  -  2*dot(n_in_effective, v.d_in) * n_in_effective;

    Vector R_check = v.R + zv * n_in_effective;
    Vector n_in_check = v.n_in  -  2*dot(n_in_effective, v.n_in) * n_in_effective;
    if (!(((r.R + R_check).length() == 0 && r.R.isZero() && R_check.isZero())
            || (r.R.isZero() && R_check.length() < 1e-10) /* ... magic number :P ... */
            || (r.R - R_check).length() / (r.R + R_check).length() < ShadowEpsilon)) {
        Log(EWarn, "virtToReal R issue: %s vs %s, diff %s, reldiff %e",
                r.R.toString().c_str(),
                R_check.toString().c_str(),
                (r.R - R_check).toString().c_str(),
                (r.R - R_check).length() / (r.R + R_check).length());
    }
    AssertWarn((r.n_in - n_in_check).length() / (r.n_in + n_in_check).length() < Epsilon);

    return true;
}


FINLINE Float DipoleModel::evalDipole(
        const Normal &n0, const Vector &u0_external,
        const Normal &nL, const Vector &uL_external,
        const Vector &R, const void *extraParams,
        const DipoleConfig &dipConf) const {
    /* If reciprocal is requested, nL should be finite and uL_external should point
     * along nL. */
    Assert(!dipConf.reciprocal || nL.isFinite());
    Assert(!dipConf.reciprocal || dot(uL_external,nL) >= -Epsilon); // positive with small margin for roundoff errors
    if (nL.isFinite() && dot(uL_external,nL) <= 0) // clamp to protect against roundoff errors
        return 0.0f;

    Float fwd = evalDipole_internal(
            n0, u0_external, nL, uL_external, R, extraParams, dipConf);

    if (!dipConf.reciprocal)
        return fwd;

    Float rev = evalDipole_internal(
            nL, -uL_external, n0, -u0_external, -R, extraParams, dipConf);
    return 0.5 * (fwd + rev);
}

/// Always returns non-reciprocal transport, regardless of dipConf.reciprocal
FINLINE Float DipoleModel::evalDipole_internal(
        const Normal &n0, const Vector &u0_external,
        const Normal &nL, const Vector &uL_external,
        const Vector &R, const void *extraParams,
        const DipoleConfig &dipConf) const {

#if MTS_DIPOLE_MODEL_REJECT_INCOMING_WRT_TRUE_SURFACE_NORMAL
    if (dot(u0_external, n0) >= 0)
        return 0.0f;
#endif

    /* Handle eta != 1 case by 'refracting' the 'external' directions
     * u0_external and uL_external to 'internal' directions u0 and uL. We
     * keep the directions pointing along the propagation direction of
     * light (i.e. not the typical refract as in BSDFs, for instance, which
     * flips to the other side of the boundary). */
    Float _cosThetaT, F0, FL;
    Vector u0 = refract(-u0_external, n0, m_eta, _cosThetaT, F0);
    Vector uL = -refract(uL_external, nL, m_eta, _cosThetaT, FL);
    Float fresnelTransmittance = (1-F0)*(1-FL);

    if (m_eta == 1)
        Assert(u0 == u0_external  &&  uL == uL_external);

    if (u0.isZero() || uL.isZero()) {
        if (m_eta > 1)
            Log(EWarn, "Could not refract, which is weird because we have a "
                    "higher ior! (eta=%f)", m_eta);
        return 0.0f;
    }


    Monopole r(n0, u0, nL, uL, R, extraParams, dipConf.useEffectiveBRDF);
    Monopole v;

    if (shouldRejectRealDin(dipConf, r))
        return 0;

    if (!realToVirt(dipConf, r, v) && (dipConf.dipoleMode & EVirt))
        return 0.0f;

    Float real = 0, virt = 0;
    if (dipConf.dipoleMode & EReal)
        real = evalMonopole(r);
    if (dipConf.dipoleMode & EVirt)
        virt = evalMonopole(v);
    Float transport;
    switch (dipConf.dipoleMode) {
        case ERealAndVirt: transport = real - virt; break;
        case EReal:        transport = real; break;
        case EVirt:        transport = virt; break; // note: positive sign
        default: Log(EError, "Unknown dipoleMode: %d", dipConf.dipoleMode); return 0;
    }

    return transport * fresnelTransmittance;
}



FINLINE Float DipoleModel::sampleExtraParamsDipole(
        const Normal &n_in,  const Vector *d_in,
        const Normal &n_out, const Vector &d_out,
        const Vector &R, void *extraParams,
        const DipoleConfig &dipConf, Sampler *sampler) const {

    Monopole real(n_in, d_in, n_out, d_out, R, extraParams, dipConf.useEffectiveBRDF);
    if (shouldRejectRealDin(dipConf, real))
        return 0;

    if (dipConf.dipoleMode == EReal) {
        return sampleExtraParamsMonopole(real, extraParams, sampler);
    }

    Monopole virt;
    if (!realToVirt(dipConf, real, virt))
        return 0.0f;

    if (dipConf.dipoleMode == EVirt) {
        return sampleExtraParamsMonopole(virt, extraParams, sampler);
    }

    Float realWeight;
    if (d_in == NULL) {
        realWeight = realSourceWeight_margOverParamsAndDirections(real, virt);
    } else {
        realWeight = realSourceWeight_margOverParams(real, virt);
    }

    Float realPdf, virtPdf;
    Assert(real.extraParams == virt.extraParams);
    if (realWeight == 1 || sampler->next1D() < realWeight) {
        realPdf = sampleExtraParamsMonopole(real, extraParams, sampler);
        virtPdf = pdfExtraParamsMonopole(virt);
    } else {
        virtPdf = sampleExtraParamsMonopole(virt, extraParams, sampler);
        realPdf = pdfExtraParamsMonopole(real);
    }
    Float pdf = realWeight * realPdf + (1 - realWeight) * virtPdf;
#ifdef MTS_DIPOLE_DEBUG
    if (pdf == 0)
        return 0;
    Float pdf2 = pdfExtraParamsDipole(n_in, d_in, n_out, d_out, R, 
            extraParams, dipConf);
    if (!(math::abs((pdf - pdf2)/(pdf + pdf2)) < 1e-3))
        SLog(EWarn, "Inconsistent pdfs: %e vs %e, rel diff %e",
                pdf, pdf2, (pdf-pdf2)/(pdf+pdf2));
#endif
    return pdf;
}

FINLINE Float DipoleModel::pdfExtraParamsDipole(
            const Normal &n_in,  const Vector *d_in,
            const Normal &n_out, const Vector &d_out,
            const Vector &R, const void *extraParams,
            const DipoleConfig &dipConf) const {

    Monopole real(n_in, d_in, n_out, d_out, R, extraParams, dipConf.useEffectiveBRDF);
    if (shouldRejectRealDin(dipConf, real))
        return 0;

    if (dipConf.dipoleMode == EReal) {
        return pdfExtraParamsMonopole(real);
    }

    Monopole virt;
    if (!realToVirt(dipConf, real, virt))
        return 0.0f;

    if (dipConf.dipoleMode == EVirt) {
        return pdfExtraParamsMonopole(virt);
    }

    Float realWeight;
    if (d_in == NULL) {
        realWeight = realSourceWeight_margOverParamsAndDirections(real, virt);
    } else {
        realWeight = realSourceWeight_margOverParams(real, virt);
    }

    Float realPdf, virtPdf;
    Assert(real.extraParams == virt.extraParams);
    realPdf = pdfExtraParamsMonopole(real);
    virtPdf = pdfExtraParamsMonopole(virt);
    return realWeight * realPdf + (1 - realWeight) * virtPdf;
}

FINLINE Float DipoleModel::sampleDirectionDipole(
        const Normal &n_in,  Vector       &d_in,
        const Normal &n_out, const Vector &d_out,
        const Vector &R, const void *extraParams,
        const DipoleConfig &dipConf, Sampler *sampler) const {
#ifdef MTS_DIPOLE_DEBUG
    d_in = Vector(0./0.);
#endif

    Monopole real(n_in, NULL, n_out, d_out, R, extraParams, dipConf.useEffectiveBRDF);

    if (dipConf.dipoleMode == EReal) {
        Float pdf = sampleDirectionMonopole(real, sampler);
        if (!(dot(real.d_in, real.n_in) <= 0)) { // protect against roundoff
            Assert(dot(real.d_in, real.n_in) <= Epsilon);
            return 0;
        }
        if (shouldRejectRealDin(dipConf, real))
            return 0;
        d_in = real.d_in;
        return pdf;
    }

    Monopole virt;
    if (!realToVirt(dipConf, real, virt))
        return 0;

    if (dipConf.dipoleMode == EVirt) {
        Float pdf = sampleDirectionMonopole(virt, sampler);
        Assert(dot(virt.d_in, virt.n_in) <= Epsilon);

        /* Don't forget, we have to return the 'real source' direction! */
        if (!virtToReal(R, n_in, dipConf, virt, real))
            return 0;
        if (!(dot(real.d_in, real.n_in) <= 0)) { // protect against roundoff
            Assert(dot(real.d_in, real.n_in) <= Epsilon);
            return 0;
        }

        if (shouldRejectRealDin(dipConf, real))
            return 0;
        d_in = real.d_in;
        return pdf;
    }

    Float realWeight = realSourceWeight_margOverDirections(real, virt);

    Float realPdf, virtPdf;
    Assert(real.extraParams == virt.extraParams);
    if (realWeight == 1 || sampler->next1D() < realWeight) {
        realPdf = sampleDirectionMonopole(real, sampler);
        if (!(dot(real.d_in, real.n_in) <= 0)) { // protect against roundoff
            Assert(dot(real.d_in, real.n_in) <= Epsilon);
            return 0;
        }
        if (shouldRejectRealDin(dipConf, real))
            return 0;
        if (!realToVirt(dipConf, real, virt))
            return 0;
        Assert(dot(virt.d_in, virt.n_in) <= Epsilon);
        virtPdf = pdfDirectionMonopole(virt);
    } else {
        virtPdf = sampleDirectionMonopole(virt, sampler);
        Assert(dot(virt.d_in, virt.n_in) <= Epsilon);
        if (!virtToReal(R, n_in, dipConf, virt, real))
            return 0;
        if (!(dot(real.d_in, real.n_in) <= 0)) { // protect against roundoff
            Assert(dot(real.d_in, real.n_in) <= Epsilon);
            return 0;
        }
        if (shouldRejectRealDin(dipConf, real))
            return 0;
        realPdf = pdfDirectionMonopole(real);
    }
    d_in = real.d_in;
    Float pdf = realWeight * realPdf + (1 - realWeight) * virtPdf;
#ifdef MTS_DIPOLE_DEBUG
    if (pdf == 0)
        return 0;
    Float pdf2 = pdfDirectionDipole(n_in, d_in, n_out, d_out, R, 
            extraParams, dipConf);
    if (!(math::abs((pdf - pdf2)/(pdf + pdf2)) < 1e-3))
        SLog(EWarn, "Inconsistent pdfs: %e vs %e, rel diff %e",
                pdf, pdf2, (pdf-pdf2)/(pdf+pdf2));
#endif
    return pdf;
}

FINLINE Float DipoleModel::pdfDirectionDipole(
        const Normal &n_in,  const Vector &d_in,
        const Normal &n_out, const Vector &d_out,
        const Vector &R, const void *extraParams,
        const DipoleConfig &dipConf) const {

    Monopole real(n_in, d_in, n_out, d_out, R, extraParams, dipConf.useEffectiveBRDF);
    Assert(dot(real.d_in, real.n_in) <= Epsilon);

    if (shouldRejectRealDin(dipConf, real))
        return 0;

    if (dipConf.dipoleMode == EReal)
        return pdfDirectionMonopole(real);

    Monopole virt;
    if (!realToVirt(dipConf, real, virt))
        return 0;
    Assert(dot(virt.d_in, virt.n_in) <= Epsilon);

    if (dipConf.dipoleMode == EVirt)
        return pdfDirectionMonopole(virt);

    Float realWeight = realSourceWeight_margOverDirections(real, virt);

    Float realPdf = pdfDirectionMonopole(real);
    Float virtPdf = pdfDirectionMonopole(virt);
    return realWeight * realPdf + (1 - realWeight) * virtPdf;
}

MTS_NAMESPACE_END

#endif /* __DIPOLE_MODEL_H */

