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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/dss.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/truncnorm.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/vmf.h>
#include <gsl/gsl_sf_lambert.h>
#include "../medium/materials.h"
#include "fwdscat.h"

/* For the soft inverse sample ray cosine sampling */
#define MARG_INVSAMPLER_EPSILON 1e-5f

#define MTS_FWDDIP_DEBUG_USE_SIMPLE_JENSEN_SURFACE_SAMPLER false

MTS_NAMESPACE_BEGIN


/// Helper functions to sample proportinal to 1/(xEpsilon + x) for x on [0..xMax]
static inline Float inverseSampler_sample(Float xEps, Float xMax, Float u) {
    SAssert(u >= 0 && u <= 1);
    return -xEps  -  (xMax+xEps) * gsl_sf_lambert_W0(
                -exp((-u*xMax - xEps)/(xEps + xMax))
                    * pow(xEps/(xEps + xMax), (Float)1.-u));
}
static inline Float inverseSampler_pdf(Float xEps, Float xMax, Float x) {
    if (x <= 0 || x >= xMax)
        return 0;
    return (xMax - x) / (
            (xEps + x) * ((xEps + xMax)*log((xEps+xMax)/xEps)  -  xMax));
}

/// Helper functions to sample according to pdf(x) = -log(x) for x on [0..1]
static inline Float logDivergenceSampler_sample(Float u) {
    SAssert(u >= 0 && u <= 1);
    return -u/gsl_sf_lambert_Wm1(-u/M_E);
}
static inline Float logDivergenceSampler_pdf(Float x) {
    if (x <= 0 || x >= 1)
        return 0;
    return -log(x);
}


class MTS_EXPORT_RENDER FwdDipSmallLengthRadialSampler2D final : public Sampler1D {
public:
    FwdDipSmallLengthRadialSampler2D(Spectrum sigma_s, Spectrum g) :
        m_p(0.5*sigma_s*(Spectrum(1.0f) - g)) { }

    virtual bool sample(int channel, Float &r,
            Sampler *sampler, Float *thePdf = NULL) const {
        Float p = m_p[channel];
        Float u = sampler->next1D();
        r = sqrt(M_PI/6) * (1 - sqrt((1-u))) / p;
        if (thePdf)
            *thePdf = pdf(channel, r);
        return true;
    }
    virtual Float pdf(int channel, Float r) const {
        // pdf in the plane! -> includes 1/(2*pi*r) factor
        Float p = m_p[channel];
        Float r_p1 = r*p;
        Float max_r_p1 = sqrt(M_PI/6);
        if (r_p1 > max_r_p1)
            return 0;
        return p*p * sqrt(6.)*(sqrt(M_PI) - sqrt(6.)*r_p1) / (M_PI*M_PI*r_p1);
    }

    MTS_DECLARE_CLASS();
protected:
    const Spectrum m_p;
};

/* Higher than 1: confines pdf more so that its peak gets higher but we can
 * cut off too soon -- lower than 1; broader pdf, but we might be
 * undersampling the extreme peak then */
#define DEFAULT_P_SAFETY_FACTOR 1

/**
 * \brief Sampler in the 2D plane perpendicular to the outgoing direction,
 * i.e. with the projection along the outgoing direction.
 */
class MTS_EXPORT_RENDER FwdDipSmallLengthSamplerPerpToDir final : public TangentSampler2D {
public:
    FwdDipSmallLengthSamplerPerpToDir(Spectrum sigma_s, Spectrum g,
            Float pSafetyFactor = DEFAULT_P_SAFETY_FACTOR) :
                m_p(pSafetyFactor * 0.5*sigma_s*(Spectrum(1.0f) - g)) {
        Assert(sigma_s.isFinite());
        Assert(g.isFinite());
        Assert(m_p.isFinite());
    }

    /**
     * \brief Sample a point x in the 2D plane.
     *
     * \param x The point to sample in the (tangent) plane, with respect to
     *          the intersection point. The first coordinate is along the
     *          'forward' direction, i.e. along the tangential component of
     *          the queried radiance direction.
     *          Due to symmetry, the sign of the second coordinate should be
     *          irrelevant and during sampling the sign will be randomized.
     * \param cosTheta The cosine between the query direction and the
     *                 plane normal. Taken to be positive.
     */
    virtual bool sample(int channel, Vector2 &x, Float cosTheta,
            const Vector2 &xLo, const Vector2 &xHi,
            Sampler *sampler, Float *thePdf = NULL) const {
        Assert(channel>=0);
        Float p = m_p[channel];
        if (p == 0)
            return false;
        Float phi = sampler->next1D() * TWO_PI;
        Float sinPhi,cosPhi;
        math::sincos(phi, &sinPhi, &cosPhi);
        Float u = sampler->next1D();
        Float s = sin(1./3.*atan2(2*sqrt(u*(1-u)), 1-2*u));
        Float c = math::safe_sqrt(1 - s*s);
        Float maxBoundDist = std::max(xLo.length(), xHi.length());
        Float cutoff = std::min(m_rCutoff, maxBoundDist);
        Float r = cutoff * pow(
                (Float) 0.5 * (1 - c) * (2 + c + sqrt(3.)*s),
                (Float) 3./2.);
        x[0] = sinPhi * r / p;
        x[1] = cosPhi * r / p;
        if (thePdf)
            *thePdf = pdf(channel, x, cosTheta, xLo, xHi);
        return true;
    }

    virtual Float pdf(int channel, Vector2 x, Float cosTheta,
            const Vector2 &xLo, const Vector2 &xHi) const {
        Assert(channel>=0);
        Float p = m_p[channel];
        Float r = x.length() * p;
        Float maxBoundDist = std::max(xLo.length(), xHi.length());
        Float cutoff = std::min(m_rCutoff, maxBoundDist);
        if (p == 0 || r >= cutoff)
            return 0;
        Float pdf_r = 2*(pow(cutoff*cutoff * r, (Float) -1./3.) - 1/cutoff);
        return pdf_r * p*p / (TWO_PI*r);
    }

    MTS_DECLARE_CLASS();

protected:
    virtual ~FwdDipSmallLengthSamplerPerpToDir() { }
    const Spectrum m_p;
    const Float m_rCutoff = 1; // max radial distance
};

/*
 * Sampler in a plane that contains the outgoing direction
 */
class MTS_EXPORT_RENDER FwdDipSmallLengthSamplerAlongDir final : public TangentSampler2D {
public:
    FwdDipSmallLengthSamplerAlongDir(Spectrum sigma_s, Spectrum g,
            Float pSafetyFactor = DEFAULT_P_SAFETY_FACTOR) :
                m_p(pSafetyFactor * 0.5*sigma_s*(Spectrum(1.0f) - g)) {
        Assert(sigma_s.isFinite());
        Assert(g.isFinite());
        Assert(m_p.isFinite());
    }

    /**
     * \brief Sample a point x in the 2D plane.
     *
     * \param x The point to sample in the (tangent) plane, with respect to
     *          the intersection point. The first coordinate is along the
     *          'forward' direction, i.e. along the tangential component of
     *          the queried radiance direction.
     *          Due to symmetry, the sign of the second coordinate should be
     *          irrelevant and during sampling the sign will be randomized.
     * \param cosTheta The cosine between the query direction and the
     *                 plane normal. Taken to be positive.
     */
    virtual bool sample(int channel, Vector2 &x, Float cosTheta,
            const Vector2 &xLo, const Vector2 &xHi,
            Sampler *sampler, Float *thePdf = NULL) const {
        Assert(channel>=0);
        Float p = m_p[channel];
        if (p == 0)
            return false;

        // We work in p=1

        Float Rmax, dMin, dMax;
        getSamplingBounds(channel, xLo, xHi, Rmax, dMin, dMax);

        /* Sampling R: MIS between uniform back-up and importance sampled 
         * weight. */
        Float R;
        Float u = sampler->next1D();
        if (u < RUniformWeight) {
            u = u / RUniformWeight;
            R = Rmax * u;
        } else {
            u = (u - RUniformWeight) / (1 - RUniformWeight);
            /* R is sampled according to "1/(Rmin + R) with a cutoff at Rmax" 
             * (see pdf() for the exact pdf) */
            R = inverseSampler_sample(Rmin, Rmax, u);
        }
        Assert(std::isfinite(R));
        AssertWarn(R >= 0 && R <= Rmax);

        /* Sampling 'sideways' displacement d */
        Float stddev = dMaxSafetyScale * sqrt(R*R*R/6);
        Float d = truncnorm(0, stddev, dMin, dMax, sampler);

        /* We should displace ourselves backwards along the outgoing
         * direction to sample the query point! */
        x[0] = -R/p; // 'backwards'
        x[1] = d/p;

        if (thePdf)
            *thePdf = pdf(channel, x, cosTheta, xLo, xHi);

        return true;
    }

    virtual Float pdf(int channel, Vector2 x, Float cosTheta,
            const Vector2 &xLo, const Vector2 &xHi) const {
        Assert(channel>=0);
        Float p = m_p[channel];
        Float R = -x[0]*p; // displacement is backwards along the query point!
        Float d = x[1]*p;
        Float Rmax, dMin, dMax;
        getSamplingBounds(channel, xLo, xHi, Rmax, dMin, dMax);
        if (p == 0 || R < 0 || R > Rmax)
            return 0;

        // check consistency of bounds (note: dMin < 0, hence the +Epsilon!):
        AssertWarn(d >= dMin*(1+Epsilon) && d <= dMax*(1+Epsilon));

        Float stddev = dMaxSafetyScale * sqrt(R*R*R/6);
        Float dPdf = truncnormPdf(0, stddev, dMin, dMax, d);
        Float RpdfImp = inverseSampler_pdf(Rmin, Rmax, R);
        Float RpdfUnif = 1. / Rmax;
        Float Rpdf = RpdfUnif * RUniformWeight  +  RpdfImp * (1 - RUniformWeight);
        Assert(std::isfinite(Rpdf) && Rpdf >= 0);
        return Rpdf * dPdf * p*p;
    }

    inline void getSamplingBounds(int channel,
            const Vector2 &xLo, const Vector2 &xHi,
            Float &Rmax, Float &dMin, Float &dMax) const {
        Float p = m_p[channel];
        // First dimension (u or d_out):
        /* Remember: d_out is opposite direction than our 
         * projection/sampling direction, so the Rmax bound is given by xLo 
         * (which is negative) */
        Rmax = std::min(RmaxDefault, -xLo[0] * p);
        // Second dimension (v or 'sideways'):
        dMin = p * xLo[1];
        dMax = p * xHi[1];
    }

    MTS_DECLARE_CLASS();

protected:
    virtual ~FwdDipSmallLengthSamplerAlongDir() { }

    const Spectrum m_p;

    /* How much wider we make the sample area, to make sure we have covered
     * the peak properly. */
    const Float dMaxSafetyScale = 2;

    /* Boundaries for sampling R according to "~1/(Rmin + R)" with cut-off 
     * above Rmax and with Rmin denoting the inflection point towards an 
     * asymptotically uniform distribution for R < Rmin (i.e. R -> 0). */
    const Float Rmin = 1e-10;
    const Float RmaxDefault = 0.2; // gets clipped if bounding box is tighter

    /* Sample R uniformly between 0 and Rmax with this weight: */
    const Float RUniformWeight = 0.4;
};

inline IntersectionWeightFunc fwdDipSmallLengthWeightFunc(
        Spectrum sigma_s, Spectrum sigma_a, Spectrum g, bool useR3) {
    Spectrum sigmaSPrime = sigma_s * (Spectrum(1.0f) - g);
    Spectrum sigmaTPrime = sigmaSPrime + sigma_a;
    Spectrum sigmaTr = (sigma_a * sigmaTPrime * 3.0f).sqrt();
    Spectrum pSpec = 0.5 * sigmaSPrime;

    return [=] (const Intersection &its_in, const Intersection &its_out,
            const Vector &d_out, int i) {
        SAssert(i>=0);
        Vector rVec = its_out.p - its_in.p; // dimensionfull
        if (rVec.length() == 0)
            return (Float) 1.0f;
        Float cosTheta = dot(d_out, normalize(rVec));
        Float p = pSpec[i];
        SAssert(p >= 0);
        if (p == 0)
            return (Float) 0.0f;
        Float R = rVec.length() * p; // dimensionless
        Float result;
        if (useR3) {
            // Based on revised derivation, but seems to fall off too quickly for large R?
            result = exp(2.5f*(cosTheta - 1.f)/R)/(R*R*R);
        } else {
            // Based on initial derivation
            //result = exp(3.f*(cosTheta - 1.f)/R)/R;
            result = exp(3.f*(cosTheta - 1.f)/R)/(R*R);
        }

        /* Modulate with long-length asymptotic factor in case we get 
         * intersections that are far away (in the diffusive range), where 
         * our own 1/R^n fall-off isn't fast enough (and hence we give too 
         * much weight to intersections where our approximation isn't 
         * valid!) */
        if (rVec.length() > 1/p) 
            result *= exp(-(rVec.length() - 1/p) * sigmaTr[i]);

        SAssert(std::isfinite(result));
        SAssert(result >= 0);
        return result;
    };
}

class MTS_EXPORT_RENDER RayDirectionSurfaceSampler final : public SurfaceSampler {
public:
    RayDirectionSurfaceSampler(const Spectrum &sigmaTr, const Spectrum &p,
            Float retreatFactor, bool bidirectional,
            const std::string &sampleStrategy, const IntersectionSampler *itsSampler)
            : m_sigmaTr(sigmaTr), m_p(p), m_retreatFactor(retreatFactor),
              m_bidirectional(bidirectional), m_itsSampler(itsSampler) {
        if (sampleStrategy == "pt") {
            // sense also with perturbed tetrahedron
            m_sensingWeight = 0.5;
            m_perturbedWeight = 0.5;
        } else if (sampleStrategy == "ss") {
            // single sense ray
            m_sensingWeight = 0.5;
            m_perturbedWeight = 0;
        } else if (sampleStrategy == "om") {
            // only marginalized
            m_sensingWeight = 0;
            m_perturbedWeight = 0;
        } else {
            Log(EError, "Unknown sampling strategy '%s' for RayDirectionSurfaceSampler", sampleStrategy.c_str());
        }
    }

    virtual Float sample(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes,
            Intersection &newIts, const Spectrum &throughput,
            Sampler *sampler) const {
        Vector rayDir;
        Float dirPdf = sampleDir(its, d_out, scene, shapes, throughput, sampler, rayDir);

        // Intersection prob
        Point startPoint = getStartPoint(its, d_out);
        IntersectionSampler::IntersectionCollectionMode itsCollMode = m_bidirectional
                ? IntersectionSampler::EBidirectionalWithEpsilon
                : IntersectionSampler::EUnidirectionalWithEpsilon;
        Float intersectionProb = m_itsSampler->sample(newIts, scene,
                startPoint, rayDir, its.time, itsCollMode,
                shapes, its, d_out, throughput, sampler);
        if (intersectionProb == 0)
            return 0.0f;

        // Change of var
        Float r = distance(startPoint, newIts.p);
        AssertWarn(r <= m_itsSampler->getItsDistanceCutoff() * 1.01);
        Float surfaceCosine = dot(rayDir, newIts.geoFrame.n);
        Float absSurfCosine = math::abs(surfaceCosine);
        if (absSurfCosine <= MTS_DSS_COSINE_CUTOFF) {
            /* Sampling this surface was badly conditioned -- abandon ship here
             * and let the other sampling strategies save us from bias.
             * TODO: NOTE: we are still slightly biased in the way that
             * our pdf is now not properly normalized anymore! But these
             * occurences should only account for a tiny fraction of
             * 'integrated weight' precisely because the cosine is so
             * small. */
            return 0.0f;
        }
        Float solidAngleToAreaWeight = absSurfCosine/(r*r);
        Float pdf = solidAngleToAreaWeight * dirPdf * intersectionProb;
        if (!std::isfinite(pdf)) {
            Log(EWarn, "Invalid pdf in samplePointFromPhaseFunction: %f, "
                    "saToWeight %e, phase %e, intersect %e",
                    pdf, solidAngleToAreaWeight, dirPdf,
                    intersectionProb);
            return 0;
        }
#ifdef MTS_FWDDIP_DEBUG
        Float pdfCheck = this->pdf(its, d_out, scene, shapes, newIts, throughput);
        if (fabs((pdf - pdfCheck)/pdf) > 1e-4) {
            Log(EWarn, "Inconsistent pdfs: %e vs %e, rel %e",
                    pdf,pdfCheck, (pdf-pdfCheck)/pdf);
        }
#endif
        return pdf;
    }

    inline Float sampleDir(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes, const Spectrum &throughput,
            Sampler *sampler, Vector &rayDir) const {
        /* MIS combo of marginalized and depthSensing, but if 
         * depthSensing fails, then 100% marginalized */
        Float dirPdf;
        if (sampler->next1D() < m_sensingWeight) {
            Float sensePdf;
            if (sampleWithSensing(its, d_out, scene, shapes,
                    throughput, sampler, rayDir, sensePdf)) {
                dirPdf = m_sensingWeight * sensePdf
                        + (1 - m_sensingWeight) * pdfMarginalized(d_out, rayDir);
            } else {
                sampleMarginalized(d_out, sampler, rayDir, dirPdf);
            }
        } else {
            Float margPdf, sensePdf;
            sampleMarginalized(d_out, sampler, rayDir, margPdf);
            if (m_sensingWeight != 0 && sampleWithSensing(
                    its, d_out, scene, shapes, throughput, NULL, rayDir, sensePdf)) {
                dirPdf = m_sensingWeight * sensePdf
                        + (1 - m_sensingWeight) * margPdf;
            } else {
                dirPdf = margPdf;
            }
        }
        return dirPdf;
    }

    /** Sample direction with the length of R marginalized over 0 to 1 
     *  reduced scattering length. */
    inline void sampleMarginalized(const Vector &d_out,
            Sampler *sampler, Vector &rayDir, Float &dirPdf) const {
        /* Note: we are marginalized up to 1 reduced scattering length, 
         * which means that we are no longer dependent on the actual medium 
         * parameters here. */

        Float cosTheta; // dot(rayDir, d_out)
        /* Note: d_out == -w_o^i in text, and \hat{r} == -rayDir -> 
         * cosTheta stays within [-1..0] and is close to -1. */
        Float u = sampler->next1D();
        if (u < 0.5) {
            u = 2*u;
            // inverse logarithm divergence
            cosTheta = logDivergenceSampler_sample(u) - 1;
        } else {
            u = 2*(u-0.5);
            // soft 1/x divergence
            cosTheta = inverseSampler_sample(MARG_INVSAMPLER_EPSILON, 1.0f, u) - 1;
        }
        Float sinTheta = math::safe_sqrt(1 - cosTheta*cosTheta);

        Float phi = sampler->next1D() * TWO_PI;
        Float sinPhi,cosPhi;
        math::sincos(phi, &sinPhi, &cosPhi);

        Frame outFrame(d_out);
        rayDir = outFrame.toWorld(Vector(cosPhi*sinTheta, sinPhi*sinTheta, cosTheta));

        dirPdf = pdfMarginalized(d_out, rayDir);
    }


    /** 
     * \brief Trace sensing rays along -d_out to query the depth until we find 
     * geometry.
     *
     * This is a MIS combo of original and perturbed sampling (= 
     * small perturbations on d_out for robustness). If the initially 
     * chosen sampling method (perturbed vs nonperturbed) fails, then try 
     * the other one as well. 
     *
     * If sampler != NULL, then sample rayDir and set the dirPdf.
     * If sampler == NULL, then only set the dirPdf from the given rayDir.
     * Returns false if we couldn't sample or compute pdf. */
    inline bool sampleWithSensing(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes, const Spectrum &throughput,
            Sampler *sampler, Vector &rayDir, Float &dirPdf) const {

        Float pertPdf, unpertPdf;
        if (sampler) {
            // Sampling was requested! (BEWARE, THIS IS A BIT OF A LOGICAL JUNGLE :P)
            if (m_perturbedWeight != 0 && sampler->next1D() < m_perturbedWeight) {
                // Try to sample perturbed
                if (sampleWithSensing_perturbed(its, d_out, scene, shapes,
                        throughput, sampler, rayDir, pertPdf)) {
                    if (sampleWithSensing_internal(its, d_out, scene, shapes,
                            throughput, NULL, rayDir, unpertPdf)) {
                        // could sample pert and unpert
                        dirPdf = m_perturbedWeight * pertPdf
                                + (1 - m_perturbedWeight) * unpertPdf;
                        return true;
                    } else {
                        // could only sample pert
                        dirPdf = pertPdf;
                        return true;
                    }
                } else {
                    // could not sample perturbed, fall back to unperturbed
                    if (sampleWithSensing_internal(its, d_out, scene, shapes,
                            throughput, sampler, rayDir, unpertPdf)) {
                        // could only sample unpert
                        dirPdf = unpertPdf;
                        return true;
                    } else {
                        // could not sample pert nor unpert!
                        return false;
                    }
                }
            } else {
                // Try to sample unperturbed
                if (sampleWithSensing_internal(its, d_out, scene, shapes,
                        throughput, sampler, rayDir, unpertPdf)) {
                    if (m_perturbedWeight != 0 && sampleWithSensing_perturbed(
                            its, d_out, scene, shapes, throughput, NULL, rayDir, pertPdf)) {
                        // could sample pert and unpert
                        dirPdf = m_perturbedWeight * pertPdf
                                + (1 - m_perturbedWeight) * unpertPdf;
                        return true;
                    } else {
                        // could only sample unpert
                        dirPdf = unpertPdf;
                        return true;
                    }
                } else {
                    // could not sample unperturbed, fall back to perturbed
                    if (sampleWithSensing_perturbed(its, d_out, scene, shapes,
                            throughput, sampler, rayDir, pertPdf)) {
                        // could only sample pert
                        dirPdf = pertPdf;
                        return true;
                    } else {
                        // could not sample pert nor unpert!
                        return false;
                    }
                }
            }
        } else {
            // only pdf is requested (this makes the logic a bit cleaner ;-) )
            bool pertWorks   = m_perturbedWeight != 0 && sampleWithSensing_perturbed(
                    its, d_out, scene, shapes, throughput, NULL, rayDir, pertPdf);
            bool unpertWorks = sampleWithSensing_internal(
                    its, d_out, scene, shapes, throughput, NULL, rayDir, unpertPdf);
            if (!pertWorks && !unpertWorks) {
                // could not sample pert nor unpert!
                return false;
            } else if (pertWorks && unpertWorks) {
                // could sample pert and unpert
                dirPdf = m_perturbedWeight * pertPdf
                        + (1 - m_perturbedWeight) * unpertPdf;
            } else if (pertWorks) {
                dirPdf = pertPdf;
            } else {
                dirPdf = unpertPdf;
            }
            return true;
        }
    }

    /** MIS of small perturbations on d_out for robustness. This sampler 
     * 'fails' if any of its perturbations fail. (Note: Might 'rebalance' 
     * to reweight MIS based on those perturbations that don't fail...)*/
    inline bool sampleWithSensing_perturbed(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes, const Spectrum &throughput,
            Sampler *sampler, Vector &rayDir, Float &dirPdf) const {
        const Float perturbSize = 0.1;
        const std::vector<Vector> perturbs = {
            // Vertices of a tetrahedron
            Vector(1, 1, 1) / SQRT_TWO,
            Vector(-1,-1,1) / SQRT_TWO,
            Vector(-1,1,-1) / SQRT_TWO,
            Vector(1,-1,-1) / SQRT_TWO,
        };

        // Sample if that is requested
        int sampledIdx;
        Float sampledDirPdf = 0./0.;
        if (sampler) {
            Float u = sampler->next1D();
            sampledIdx = std::min((int)(u * perturbs.size()), (int)(perturbs.size() - 1));
            Vector pert = perturbSize * perturbs[sampledIdx];
            if (!sampleWithSensing_internal(its, d_out, scene, shapes,
                    throughput, sampler, rayDir, sampledDirPdf, &pert))
                return false;
        }

        // Compute MIS pdf combination
        dirPdf = 0;
        for (int i = 0; i < (int)perturbs.size(); i++) {
            if (sampler && i == sampledIdx) {
                dirPdf += sampledDirPdf / perturbs.size();
            } else {
                Float thisPertPdf;
                Vector pert = perturbSize * perturbs[i];
                if (!sampleWithSensing_internal(its, d_out, scene, shapes, throughput, 
                        NULL, rayDir, thisPertPdf, &pert))
                    return false;
                dirPdf += thisPertPdf / perturbs.size();
            }
        }
        return true;
    }

    /** If sampler != NULL, then sample rayDir and set the dirPdf.
     *  If sampler == NULL, then only set the dirPdf from the given rayDir.
     *  Returns false if we couldn't sample or compute pdf. */
    inline bool sampleWithSensing_internal(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes,
            const Spectrum &throughput, Sampler *sampler,
            Vector &rayDir, Float &dirPdf,
            const Vector *sensePerturb = NULL) const {

        /* Collect all intersections up to a distance of 1/p */
        int chan = throughput.getNonZeroChannel();
        Float tMax = chan == -1 ?
                1.0/m_p.min() : 1.0/m_p[chan];
        std::vector<Intersection> senseIntersections;
        Vector senseDir = sensePerturb ? 
                normalize(-d_out + *sensePerturb) : -d_out;
        scene->rayIntersectFully(Ray(its.p, senseDir, Epsilon, tMax, its.time),
                senseIntersections, &shapes);
        if (senseIntersections.size() == 0)
            return false;

        /* Sample (or set localDir for the pdf if sampling isn't requested) */
        Frame dirFrame(-d_out); // point inwards; construction is deterministic
        Vector localDir;
        if (sampler) {
            int chosenChannel = throughput.sampleWeightedChannel(sampler);
            Float chosenDistance = senseIntersections[
                    (size_t)(sampler->next1D() * senseIntersections.size())].t;
            VonMisesFisherDistr vmf;
            if (!sampleWithSensing_internalVmf(chosenChannel, chosenDistance, vmf))
                return false;
            localDir = vmf.sample(sampler->next2D());
            rayDir = dirFrame.toWorld(localDir);
        } else {
            localDir = dirFrame.toLocal(rayDir);
        }

        /* PDF is MIS combination over the distances of all intersections 
         * and over all channels if that was requested */
        dirPdf = 0;
        for (const Intersection &senseIts : senseIntersections) {
            dirPdf += channelWeightedMean(throughput,
                        [=] (int chan) { return sampleWithSensing_internalPdf(
                                           chan, senseIts.t, localDir.z); });
        }
        dirPdf /= senseIntersections.size();

        return true;
    }

    inline bool sampleWithSensing_internalVmf(
            int channel, Float distance, VonMisesFisherDistr &vmf) const {
        Assert(channel>=0);
        Float r = distance * m_p[channel];
        Assert(r >= 0);
        Float kappa = m_sampleWithSensingSafetyFactor * 2.5 / r;
        if (!std::isfinite(kappa))
            return false;
        Assert(kappa >= 0);
        //kappa = std::min(kappa, 1.0 / Epsilon); // clamp cosTheta in the range of Epsilon
        vmf = VonMisesFisherDistr(kappa);
        return true;
    }

    inline Float sampleWithSensing_internalPdf(
            int channel, Float distance, Float cosTheta) const {
        VonMisesFisherDistr vmf;
        if (!sampleWithSensing_internalVmf(channel, distance, vmf))
            return 0.0f;
        return vmf.eval(cosTheta);
    }

    virtual Float pdf(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes,
            const Intersection &newIts, const Spectrum &throughput) const {
        Point startPoint = getStartPoint(its, d_out);
        Vector rayDirUnnorm = newIts.p - startPoint;
        Float t = rayDirUnnorm.length();
        if (t <= RCPOVERFLOW)
            return 0.0f;
        Vector rayDir = rayDirUnnorm / t;

        /* If we are calculating a pdf that was sampled by another sampling 
         * strategy, then it is possible that we could not have possible 
         * sampled this point ourselves due to the Epsilon when 
         * intersecting rays with geometry to avoid self-intersections! */
        Float rayMinT = ShapeKDTree::getAdaptiveRayMinT(Ray(startPoint, rayDir, its.time));
        if (t < rayMinT)
            return 0.0f;

        Float surfaceCosine = dot(rayDir, newIts.geoFrame.n);
        Float absSurfCosine = math::abs(surfaceCosine);
        if (absSurfCosine <= MTS_DSS_COSINE_CUTOFF)
            return 0.0f; // sampling & pdf are badly conditioned -> bail out

        Float dirPdf = pdfDir(its, d_out, scene, shapes, throughput, rayDir);

        // Intersection prob
        IntersectionSampler::IntersectionCollectionMode itsCollMode = m_bidirectional
                ? IntersectionSampler::EBidirectionalWithEpsilon
                : IntersectionSampler::EUnidirectionalWithEpsilon;
        Float intersectionProb = m_itsSampler->pdf(newIts, scene,
                startPoint, rayDir, its.time, itsCollMode,
                shapes, its, d_out, throughput);
        if (intersectionProb == 0)
            return 0.0f;

        // Change of var
        Float r = distance(startPoint, newIts.p);
        AssertWarn(r <= m_itsSampler->getItsDistanceCutoff() * 1.01);
        Float solidAngleToAreaWeight = absSurfCosine/(r*r);
        Float pdf = solidAngleToAreaWeight * dirPdf * intersectionProb;
        if (!std::isfinite(pdf))
            return 0;

        return pdf;
    }

    inline Float pdfDir(const Intersection &its,
            const Vector &d_out, const Scene *scene,
            const std::vector<const Shape *> &shapes, const Spectrum &throughput,
            const Vector &rayDir) const {
        // MIS combo between marginalized and 'sensing' for direction sampling
        Float dirPdf;
        Float margPdf = pdfMarginalized(d_out, rayDir);
        Float sensePdf;
        Vector nonConstRayDir(rayDir); // to make compiler happy
        if (m_sensingWeight != 0 && sampleWithSensing(
                its, d_out, scene, shapes, throughput, NULL, nonConstRayDir, sensePdf)) {
            dirPdf = m_sensingWeight * sensePdf
                    + (1 - m_sensingWeight) * margPdf;
        } else {
            dirPdf = margPdf;
        }
        return dirPdf;
    }

    inline Float pdfMarginalized(const Vector &d_out,
            const Vector &rayDir) const {
        /* Note: we are marginalized up to 1 reduced scattering length, 
         * which means that we are no longer dependent on the actual medium 
         * parameters here. */
        Float cosTheta = dot(d_out, rayDir);
        if (cosTheta >= 0 || cosTheta <= -1 /* <- for roundoff issues */)
            return 0.0f;
        Float x = cosTheta + 1; // back to [0..1] with peak at 0
        Float logCosPdf = logDivergenceSampler_pdf(x);
        Float invCosPdf = inverseSampler_pdf(MARG_INVSAMPLER_EPSILON, 1.0f, x);
        Float dirPdf = 0.5 * (logCosPdf + invCosPdf) * INV_TWOPI;
        if (dirPdf <= RCPOVERFLOW)
            return 0.0f;
        return dirPdf;
    }


    MTS_DECLARE_CLASS();
protected:
    Point getStartPoint(const Intersection &its, const Vector &d_out) const {
        /* So we can actually intersect ourself locally in case of planar
         * geometry: retreat a bit below the surface in the (opposite)
         * direction of the outgoing query -- d_out points outwards */
        return its.p - m_retreatFactor/m_sigmaTr.max() * d_out;
    }

    /// Virtual destructor
    virtual ~RayDirectionSurfaceSampler() { }

    const Spectrum m_sigmaTr;

    const Spectrum m_p;

    /// In fractions of the (channel-wise minimal) reduced scattering length
    const Float m_retreatFactor;
    const bool m_bidirectional;

    /** When sampling with 'depth sensing': scale the VMF kappa by this 
     * amount (<1 to widen peak) */
    const Float m_sampleWithSensingSafetyFactor = 1;

    /** Fraction of times to attempt 'depth sensing' sampling as opposed to 
     * marginalized sampling */
    Float m_sensingWeight;

    /** If we are 'depth sensing', fraction of times to use the ensemble of
     * perturbed rays */
    Float m_perturbedWeight;

    const ref<const IntersectionSampler> m_itsSampler;
};





/*!\plugin{fwddip}{Forward Scattering Dipole subsurface scattering model}
 * \parameters{
 *     \parameter{material}{\String}{
 *         Name of a material preset, see
 *         \tblref{medium-coefficients}. \default{\texttt{skin1}}
 *     }
 *     \parameter{sigmaA, sigmaS}{\Spectrum}{
 *         Absorption and scattering
 *         coefficients of the medium in inverse scene units.
 *         These parameters are mutually exclusive with \code{sigmaT} and \code{albedo}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{sigmaT, albedo}{\Spectrum}{
 *         Extinction coefficient in inverse scene units
 *         and a (unitless) single-scattering albedo.
 *         These parameters are mutually exclusive with \code{sigmaA} and \code{sigmaS}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{g}{\Float}{
 *         Anisotropy of the phase function. The Forward Scattering Dipole
 *         model works best for strongly forward scattering media, i.e.
 *         \code{g} close to 1.
 *     }
 *     \parameter{scale}{\Float}{
 *         Optional scale factor that will be applied to the \code{sigma*} parameters.
 *         It is provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}
 *     }
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *         numerically or using a known material name.
 *         For index-mismatched media, it is advised to keep this Forward
 *         Scattering Dipole subsurface scattering model index-matched
 *         internally (i.e. \code{intIOR} = \code{extIOR}) and to couple this
 *         to a proper refractive BSDF that will explicitly handle the
 *         boundary conditions.
 *         \default{based on \code{material}}
 *     }
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *         numerically or using a known material name.
 *         For index-mismatched media, it is advised to keep this Forward
 *         Scattering Dipole subsurface scattering model index-matched
 *         internally (i.e. \code{intIOR} = \code{extIOR}) and to couple this
 *         to a proper refractive BSDF that will explicitly handle the
 *         boundary conditions.
 *         \default{based on \code{material}}
 *     }
 *     \parameter{maxInternalReflections}{\Integer}{
 *         Maximum number of internal reflections within a single subsurface
 *         scattering model evaluation. Limit this to prevent issues with
 *         energy violation and excessive fireflies. This is only relevant
 *         when we are index-matched internally, and coupled to an explicit
 *         (index-mismatched) boundary BSDF (also see \code{intIOR} and
 *         \code{extIOR}).
 *         \default{2}
 *     }
 *     \parameter{numSIRsurface}{\Integer}{
 *         Number of tentative Resampled Importance Sampling samples to generate
 *         when sampling incoming query locations and all subsequent parameters
 *         (directions and path length). This can temper outliers, especially
 *         for self-illuminating media and/or internal reflections, in return
 *         for an increase in computation time (which is typically worth it).
 *         A good starting value for the Forward Scattering Dipole Model is 4.
 *         See also \code{SIRnonSurfaceOversamplingFactor}.
 *     }
 *     \parameter{SIRnonSurfaceOversamplingFactor}{\Integer}{
 *         Additional factor of Resampled Importance Sampling: number of
 *         proposals to generate for the directions and path length for each
 *         tentatively sampled incoming query location. This can be tuned
 *         separately because those sampling steps are typically cheaper
 *         than sampling an incoming query location, which uses ray tracing
 *         to project onto the geometry.
 *         See also \code{numSIRsurface}.
 *         \default{1}
 *     }
 *     \parameter{directSampling}{\Boolean}{
 *         Use direct sampling strategies?
 *         \default{\code{true}}
 *     }
 *     \parameter{directSamplingMIS}{\Boolean}{
 *         Use Multiple Importance Sampling to combine the direct sampling
 *         strategy with indirect sampling? (Only has an effect when
 *         \code{directSampling} is \code{true})
 *         \default{\code{true}}
 *     }
 *     \parameter{singleChannel}{\Boolean}{
 *         Force the BSSRDF to only select a single colour channel?
 *         \default{\code{false}}
 *     }
 *     \parameter{cutoffNumEffectiveTransportLengths}{\Float}{
 *         When sampling an incoming query point on the surface of the
 *         associated geometry, don't consider points that are further away
 *         from the outgoing point than this number of effective transport
 *         lengths. For `coloured' absorption, this corresponds to the longest
 *         effective transport length over all channels.
 *         \default{10}
 *     }
 *     \parameter{direction_min_mu}{\Float}{
 *         Regularisation for extremely small overall direction deviations, 
 *         with overall deviation cosine equal to 1-direction_min_mu.
 *         For experiments with extremely sparse media (nearly transparent),
 *         this should be reduced to a smaller value, e.g. 1e-11
 *         (when compiled for using double precision!).
 *     }
 *     \parameter{reciprocal}{\Boolean}{
 *         Force reciprocity of the model?
 *         \default{\code{false}}
 *     }
 *     \parameter{zvMode}{\String}{
 *         The method used for determining the $z_\mathrm{v}$ displacement
 *         of the virtual source.
 *         \begin{enumerate}[(i)]
 *             \item \code{diff}:    Like the classical Jensen et al. dipole.
 *             \item \code{better}:  Like the better dipole of d'Eon.
 *             \item \code{Frisvad}: Like the directional dipole of
 *                                   Frisvad et al.
 *         \end{enumerate}
 *         \default{\texttt{diff}}
 *     }
 *     \parameter{tangentMode}{\String}{
 *         The method used for determining the normal of the tangent plane for the dipole configuration.
 *         \begin{enumerate}[(i)]
 *             \item \code{incoming}: Use the normal of the incoming point.
 *             \item \code{outgoing}: Use the normal of the outgoing point.
 *             \item \code{Frisvad}:  Use the modified tangent plane of
 *                                    the directinal dipole of Frisvad et al.
 *             \item \code{FrisvadMean}: Like the modified tangent plane of
 *                                    the directinal dipole of Frisvad et al.,
 *                                    but based on the 'average' normal at
 *                                    the incoming and outgoing points
 *                                    instead of on the incoming normal.
 *         \end{enumerate}
 *         \default{\texttt{Frisvad}}
 *     }
 *     \parameter{useEffectiveBRDF}{\Boolean}{
 *         Instead of using a full BSSRDF, use the associated effective
 *         BRDF approximation. This can be useful for objects that are
 *         densely scattering or far away from the camera, as the effective
 *         BRDF is much more easy to compute and its estimator has a
 *         greatly diminished variance.
 *         \default{\code{false}}
 *     }
 *     \parameter{onlyReal}{\Boolean}{
 *         Only show the contributions of the real source of the dipole.
 *         This is mainly for testing/debugging purposes.
 *         \default{\code{false}}
 *     }
 *     \parameter{onlyVirt}{\Boolean}{
 *         Only show the contributions of the (positive) virtual source of
 *         the dipole. This is mainly for testing/debugging purposes.
 *         \default{\code{false}}
 *     }
 *     \parameter{allowIncomingOutgoingDirections}{\Boolean}{
 *         Do we allow outgoing directions of the dipole model to actually
 *         point back towards the inside? This option is mostly only useful
 *         for debugging. For instance, it can be used to check the
 *         validity of the boundary condition assumptions, which state that
 *         (for an index-matched boundary) \emph{no} radiance should be
 *         emitted back to the inside of the boundary.
 *         \default{\code{false}}
 *     }
 * }
 *
 * This plugin implements the forward scattering dipole model from
 * Frederickx and Dutr\'e (`A Forward Scattering Dipole Model from a
 * Functional Integral Approximation', SIGGRAPH2017).
 *
 * This subsurface scattering model currently only works with the
 * \c volpath integrator. It is recommended to use path splitting (combined
 * with Russian roulette) and ideally with a statistically robust rendering
 * method such as provided by the \c RobustAdaptiveMC integrator, which can
 * wrap around \c volpath. Doing so should drastically increase the
 * convergence speed by lessening the effect of fireflies.
 *
 * For scenes where the geometry associated with this subsurface model is
 * less complex than the geometry of the rest of the scene, it might also
 * be worthwile to increase \c numSIR to obtain better estimates of the
 * subsurface contribution per ray.
 */
typedef DipoleDSS<FwdScat> FwdDip; // TODO can the documentation generation handle/find this?

template <>
void FwdDip::configure() {

    if (getNumSamplers() > 0)
        return; // Already configured!

    Spectrum sigmaSPrime = m_sigmaS * (Spectrum(1.0f) - m_g);
    Spectrum sigmaTPrime = sigmaSPrime + m_sigmaA;
    /* Effective transport extinction coefficient */
    Spectrum sigmaTr = (3 * m_sigmaA * sigmaTPrime).sqrt();

    //Float mu = 1 - m_g.average();

    Spectrum p_spectrum = (0.5*sigmaSPrime);


    // No need for fancy planar samplers if we are using the effective BRDF!
    if (m_dipConf.useEffectiveBRDF) {
        registerSampler(1.0, new BRDFDeltaSurfaceSampler());
    } else {
        /* MIS weights for surface sampler: */
        /* Classical dipole for large lengths, weight 1/3 because this 
         * should work almost equally wel in any of the three 
         * projection directions */
        const Float jensenWeight = 1./3.;
        /* Dedicated sampler for small lengths, much more sensitive to 
         * the projection direction, so give unit weight to each one. */
        const Float smallLengthWeight = 1;

        ref<InstanceManager> manager = new InstanceManager();

        /* Intersection sampler */
        ref<IntersectionSampler> itsSamplerExactJensenDipole =
                new WeightIntersectionSampler(distanceWeightWrapper(
                        makeExactDiffusionDipoleDistanceWeight(
                        m_sigmaA, m_sigmaS, m_g, m_eta)),
                m_itsDistanceCutoff, m_numItsLayers);

        ref<IntersectionSampler> itsSamplerEffectiveExtinction =
                new WeightIntersectionSampler(distanceWeightWrapper(
                        makeExponentialDistanceWeight(sigmaTr)),
                m_itsDistanceCutoff, m_numItsLayers);

        ref<IntersectionSampler> itsSamplerFwdDipSmallLenR2 =
                new WeightIntersectionSampler(
                        fwdDipSmallLengthWeightFunc(
                            m_sigmaS, m_sigmaA, m_g, false),
                m_itsDistanceCutoff, m_numItsLayers);
        ref<IntersectionSampler> itsSamplerFwdDipSmallLenR3 =
                new WeightIntersectionSampler(
                        fwdDipSmallLengthWeightFunc(
                            m_sigmaS, m_sigmaA, m_g, true),
                m_itsDistanceCutoff, m_numItsLayers);

        std::vector<std::pair<Float, const IntersectionSampler*> > is;
        is.push_back(std::make_pair(0.1, itsSamplerEffectiveExtinction.get()));
        is.push_back(std::make_pair(1.0, itsSamplerExactJensenDipole.get()));
        is.push_back(std::make_pair(0.5, itsSamplerFwdDipSmallLenR2.get()));
        is.push_back(std::make_pair(0.5, itsSamplerFwdDipSmallLenR3.get()));
        ref<IntersectionSampler> itsSampler = new MISIntersectionSampler(is);


        /* Plane projection sampler */
        ref<TangentSampler2D> exactJensenDipoleSampler = new RadialSampler2D(
                new RadialExactDipoleSampler2D(m_sigmaA, m_sigmaS, m_g, m_eta));

#if MTS_FWDDIP_DEBUG_USE_SIMPLE_JENSEN_SURFACE_SAMPLER
        /* Double the weight because this is the most 'natural' sampler 
         * for truly near-planar boundaries */
        registerSampler(2.0, new ProjSurfaceSampler(
                DSSProjFrame::EDirectionOut,
                exactJensenDipoleSampler, itsSampler.get()));
        registerSampler(1.0, new ProjSurfaceSampler(
                DSSProjFrame::EDirectionDirection,
                exactJensenDipoleSampler, itsSampler.get()));
        registerSampler(1.0, new ProjSurfaceSampler(
                DSSProjFrame::EDirectionSide,
                exactJensenDipoleSampler, itsSampler.get()));
#else
        std::vector<std::pair<Float, const TangentSampler2D*> > perp;
        perp.push_back(std::make_pair(smallLengthWeight,
                new FwdDipSmallLengthSamplerPerpToDir(m_sigmaS, m_g, 1)));
        perp.push_back(std::make_pair(jensenWeight,
                exactJensenDipoleSampler));

        std::vector<std::pair<Float, const TangentSampler2D*> > along;
        along.push_back(std::make_pair(smallLengthWeight,
                new FwdDipSmallLengthSamplerAlongDir(m_sigmaS, m_g, 1)));
        along.push_back(std::make_pair(jensenWeight,
                exactJensenDipoleSampler));

        ref<TangentSampler2D> smallLengthSampler_perp =
                new MISTangentSampler2D(perp);
        ref<TangentSampler2D> smallLengthSampler_along =
                new MISTangentSampler2D(along);



# if 1
        registerSampler(1.0, new ProjSurfaceSampler(
                DSSProjFrame::EDirectionDirection,
                smallLengthSampler_perp, itsSampler.get()));
# endif
# if 1
        /* Double the weight because this is the most 'natural' sampler 
         * for truly near-planar boundaries */
        registerSampler(2.0, new ProjSurfaceSampler(
                DSSProjFrame::EDirectionOut,
                smallLengthSampler_along, itsSampler.get()));
# endif
# if 1
        registerSampler(1.0, new ProjSurfaceSampler(
                DSSProjFrame::EDirectionSide,
                smallLengthSampler_along, itsSampler.get()));
# endif

        if (m_dipoles[0]->useRayDirSurfSampler()) { // quick and dirty, not the best place to store this
            bool bidirectional = m_dipoles[0]->useBidirectionalrayDirSurfSampler(); // quick and dirty, not the best place to store this
            const std::string &rdsSampStrategy = m_dipoles[0]->getRayDirSurfSamplerStrategy(); // quick and dirty, not the best place to store this
            registerSampler(1.0, new RayDirectionSurfaceSampler(sigmaTr, p_spectrum, 0, bidirectional, rdsSampStrategy, itsSampler));
            //registerSampler(1.0, new RayDirectionSurfaceSampler(sigmaTr, p_spectrum, ShadowEpsilon, bidirectional, rdsSampStrategy, itsSampler));
        }
#endif /* MTS_FWDDIP_DEBUG_USE_SIMPLE_JENSEN_SURFACE_SAMPLER*/
    }

    normalizeSamplers();
}

MTS_IMPLEMENT_CLASS_TS(FwdDip, false, DirectSamplingSubsurface);
MTS_IMPLEMENT_CLASS(FwdDipSmallLengthSamplerPerpToDir, false, TangentSampler2D);
MTS_IMPLEMENT_CLASS(FwdDipSmallLengthSamplerAlongDir, false, TangentSampler2D);
MTS_IMPLEMENT_CLASS(FwdDipSmallLengthRadialSampler2D, false, Sampler1D);
MTS_IMPLEMENT_CLASS(RayDirectionSurfaceSampler, false, SurfaceSampler);
MTS_EXPORT_PLUGIN(FwdDip, "Forward scattering dipole model based on a "
        "functional integral approximation");
MTS_NAMESPACE_END
