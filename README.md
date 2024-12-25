Mitsuba — Forward Scattering Dipole Implementation
==================================================

**NOTE: If you stumbled upon this repo looking for an implementation of
Adaptive Lightslice for Virtual Ray Lights, checkout the ALVRL branch!**


## About the Forward Scattering Dipole

This repository implements the Forward Scattering Dipole model within the
Mitsuba renderer. This model was first published in the SIGGRAPH 2017 paper
'A Forward Scattering Dipole Model from a Functional Integral Approximation'
(Frederickx and Dutré). The code in this repository has subsequently been
improved to reflect the updated model as described in Roald Frederickx's
PhD thesis 'Functional Integrals for a Forward Scattering Dipole Model'
(KU Leuven, 2025).

Project page: http://graphics.cs.kuleuven.be/publications/phdFrederickx/

Citation:

    @phdthesis{FrederickxPhD,
        title   = {Functional Integrals for a Forward Scattering Dipole Model},
        school  = {KU Leuven},
        author  = {Frederickx, Roald},
        year    = {2025},
        month   = {January}
    }


## Dependencies

All the usual dependencies of Mitsuba with an additional dependency on
GSL (the Gnu Scientific Library).
See the [official Mitsuba documentation](http://mitsuba-renderer.org/docs.html)
for more information.

On an Ubuntu machine (e.g. 20.04 LTS), this should get you close:

     sudo apt install build-essential g++ cmake scons libgsl-dev libopenexr-dev libboost-all-dev libeigen3-dev libjpeg-dev libpng-dev libxerces-c-dev libfftw3-dev mesa-common-dev libgl1-mesa-dev libglewmx-dev xmlstarlet libnlopt-dev libarmadillo-dev libxml2-dev python3-collada libcollada-dom2.4-dp-dev qtbase5-dev qttools5-dev libxxf86vm-dev

If there are compile errors about not finding the correct Qt header files and
you can live without the GUI, then you can simply force `hasQt = False` in
build/SConscript.configure to avoid building the GUI. At the time of writing,
the code builds correctly on Ubuntu 20.04.6 LTS, when compilation of the GUI
has been disabled in this manner.


## How to Build

Get the repository:

    git clone https://github.com/roaldfre/mitsuba-ALVRL-fwddip

The forward scattering dipole code lives in the fwddip branch.

    git checkout fwddip

Select a build configuration (currently, only the gcc profile for double
precision is fully supported)
   
    cp build/config-linux-gcc-double.py config.py

Alternatively, if you want spectral rendering (e.g. to render a realistic skin
material), you can choose a spectral build profile:

    cp build/config-linux-gcc-spectral-double.py config.py

Once you have chosen a build profile and the dependencies are installed
(see above), start the build process (e.g. assuming eight cores)

    scons -j 8

If the build finished successfully, you can add the proper directories to the
`$PATH` (e.g. assuming sh/bash shell) by sourcing setpath.sh as follows (note
the leading period)

    . setpath.sh

You should now be able to run the `mitsuba` command (and several others,
such as `mtsutil`, and `mtsgui` if the GUI was built as well).


## Using the Forward Scattering Dipole Model

You can find documented example scenes on the project page
http://graphics.cs.kuleuven.be/publications/phdFrederickx/
or at https://doi.org/10.48804/8AX5X5

The forward scattering dipole model is a subsurface scattering model living
under the `fwddip` name. A typical configuration snippet will look like this
(see the example scenes above for more documented options):

```xml
<!-- Define the subsurface scattering model -->
<subsurface type="fwddip" id="medium_fwddip">
    <rgb name="sigmaS" value="YOUR_VALUE_HERE"/>
    <rgb name="sigmaA" value="YOUR_VALUE_HERE"/>
    <float name="g" value="YOUR_VALUE_HERE"/>

    <!-- The forward scattering dipole model itself is index-matched and
         it gets coupled to an explicit BSDF on the geometry below. -->
    <float name="intIOR" value="1"/>
    <float name="extIOR" value="1"/>

    <!-- Maximum length of an internal reflection chain within the
         medium. Setting this to 0 avoids internal reflection completely.
         Setting this too high can increase the variance and create
         energy conservation issues. -->
    <integer name="maxInternalReflections" value="2"/>

    <!-- Number of tentative incoming surface positions to generate for
         Resampled Importance Sampling ('RIS'; or 'SIR': Sample Importance
         Resampling). For the forward scattering dipole model, a value of 4
         is a good overall choice, especially when coupled with a robust
         pixel estimator or post-processing denoising step to handle
         occasional fireflies due to outlier samples at relatively low sample
         counts. See also SIRnonSurfaceOversamplingFactor. -->
    <integer name="numSIRsurface" value="4"/>

    <!-- Additional factor of Resampled Importance Sampling: number of
         proposals to generate for the directions and path length for each
         tentatively sampled incoming surface position. This can be tuned
         separately because those sampling steps are typically cheaper
         than sampling an incoming surface position, which uses ray tracing
         to project onto the geometry. See also numSIRsurface. -->
    <integer name="SIRnonSurfaceOversamplingFactor" value="1"/>
</subsurface>

<!-- Define the geometry that uses the subsurface scattering model -->
<shape type="obj">
    <string name="filename" value="MY_AWESOME_OBJECT.obj"/>

    <!-- Boundary conditions for the forward scattering dipole are given
         through an explicit BSDF, in contrast to the standard subsurface
         scattering model in Mitsuba. -->
    <bsdf type="dielectric">
        <float name="intIOR" value="1.6"/>
        <float name="extIOR" value="1.0"/>
    </bsdf>

    <!-- Use the forward scattering dipole subsurface model defined above -->
    <ref id="medium_fwddip"/>
</shape>

<!-- Define the overall integrator -->
<integrator type="volpath">
    <!-- Maximum number of subsurface scattering model interactions in a
         path. Setting this to 1 avoids self-illumination (or
         inter-illumination) of subsurface scattering materials. Setting this
         too high can increase the variance and create energy conservation
         issues. -->
    <integer name="maxSubsurfInteractions" value="2"/>
</integrator>
```



## Simple Robust Pixel Estimator

The tonemap utility (`mtsutil tonemap` command) can be used to merge multiple
noisy renders of the same scene (with independent random number seeds) into a
single, robust estimate.
For example:

    mitsuba scene.xml -o noisy1.exr
    mitsuba scene.xml -o noisy2.exr
    mitsuba scene.xml -o noisy3.exr
    mtsutil tonemap -M -R 3 -o robust.png noisy1.exr noisy2.exr noisy3.exr

This produces a denoised `robust.exr` and tonemapped `robust.png` image. You
can use more bins (e.g. 10 with -R 10) and input arbitrarily more noisy renders
to the tonemap command (ideally a multiple of the number of bins).

Note that this robust averaging requires *independent* random numbers for each
individual noisy render, which currently is not the case on Windows machines
(which use a fixed seed and will generate identical noisy renders). See
src/libcore/random.cpp in the Mitsuba source code and make the appropriate
changes when using Windows.



## Additional Subsurface Scattering Implementations

This branch also contains implementations for:
* a classic Jensen dipole without irradiance caching (`uncacheddipole`)
* the Dual-Beam 3D Searchlight BSSRDF of d'Eon (`dEonDualBeam`)
* the Directional Dipole Model for Subsurface Scattering of Frisvad et al. (`dirpole`)

