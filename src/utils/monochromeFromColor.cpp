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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/util.h>

MTS_NAMESPACE_BEGIN

class AddImages : public Utility {
public:
    int run(int argc, char **argv) {
        if (argc != 5) {
            cout << "Create a colored image from monochrome images that represent the R,G,B channels" << endl;
            cout << "Syntax: mtsutil monochromeFromColor <input.exr> <redChannelOutput.exr> <greenChannelOutput.exr> <blueChannelOutput.exr>" << endl;
            return -1;
        }
        ref<FileStream> inFile = new FileStream(argv[1], FileStream::EReadOnly);
        ref<FileStream> rFile  = new FileStream(argv[2], FileStream::ETruncReadWrite);
        ref<FileStream> gFile  = new FileStream(argv[3], FileStream::ETruncReadWrite);
        ref<FileStream> bFile  = new FileStream(argv[4], FileStream::ETruncReadWrite);

        ref<Bitmap> inBitmap = new Bitmap(Bitmap::EOpenEXR, inFile);

        ref<Bitmap> rBitmap = new Bitmap(inBitmap->getPixelFormat(), inBitmap->getComponentFormat(), inBitmap->getSize());
        ref<Bitmap> gBitmap = new Bitmap(inBitmap->getPixelFormat(), inBitmap->getComponentFormat(), inBitmap->getSize());
        ref<Bitmap> bBitmap = new Bitmap(inBitmap->getPixelFormat(), inBitmap->getComponentFormat(), inBitmap->getSize());

        size_t nPixels =
            (size_t) inBitmap->getSize().x *
            (size_t) inBitmap->getSize().y;
        size_t nChannels = inBitmap->getChannelCount();
        if (nChannels != 3)
            Log(EError, "Error: I currently only work with R,G,B images, but got an image with %zu channels!", nChannels);

        switch (rBitmap->getComponentFormat()) {
            case Bitmap::EFloat16: {
                    half *inData = inBitmap->getFloat16Data();
                    half *rData = rBitmap->getFloat16Data();
                    half *gData = gBitmap->getFloat16Data();
                    half *bData = bBitmap->getFloat16Data();
                    for (size_t i=0; i<nPixels; ++i) {
                        half r = *inData++;
                        half g = *inData++;
                        half b = *inData++;
                        for (int j = 0; j < 3; ++j) {
                            *rData++ = r;
                            *gData++ = g;
                            *bData++ = b;
                        }
                    }
                }
                break;

            case Bitmap::EFloat32: {
                    float *inData = inBitmap->getFloat32Data();
                    float *rData = rBitmap->getFloat32Data();
                    float *gData = gBitmap->getFloat32Data();
                    float *bData = bBitmap->getFloat32Data();
                    for (size_t i=0; i<nPixels; ++i) {
                        float r = *inData++;
                        float g = *inData++;
                        float b = *inData++;
                        for (int j = 0; j < 3; ++j) {
                            *rData++ = r;
                            *gData++ = g;
                            *bData++ = b;
                        }
                    }
                }
                break;

            case Bitmap::EUInt32: {
                    uint32_t *inData = inBitmap->getUInt32Data();
                    uint32_t *rData = rBitmap->getUInt32Data();
                    uint32_t *gData = gBitmap->getUInt32Data();
                    uint32_t *bData = bBitmap->getUInt32Data();
                    for (size_t i=0; i<nPixels; ++i) {
                        uint32_t r = *inData++;
                        uint32_t g = *inData++;
                        uint32_t b = *inData++;
                        for (int j = 0; j < 3; ++j) {
                            *rData++ = r;
                            *gData++ = g;
                            *bData++ = b;
                        }
                    }
                }
                break;

            default:
                Log(EError, "Unsupported component format!");
        }

        rBitmap->write(Bitmap::EOpenEXR, rFile);
        gBitmap->write(Bitmap::EOpenEXR, gFile);
        bBitmap->write(Bitmap::EOpenEXR, bFile);
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AddImages, "Generate monochrome channel images from a color image")
MTS_NAMESPACE_END
