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
            cout << "Syntax: mtsutil colorFromMonochrome <redChannel.exr> <greenChannel.exr> <blueChannel.exr> <output.exr>" << endl;
            return -1;
        }
        ref<FileStream> rFile   = new FileStream(argv[1], FileStream::EReadOnly);
        ref<FileStream> gFile   = new FileStream(argv[2], FileStream::EReadOnly);
        ref<FileStream> bFile   = new FileStream(argv[3], FileStream::EReadOnly);
        ref<FileStream> outFile = new FileStream(argv[4], FileStream::ETruncReadWrite);

        ref<Bitmap> rBitmap = new Bitmap(Bitmap::EOpenEXR, rFile);
        ref<Bitmap> gBitmap = new Bitmap(Bitmap::EOpenEXR, gFile);
        ref<Bitmap> bBitmap = new Bitmap(Bitmap::EOpenEXR, bFile);

        /* A few sanity checks */
        if (rBitmap->getPixelFormat() != gBitmap->getPixelFormat() || gBitmap->getPixelFormat() != bBitmap->getPixelFormat())
            Log(EError, "Error: Input bitmaps have a different pixel format!");
        if (rBitmap->getComponentFormat() != gBitmap->getComponentFormat() || gBitmap->getComponentFormat() != bBitmap->getComponentFormat())
            Log(EError, "Error: Input bitmaps have a different component format!");
        if (rBitmap->getSize() != gBitmap->getSize() || gBitmap->getSize() != bBitmap->getSize())
            Log(EError, "Error: Input bitmaps have a different size!");

        ref<Bitmap> outBitmap = new Bitmap(rBitmap->getPixelFormat(),
                rBitmap->getComponentFormat(), rBitmap->getSize());

        size_t nPixels =
            (size_t) rBitmap->getSize().x *
            (size_t) rBitmap->getSize().y;
        size_t nChannels = rBitmap->getChannelCount();
        if (nChannels != 3)
            Log(EError, "Error: I currently only work with R,G,B images, but got an image with %zu channels!", nChannels);

        switch (rBitmap->getComponentFormat()) {
            case Bitmap::EFloat16: {
                    half *rData = rBitmap->getFloat16Data();
                    half *gData = gBitmap->getFloat16Data();
                    half *bData = bBitmap->getFloat16Data();
                    half *outData = outBitmap->getFloat16Data();
                    for (size_t i=0; i<nPixels; ++i) {
                        *outData++ = *(rData+0);
                        *outData++ = *(gData+1);
                        *outData++ = *(bData+2);
                        rData += 3;
                        gData += 3;
                        bData += 3;
                    }
                }
                break;

            case Bitmap::EFloat32: {
                    float *rData = rBitmap->getFloat32Data();
                    float *gData = gBitmap->getFloat32Data();
                    float *bData = bBitmap->getFloat32Data();
                    float *outData = outBitmap->getFloat32Data();
                    for (size_t i=0; i<nPixels; ++i) {
                        *outData++ = *(rData+0);
                        *outData++ = *(gData+1);
                        *outData++ = *(bData+2);
                        rData += 3;
                        gData += 3;
                        bData += 3;
                    }
                }
                break;

            case Bitmap::EUInt32: {
                    uint32_t *rData = rBitmap->getUInt32Data();
                    uint32_t *gData = gBitmap->getUInt32Data();
                    uint32_t *bData = bBitmap->getUInt32Data();
                    uint32_t *outData = outBitmap->getUInt32Data();
                    for (size_t i=0; i<nPixels; ++i) {
                        *outData++ = *(rData+0);
                        *outData++ = *(gData+1);
                        *outData++ = *(bData+2);
                        rData += 3;
                        gData += 3;
                        bData += 3;
                    }
                }
                break;

            default:
                Log(EError, "Unsupported component format!");
        }

        outBitmap->write(Bitmap::EOpenEXR, outFile);
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AddImages, "Generate a color image from monochrome channel images")
MTS_NAMESPACE_END
