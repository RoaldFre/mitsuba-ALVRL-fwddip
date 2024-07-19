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

class Downsample : public Utility {
public:
    int run(int argc, char **argv) {
        if (argc != 4) {
            cout << "Create an NxN-downsampled image" << endl;
            cout << "Syntax: mtsutil downsample <input.exr> <output.exr> <N>" << endl;
            return -1;
        }
        ref<FileStream> inFile = new FileStream(argv[1], FileStream::EReadOnly);
        ref<FileStream> outFile  = new FileStream(argv[2], FileStream::ETruncReadWrite);

        char *end_ptr = NULL;
        int N = strtol(argv[3], &end_ptr, 10);
        if (*end_ptr != '\0')
            SLog(EError, "Could not parse N from '%s'", argv[3]);

        ref<Bitmap> inBitmap = new Bitmap(Bitmap::EOpenEXR, inFile);

        Vector2i inSize = inBitmap->getSize();
        size_t inWidth = inSize.x;
        size_t outWidth = inSize.x / N;
        size_t outHeight = inSize.y / N;
        size_t nChannels = inBitmap->getChannelCount();

        ref<Bitmap> outBitmap = new Bitmap(inBitmap->getPixelFormat(), inBitmap->getComponentFormat(), Vector2i(outWidth, outHeight));

        switch (outBitmap->getComponentFormat()) {
            case Bitmap::EFloat16: {
                    half *inData = inBitmap->getFloat16Data();
                    half *outData = outBitmap->getFloat16Data();
                    for (size_t out_i=0; out_i<outHeight; ++out_i) {
                        for (size_t out_j=0; out_j<outWidth; ++out_j) {
                            for (size_t c = 0; c < nChannels; ++c) {
                                double sum = 0; // compute average in high precision
                                for (int u = 0; u < N; ++u) {
                                    for (int v = 0; v < N; ++v) {
                                        sum += inData[nChannels*(inWidth*(N*out_i+u) + N*out_j+v) + c];
                                    }
                                }
                                outData[nChannels*(outWidth*out_i + out_j) + c] = sum / (N*N);
                            }
                        }
                    }
                }
                break;

            case Bitmap::EFloat32: {
                    float *inData = inBitmap->getFloat32Data();
                    float *outData = outBitmap->getFloat32Data();
                    for (size_t out_i=0; out_i<outHeight; ++out_i) {
                        for (size_t out_j=0; out_j<outWidth; ++out_j) {
                            for (size_t c = 0; c < nChannels; ++c) {
                                double sum = 0; // compute average in high precision
                                for (int u = 0; u < N; ++u) {
                                    for (int v = 0; v < N; ++v) {
                                        sum += inData[nChannels*(inWidth*(N*out_i+u) + N*out_j+v) + c];
                                    }
                                }
                                outData[nChannels*(outWidth*out_i + out_j) + c] = sum / (N*N);
                            }
                        }
                    }
                }
                break;

            case Bitmap::EUInt32:
                Log(EError, "Cannot average integer bitmaps!");
                break;

            default:
                Log(EError, "Unsupported component format!");
        }

        outBitmap->write(Bitmap::EOpenEXR, outFile);
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(Downsample, "Create an NxN-downsampled image")
MTS_NAMESPACE_END
