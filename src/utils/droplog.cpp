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

class DropLog : public Utility {
public:
    int run(int argc, char **argv) {
        if (argc != 2) {
            cout << "Drop the logfile from the EXIF of a given EXR image" << endl;
            cout << "Syntax: mtsutil droplog <target.exr>" << endl;
            return -1;
        }
        std::string outFileNameBase = boost::filesystem::change_extension(argv[1], "").string();
        std::ostringstream oss;
        oss << outFileNameBase << "_droplog.exr";
        ref<FileStream> inFile  = new FileStream(argv[1], FileStream::EReadOnly);
        ref<FileStream> outFile = new FileStream(oss.str().c_str(), FileStream::ETruncReadWrite);

        Log(EInfo, "Infile:  %s", argv[1]);
        Log(EInfo, "Outfile: %s", oss.str().c_str());

        ref<Bitmap> inBitmap = new Bitmap(Bitmap::EOpenEXR, inFile);
        ref<Bitmap> outBitmap = new Bitmap(inBitmap->getPixelFormat(),
                inBitmap->getComponentFormat(), inBitmap->getSize());

        size_t nEntries =
            (size_t) inBitmap->getSize().x *
            (size_t) inBitmap->getSize().y *
            inBitmap->getChannelCount();

        switch (inBitmap->getComponentFormat()) {
            case Bitmap::EFloat16: {
                    half *inData = inBitmap->getFloat16Data();
                    half *outData = outBitmap->getFloat16Data();
                    for (size_t i=0; i<nEntries; ++i)
                        *outData++ = (half) (*inData++);
                }
                break;

            case Bitmap::EFloat32: {
                    float *inData = inBitmap->getFloat32Data();
                    float *outData = outBitmap->getFloat32Data();
                    for (size_t i=0; i<nEntries; ++i)
                        *outData++ = (float) (Float) (*inData++);
                }
                break;

            case Bitmap::EUInt32: {
                    uint32_t *inData = inBitmap->getUInt32Data();
                    uint32_t *outData = outBitmap->getUInt32Data();
                    for (size_t i=0; i<nEntries; ++i)
                        *outData++ = (uint32_t) (*inData++);
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

MTS_EXPORT_UTILITY(DropLog, "Drop the logfile from the EXIF of a given EXR image")
MTS_NAMESPACE_END
