using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Sounds
{
    public class WAVStream
    {
        public int chunkID, fileSize, riffType, fmtID, fmtSize, fmtCode, channels, sampleRate, fmtAvgBPS, fmtBlockAlign, bitDepth, dataID, dataSize;
        public string path;
        public FileMode mode;
        public WAVStream(string path, FileMode mode)
        {
            this.path = path;
            this.mode = mode;
        }
        public Vector Load()
        {

            Stream waveFileStream = File.OpenRead(path);
            BinaryReader reader = new BinaryReader(waveFileStream);

            chunkID = reader.ReadInt32();
            fileSize = reader.ReadInt32();
            riffType = reader.ReadInt32();
            fmtID = reader.ReadInt32();
            fmtSize = reader.ReadInt32();
            fmtCode = reader.ReadInt16();
            channels = reader.ReadInt16();
            sampleRate = reader.ReadInt32();
            fmtAvgBPS = reader.ReadInt32();
            fmtBlockAlign = reader.ReadInt16();
            bitDepth = reader.ReadInt16();

            if (fmtSize == 18)
            {
                // Read any extra values 
                int fmtExtraSize = reader.ReadInt16();
                reader.ReadBytes(fmtExtraSize);
            }

            dataID = reader.ReadInt32();
            dataSize = reader.ReadInt32();

            List<double> fl = new List<double>();

            while (true)
            {
                try
                {
                    fl.Add(reader.ReadInt16() / 32000.0);
                }
                catch
                {
                    break;
                }
            }

            return new Vector(fl);
        }
    }
}
