using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MatLib;

namespace MatLib.Linguistics
{
    [Serializable]
    public class Word
    {
        public string word;
        public int index;
        public double freq;
        public Word(string word, int index = 0)
        {
            this.word = word;
            this.index = index;
        }
        public Vector ToIndexes()
        {
            return Linguistic.Convertor.TextToVectorIndexes(word);
        }
    }
}
