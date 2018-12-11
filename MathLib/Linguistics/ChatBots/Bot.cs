using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using NeuralNetwork;
using MatLib;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace MatLib.Linguistics.ChatBots
{
    public class Bot
    {
        NeuralNet net;
        string inputText;
        Matrix inputMatr;
        List<Word> specificWords;

        public void LoadText(string path)
        {
            string text;
            using (StreamReader stream = new StreamReader(path, Encoding.Default))
            {
                text = stream.ReadToEnd();
            }
            specificWords = Linguistic.Convertor.GetSpecificWords(PreprocText(text));
            specificWords.Add(new Word(" ", specificWords.Count));
        }

        public void SaveDictionary(string path)
        {
            using (FileStream stream = new FileStream(path, FileMode.Create))
            {
                BinaryFormatter serialize = new BinaryFormatter();
                serialize.Serialize(stream, specificWords);
            }
        }
        public void LoadDictionary(string path)
        {
            using (FileStream stream = new FileStream(path, FileMode.Open))
            {
                BinaryFormatter deserialize = new BinaryFormatter();
                specificWords = (List<Word>)deserialize.Deserialize(stream);
            }
        }
        public Bot(string[] net)
        {
            this.net = new NeuralNet(net);
        }
        public Bot()
        {
            
        }
        public void AutoArchitecture()
        {
            GetTrainSample();
            
            net = new NeuralNet(new string[] { 
                "input:" + specificWords.Count + ":" + x[0].height + ":1:1",
                 
                "conv:" + specificWords.Count + ":1:10",
                "sin",

                "conv:1:1:878",
                //"Tensor3ToTensor3:" + specificWords.Count + ":" + x[0].height + ":1"

               "deepToMatrix"
                
            
            });
            // net.LoadWeights("weights");
            //string s = net.GetInfoAboutOutputs();
            /*
            List<Word>[] words1 = new List<Word>[x.Length];
            List<Word>[] words2 = new List<Word>[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
              
            words1[i] = Linguistic.Convertor.GetWords(x[i]);
            words2[i] = Linguistic.Convertor.GetWords(y[i]);
            }
            Matrix[] inp = new Matrix[x.Length];
            Matrix[] outp = new Matrix[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                inp[i] = Linguistic.Convertor.WordsToMatrix(words1[i], specificWords);
                outp[i] = Linguistic.Convertor.WordsToMatrix(words2[i], specificWords);
            }
            for (int i = 0; i < 1000; i++)
            {
                for (int j = 0; j < x.Length; j++)
                {
                    //Tensor4 answ = new Tensor4()
                    net.TrainWithTeach(inp[j].ToTensor4(), outp[j].ToTensor4(), 0.001, 0.5, 1e-4, NeuralNet.Optimizer.Adam);
                }
            }
            //*/
        }
        void GetTrainSample(string path = "dial.txt")
        {
            List<string> x = new List<string>();
            List<string> y = new List<string>();
            string text;
            using (StreamReader stream = new StreamReader(path, Encoding.Default))
            {
                text = stream.ReadToEnd();
            }
            var answ = text.Split(new string[]{"\\"}, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 0; i < answ.Length - 1; i+=2)
            {
                x.Add(PreprocText(answ[i]));
                y.Add(PreprocText(answ[i + 1]));
            }
            List<Word>[] words1 = new List<Word>[x.Count];
            List<Word>[] words2 = new List<Word>[x.Count];
            for (int i = 0; i < x.Count; i++)
            {
                words1[i] = Linguistic.Convertor.GetWords(x[i]);
                words2[i] = Linguistic.Convertor.GetWords(y[i]);
            }

            int maxWords = words1[0].Count;
            for (int i = 1; i < words1.Length; i++)
                if (maxWords < words1[i].Count) maxWords = words1[i].Count;
            for (int i = 0; i < words2.Length; i++)
                if (maxWords < words2[i].Count) maxWords = words2[i].Count;
            maxWords++;
            for (int i = 0; i < x.Count; i++)
            {
                if(words1[i].Count != maxWords)
                {
                    int length = maxWords - words1[i].Count;
                    for (int j = 0; j < length; j++)
                    {
                        words1[i].Add(new Word("", -1));
                    }
                }
                if (words2[i].Count != maxWords)
                {
                    int length = maxWords - words2[i].Count;
                    for (int j = 0; j < length; j++)
                    {
                        words2[i].Add(new Word("", -1));
                    }
                }
                #region
                /*
                if(words1[i].Count > words2[i].Count)
                {
                    int l = words1[i].Count - words2[i].Count;
                    for (int j = 0; j < l; j++)
                        words1[i].RemoveAt(words1[i].Count - 1);
                }
                else if(words1[i].Count < words2[i].Count)
                {
                    int l = -words1[i].Count + words2[i].Count;
                    for (int j = 0; j < l; j++)
                        words2[i].RemoveAt(words2[i].Count - 1);
                }
                //*/
                #endregion
            }
            for (int i = 0; i < x.Count; i++)
            {
               this.x.Add(Linguistic.Convertor.WordsToMatrix(words1[i], specificWords).ToTensor4());
               this.y.Add(Linguistic.Convertor.WordsToMatrix(words2[i], specificWords).ToTensor4());
            }
        }
        List<Tensor4> x = new List<Tensor4>();
        List<Tensor4> y = new List<Tensor4>();
        bool ChecksymbolsSpecific(char sym)
        {
            for (int k = 0; k < Linguistic.Dictionary.symbolsSpecific.Length; k++)
                if (sym == Linguistic.Dictionary.symbolsSpecific[k])
                    return true;
            return false;
        }
        string PreprocText(string text)
        {
            string res = text.ToLower();
            res = res.Replace("\r\n", " ");
            for(int i = 0; i < res.Length - 1; i++)
            {
                if(ChecksymbolsSpecific(res[i]))
                {
                    for(int j = i + 1; j < res.Length; j++)
                        if(!ChecksymbolsSpecific(res[j]) || j == res.Length - 1)
                        {
                            if (j == res.Length - 1)
                                res = res.Remove(i, j - i + 1);
                            else if (res[j] != ' ') res = res.Remove(i, j - i).Insert(i, " ");
                            else res = res.Remove(i, j - i + 1).Insert(i, " ");
                            break;
                        }
                }
            }
            if (res.Length == 0) return null;
            if (ChecksymbolsSpecific(res[res.Length - 1]))
                res = res.Remove(res.Length - 1, 1);
            return res;
               // text.ToLower()..Trim(Linguistic.Dictionary.symbolsSpecific.ToCharArray());
        }
        public string IndexToWord(int index)
        {
            return specificWords.ElementAt(index).word;
        }
        public void Teach()
        {
            for (int i = 0; i < 10; i++)
            {
                for (int j = 0; j < x.Count; j++)
                {
                    //Tensor4 answ = new Tensor4()
                    net.TrainWithTeach(x[j], y[j], 0.0014, 0.5, 1e-4, NeuralNet.Optimizer.Adam);
                }
            }
            //net.SaveWeights("weights");
        }
        public string GetAnswer(string text)
        {
            string output = "";

            this.inputText = PreprocText(text);
            var a = Linguistic.Convertor.GetWords(this.inputText);
            if(a.Count != net.layers[0].input.height)
            {
                int length = net.layers[0].input.height - a.Count;
                for (int i = 0; i < length; i++)
                    a.Add(new Word("", -1));    
            }

            inputMatr = Linguistic.Convertor.WordsToMatrix(a, specificWords);

            net.Calculation(inputMatr);

            var y = net.output.ToMatrix(0, 0);
            Vector indexes = new Vector(y.height);

            for (int i = 0; i < y.height; i++)
            {
                var vect = y.GetVector(i, 1);
                if (vect.elements.Max() > 0.3)
                    indexes[i] = vect.IndexMax();
                else indexes[i] = -1.0;
            }

            for (int i = 0; i < indexes.Length; i++)
            {
                if (indexes[i] != -1)
                {
                    if (i != indexes.Length - 1)
                        output += IndexToWord((int)indexes[i]) + " ";
                    else output += IndexToWord((int)indexes[i]);
                }
            }

            return output;
        }
        public void Load(string settings)
        {

        }
        public void Save(string path)
        {

        }
    }
}
