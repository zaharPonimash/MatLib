using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MatLib;
using MatLib.Linguistics;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
namespace MatLib.Linguistics.ChatBots.Database
{
    public class Data
    {
        List<Word[]> questions = new List<Word[]>(), answers = new List<Word[]>();
        public Data()
        {

        }
        public Word[] GetAnswer(string text)
        {
            var words = Linguistic.Convertor.GetWords(text);
            List<int> indexes = new List<int>();
            int max = 0;
            for (int i = 0; i < questions.Count; i++)
            {
                int ctr = 0;
                int l = (words.Count < questions[i].Length) ? words.Count : questions[i].Length;
                for (int j = 0; j < l; j++)
                {
                        if (words[j].word == questions[i][j].word)
                            ctr++;
                }
                if (ctr == words.Count && words.Count == questions[i].Length) indexes.Add(i);
                if(max < ctr) max = ctr;
            }
            if(indexes.Count == 0 || max != 1)
            {
                //throw new Exception("Просто, по приколу вызвал ошибку");
                for (int i = 0; i < questions.Count; i++)
                {
                    int ctr = 0;
                    int l = (words.Count < questions[i].Length) ? words.Count : questions[i].Length;
                    for (int j = 0; j < l; j++)
                    {
                        for (int k = 0; k < l; k++)
                        {
                            if (words[j].word == questions[i][k].word) ctr++;
                        }
                    }
                    if (ctr == max) indexes.Add(i);
                }
            }
            int index = (int)(matlib.random.NextDouble() * indexes.Count);
            try
            {
                return answers[indexes[index]];
            }
            catch { return new Word[] { new Word("") }; }
        }
        public void Load(string path)
        {
            using(FileStream stream = new FileStream(path, FileMode.Open))
            {
                BinaryFormatter bin = new BinaryFormatter();
                questions = (List<Word[]>)bin.Deserialize(stream);
                answers = (List<Word[]>)bin.Deserialize(stream);
            }
        }
        public void Save(string path)
        {
            using(FileStream stream = new FileStream(path, FileMode.Create))
            {
                BinaryFormatter bin = new BinaryFormatter();
                bin.Serialize(stream, questions);
                bin.Serialize(stream, answers);
            }
        }
        public void GetDatabase(string text)
        {
            var pairs = text.Split(new string[]{@"\"}, StringSplitOptions.None);
            if(pairs.Length % 2 != 0) throw new Exception("Error. ");

            for (int i = 0; i < pairs.Length; i+=2)
			{
                questions.Add(Linguistic.Convertor.GetWords(pairs[i]).ToArray());
                answers.Add(Linguistic.Convertor.GetWords(pairs[i + 1]).ToArray());
			}
        }
    }
}
