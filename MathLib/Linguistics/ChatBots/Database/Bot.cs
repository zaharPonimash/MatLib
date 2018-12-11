using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MatLib;
namespace MatLib.Linguistics.ChatBots.Database
{
    public class Bot
    {
        public string input;
        public Data data;
        public Bot()
        {
            data = new Data();
        }

        bool ChecksymbolsSpecific(char sym)
        {
            for (int k = 0; k < Linguistic.Dictionary.symbolsSpecific.Length; k++)
                if (sym == Linguistic.Dictionary.symbolsSpecific[k] && sym != (@"\")[0])
                    return true;
            return false;
        }
        string PreprocText(string text)
        {
            string res = text.ToLower();
            res = res.Replace("\r\n", " ");
            for (int i = 0; i < res.Length - 1; i++)
            {
                if (ChecksymbolsSpecific(res[i]))
                {
                    for (int j = i + 1; j < res.Length; j++)
                        if (!ChecksymbolsSpecific(res[j]) || j == res.Length - 1)
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
        string ToText(Word[] words)
        {
            string text = "";
            for (int i = 0; i < words.Length; i++)
                text += words[i].word + " ";

            return text;
            
        }
        public string GetAnswer(string text)
        {
            input = text;
            input = PreprocText(input);

            return ToText(data.GetAnswer(input));
        }
        public void Load(string path)
        {
            data.Load(path);
        }
        public void Save(string path)
        {
            data.Save(path);
        }
        public void GetDatabase(string text)
        {
            text = PreprocText(text);
            data.GetDatabase(text);
        }
    }
}
