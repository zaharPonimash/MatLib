using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Linguistics
{
    public static class Linguistic
    {
            public class Convertor
            {
                public static Vector TextToVectorIndexes(string text)
                {
                    Vector res = new Vector(text.Length);
                    for (int i = 0; i < text.Length; i++)
                        res[i] = Dictionary.getIndex(text[i]);
                    return res;
                }
                public static string VectorIndexesToText(Vector indexes)
                {
                    string res = "";
                    for (int i = 0; i < indexes.Length; i++)
                    {
                        var a = matlib.Round(indexes[i]);
                        if (a < 0) a = 0;
                        if (a >= Dictionary.symbols.Length) a = Dictionary.symbols.Length - 1;
                        res += Dictionary.symbols[a];

                    }
                    return res;
                }
                public static List<Word> GetSpecificWords(string text)
                {
                    List<Word> result = new List<Word>();
                    List<string> words = text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).ToList();
                    for (int i = 0; i < words.Count - 1; i++)
                        for (int j = i + 1; j < words.Count; j++)
                            if (words[i] == words[j]) { words.RemoveAt(j); j--; }

                    for (int i = 0; i < words.Count; i++)
                        result.Add(new Word(words[i], i));

                    return result;

                }
                public static List<Word> GetWords(string text)
                {
                    List<Word> result = new List<Word>();
                    List<string> words = text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).ToList();

                    for (int i = 0; i < words.Count; i++)
                        result.Add(new Word(words[i], i));

                    return result;

                }
                public static int SearchWordIndex(List<Word> words, Word word)
                {
                    for (int i = 0; i < words.Count; i++)
                        if (words[i].word == word.word) return i;

                    return -1;
                }
                public static Matrix WordsToMatrix(List<Word> words, List<Word> wordsSpec)
                {

                    Matrix result = new Matrix(wordsSpec.Count, words.Count);
                    for (int i = 0; i < words.Count; i++)
                    {
                        if (words[i].index != -1)
                        {
                            var a = SearchWordIndex(wordsSpec, words[i]);

                            if (a == -1)
                                return null;

                            result[i, wordsSpec[a].index] = 1.0;
                        }
                    }
                    return result;          
                }
                public static Vector[] WordsToVectors(List<Word> words, int numWords)
                {

                    Vector[] result = new Vector[words.Count];
                    for (int i = 0; i < result.Length; i++)
                        result[i] = new Vector(numWords);
                    
                    for (int i = 0; i < words.Count; i++)
                        result[i][words[i].index] = 1.0;

                    return result;
                }
            }
            public class Dictionary
            {
                public const string symbols = " абвгдеёжзийклмнопрстуфхцчшщъыьэюя)(.,?!-;:";
                public const string symbolsSpecific = ")(.,?!-;:/|\n\t\r+_><" + @"\";
                public static int getIndex(char symbol)
                {
                    return symbols.IndexOf(symbol);
                }

            }
            public static int CountAllWords(string text)
            {
                int ctr = 0;
                for (int i = 0; i < text.Length; i++)
                    if (text[i] == ' ') ctr++;
                
                return ctr;
            }
            public static int CountSpecificWords(string text)
            {
                List<string> words = text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).ToList();
                for (int i = 0; i < words.Count - 1; i++)
                    for (int j = i + 1; j < words.Count; j++)
                        if (words[i] == words[j]) { words.RemoveAt(j); j--; }

                return words.Count;
            }
    }
}
