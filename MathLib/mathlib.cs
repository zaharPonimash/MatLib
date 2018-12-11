using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Drawing;
using ELW.Library.Math;
using ELW.Library.Math.Expressions;
using ELW.Library.Math.Tools;
using ELW.Library.Math.Calculators;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using ZedGraph;

using AI.MathMod.ML.Regression;
namespace MatLib
{
    [Serializable]
    public class matlib
    {
    	static public Random random = new Random();
        public class Function
        {
                public class Regression
                {
                    public static Vector Polym(Vector x, Vector y, int n = 5)
                    {
                        AI.MathMod.Vector vecX = new AI.MathMod.Vector(new Vector(x).elements);
                        AI.MathMod.Vector vecY = new AI.MathMod.Vector(new Vector(y).elements);
                        var lag = new RBFGauss(vecX, vecY, n);
                        
                        return new Vector(lag.param.Vecktor);
                        

                    }
                }
                public class AutoRegression
                {
                    public static Vector ARFIMA(Vector x, Vector y, Vector X, int K = 5, int D = 4)
                    {
                        Vector Y = new Vector(X.Length);
                        for (int i = 0; i < X.Length; i++)
                        {

                            for (int k = 0; k < K; k++)
                            {
                                double multip = 1.0;
                                for (int d = 0; d < k - 1; d++)
                                {
                                  //  multip *= (D - d) * ((k % 2 == 0) ? 0 : 1) * matlib.Function.Interpolation.Lagrange()
                                }
                            }

                        }

                        return Y;

                    }
                }
                static int indexOf(Vector vals, double val)
                {
                    int res = -1;
                    for (int i = 0; i < vals.Length; i++)
                    {
                        if (vals[i] == val) return i;
                    }
                    return res;
                }
                static public double Taylor(double x, double a, Vector derivs)
                {
                    if (derivs.Length < 1) return double.NaN;
                    double res = 0.0;

                    for (int i = 0; i < derivs.Length; i++)
                        res += derivs[i] * Math.Pow(x - a, i) / (double)matlib.Factorial(i);
                    
                    return res;
                }
            
            public class Func
            {
                public class SigmoidalFunc
                {
                    static public double Sigmoida(double x)
                    {
                        return 1.0 / (1.0 + Math.Exp(-x));
                    }
                    static public double HyperbolicTangent(double x)
                    {
                        return (Math.Exp(2.0 * x) - 1.0) / (Math.Exp(2.0 * x) + 1.0);
                    }
                }
                public class PiecewiseDefined
                {
                    public static int Threshold(double x, double a)
                    {
                        if (x >= a) return 1;
                        else return 0;
                    }
                    static public double Relu(double x)
                    {
                        if (x > 0.0) return x;
                        else return 0.0;
                    }
                }
                static public double f(string func, double x)
                {
                    return ToolsHelper.Calculator.Calculate(ToolsHelper.Compiler.Compile(ToolsHelper.Parser.Parse(func)), new List<VariableValue>() { new VariableValue(x, "x") });
                }
                public static Vector FuncToVector(string func, int a, int b, int N = 10000)
                {
                    Vector y = new Vector(N);
                    PreparedExpression PE = ToolsHelper.Parser.Parse(func);
                    CompiledExpression CE = ToolsHelper.Compiler.Compile(PE);
                    List<VariableValue> variables = new List<VariableValue>();
                    double step = (b - a) / (double)N;
                    for (int i = 0; i < N; i++)
                    {
                        variables.Clear();
                        variables.Add(new VariableValue(a + step * i, "x"));
                        y[i] = ToolsHelper.Calculator.Calculate(CE, variables);
                    }
                    return y;
                }
                public static Vector Sin(int counts, double frequency, double fi, double noisy = 0.0)
                {
                    Vector y = new Vector(counts);
                    for (int i = 0; i < counts; i++)
                        y[i] = Math.Sin(2.0 * Math.PI * i / (double)counts * frequency + fi) + noisy * Statistic.Random.Gauss();

                    return y;
                }
                public static Vector Rect(int counts, double frequency, double fi, double noisy)
                {
                    Vector y = new Vector(counts);
                    double T = 1.0 / frequency;
                    int k = 0;
                    for (int i = 0; i < counts; i++)
                    {
                        double x = (double)i / counts + fi;
                        if (x >= k * T)
                            k++;
                        if (x < T * (k - 1) + T / 2.0)
                            y[i] = 1.0;
                        else
                            y[i] = 0.0;
                        y[i] += noisy * Statistic.Random.Gauss();
                    }
                    return y;
                }
            }
            public class FuncDerivative
            {
                static public double Sigmoida(double x)
                {
                    return Math.Exp(-x) / Math.Pow(1.0 + Math.Exp(-x), 2.0);
                }
                static public double HyperbolicTangent(double x)
                {
                    return (4.0 * Math.Exp(2.0 * x)) / Math.Pow(Math.Exp(2.0 * x) + 1.0, 2.0);
                }
                static public double Func(string func, double x)
                {
                    double increment = 1e-8;
                    return (Function.Func.f(func, x + increment) - Function.Func.f(func, x)) / increment;
                }
            }
            public class Optimize
            {
                public class MinimizationObjFunc
                {
                    static CompiledExpression CE;
                    static List<VariableValue> variables;
                    static int LengthSpecifVars;
                    static int LengthMinimizVars;
                    static Vector GRAD;
                    static double incriment;
                    static Vector result;
                    static string[] MinimizeVars;
                    public static Vector MethodConjGrad(string ErrFunc, Vector[] x, Vector y, string[] MinimizVars, string[] SpecifVars, double Err, Vector startValues)
                    {
                        LengthSpecifVars = SpecifVars.Length;
                        LengthMinimizVars = MinimizVars.Length;
                        MinimizeVars = MinimizVars;

                        incriment = 10e-9;

                        result = new Vector(MinimizVars.Length);

                        GRAD = new Vector(MinimizVars.Length);

                        Vector LastGRAD = new Vector(MinimizVars.Length);

                        Vector B = new Vector(MinimizVars.Length);

                        Vector LastB = new Vector(MinimizVars.Length);

                        variables = new List<VariableValue>();

                        PreparedExpression PE = ToolsHelper.Parser.Parse(ErrFunc);

                        CE = ToolsHelper.Compiler.Compile(PE);

                        result = new Vector(startValues);

                        int IterNum = 0;

                        double w = 0.0;

                        do
                        {
                            for (int index = 0; index < x.Length; index++)
                            {
                                variables.Clear();
                                for (int i = 0; i < SpecifVars.Length; i++)
                                    variables.Add(new VariableValue(x[index][i], SpecifVars[i]));

                                variables.Add(new VariableValue(y[index], "y"));

                                calcGRAD();

                                if (IterNum >= 1)
                                {
                                    w = GRAD.EuclidNorm() / LastGRAD.EuclidNorm();
                                    B = GRAD + (LastB * w);
                                }
                                else B = GRAD;
                                result = result - (B * 0.0001);
                                LastB = new Vector(B);
                                LastGRAD = new Vector(GRAD);
                                calcGRAD();
                                variables.RemoveAt(variables.Count - 1);
                                //
                                variables.Clear();
                                for (int i = 0; i < SpecifVars.Length; i++)
                                    variables.Add(new VariableValue(x[index][i], SpecifVars[i]));
                                variables.Add(new VariableValue(y[index], "y"));
                                for (int i = 0; i < MinimizeVars.Length; i++)
                                    variables.Add(new VariableValue(result[i], MinimizeVars[i]));
                            }
                            if (result.IsNaN() || ++IterNum == 100000)
                            {
                                return result;
                            }
                        } while (CalcErr(ErrFunc) > Err);
                        return result;
                    }
                    public static void GeneticAlgorithm(Delegate FitnessFunc)
                    {

                    }
                    private static double CalcErr(string Func)
                    {
                        double res = Math.Abs(ToolsHelper.Calculator.Calculate(ToolsHelper.Compiler.Compile(ToolsHelper.Parser.Parse(Func)), variables));
                        if (!Double.IsNaN(res))
                            return res;
                        else return 100000000;
                        //return Math.Abs(ToolsHelper.Calculator.Calculate(ToolsHelper.Compiler.Compile(ToolsHelper.Parser.Parse(Func)), variables));
                    }
                    static double f(double[] x, string[] name)
                    {
                        int l = variables.Count;
                        for (int i = 0; i < x.Length; i++)
                            variables.Add(new VariableValue(x[i], name[i]));
                        double result = ToolsHelper.Calculator.Calculate(CE, variables);
                        for (int i = 0; i < x.Length; i++)
                            variables.RemoveAt(l);
                        //clear();
                        return result;
                    }
                    /*
                    static void clear()
                    {
                        for (int i = 0; i < LengthMinimizVars; i++)
                            variables.RemoveAt(LengthSpecifVars + 1 + i);
                    }
                    //*/
                    static void calcGRAD()
                    {
                        for (int i = 0; i < GRAD.Length; i++)
                        {
                            result[i] += incriment;
                            double value1 = f(result.elements, MinimizeVars);
                            result[i] -= incriment;
                            double value2 = f(result.elements, MinimizeVars);
                            GRAD[i] = (value1 - value2) / incriment;
                        }
                    }

                }
            }
        }
        public class Statistic
        {
            public class Random
            {
                static public double Gauss()
                {
                    double u, v, s;
                    do
                    {
                        u = 2.0 * random.NextDouble() - 1.0;
                        v = 2.0 * random.NextDouble() - 1.0;
                        s = u * u + v * v;
                    } while (s == 0);
                    return v * Math.Sqrt(Math.Abs(2.0 * Math.Log(s) / s));
                }
            }
        }
        public static double Sqrt(double x, double num = 20)
        {
            double a = 0.0, result = 0.0, n = 1.0 / 2.0;
            double[] array = new double[501];
            for (int i = 0; i <= 500; i++)
                array[i] = Math.Abs(x - Math.Pow(i, 2.0));

            var min = array.Min();
            for (int i = 0; i < array.Length; i++)
                if (array[i] == min)
                {
                    a = Math.Pow(i, 2.0);
                    break;
                }
            for (int k = 0; k < num; k++)
            {
                double m = 1.0;
                for(int j = 0; j < k; j++)
                    m *= n - (double)j;
                result += m * Math.Pow(a, n - k) * Math.Pow(x - a, k) / (double)Factorial(k);
            }
            return result;
        }
        /// <summary>
        /// x!
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static public decimal Factorial(int x)
        {
            if (x >= 0)
            {
                decimal result = 1;
                for (int i = 2; i <= x; i++)
                    result *= i;
                return result;
            }
            else { throw new InvalidOperationException(); }
        }
        static public int gcd(int a, int b)
        {
            if (a < b) { int c = a; a = b; b = c; } 
            return (a % b == 0) ? b : gcd(a % b, b);
        }

        public class Transofm
        {
            public class FFT
            {
                /// <summary>
                /// Вычисление поворачивающего модуля e^(-i*2*PI*k/N)
                /// </summary>
                /// <param name="k"></param>
                /// <param name="N"></param>
                /// <returns></returns>
                private static Complex w1(int k, int N)
                {
                    if (k % N == 0) return 1;
                    double arg = -2 * Math.PI * k / N;
                    return new Complex(Math.Cos(arg), Math.Sin(arg));
                }
                private static Complex w2(int k, int N)
                {
                    if (k % N == 0) return 1;
                    double arg = 2 * Math.PI * k / N;
                    return new Complex(Math.Cos(arg), Math.Sin(arg));
                }
                /// <summary>
                /// Возвращает спектр сигнала
                /// </summary>
                /// <param name="x">Массив значений сигнала. Количество значений должно быть степенью 2</param>
                /// <returns>Массив со значениями спектра сигнала</returns>
                public static Complex[] fft(Complex[] x)
                {
                    Complex[] X;
                    int N = x.Length;
                    if (N == 2)
                    {
                        X = new Complex[2];
                        X[0] = x[0] + x[1];
                        X[1] = x[0] - x[1];
                    }
                    else
                    {
                        Complex[] x_even = new Complex[N / 2];
                        Complex[] x_odd = new Complex[N / 2];
                        for (int i = 0; i < N / 2; i++)
                        {
                            x_even[i] = x[2 * i];
                            x_odd[i] = x[2 * i + 1];
                        }
                        Complex[] X_even = fft(x_even);
                        Complex[] X_odd = fft(x_odd);
                        X = new Complex[N];
                        for (int i = 0; i < N / 2; i++)
                        {
                            X[i] = X_even[i] + w1(i, N) * X_odd[i];
                            X[i + N / 2] = X_even[i] - w1(i, N) * X_odd[i];
                        }
                    }
                    return X;
                }
                public static Complex[] fft(Vector x)
                {
                    Complex[] X;
                    int N = x.Length;
                    if (N == 2)
                    {
                        X = new Complex[2];
                        X[0] = x[0] + x[1];
                        X[1] = x[0] - x[1];
                    }
                    else
                    {
                        Complex[] x_even = new Complex[N / 2];
                        Complex[] x_odd = new Complex[N / 2];
                        for (int i = 0; i < N / 2; i++)
                        {
                            x_even[i] = x[2 * i];
                            x_odd[i] = x[2 * i + 1];
                        }
                        Complex[] X_even = fft(x_even);
                        Complex[] X_odd = fft(x_odd);
                        X = new Complex[N];
                        for (int i = 0; i < N / 2; i++)
                        {
                            X[i] = X_even[i] + w1(i, N) * X_odd[i];
                            X[i + N / 2] = X_even[i] - w1(i, N) * X_odd[i];
                        }
                    }
                    return X;
                }
                public static Complex[] ifft(Complex[] x)
                {
                    Complex[] X;
                    int N = x.Length;
                    if (N == 2)
                    {
                        X = new Complex[2];
                        X[0] = x[0] + x[1];
                        X[1] = x[0] - x[1];
                    }
                    else
                    {
                        Complex[] x_even = new Complex[N / 2];
                        Complex[] x_odd = new Complex[N / 2];
                        for (int i = 0; i < N / 2; i++)
                        {
                            x_even[i] = x[2 * i];
                            x_odd[i] = x[2 * i + 1];
                        }
                        Complex[] X_even = ifft(x_even);
                        Complex[] X_odd = ifft(x_odd);
                        X = new Complex[N];
                        for (int i = 0; i < N / 2; i++)
                        {
                            X[i] = X_even[i] + w2(i, N) * X_odd[i];
                            X[i + N / 2] = X_even[i] - w2(i, N) * X_odd[i];
                        }
                    }
                    return X;
                }
                public static Complex[] ifft(Vector x)
                {
                    Complex[] X;
                    int N = x.Length;
                    if (N == 2)
                    {
                        X = new Complex[2];
                        X[0] = x[0] + x[1];
                        X[1] = x[0] - x[1];
                    }
                    else
                    {
                        Complex[] x_even = new Complex[N / 2];
                        Complex[] x_odd = new Complex[N / 2];
                        for (int i = 0; i < N / 2; i++)
                        {
                            x_even[i] = x[2 * i];
                            x_odd[i] = x[2 * i + 1];
                        }
                        Complex[] X_even = ifft(x_even);
                        Complex[] X_odd = ifft(x_odd);
                        X = new Complex[N];
                        for (int i = 0; i < N / 2; i++)
                        {
                            X[i] = X_even[i] + w2(i, N) * X_odd[i];
                            X[i + N / 2] = X_even[i] - w2(i, N) * X_odd[i];
                        }
                    }
                    return X;
                }
                /// <summary>
                /// Центровка массива значений полученных в fft (спектральная составляющая при нулевой частоте будет в центре массива)
                /// </summary>
                /// <param name="X">Массив значений полученный в fft</param>
                /// <returns></returns>
                public static Complex[] nfft(Complex[] X)
                {
                    int N = X.Length;
                    Complex[] X_n = new Complex[N];
                    for (int i = 0; i < N / 2; i++)
                    {
                        X_n[i] = X[N / 2 + i];
                        X_n[N / 2 + i] = X[i];
                    }
                    return X_n;
                }
            }
            static public Vector Fourier(string function, Vector x, double step, int a, int b)
            {
                try
                {
                    Complex[] y = new Complex[x.Length];
                    for (int i = 0; i < y.Length; i++) y[i] = new Complex();
                    Vector res = new Vector(x.Length);
                    Complex j = new Complex(0.0, -1.0);
                    double X;
                    int N = (int)((b - a) / step);
                    for (int k = 0; k < x.Length; k++)
                    {
                        for (int i = 0; i <= N; i++)
                        {
                            X = step * i + a;
                            y[k] += Function.Func.f(function, X) * Complex.Exp(-j * X * x.elements[k]);
                        }
                        y[k] *= step;
                    }
                    double w = 1.0 / Math.Sqrt(2.0 * Math.PI);
                    for (int i = 0; i < y.Length; i++)
                    {
                        y[i] *= w;
                        res.elements[i] = y[i].Magnitude;
                    }
                    return res;
                }
                catch { return null; }
            }
            static public Vector Kepstral(Complex[] x)
            {
                Vector val = Convertor.Complex_.ComplexInMagnitude(FFT.fft(x)).ElementsPow(2.0);
                for (int i = 0; i < val.Length; i++)
                    val[i] = Math.Log(val[i]);
                return Convertor.Complex_.ComplexInReal(FFT.fft(val));
            }
            static public Vector Kepstral(Vector x)
            {
                Vector val = Convertor.Complex_.ComplexInMagnitude(FFT.fft(x)).ElementsPow(2.0);
                for (int i = 0; i < val.Length; i++)
                    val[i] = Math.Log(val[i]);
                return Convertor.Complex_.ComplexInReal(FFT.fft(val));
            }
        }
        public class ZedGraph
        {
            public static void DrawZedGraphCurve(ZedGraphControl zedGraphControl, double[] y, double a = -1, double b = -1, Color? col = null)
            {
                Color color = col ?? Color.Black;
                if (a == -1) a = 1;
                if (b == -1) b = y.Length;
                double[] x = new double[y.Length];
                double step = Math.Abs(b - a + 0.0) / (double)x.Length;
                for(int i = 0; i < x.Length; i++)
                {
                    x[i] = a + step * i;
                }
                zedGraphControl.GraphPane.CurveList.Clear();
                zedGraphControl.GraphPane.AddCurve("", x, y, color, SymbolType.None);
                zedGraphControl.Invalidate();
                zedGraphControl.AxisChange();
            }

            public static void DrawZedGraphCurveAdd(ZedGraphControl zedGraphControl, double[] y, double a = -1, double b = -1, Color? col = null)
            {
                Color color = col ?? Color.Black;
                if (a == -1) a = 1;
                if (b == -1) b = y.Length;
                double[] x = new double[y.Length];
                double step = Math.Abs(b - a + 0.0) / (double)x.Length;
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] = a + step * i;
                }
                //zedGraphControl.GraphPane.CurveList.Clear();
                zedGraphControl.GraphPane.AddCurve("", x, y, color, SymbolType.None);
                zedGraphControl.Invalidate();
                zedGraphControl.AxisChange();
            }
            public static void DrawZedGraphClear(ZedGraphControl zedGraphControl)
            {
                zedGraphControl.GraphPane.CurveList.Clear();
            }
            public static void DrawZedGraphBar(ZedGraphControl zedGraphControl, double[] y, int a, int b, Color color)
            {
                double[] x = new double[y.Length];
                double step = Math.Abs(b - a + 0.0) / (double)x.Length;
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] = a + step * i;
                }
                zedGraphControl.GraphPane.CurveList.Clear();
                zedGraphControl.GraphPane.AddBar("", x, y, color);
                zedGraphControl.Invalidate();
                zedGraphControl.AxisChange();
            }
        }
        public class Convertor
        {
            public class Picture
            {
                unsafe static public Bitmap MatrixInBitmap(Matrix x)
                {
                    Bitmap pict = new Bitmap(x.width, x.height);
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    try
                    {
                        byte* curpos;
                        fixed (double* res_ = x.elements)
                        {
                            double* elem = res_;
                            for (int i = 0; i < x.height; i++)
                            {
                                curpos = ((byte*)pictData.Scan0) + i * pictData.Stride;
                                for (int j = 0; j < x.width; j++)
                                {
                                    int val = Round(*elem++);
                                    if (val > 255) val = 255;
                                    else if (val < 0) val = 0;
                                    byte value = Convert.ToByte(val);
                                    *(curpos++) = value;
                                    *(curpos++) = value;
                                    *(curpos++) = value;
                                }
                            }
                        }
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return pict;
                }
                unsafe static public Bitmap Tensor3InBitmap(Tensor3 x)
                {
                    if (x.deep != 3) throw new Exception();

                    //double[, ,] res = x.GetElem();
                    Bitmap pict = new Bitmap(x.width, x.height);
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    int Width = x.width;
                    int Height = x.height;
                    try
                    {
                        byte* curpos;
                        //fixed (double* res_ = res)
                        {
                            //double* r = res_, g = res_ + Width * Height, b = res_ + 2 * Width * Height;
                            for (int i = 0; i < Height; i++)
                            {
                                curpos = ((byte*)pictData.Scan0) + i * pictData.Stride;
                                for (int j = 0; j < Width; j++)
                                {
                                    //int r_ = Round(*r++), g_ = Round(*g++), b_ = Round(*b++);
                                    int r_ = Round(x[0, i, j]), g_ = Round(x[1, i, j]), b_ = Round(x[2, i, j]);

                                    if (r_ < 0) r_ = 0;
                                    else if (r_ > 255) r_ = 255;

                                    if (g_ < 0) g_ = 0;
                                    else if (g_ > 255) g_ = 255;

                                    if (b_ < 0) b_ = 0;
                                    else if (b_ > 255) b_ = 255;

                                    byte R = Convert.ToByte(r_), G = Convert.ToByte(g_), B = Convert.ToByte(b_);

                                    *(curpos++) = R;
                                    *(curpos++) = G;
                                    *(curpos++) = B;
                                }
                            }
                        }
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return pict;
                }
                unsafe static public Tensor3 BitmapInTensor3(Bitmap pict)
                {
                    Tensor3 res = new Tensor3(pict.Width, pict.Height, 3);
                    //double[,,] res = new double[3, pict.Height, pict.Width];
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    int Width = res.width;
                    int Height = res.height;
                    try
                    {
                        byte* curpos;
                        //fixed (double* res_ = res)
                        {
                            //double* r = res_, g = res_ + Width * Height, b = res_ + 2 * Width * Height;
                            for (int i = 0; i < Height; i++)
                            {
                                curpos = ((byte*)pictData.Scan0) + i * pictData.Stride;
                                for (int j = 0; j < Width; j++)
                                {
                                    res[0, i, j] = *(curpos++);
                                    res[1, i, j] = *(curpos++);
                                    res[2, i, j] = *(curpos++);
                                    //*(b++) = *(curpos++);
                                    //*(g++) = *(curpos++);
                                    //*(r++) = *(curpos++);
                                }
                            }
                        }
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return res;
                }
                unsafe static public Tensor4 BitmapInTensor4(Bitmap pict)
                {
                    Tensor4 res = new Tensor4(pict.Width, pict.Height, 3, 1);
                    //double[,,] res = new double[3, pict.Height, pict.Width];
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    int Width = res.width;
                    int Height = res.height;
                    try
                    {
                        byte* curpos;
                        //fixed (double* res_ = res)
                        {
                            //double* r = res_, g = res_ + Width * Height, b = res_ + 2 * Width * Height;
                            for (int i = 0; i < Height; i++)
                            {
                                curpos = ((byte*)pictData.Scan0) + i * pictData.Stride;
                                for (int j = 0; j < Width; j++)
                                {
                                    res[0, 0, i, j] = *(curpos++);
                                    res[0, 1, i, j] = *(curpos++);
                                    res[0, 2, i, j] = *(curpos++);
                                    //*(b++) = *(curpos++);
                                    //*(g++) = *(curpos++);
                                    //*(r++) = *(curpos++);
                                }
                            }
                        }
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return res;
                }
                public unsafe static Matrix BitmapInMatrix(Bitmap pict)
                {
                    Matrix res = new Matrix(pict.Width, pict.Height);
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    try
                    {
                        byte* curpos;
                        fixed (double* res_ = res.elements)
                        {
                            double* elem = res_;
                            for (int i = 0; i < res.height; i++)
                            {
                                curpos = ((byte*)pictData.Scan0) + i * pictData.Stride;
                                for (int j = 0; j < res.width; j++)
                                {
                                    int value = Round((*(curpos++) + *(curpos++) + *(curpos++)) / 3.0);
                                    *elem++ = value;
                                }
                            }
                        }
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return res;
                }
                public static unsafe Tensor4 BitmapInBlocks(Bitmap pict, Size windowSize)
                {
                    int numH = pict.Height / windowSize.Height;
                    int numW = pict.Width / windowSize.Width;

                    var tensor4 = new Tensor4(windowSize.Width, windowSize.Height, 3, numH * numW);

                    int pictWidth = pict.Width;
                    int pictHeight = pict.Height;
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    try
                    {
                        Parallel.For(0, numH, y =>
                        {//for(int y = 0; y < numH; y++)
                            for (int x = 0; x < numW; x++)
                            {
                                // double[, ,] res = new double[3, windowSize.Height, windowSize.Width];
                                int Width = tensor4.width;
                                int Height = tensor4.height;
                                byte* curpos;
                                //fixed (double* res_ = res)
                                {
                                    //double* r = res_, g = res_ + Width * Height, b = res_ + 2 * Width * Height;
                                    for (int i = 0; i < Height; i++)
                                    {
                                        curpos = ((byte*)pictData.Scan0) + i * pictData.Stride + y * Height * pictData.Stride + x * Width * 3;
                                        for (int j = 0; j < Width; j++)
                                        {
                                            tensor4[y * numW + x, 0, i, j] = *(curpos++);
                                            tensor4[y * numW + x, 1, i, j] = *(curpos++);
                                            tensor4[y * numW + x, 2, i, j] = *(curpos++);
                                            //*(b++) = *(curpos++);
                                            //*(g++) = *(curpos++);
                                            //*(r++) = *(curpos++);
                                        }
                                    }
                                }
                                //tensors4[y, x] = new Tensor3(res).ToTensor4();
                            }
                        });
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return tensor4;
                }
                public static unsafe Tensor4 BitmapInBlocksMonochrome(Bitmap pict, Size windowSize)
                {
                    int numH = pict.Height / windowSize.Height;
                    int numW = pict.Width / windowSize.Width;

                    var tensor4 = new Tensor4(windowSize.Width, windowSize.Height, 1, numH * numW);

                    int pictWidth = pict.Width;
                    int pictHeight = pict.Height;
                    BitmapData pictData = pict.LockBits(new Rectangle(0, 0, pict.Width, pict.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    try
                    {
                        Parallel.For(0, numH, y =>
                        {//for(int y = 0; y < numH; y++)
                            for (int x = 0; x < numW; x++)
                            {
                                // double[, ,] res = new double[3, windowSize.Height, windowSize.Width];
                                int Width = tensor4.width;
                                int Height = tensor4.height;
                                byte* curpos;
                                //fixed (double* res_ = res)
                                {
                                    //double* r = res_, g = res_ + Width * Height, b = res_ + 2 * Width * Height;
                                    for (int i = 0; i < Height; i++)
                                    {
                                        curpos = ((byte*)pictData.Scan0) + i * pictData.Stride + y * Height * pictData.Stride + x * Width * 3;
                                        for (int j = 0; j < Width; j++)
                                            tensor4[y * numW + x, 0, i, j] = (*(curpos++) + *(curpos++) + *(curpos++)) / 3;
                                        
                                    }
                                }
                                //tensors4[y, x] = new Tensor3(res).ToTensor4();
                            }
                        });
                    }
                    finally
                    {
                        pict.UnlockBits(pictData);
                    }
                    return tensor4;
                }
            }
            public class Complex_
            {
                public static Vector ComplexInMagnitude(Complex[] complex)
                {
                    Vector res = new Vector(complex.Length);
                    for (int i = 0; i < res.Length; i++)
                        res[i] = Math.Sqrt(complex[i].Real * complex[i].Real + complex[i].Imaginary * complex[i].Imaginary);
                    return res;
                }
                public static Vector ComplexInReal(Complex[] complex)
                {
                    Vector res = new Vector(complex.Length);
                    for (int i = 0; i < res.Length; i++)
                        res[i] = complex[i].Real;
                    return res;
                }
            }
            public class Notation
            {
                //
                public static double BinToDouble(byte[] bytes) // max value: 65 535 (2^15+2^14+...+2^0)
                {

                    if (bytes.Length != 37) throw new Exception();
                    // 1 байт с конца обозначает знак. (0 = -; 1 = +)
                    // 16 разрядов отведено под целую часть.
                    // 20 разрядов отведено под мантиссу. (6 знаков)
                    double res = 0.0;
                    int n = 1;
                    for (int i = 35; i >= 20; i--)
                    {
                        res += n * bytes[i];
                        n = n << 1;
                    }
                    n = 1;
                    double mantissa = 0.0;
                    for (int i = 0; i < 20; i++)
                    {
                        mantissa += n * bytes[19 - i];
                        n = n << 1;
                    }
                    if (mantissa > 999999) mantissa = 999999;
                    res += mantissa / 1000000.0;

                    if (bytes[bytes.Length - 1] == 0) res *= -1.0;

                    return res;
                }
                public static byte[] DoubleToBin(double value) // max value: 65 535
                {
                    if(value > 65535) throw new Exception();
                    byte[] bytes = new byte[37];
                    bytes[bytes.Length - 1] = (value > 0) ? (byte)1 : (byte)0;
                    var val = Math.Round(value, 6);
                    int z = (int)val;

                    for (int i = 0; z != 0; i++)
                    {
                        bytes[bytes.Length - 2 - i] = Convert.ToByte((z - 2*(z = z >> 1)));
                    }

                    int m = (int)((val - (int)val) * 1000000.0);

                    for (int i = 0; m != 0; i++)
                    {
                        bytes[19 - i] = Convert.ToByte((m - 2*(m = m >> 1)));
                        //m /= 2;
                    }

                    return bytes;
                }
            }
        }
        public class Preprocessing
        {
            public class Filter
            {
                /// <summary>
                /// Медианный фильтр для цветного изображения.
                /// </summary>
                /// <param name="pict">входное изображение</param>
                /// <param name="Width">ширина фильтра</param>
                /// <param name="Height">высота фильтра</param>
                /// <returns></returns>
                public static Tensor3 Median(Tensor3 inp, int Width, int Height)
                {
                    int H = inp.height - Height + 1, W = inp.width - Width + 1;
                    Parallel.For(0, inp.deep, z =>
                    {
                        for (int y = 0; y < H; y += Height)
                            for (int x = 0; x < W; x += Width)
                            {
                                Vector a = new Vector(Height * Width);
                                for (int dy = 0; dy < Height; dy++)
                                {
                                    for (int dx = 0; dx < Width; dx++)
                                    {
                                        a[dy * Width + dx] = inp[z, y + dy, x + dx];
                                    }
                                }
                                a.SortAscending();
                                double value = a[a.Length / 2];
                                for (int dy = 0; dy < Height; dy++)
                                {
                                    for (int dx = 0; dx < Width; dx++)
                                    {
                                        inp[z, y + dy, x + dx] = value;
                                    }
                                }
                            }
                    });
                    return inp;
                }
                public static Vector MedianDecrease(Vector inp, int Width)
                {
                    Vector res = new Vector(inp.Length / Width);
                    int iter = 0;
                    int W = inp.Length - Width + 1;

                            for (int x = 0; x < W; x += Width)
                            {
                                Vector a = new Vector(Width);
                                {
                                    for (int dx = 0; dx < Width; dx++)
                                    {
                                        a[dx] = inp[x + dx];
                                    }
                                }
                                a.SortAscending();
                                double value = a[a.Length / 2];
                                res[iter++] = value;
                            }
                    return res;
                }
                public static Bitmap Gabor(Bitmap pict, int n)
                {
                    Matrix matr = Convertor.Picture.BitmapInMatrix(pict);
                    Matrix res = new Matrix(matr.width, matr.height);
                    double disp = matr.ToVector().Dispersion();
                    double dispX = 3 * disp;
                    double dispY = 2 * dispX + 1;
                    for (int y = 0; y < matr.height; y++)
                        for (int x = 0; x < matr.width; x++)
                        {
                            res[y, x] = 1.0 / (n * n);
                            for (int dy = 0; dy < n; dy++)
                            {
                                for (int dx = 0; dx < n; dx++)
                                {
                                    res[y, x] += matr[y - n / 2 + dy, x - n / 2 + dx] * G(dx, dy, 1, dispX, dispY);

                                }
                            }
                        }
                    return res.Visualize();
                }
                private static double G(int x, int y, double theta, double dispX, double dispY)
                {
                    double x_theta = x * Math.Cos(theta) + y * Math.Sin(theta);
                    double y_theta = -x * Math.Sin(theta) + y * Math.Cos(theta);
                    return Math.Exp(-0.5 * ((x_theta * x_theta) / dispX + (y_theta * y_theta) / dispY)) * Math.Cos(2.0 * Math.PI * theta * x_theta);
                }
            }
            public class Picture
            {
                unsafe static public Bitmap ConvertInMonochromePicture(Bitmap picture)
                {
                    Bitmap result = new Bitmap(picture);
                    BitmapData btmData = result.LockBits(new Rectangle(0, 0, picture.Width, picture.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
                    IntPtr ptr = btmData.Scan0;
                    int N = picture.Width * picture.Height * 3;
                    int Width = picture.Width;
                    int Height = picture.Height;
                    try
                    {
                        byte* curpos;
                        for(int i = 0; i < Height; i++)
                        {
                            curpos = ((byte*)btmData.Scan0) + i * btmData.Stride;
                            for (int j = 0; j < Width; j++)
                            {
                                byte value = Convert.ToByte(Round((*(curpos) + *(curpos+1) + *(curpos+2)) / 3.0));
                                *curpos++ = value;
                                *curpos++ = value;
                                *curpos++ = value;
                            }
                        }
                    }
                    finally
                    {
                        result.UnlockBits(btmData);
                    }
                    return result;
                }
                public class Filter
                {
                    /// <summary>
                    /// Медианный фильтр для цветного изображения.
                    /// </summary>
                    /// <param name="pict">входное изображение</param>
                    /// <param name="Width">ширина фильтра</param>
                    /// <param name="Height">высота фильтра</param>
                    /// <returns></returns>
                    public static Bitmap Median(Bitmap pict, int Width, int Height)
                    {
                        Tensor3 tensor3 = Convertor.Picture.BitmapInTensor3(pict);
                        int H = tensor3.height - Height + 1, W = tensor3.width - Width + 1;
                        Parallel.For(0, tensor3.deep, z => {
                            for(int y = 0; y < H; y+=Height)
                                for (int x = 0; x < W; x+=Width)
                                {
                                    Vector a = new Vector(Height * Width);
                                    for (int dy = 0; dy < Height; dy++)
                                    {
                                        for (int dx = 0; dx < Width; dx++)
                                        {
                                            a[dy * Width + dx] = tensor3[z, y + dy, x + dx];
                                        }
                                    }
                                    a.SortAscending();
                                    double value = a[a.Length / 2];
                                    for (int dy = 0; dy < Height; dy++)
                                    {
                                        for (int dx = 0; dx < Width; dx++)
                                        {
                                            tensor3[z, y + dy, x + dx] = value;
                                        }
                                    }
                                }
                        });
                        return Convertor.Picture.Tensor3InBitmap(tensor3);
                    }
                    public static Bitmap Gabor(Bitmap pict, int n)
                    {
                        Matrix matr = Convertor.Picture.BitmapInMatrix(pict);
                        Matrix res = new Matrix(matr.width, matr.height);
                        double disp = matr.ToVector().Dispersion();
                        double dispX = 3 * disp;
                        double dispY =2 * dispX + 1;
                        for(int y = 0; y < matr.height; y++)
                            for(int x= 0; x < matr.width; x++)
                            {
                                res[y, x] = 1.0 / (n * n);
                                for (int dy = 0; dy < n; dy++)
                                {
                                    for (int dx = 0; dx < n; dx++)
                                    {
                                        res[y, x] += matr[y - n / 2 + dy, x - n / 2 + dx] * G(dx, dy, 1, dispX, dispY);

                                    }
                                }
                            }
                        return res.Visualize();
                    }
                    private static double G(int x, int y, double theta, double dispX, double dispY)
                    {
                        double x_theta = x * Math.Cos(theta) + y * Math.Sin(theta);
                        double y_theta = -x * Math.Sin(theta) + y * Math.Cos(theta);
                        return Math.Exp(-0.5 * ((x_theta * x_theta) / dispX + (y_theta * y_theta) / dispY)) * Math.Cos(2.0 * Math.PI * theta * x_theta);
                    }
                }
                public class Noise
                {
                    public class Add
                    {
                        /// <summary>
                        /// 
                        /// </summary>
                        /// <param name="pict">входное изображение</param>
                        /// <param name="Width">макс. ширина окрестности одной шумовой точки</param>
                        /// <param name="Height">макс. высота окрестности одной шумовой точки</param>
                        /// <returns></returns>
                        public static Bitmap Impulse(Bitmap pict, int Width, int Height, double threshold)
                        {
                            Tensor3 tensor3 = Convertor.Picture.BitmapInTensor3(pict);
                            System.Random rand = new System.Random();
                            int H = Round(rand.NextDouble() * Height);
                            int W = Round(rand.NextDouble() * Width);
                            if (H == 0) H = 1;
                            if (W == 0) W = 1;
                            int Y = tensor3.height - H + 1, X = tensor3.width - W + 1;
                            for(int y = 0; y < Y; y+= H)
                                for (int x = 0; x < X; x+= W)
                                {
                                    if(Statistic.Random.Gauss() > threshold)
                                    {
                                        int value = (rand.NextDouble() > 0.5) ? 255 : 0;
                                        for (int dy = 0; dy < H; dy++)
                                        {
                                            for (int dx = 0; dx < W; dx++)
                                            {
                                                tensor3[0, y + dy, x + dx] = value;
                                                tensor3[1, y + dy, x + dx] = value;
                                                tensor3[2, y + dy, x + dx] = value;
                                            }
                                        }
                                    }
                                }
                            return Convertor.Picture.Tensor3InBitmap(tensor3);
                        }
                    }
                }
            }
        }
        public class Combinatorics
        {
            public static Vector[] GetAllPermutations(int N)
            {
                Vector[] result = new Vector[(int)Factorial(N)];
                for (int i = 0; i < result.Length; i++)
                    result[i] = new Vector(N);
                for (int i = 0; i < N; i++)
                    result[0][i] = i + 1;
                for (int i = 1; i < result.Length; i++ )
                {
                    result[i] = GetNextPermutation(result[i - 1]);
                }
                    return result;
            }
            private static Vector GetNextPermutation(Vector x)
            {
                Vector result = new Vector(x);
                int k = x.Length - 2;
                while (k != -1 && x[k] >= x[k + 1]) k--;
                if (k == -1)
                    return null;
                int j = x.Length - 1;
                while (x[j] <= x[k]) j--;
                result[j] = x[k];
                result[k] = x[j];
                Swap(k + 1, result);
                return result;
            }
            private static void Swap(int IndexMin, Vector x)
            {
                int a = IndexMin;
                int b = x.Length - 1;
                while(a < b)
                {
                    double c = x[a];
                    x[a] = x[b];
                    x[b] = c;
                    a++; b--;
                }
            }
            public static int InverseNum(Vector permutation)
            {
                int result = 0;
                for (int i = 0; i < permutation.Length; i++)
                    for (int j = i + 1; j < permutation.Length; j++)
                        if (permutation[i] > permutation[j])
                            result++;
                return result;
            }
        }
        public class SolvingSystems
        {
            public static Vector SolvingSLAY(Matrix A, Vector b)
            {
                Vector x = new Vector(A.width);
                Matrix matrix = new Matrix(A);
                Vector vector = new Vector(b);
                for (int i = 0; i < matrix.width; i++)
                {
                    for (int j = i; j < matrix.height; j++)
                    {
                        if (matrix[j, i] != 0)
                        {
                            matrix.Swap(i, j, 1);
                            vector.Swap(i, j);
                            break;
                        }
                        else if (j + 1 == matrix.height)
                        {
                            i++;
                            goto go;
                        }
                    }
                    Vector a = matrix.GetVector(i, 1);
                    a *= 1.0 / a[i];
                    double value = matrix[i, i];
                    for (int j = i + 1; j < matrix.height; j++)
                    {
                        double c = matrix[j, i];
                        for (int k = i; k < matrix.width; k++)
                            matrix[j, k] -= a[k] * c;
                        vector[j] -= vector[i] * c / value;
                    }
                go: { }
                }
                double determinant = 1.0;
                for (int i = 0; i < matrix.width; i++)
                    determinant *= matrix[i, i];
                if (determinant == 0)
                    throw new Exception("Error: determinant = 0");
                for (int i = matrix.height - 1; i >= 0; i--)
                {
                    for (int j = matrix.width - 1; j > i; j--)
                        x[i] += x[j] * matrix[i, j];
                    x[i] = (vector[i] - x[i]) / matrix[i, i];
                }
                return x;
            }
        }
        public class Text
        {
            public class EditorialDistance
            {
                static public int LevenshteinDistance(string text1, string text2)
                {
                    string[] texts = PreparingStrs(text1, text2);
                    int[,] D = new int[text1.Length, text2.Length];
                    for (int i = 0; i < text1.Length; i++)
                        for (int j = 0; j < text2.Length; j++)
                        {
                            if (i == 0 && j == 0)
                                D[i, j] = 0;
                            else if (j == 0 && i > 0)
                                D[i, j] = i;
                            else if (i == 0 && j > 0)
                                D[i, j] = j;
                            else
                            {
                                int min = D[i, j - 1] + 1;
                                if (min > D[i - 1, j] + 1)
                                    min = D[i - 1, j] + 1;
                                int m = 0;
                                if (!(texts[0][i] == texts[1][j]))
                                    m = 1;
                                if (min > D[i - 1, j - 1] + m)
                                    min = D[i - 1, j - 1] + m;
                                D[i, j] = min;
                            }
                        }
                    return D[text1.Length - 1, text2.Length - 1];
                }
                static private string[] PreparingStrs(string text1, string text2)
                {
                    string[] strs = new string[2];
                    if (!(text1.Length == text2.Length))
                    {
                        if (text1.Length > text2.Length)
                            for (int i = text2.Length; i < text1.Length; i++)
                                text2 += " ";
                        else
                            for (int i = text1.Length; i < text2.Length; i++)
                                text1 += " ";
                    }
                    strs[0] = " " + text1;
                    strs[1] = " " + text2;
                    return strs;
                }
            }
        }
        public static int Round(double x)
        {
            double mantissa = x - (int)x;
            if (mantissa < 0.48) return (int)x;
            else return (int)(x + 1.0);
        }
        public static decimal Round(decimal x)
        {
            decimal mantissa = x - (int)x;
            if (mantissa < (decimal)0.48) return (int)x;
            else return (decimal)(x + (decimal)1.0);
        }
        public class Sort
        {
            static int partition(Vector array, int start, int end)
            {
                double temp;//swap helper
                int marker = start;//divides left and right subarrays
                for (int i = start; i <= end; i++)
                {
                    if (array[i] < array[end]) //array[end] is pivot
                    {
                        temp = array[marker]; // swap
                        array[marker] = array[i];
                        array[i] = temp;
                        marker += 1;
                    }
                }
                //put pivot(array[end]) between left and right subarrays
                temp = array[marker];
                array[marker] = array[end];
                array[end] = temp;
                return marker;
            }

            static void quicksort(Vector array, int start, int end)
            {
                if (start >= end)
                {
                    return;
                }
                int pivot = partition(array, start, end);
                quicksort(array, start, pivot - 1);
                quicksort(array, pivot + 1, end);
            }
            public static Vector qsort(Vector x)
            {
                Vector X = new Vector(x);
                quicksort(X, 0, X.Length - 1);
                return X;
            }



            static int partitionMulty(Vector array, Vector indexes, int start, int end)
            {
                double temp;//swap helper
                double temp2;
                int marker = start;//divides left and right subarrays
                for (int i = start; i <= end; i++)
                {
                    if (array[i] < array[end]) //array[end] is pivot
                    {
                        temp = array[marker]; // swap
                        temp2 = indexes[marker];

                        array[marker] = array[i];
                        array[i] = temp;

                        indexes[marker] = indexes[i];
                        indexes[i] = temp2;

                        marker += 1;
                    }
                }
                //put pivot(array[end]) between left and right subarrays
                temp = array[marker];
                temp2 = indexes[marker];

                array[marker] = array[end];
                array[end] = temp;

                indexes[marker] = indexes[end];
                indexes[end] = temp2;

                return marker;
            }
            static void quicksortMulty(Vector array, Vector indexes, int start, int end)
            {
                if (start >= end)
                {
                    return;
                }
                int pivot = partitionMulty(array, indexes, start, end);
                quicksortMulty(array, indexes, start, pivot - 1);
                quicksortMulty(array, indexes, pivot + 1, end);
            }
            public static Vector[] qsort(Vector indexes, Vector x)
            {
                Vector X = new Vector(x);
                Vector Ind = new Vector(indexes);
                quicksortMulty(X, Ind, 0, x.Length - 1);

                return new Vector[] { Ind, X };
            }
        }
    }
}
