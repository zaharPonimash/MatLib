using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Cloo;
using MatLib;
using OpenCLTemplate;
using System.IO;
namespace MatLib
{
    /* Части кода
      for (int d = 0; d < bs; d++)
        for (int z = 0; z < deep; z++)
         for (int y = 0; y < height; y++)
          for (int x = 0; x < width; x++)
           {
     
           }
    //*/
    [Serializable]
    public class Tensor4
    {
        static List<CLCalc.Program.Kernel> kernels;

        public double[] elements;
        public int width, height, deep, bs;
        public int dhw, hw;
        //CLCalc.CLPrograms.doubleLinearAlgebra program = new CLCalc.CLPrograms.doubleLinearAlgebra();
        //float[,] a = new float[2000, 2000];
        //float[,] b = new float[2000, 2000];
        //float[,] syka = new float[2000, 2000];
        //CLCalc.CLPrograms.floatLinearAlgebra la = new CLCalc.CLPrograms.floatLinearAlgebra();
        int getIndex(int bs, int z, int y, int x)
        {
            return bs * dhw + z * hw + y * width + x;
        }
        public double this[int bs, int z, int y, int x]
        {
            get { return elements[bs * dhw + z * hw + y * width + x]; }
            set { elements[bs * dhw + z * hw + y * width + x] = value; }
        }
        public Tensor4(int width, int height, int deep, int bs)
        {
            this.width = width; this.height = height; this.deep = deep; this.bs = bs;
            hw = height * width;
            dhw = hw * deep;
            elements = new double[bs * dhw];
        }
        public Tensor4(double[,,,] tensor4)
        {
            this.width = tensor4.GetLength(3);
            this.height = tensor4.GetLength(2);
            this.deep = tensor4.GetLength(1);
            this.bs = tensor4.GetLength(0);

            hw = height * width;
            dhw = hw * deep;

            elements = new double[bs * dhw];

            for (int d = 0; d < bs; d++)
                for (int z = 0; z < deep; z++)
                    for (int y = 0; y < height; y++)
                        for (int x = 0; x < width; x++)
                            elements[d * dhw + z * hw + y * width + x] = tensor4[d, z, y, x];
        }
        public Tensor4(double[] vector, int width, int height, int deep, int bs)
        {
            this.width = width;
            this.height = height;
            this.deep = deep;
            this.bs = bs;
            hw = height * width;
            dhw = hw * deep;
            elements = new double[vector.Length];

            vector.CopyTo(elements, 0);
        }
        public Tensor4(Tensor3[] tensor3)
        {
            this.width = tensor3[0].width; this.height = tensor3[0].height; this.deep = tensor3[0].deep; this.bs = tensor3.Length;
            hw = height * width;
            dhw = hw * deep;

            elements = new double[bs * dhw];

            for (int d = 0; d < bs; d++)
                for (int z = 0; z < deep; z++)
                    for (int y = 0; y < height; y++)
                        for (int x = 0; x < width; x++)
                            elements[d * dhw + z * hw + y * width + x] = tensor3[d][z, y, x];
        }

        public void InitGPU()
        {
            kernels = new List<CLCalc.Program.Kernel>();

            //OpenCLTemplate.CLCalc.InitCL(Cloo.ComputeDeviceTypes.All);
            CLCalc.InitCL();
            CLCalc.Program.DefaultCQ = 0;
            string text = "";
            using (StreamReader stream = new StreamReader(@"C:\Users\Marat\Documents\Visual Studio 2013\Projects\MathLib\MathLib\Programs.txt"))
            {
                text = stream.ReadToEnd();
            }
            CLCalc.Program.Compile(new string[] { text });

            kernels.Add(new CLCalc.Program.Kernel("vecMul")); // 0 
            kernels.Add(new CLCalc.Program.Kernel("vecElemMul")); // 1
            kernels.Add(new CLCalc.Program.Kernel("vecSum")); // 2
            kernels.Add(new CLCalc.Program.Kernel("vecDif")); // 3

            kernels.Add(new CLCalc.Program.Kernel("matrMul")); // 4
        }
        public void Fill(double value)
        {
            for (int x = 0; x < elements.Length; x++)
                elements[x] = value;
        }
        public static Tensor4 operator +(Tensor4 a, Tensor4 b)
        {
            if (a.bs != b.bs || a.deep != b.deep || a.height != b.height || a.width != b.width) throw new Exception("Размерности не совпадают");
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] + b.elements[x];

            return res;
        }
        public static Tensor4 operator -(Tensor4 a, Tensor4 b)
        {
            if (a.bs != b.bs || a.deep != b.deep || a.height != b.height || a.width != b.width) throw new Exception("Размерности не совпадают");
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] - b.elements[x];

            return res;
        }
        public static Tensor4 operator *(Tensor4 a, Tensor4 b)
        {
            if (a.bs != b.bs || a.deep != b.deep || a.height != b.height || a.width != b.width) throw new Exception("Размерности не совпадают");
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] * b.elements[x];

            return res;
        }
        public static Tensor4 operator /(Tensor4 a, Tensor4 b)
        {
            if (a.bs != b.bs || a.deep != b.deep || a.height != b.height || a.width != b.width) throw new Exception("Размерности не совпадают");
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] / b.elements[x];

            return res;
        }
        public static Tensor4 operator +(Tensor4 a, double b)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] + b;

            return res;
        }
        public static Tensor4 operator +(double b, Tensor4 a)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] + b;

            return res;
        }
        public static Tensor4 operator -(Tensor4 a, double b)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] - b;
            return res;
        }
        public static Tensor4 operator -(double b, Tensor4 a)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = b - a.elements[x];
            return res;
        }
        public static Tensor4 operator *(Tensor4 a, double b)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] * b;

            return res;
        }
        public static Tensor4 operator *(double b, Tensor4 a)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] * b;
            return res;
        }
        public static Tensor4 operator /(Tensor4 a, double b)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);


            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = a.elements[x] / b;

            return res;
        }
        public static Tensor4 operator /(double b, Tensor4 a)
        {
            Tensor4 res = new Tensor4(a.width, a.height, a.deep, a.bs);

            for (int x = 0; x < a.elements.Length; x++)
                res.elements[x] = b / a.elements[x];

            return res;
        }
        public Tensor4 VectorOnMatrix(Tensor4 matr)
        {

            //CLCalc.Program.Variable varA = new CLCalc.Program.Variable(a);
            //CLCalc.Program.Variable varB = new CLCalc.Program.Variable(b);
            //CLCalc.Program.Variable varC = new CLCalc.Program.Variable(c);
            //la.MatrixSum(a, b);
            //la.MatrixMultiply(a, b);
            //syka = program.MatrixMultiply(a, b);
            //CLCalc.CLPrograms.doubleLinearAlgebra list = new CLCalc.Program.MemoryObject();
           
            //*
            if(matr.deep != 1) throw new Exception("I don't no.");
            if (dhw != matr.height) throw new Exception("Размерности не совпадают");

            Tensor4 res = new Tensor4(matr.width, 1, 1, bs);
            
            if (bs == 1)
            {
                //for (int d = 0; d < bs; d++)
                Parallel.For(0, matr.width, c =>
                {//for (int c = 0; c < matr.width; c++)
                    for (int x = 0; x < matr.height; x++)
                        res[0, 0, 0, c] += elements[x] * matr[0, 0, x, c];
                });
            }
            else
            {
                Parallel.For(0, bs, d =>
                {//for (int d = 0; d < bs; d++)
                    for (int c = 0; c < matr.width; c++)
                        for (int x = 0; x < matr.height; x++)
                            res[d, 0, 0, c] += elements[d * dhw + x] * matr[0, 0, x, c];
                });
            }
                  //*/
            return res;
        }
        public Tensor4 MatrixOnVector(Tensor4 vector)
        {
            if (deep != 1) throw new Exception("I don't no.");

            if (vector.dhw != width) throw new Exception("Размерности не совпадают");

            Tensor4 res = new Tensor4(height, 1, 1, vector.bs);

            if (vector.bs == 1)
            {
                //for (int d = 0; d < vector.bs; d++)
                Parallel.For(0, height, c =>
                {//for (int c = 0; c < height; c++)
                    for (int x = 0; x < width; x++)
                                res[0, 0, 0, c] += elements[c * width + x] * vector.elements[x];
                });
            }
            else
            {
                Parallel.For(0, vector.bs, d =>
                {//for (int d = 0; d < vector.bs; d++)
                    for (int c = 0; c < height; c++)
                        for (int x = 0; x < width; x++)
                            res[d, 0, 0, c] += elements[c * width + x] * vector.elements[d * vector.dhw + x];
                });
            }
            return res;
        }
        public double EuclidNorm()
        {
            double result = 0.0;

            for (int x = 0; x < elements.Length; x++)
                result += elements[x] * elements[x];

            result = Math.Sqrt(result);
            return result;
        }
        public double EuclideanDistance(Tensor4 b)
        {
            double result = 0.0;

            for (int x = 0; x < elements.Length; x++)
                result += Math.Pow( (elements[x] - b.elements[x]), 2.0);

            result = Math.Sqrt(result);
            return result;
        }
        public Vector Distance(Tensor4 b)
        {
            if (bs != 1 || b.bs != 1 || deep != b.deep || width != b.width || height != b.height) throw new Exception("Чет браток не то у тебя.");
            Vector result = new Vector(deep);
            for (int i = 0; i < deep; i++)
            {
                for (int j = 0; j < height; j++)
                    for (int k = 0; k < width; k++)
                        result[i] += Math.Abs(elements[getIndex(0, i, j, k)] - b[0, i, j, k]);

                result[i] /= (double)hw;
            }

            return result;
        }
        public Tensor4 DistanceSquare(Tensor4 b)
        {
            if (deep != b.deep || width != b.width || height != b.height) throw new Exception("Чет браток не то у тебя.");
            Tensor4 result = new Tensor4(b.bs, 1, 1, bs);
            for (int d = 0; d < bs; d++)
            {
                for (int s = 0; s < b.bs; s++)
                {
                    for (int i = 0; i < deep; i++)
                        for (int j = 0; j < height; j++)
                            for (int k = 0; k < width; k++)
                                result[d, 0, 0, s] += Math.Pow(elements[getIndex(d, i, j, k)] - b[s, i, j, k], 2.0);

                    result[d, 0, 0, s] /= (double)dhw;
                }
                
            }
            return result;
        }
        public Tensor4 ConvEuclid(Tensor4 b)
        {
            if (deep != b.deep) throw new Exception("Глубины не совпадают.");
            Tensor4 res = new Tensor4(width - b.width + 1, height - b.height + 1, b.bs, bs);
            if (bs == 1)
            {
                Parallel.For(0, b.bs, s =>
                {//for (int s = 0; s < b.bs; s++)
                    for (int y = 0; y < res.height; y++)
                        for (int x = 0; x < res.width; x++)
                        {
                            for (int z = 0; z < deep; z++)
                                for (int dy = 0; dy < b.height; dy++)
                                    for (int dx = 0; dx < b.width; dx++)
                                        res[0, s, y, x] += Math.Pow(elements[z * hw + (y + dy) * width + x + dx] - b[s, z, dy, dx], 2.0);

                            res[0, s, y, x] /= b.hw;
                        }
                });
            }
            else
            {
                Parallel.For(0, bs, d =>
                {//for (int d = 0; d < bs; d++)
                    for (int s = 0; s < b.bs; s++)
                        for (int y = 0; y < res.height; y++)
                            for (int x = 0; x < res.width; x++)
                            {
                                for (int z = 0; z < deep; z++)
                                    for (int dy = 0; dy < b.height; dy++)
                                        for (int dx = 0; dx < b.width; dx++)
                                            res[0, s, y, x] += Math.Pow(elements[z * hw + (y + dy) * width + x + dx] - b[s, z, dy, dx], 2.0);

                                res[0, s, y, x] /= b.hw;
                            }
                });
            }
            return res;
        }
        public Tensor4 BatchNormalization()
        {
            if (bs > 1)
            {
                Tensor4 average = AverageOnBatch();
                Tensor4 res = new Tensor4(width, height, deep, bs);
                for (int x = 0; x < dhw; x++)
                {
                    Tensor4 dispersion = new Tensor4(width, height, deep, 1);

                    for (int d = 0; d < bs; d++)
                        dispersion.elements[x] += (elements[d * dhw + x] - average.elements[x]) * (elements[d * dhw + x] - average.elements[x]);

                    dispersion.elements[x] /= (double)bs - 1.0;
                    dispersion.elements[x] = Math.Sqrt(dispersion.elements[x]);
                    for (int d = 0; d < bs; d++)
                        res.elements[d * dhw + x] = elements[d * dhw + x] / dispersion.elements[x];
                }
                return res;
            }
            else return new Tensor4(elements, width, height, deep, bs);

        }
        public Tensor4 Conv(Tensor4 b)
        {
            if (deep != b.deep) throw new Exception("Глубины не совпадают.");
            Tensor4 res = new Tensor4(width - b.width + 1, height - b.height + 1, b.bs, bs);
            if (bs == 1)
            {
                Parallel.For(0, b.bs, s =>
                {//for (int s = 0; s < b.bs; s++)
                    for (int y = 0; y < res.height; y++)
                        for (int x = 0; x < res.width; x++)
                            for (int z = 0; z < deep; z++)
                                for (int dy = 0; dy < b.height; dy++)
                                    for (int dx = 0; dx < b.width; dx++)
                                        res[0, s, y, x] += elements[z * hw + (y + dy) * width + x + dx] * b[s, z, dy, dx];     
                });
            }
            else
            {
                Parallel.For(0, bs, d =>
                {//for (int d = 0; d < bs; d++)
                    for (int s = 0; s < b.bs; s++)
                        for (int y = 0; y < res.height; y++)
                            for (int x = 0; x < res.width; x++)
                                for (int z = 0; z < deep; z++)
                                    for (int dy = 0; dy < b.height; dy++)
                                        for (int dx = 0; dx < b.width; dx++)
                                            res[d, s, y, x] += elements[d * dhw + z * hw + (y + dy) * width + x + dx] * b[s, z, dy, dx];   
                });
            }
            return res;
        }
        public Tensor4 AverageOnBatch()
        {
            Tensor4 res = new Tensor4(width, height, deep, 1);
            for (int x = 0; x < dhw; x++)
            {
                for (int d = 0; d < bs; d++)
                    res.elements[x] += elements[d * dhw + x];

                res.elements[x] /= (double)bs;
            }
                return res;
        }
        public void AverageNormaliz(int butch = 0)
        {
            double aver = 0;
            for (int x = 0; x < dhw; x++) {
                    aver += elements[butch * dhw + x];
            }
            aver /= (double)dhw;
            for (int x = 0; x < dhw; x++) {
                elements[butch * dhw + x] -= aver;
            }
        }
        public void DispersionNormaliz(int butch = 0)
        {
            Vector vec = new Vector(dhw);
            for (int i = 0; i < dhw; i++){
                vec[i] = elements[butch * dhw + i];
            }

            double disp = Math.Sqrt(vec.Dispersion());

            for (int x = 0; x < dhw; x++)
            {
                elements[butch * dhw + x] /= disp;
            }
        }
        public void Random()
        {
            for (int x = 0; x < elements.Length; x++)
                elements[x] = matlib.Statistic.Random.Gauss();
        }
        public Vector ToVector(int butch)
        {
            Vector res = new Vector(dhw);
            int a = butch * dhw;
            for (int x = 0; x < dhw; x++)
                res[x] = elements[x + a];

            return res;
        }
        public Vector ToVector()
        {
            Vector res = new Vector(elements.Length);
            res.elements = elements;
            return res;
        }
        public Matrix ToMatrix(int butch, int deep)
        {
            Matrix res = new Matrix(width, height);
            int a = butch * dhw + deep * hw;

            for (int y = 0; y < height; y++)
                for (int x = 0; x < width; x++)
                    res[y, x] = elements[y * width + x + a];

            return res;
        }
        public double Sum(int butch)
        {
            double sum = 0.0;
            int a = butch * dhw;
            for (int x = 0; x < dhw; x++)
                sum += elements[a + x];

            return sum;
        }
        public double Sum()
        {
            return elements.Sum();
        }
        public double MaxAbs(int butch)
        {
            int offset = butch * dhw;
            double Max = Math.Abs(elements[offset]);

            for (int i = offset + 1; i < offset + dhw; i++)
                if (Max < Math.Abs(elements[i])) Max = Math.Abs(elements[i]);

            return Max;
        }
        /// <summary>
        /// Нормализует тензор, деля на максимальное по модулю отклонение в каждом батче.
        /// </summary>
        /// <returns></returns>
        public Tensor4 NormalizeAbs()
        {
            Tensor4 res = new Tensor4(width, height, deep, bs);
            for (int i = 0; i < bs; i++)
            {
                var max = MaxAbs(i);
                for (int j = 0; j < dhw; j++)
                    res.elements[i * dhw + j] = elements[i * dhw + j] / max;
            }
            return res;
        }
        public Tensor3 ToTensor3(int butch)
        {
            Tensor3 res = new Tensor3(width, height, deep);
            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        res[z, y, x] = elements[butch * dhw + z * hw + y * width + x];
            return res;
        }
        public Tensor4 Copy()
        {
            Tensor4 copy = new Tensor4(width, height, deep, bs);
            elements.CopyTo(copy.elements, 0);
            return copy;
        }
        public Tensor4 ElementsPow(double degree)
        {
            Tensor4 res = new Tensor4(width, height, deep, bs);

            for (int x = 0; x < elements.Length; x++)
                res.elements[x] = Math.Pow(elements[x], degree);

            return res;
        }
        public Tensor4 NormalizeTensor3(double degree)
        {
            Tensor4 res = new Tensor4(width, height, deep, bs);
            for (int i = 0; i < bs; i++)
            {
                double sum = 0.0;
                var offset = i * dhw;
                for (int k = offset; k < offset + dhw; k++ )
                    sum += Math.Pow(elements[k], degree);

                sum = Math.Pow(sum, 1.0 / degree);
                for (int k = offset; k < offset + dhw; k++)
                    res.elements[k] = elements[k] / sum;
            }
            return res;
        }
    }
}
