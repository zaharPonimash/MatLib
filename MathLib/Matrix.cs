using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
namespace MatLib
{
    [Serializable]
    public class Matrix
    {
        public double[] elements { get; set; }
        public int width, height;
        private Matrix resultM;
        private Vector vector;
        public Complex[,] elementsF;
        public double this[int i, int j]
        {
            get { return elements[getInd(i, j)]; }
            set { elements[getInd(i, j)] = value; }
        }
        int getInd(int y, int x)
        {
            return y * width + x;
        }
        public override string ToString()
        {
            string Text = "";
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                    Text += elements[getInd(i, j)] + "  ";
                Text += "\n";
            }
            return Text;
                    //return base.ToString();
        }
        public static Matrix operator +(Matrix a, Matrix b)
        {
            Matrix res = new Matrix(a.width, a.height);
            if ((a.height == b.height) && (a.width == b.width))
                for (int i = 0; i < a.elements.Length; i++)
                    res.elements[i] = a.elements[i] + b.elements[i];

            else return null;
            return res;
        }
        public static Matrix operator +(Matrix a, double b)
        {
            Matrix res = new Matrix(a.width, a.height);

            for (int i = 0; i < a.elements.Length; i++)
                res.elements[i] = a.elements[i] + b;

            return res;
        }
        public static Matrix operator -(Matrix a, Matrix b)
        {
            Matrix res = new Matrix(a.width, a.height);
            if ((a.height == b.height) && (a.width == b.width))
                for (int i = 0; i < a.elements.Length; i++)
                    res.elements[i] = a.elements[i] - b.elements[i];
            else return null;
            return res;
        }
        public static Matrix operator -(Matrix a, double b)
        {
            Matrix res = new Matrix(a.width, a.height);

            for (int i = 0; i < a.elements.Length; i++)
                res.elements[i] = a.elements[i] - b;

            return res;
        }
        public static Matrix operator *(Matrix a, Matrix b)
        {
            Matrix res = new Matrix(b.width, a.height);
            if (a.width == b.height)
            {
                if (b.width * b.height >= 2500 && a.height >= 80)
                {

                    Parallel.For(0, a.height, i =>
                    {
                        for (int j = 0; j < b.width; j++)
                            for (int k = 0; k < a.width; k++)
                                res[i, j] += a[i, k] * b[k, j];
                    });
                }
                else
                {
                    for (int i = 0; i < a.height; i++)
                        for (int j = 0; j < b.width; j++)
                            for (int k = 0; k < a.width; k++)
                                res[i, j] += a[i, k] * b[k, j];
                }
            }
            return res;
        }
        public static Vector operator *(Matrix a, Vector b)
        {
            if (a.width == b.Length)
            {
                Vector res = new Vector(a.height);
                if (a.width * a.height > 4000)
                {
                    Parallel.For(0, a.height, i => {
                        for (int j = 0; j < a.width; j++)
                            res[i] += a[i, j] * b[j];
                    });
                }
                else
                {
                    for (int i = 0; i < a.height; i++)
                        for (int j = 0; j < a.width; j++)
                            res[i] += a[i, j] * b[j];
                }
                return res;
            }
            return null;
        }
        public static Matrix operator /(Matrix a, double b)
        {
            Matrix result = new Matrix(a.width, a.height);
            for (int i = 0; i < a.height; i++)
                for (int j = 0; j < a.width; j++)
                    result[i, j] = a[i, j] / b;
            return result;
        }
        public static Matrix operator *(Matrix a, double b)
        {
            Matrix result = new Matrix(a.width, a.height);
            for (int i = 0; i < a.height; i++)
                for (int j = 0; j < a.width; j++)
                    result[i, j] = a[i, j] * b;
            return result;
        }
        public static Matrix operator *(double b, Matrix a)
        {
            Matrix result = new Matrix(a.width, a.height);
            for (int i = 0; i < a.height; i++)
                for (int j = 0; j < a.width; j++)
                    result[i, j] = a[i, j] * b;
            return result;
        }
        public static Matrix operator ^(Matrix a, double b)
        {
            Matrix result = new Matrix(a.width, a.height);
            for (int i = 0; i < a.height; i++)
                for (int j = 0; j < a.width; j++)
                    result[i, j] = Math.Pow(a[i, j], b);
            return result;
        }
        public Matrix(int width, int height)
        {
            this.width = width;
            this.height = height;
            elements = new double[height * width];
        }
        public Matrix(double[,] x)
        {
            this.width = x.GetLength(1);
            this.height = x.GetLength(0);
            elements = new double[height * width];

            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    elements[getInd(i, j)] = x[i, j];
        }
        public Matrix(double[] x, int width, int height)
        {
            this.width = width;
            this.height = height;

            elements = new double[width * height];

            for (int j = 0; j < x.Length; j++)
                elements[j] = x[j];
        }
        public Matrix(Vector[] vectors)
        {
            this.height = vectors.Length;
            this.width = vectors[0].Length;
            elements = new double[height * width];
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    elements[getInd(i, j)] = vectors[i].elements[j];
        }
        public Matrix(Matrix x)
        {
            this.width = x.width;
            this.height = x.height;
            elements = new double[height * width];

            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    elements[getInd(i, j)] = x[i, j];
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x">изображение с оттенками серого.</param>

        public Matrix(Vector x, int dimension)
        {
            switch (dimension)
            {
                case 0:
                    width = 1;
                    height = x.Length;
                    elements = new double[x.Length];
                    for (int i = 0; i < x.Length; i++)
                        elements[getInd(i, 0)] = x[i];
                    break;
                case 1:
                    height = 1;
                    width = x.Length;
                    elements = new double[x.Length];
                    for (int i = 0; i < x.Length; i++)
                        elements[getInd(0, i)] = x[i];
                    break;
            }
        }
        public Vector GetVector(int index, int dimension)
        {
            Vector result;
            switch (dimension)
            {
                case 0:
                    result = new Vector(height);
                    for (int i = 0; i < height; i++)
                        result[i] = elements[getInd(i, index)];
                    return result;
                case 1:
                    result = new Vector(width);
                    for (int i = 0; i < width; i++)
                        result[i] = elements[getInd(index, i)];
                    return result;
            }
            return null;
        }
        public void Set(Vector x, int index, int dimension)
        {
            switch (dimension)
            {
                case 0:
                    for (int i = 0; i < height; i++)
                        elements[getInd(i, index)] = x.elements[i];
                    break;
                case 1:
                    for (int i = 0; i < width; i++)
                        elements[getInd(index, i)] = x.elements[i];
                    break;
            }
        }
        public void Set(double[] x, int index, int dimension)
        {
            switch (dimension)
            {
                case 0:
                    for (int i = 0; i < height; i++)
                        elements[getInd(i, index)] = x[i];
                    break;
                case 1:
                    for (int i = 0; i < width; i++)
                        elements[getInd(index, i)] = x[i];
                    break;
            }
        }
        public double EuclidNorm(int index, int dimension)
        {
            double result = 0.0;
            switch (dimension)
            {
                case 0:
                    for (int i = 0; i < height; i++)
                        result += Math.Pow(elements[getInd(i, index)], 2.0);
                    break;
                case 1:
                    for (int i = 0; i < width; i++)
                        result += Math.Pow(elements[getInd(index, i)], 2.0);
                    break;
            }
            result = Math.Sqrt(result);
            return result;
        }
        public double DistanceEuclidean(Matrix b)
        {
            double result = 0.0;
            for (int i = 0; i < elements.Length; i++)
                result += Math.Pow(elements[i] - b.elements[i], 2.0); 
            
            result = Math.Sqrt(result);
            return result;
        }
        public double Distance(Matrix b)
        {
            double result = 0.0;
            for (int i = 0; i < elements.Length; i++)
                result += Math.Abs(elements[i] - b.elements[i]);

            result /= (double)elements.Length;

            return result;
        }
        public double DistanceSquare(Matrix b)
        {
            double result = 0.0;
            for (int i = 0; i < elements.Length; i++)
                result += Math.Pow(elements[i] - b.elements[i], 2.0);

            result /= (double)elements.Length;

            return result;
        }
        public Matrix Normalize(double degree)
        {
            Matrix res = new Matrix(width, height);
            double sum = 0.0;
            foreach (double value in elements)
                sum += Math.Pow(value, degree);
            sum = Math.Pow(sum, 1.0 / degree);
            for (int i = 0; i < elements.Length; i++)
                res.elements[i] = elements[i] / sum;
            return res;
        }
        public void GaussRandom(System.Random random)
        {
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    elements[getInd(i, j)] = matlib.Statistic.Random.Gauss();
        }
        public Matrix Copy()
        {
            Matrix copy = new Matrix(width, height);

            elements.CopyTo(copy.elements, 0);

            return copy;
        }
        public Vector ToVector()
        {
            Vector result = new Vector(width * height);

            for (int i = 0; i < result.Length; i++)
                result[i] = elements[i];

            return result;
        }
        public Tensor3 ToTensor3()
        {
            Tensor3 answ = new Tensor3(width, height, 1);
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    answ[0, i, j] = elements[getInd(i, j)];
            return answ;
        }
        public Tensor4 ToTensor4()
        {
            Tensor4 res = new Tensor4(width, height, 1, 1);
            for (int y = 0; y < height; y++)
                for (int x = 0; x < width; x++)
                    res[0, 0, y, x] = elements[getInd(y, x)];
            return res;
        }
        public Matrix Transpose()
        {
            Matrix result = new Matrix(height, width);
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    result[j, i] = elements[getInd(i, j)];
            return result;
        }
        public Bitmap Visualize()
        {
            Matrix elem = new Matrix(elements, width, height);
            double min = Min();
            if (min < 0)
            {
                for (int i = 0; i < height; i++)
                    for (int j = 0; j < width; j++)
                        elem[i, j] -= min;
            }
            //
            double max = elem.Max();
            double c = 255.0 / max;
            elem *= c;
            return matlib.Convertor.Picture.MatrixInBitmap(elem);
        }
        public void Sigmoida()
        {
                for (int j = 0; j < elements.Length; j++)
                    elements[j] = matlib.Function.Func.SigmoidalFunc.Sigmoida(elements[j]);
        }
        public void Threshold(double a)
        {
                for (int j = 0; j < elements.Length; j++)
                    elements[j] = matlib.Function.Func.PiecewiseDefined.Threshold(elements[j], a);
        }
        public void HyperbolicTangent()
        {
                for (int j = 0; j < elements.Length; j++)
                    elements[j] = matlib.Function.Func.SigmoidalFunc.HyperbolicTangent(elements[j]);
        }
        public void Relu()
        {
                for (int j = 0; j < elements.Length; j++)
                    elements[j] = matlib.Function.Func.PiecewiseDefined.Relu(elements[j]);
        }
        public void DsigmoidaX()
        {
            for (int i = 0; i < elements.Length; i++)
                    elements[i] = matlib.Function.FuncDerivative.Sigmoida(elements[i]);
        }
        public void DsigmoidaY()
        {
            for (int i = 0; i < elements.Length; i++)
                    elements[i] *= 1.0 - elements[i]; 
        }
        public void DhyperbolicTangentX()
        {
            for (int i = 0; i < elements.Length; i++)
                    elements[i] = matlib.Function.FuncDerivative.HyperbolicTangent(elements[i]);
        }
        public void DhyperbolicTangentY()
        {
            for (int i = 0; i < elements.Length; i++)
                    elements[i] = 1.0 - Math.Pow(elements[i], 2.0);
        }
        public Matrix Difference(Vector x, int index, int dimension)
        {
            Matrix result = new Matrix(elements, width, height);
            switch(dimension)
            {
                case 0:
                    if (height == x.Length)
                        for (int i = 0; i < height; i++)
                            result[i, index] -= x[i];
                    else return null;
                    break;
                case 1:
                    if (width == x.Length)
                        for (int i = 0; i < width; i++)
                            result[index, i] -= x[i];
                    else return null;
                    break;
            }
            return result;
        }
        public double Sum()
        {
            double result = 0.0;
                for (int i = 0; i < elements.Length; i++)
                        result += elements[i];
            return result;
        }
        public Matrix ElementsPow(double a)
        {
            Matrix result = new Matrix(width, height);
            for (int i = 0; i < elements.Length; i++)
                result.elements[i] = Math.Pow(elements[i], a);
            return result;
        }
        public void Fourier()
        {
            Matrix result = new Matrix(width, height);
            elementsF = new Complex[height, width];
            Complex i = new Complex(0.0, -1.0);
            for (int u = 0; u < width; u++)
                for (int v = 0; v < height; v++)
                {
                    Complex a = new Complex(0, 0);
                        for (int x = 0; x < width; x++)
                            for (int y = 0; y < height; y++)
                            a += elements[getInd(y, x)] * Complex.Exp(-2.0 * i * Math.PI * (u * x / (double)width + v * y / (double)height));
                    result.elements[getInd(v, u)] = a.Magnitude;
                    elementsF[v, u] = a;
                }
            elements = result.elements;
        }
        public void FFT()
        {
            int N = width * height;
            for (int i = 1; i <= 30; i++) ;
        }
        public void InverseFourier()
        {
            Complex i = new Complex(0.0, -1.0);
            double S = (double)(width * height);
            for (int u = 0; u < width; u++)
                for (int v = 0; v < height; v++)
                {
                    Complex a = new Complex(0, 0);
                    for (int x = 0; x < width; x++)
                        for (int y = 0; y < height; y++)
                            a += elementsF[y, x] * Complex.Exp(2.0 * i * Math.PI * (u * x / (double)width + v * y / (double)height));
                    a /= S;
                    elements[getInd(v, u)] = a.Magnitude;
                }
        }
        public double Max()
        {
            double Max = elements[0];
            for (int i = 0; i < elements.Length; i++)
                    if (Max < elements[i]) Max = elements[i];
            return Max;
        }
        public double MaxAbs()
        {
            double Max = Math.Abs(elements[0]);
            for (int i = 0; i < elements.Length; i++)
                if (Max < Math.Abs(elements[i])) Max = Math.Abs(elements[i]);
            return Max;
        }
        public double Min()
        {
            double Min = elements[0];
            for (int i = 0; i < elements.Length; i++)
                    if (Min > elements[i]) Min = elements[i];
            return Min;
        }
        /// <summary>
        /// Заменяет строки/столбцы
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <param name="dimension"></param>
        public void Swap(int i, int j, int dimension)
        {
            if (i != j)
            {
                double c;
                switch (dimension)
                {
                    case 0:
                        for (int k = 0; k < height; k++)
                        {
                            c = elements[getInd(k, i)];
                            elements[getInd(k, i)] = elements[getInd(k, j)];
                            elements[getInd(k, j)] = c;
                        }
                        break;
                    case 1:
                        for (int k = 0; k < width; k++)
                        {
                            c = elements[getInd(i, k)];
                            elements[getInd(i, k)] = elements[getInd(j, k)];
                            elements[getInd(j, k)] = c;
                        }
                        break;
                }
            }
        }
        public double Determinant()
        {
            double result = 1.0;
            Matrix matrix = new Matrix(elements, width, height);
            for (int i = 0; i < width; i++)
            {
                for (int j = i; j < height; j++)
                {
                    if (matrix[j, i] != 0)
                    {
                        matrix.Swap(i, j, 1);
                        break;
                    }
                    else if(j + 1 == height)
                    {
                        i++;
                        goto go;
                    }
                }
                Vector a = matrix.GetVector(i, 1);
                a *= 1.0 / a[i];
                for (int j = i + 1; j < height; j++)
                {
                    double c = matrix[j, i];
                    for (int k = i; k < width; k++)
                        matrix[j, k] -= a[k] * c;
                }
            go: { }
            }
            for (int i = 0; i < width; i++)
                result *= matrix[i, i];
                return result;
        }
        public Matrix AlgebraicComplements()
        {
            Matrix result = new Matrix(width, height);
            for (int i = 0; i < height; i++)
                for (int j = 0; i < width; i++)
                    result[i, j] = MyPow1((i + 1) * (j + 1)) * elements[getInd(i, j)];
            return result;
        }
        private int MyPow1(int degree)
        {
            if (degree % 2 == 0)
                return 1;
            return 0;
        }
    }
}
