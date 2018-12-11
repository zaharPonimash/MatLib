/*
 * Создано в SharpDevelop.
 * Пользователь: Marat
 * Дата: 06.01.2018
 * Время: 22:25
 * 
 * Для изменения этого шаблона используйте меню "Инструменты | Параметры | Кодирование | Стандартные заголовки".
 */
using System;
using System.Threading.Tasks;
namespace MatLib
{
    [Serializable]
	/// <summary>
	/// Description of Tensor3.
	/// </summary>
	public class Tensor3
	{
        public Matrix[] elements { get; set; }
		public int width, height, deep;
        static public Tensor3 operator * (Tensor3 tensor3, double a)
        {
            Tensor3 res = new Tensor3(tensor3.width, tensor3.height, tensor3.deep);
            for (int z = 0; z < tensor3.deep; z++)
                for (int y = 0; y < tensor3.height; y++)
                    for (int x = 0; x < tensor3.width; x++)
                        res[z, y, x] = tensor3[z, y, x] * a;
            return res;
        }
        static public Tensor3 operator *(double a, Tensor3 tensor3)
        {
            Tensor3 res = new Tensor3(tensor3.width, tensor3.height, tensor3.deep);
            for (int z = 0; z < tensor3.deep; z++)
                for (int y = 0; y < tensor3.height; y++)
                    for (int x = 0; x < tensor3.width; x++)
                        res[z, y, x] = tensor3[z, y, x] * a;
            return res;
        }
        static public Tensor3 operator /(Tensor3 tensor3, double a)
        {
            Tensor3 res = new Tensor3(tensor3.width, tensor3.height, tensor3.deep);
            for (int z = 0; z < tensor3.deep; z++)
                for (int y = 0; y < tensor3.height; y++)
                    for (int x = 0; x < tensor3.width; x++)
                        res[z, y, x] = tensor3[z, y, x] / a;
            return res;
        }
        static public Tensor3 operator -(Tensor3 a, Tensor3 b)
        {
            if (a.width != b.width || a.height != b.height || a.deep != b.deep) throw new Exception();
            Tensor3 res = new Tensor3(a.width, a.height, a.deep);
                for (int z = 0; z < a.deep; z++)
                    for (int y = 0; y < a.height; y++)
                        for (int x = 0; x < a.width; x++)
                            res[z, y, x] = a[z, y, x] - b[z, y, x];
            
            return res;
        }
        static public Tensor3 operator +(Tensor3 a, Tensor3 b)
        {
            if (a.width != b.width || a.height != b.height || a.deep != b.deep) throw new Exception();
            Tensor3 res = new Tensor3(a.width, a.height, a.deep);
            for (int z = 0; z < a.deep; z++)
                for (int y = 0; y < a.height; y++)
                    for (int x = 0; x < a.width; x++)
                        res[z, y, x] = a[z, y, x] + b[z, y, x];
            return res;
        }
        static public Tensor3 operator +(Tensor3 a, double b)
        {
            Tensor3 res = new Tensor3(a.width, a.height, a.deep);
            for (int z = 0; z < a.deep; z++)
                for (int y = 0; y < a.height; y++)
                    for (int x = 0; x < a.width; x++)
                        res[z, y, x] = a[z, y, x] + b;
            return res;
        }
        static public Tensor3 operator +(double b, Tensor3 a)
        {
            Tensor3 res = new Tensor3(a.width, a.height, a.deep);
            for (int z = 0; z < a.deep; z++)
                for (int y = 0; y < a.height; y++)
                    for (int x = 0; x < a.width; x++)
                        res[z, y, x] = a[z, y, x] + b;
            return res;
        }
        static public Tensor3 operator -(double b, Tensor3 a)
        {
            Tensor3 res = new Tensor3(a.width, a.height, a.deep);
            for (int z = 0; z < a.deep; z++)
                for (int y = 0; y < a.height; y++)
                    for (int x = 0; x < a.width; x++)
                        res[z, y, x] = b - a[z, y, x];
            return res;
        }
        static public Tensor3 operator -(Tensor3 a, double b)
        {
            Tensor3 res = new Tensor3(a.width, a.height, a.deep);
            for (int z = 0; z < a.deep; z++)
                for (int y = 0; y < a.height; y++)
                    for (int x = 0; x < a.width; x++)
                        res[z, y, x] = a[z, y, x] - b;
            return res;
        }
		public Tensor3(int width, int height, int deep)
		{
            this.width = width;
            this.height = height;
            this.deep = deep;
            elements = new Matrix[deep];
            for(int i = 0; i < deep; i++)
            elements[i] = new Matrix(width, height);
		}
        public Tensor3(Tensor3 tensor3)
        {
            this.width = tensor3.width;
            this.height = tensor3.height;
            this.deep = tensor3.deep;

            elements = new Matrix[tensor3.deep];
            for (int i = 0; i < tensor3.deep; i++)
                elements[i] = new Matrix(tensor3.width, tensor3.height);

            Parallel.For(0, deep, z =>
            {//for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        elements[z][y, x] = tensor3[z, y, x];
            });
        }
        public Tensor3(Matrix[] tensor3)
        {
            this.width = tensor3[0].width;
            this.height = tensor3[0].height;
            this.deep = tensor3.Length;

            elements = new Matrix[tensor3.Length];
            for (int i = 0; i < tensor3.Length; i++)
            	elements[i] = new Matrix(tensor3[0].width, tensor3[0].height);
            Parallel.For(0, deep, z =>
            {//for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                		elements[z][y, x] = tensor3[z][y, x];
            });
        }
        public Tensor3(double[,,] tensor3)
        {
            this.width = tensor3.GetLength(2);
            this.height = tensor3.GetLength(1);
            this.deep = tensor3.GetLength(0);

            elements = new Matrix[tensor3.GetLength(0)];
            for (int i = 0; i < elements.Length; i++)
                elements[i] = new Matrix(tensor3.GetLength(2), tensor3.GetLength(1));
            Parallel.For(0, deep, z =>
            {//for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        elements[z][y, x] = tensor3[z, y, x];
            });
        }
		public double this[int z, int y, int x]
		{
			get{ return elements[z][y, x]; }
            set { elements[z][y, x] = value; }
		}
        public double Max()
        {
            double max = elements[0].Max();
            for(int i = 1; i < deep; i++)
            {
                double max2 = elements[i].Max();
                if (max < max2) max = max2;
            }
            return max;
        }
		public Vector ToVector()
		{
			Vector result = new Vector(width * height * deep);
			int mul = height * width;
            for(int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        result[z * mul + y * width + x] = elements[z][y, x];
			return result;
		}
        public Matrix ToMatrix()
        {
            if (deep != 1) throw new Exception();
            return new Matrix(elements[0]);
        }
        public Tensor4 ToTensor4()
        {
            Tensor4 res = new Tensor4(width, height, deep, 1);

            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        res[0, z, y, x] = elements[z][y, x];

            return res;
        }
        public Matrix GetMatrix(int deep)
        {
            return new Matrix(elements[deep]);
        }
        public void SetMatrix(Matrix x, int deep)
        {
            elements[deep] = new Matrix(x);
        }
		public Matrix Conv(Tensor3 b)
		{
			if (deep != b.deep)
				throw new Exception();

			Matrix result = new Matrix(width - b.width + 1, height - b.height + 1);

               Parallel.For(0, result.height, y =>
               //for (int y = 0; y < result.height; y++)
                {
                    //Parallel.For(0, result.width, x =>
                    for (int x = 0; x < result.width; x++)
                        for (int dy = 0; dy < b.height; dy++)
                            for (int dx = 0; dx < b.width; dx++)
                                for (int dz = 0; dz < b.deep; dz++)
                                {
                                    result[y, x] += elements[dz][y + dy, x + dx] * b[dz, dy, dx];
                                }
                    //});
                });

			return result;
		}
        public double[,,] GetElem()
        {
            double[,,] res = new double[deep, height, width];
            if (deep > 2)
                Parallel.For(0, deep, i =>
                {
                    for (int j = 0; j < height; j++)
                        for (int k = 0; k < width; k++)
                            res[i, j, k] = elements[i][j, k];
                });
            else
            {
                for (int i = 0; i < deep; i++)
                    for (int j = 0; j < height; j++)
                        for (int k = 0; k < width; k++)
                            res[i, j, k] = elements[i][j, k];
            }
            return res;
        }
        public void Normalize(double degree)
        {
            double sum = 0.0;
            for (int i = 0; i < deep; i++)
                for (int j = 0; j < height; j++)
                    for (int k = 0; k < width; k++)
                        sum += Math.Pow(elements[i][j, k], degree);
            sum = Math.Pow(sum, 1.0 / degree);
            for (int i = 0; i < deep; i++)
                for (int j = 0; j < height; j++)
                    for (int k = 0; k < width; k++)
                        elements[i][j, k] = elements[i][j, k] / sum;
        }
        public void Sigmoida()
        {
            //Parallel.For(0, deep, z =>
            {for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        elements[z][y, x] = matlib.Function.Func.SigmoidalFunc.Sigmoida(elements[z][y, x]);
            }//);
        }
        public Tensor3 ElementsPow(double pow)
        {
            Tensor3 res = new Tensor3(width, height, deep);
            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                    {
                        res[z, y, x] = Math.Pow(elements[z][y, x], pow);
                        //if (double.IsNaN(res[z, y, x])) res[z, y, x] = 0.0;
                    }
            return res;
        }
        public Tensor3 Abs()
        {
            Tensor3 res = new Tensor3(width, height, deep);
            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                    {
                        res[z, y, x] = Math.Abs(elements[z][y, x]);
                        //if (double.IsNaN(res[z, y, x])) res[z, y, x] = 0.0;
                    }
            return res;
        }
        public Tensor3 ElementsDivision(Tensor3 b)
        {
            if (width != b.width || height != b.height || deep != b.deep) throw new Exception();
            Tensor3 res = new Tensor3(width, height, deep);
            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        res[z, y, x] = elements[z][y, x] / b[z, y, x];
            return res;
        }
        public Tensor3 ElementsMultiply(Tensor3 b)
        {
            if (width != b.width || height != b.height || deep != b.deep) throw new Exception();
            Tensor3 res = new Tensor3(width, height, deep);
            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        res[z, y, x] = elements[z][y, x] * b[z, y, x];
            return res;
        }
        public void GaussRandom()
        {
            for (int z = 0; z < deep; z++)
                for (int y = 0; y < height; y++)
                    for (int x = 0; x < width; x++)
                        elements[z][y, x] = matlib.Statistic.Random.Gauss();
        }
	}
}
