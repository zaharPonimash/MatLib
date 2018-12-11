using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Interpolation
{
    class InterpolGauss
    {
               public Vector a, b;
       public int n;
       /// <summary>
       /// 
       /// </summary>
       /// <param name="n">кол-во гармонических составляющих</param>
       public InterpolGauss(int n = 2)
       {
           this.n = n;
           a = new Vector(n);
           b = new Vector(n);
       }
       void CalcA(Vector y)
       {

           for (int i = 0; i < n; i++)
           {
               double S = 0;
               for (int j = -n; j < n + 1; j++)
               {
                   S += y[j + n] * Math.Cos(2.0 * Math.PI * (double)(i * j) / (double)(2.0 * n + 1.0));
               }
               if (i == 0) S = S / (2.0*n+1.0);
               else S = 2.0 * S / (2.0 * n + 1.0);

               a[i] = S;
           }
       }
       void CalcB(Vector y)
       {
           for (int i = 0; i < n; i++)
           {
               double S = 0;
               //int ii;
               for (int j = -n; j < n + 1; j++)
               {
                   //ii = -n + j;
                   S += y[j + n] * Math.Sin(2.0 * Math.PI * (double)(i * j) / (double)(2.0 * n + 1.0));
               }
               S *= 2.0 / (2.0 * n + 1.0);

               b[i] = S;
           }
       }
       public void Learn(Vector y)
       {
           Vector y_ = new Vector(2 * n + 1);
           int step = y.Length / y_.Length;

           int j = 0;
           for (int i = 0; i < y_.Length; i++)
           {
               y_[i] = y.GetInterval(j, j + step - 1).AverageValue();
               j += step;
           }

           CalcA(y_);
           CalcB(y_);
       }
       public double Interpolation(double x)
       {
           double S = a[0];
           for (int i = 1; i < n; i++)
           {
               S += a[i] * Math.Cos((double)i * x) + b[i] * Math.Sin((double)i * x);
           }
           return S;
       }
       public Vector Interpolation(Vector x)
       {
           Vector res = new Vector(x.Length);
           for (int i = 0; i < x.Length; i++)
           {
               res[i] = Interpolation(x[i]);
           }
           return res;
       }
    }
}
