using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Interpolation
{
    public class InterpolBase
    {
        public static Vector Lagrange(Vector x, Vector y, Vector X, int n = -1)
        {
            if (n == -1) n = x.Length;

            Vector xPreproc = new Vector(n);
            Vector yPreproc = new Vector(xPreproc.Length);
            int sh = x.Length / n;
            int b = 0;
            for (int i = 0; i < n; i++)
            {
                xPreproc[i] = x.GetInterval(b, b + sh - 1).AverageValue();
                yPreproc[i] = y.GetInterval(b, b + sh - 1).AverageValue();
                b += sh;
            }



            Vector res = new Vector(X.Length);
            Vector[] l = new Vector[X.Length];
            for (int i = 0; i < l.Length; i++)
                l[i] = new Vector(xPreproc.Length) + 1;

            for (int i = 0; i < l.Length; i++)
            {
                for (int j = 0; j < l[i].Length; j++)
                {
                    for (int k = 0; k < xPreproc.Length; k++)
                        if (k != j)
                            l[i][j] *= (X[i] - xPreproc[k]) / (xPreproc[j] - xPreproc[k]);
                    //l[i][j] /= z[j];
                }
            }
            for (int i = 0; i < l.Length; i++)
                for (int j = 0; j < l[i].Length; j++)
                    res[i] += l[i][j] * yPreproc[j];

            return res;
        }

        private static fx[] SearchNearestNum(fx[] x, double val)
        {
            fx[] fxs = new fx[2];

            fxs[0] = new fx(); fxs[1] = new fx();

            #region fxs[0]
            fxs[0].x = x[0].x;
            fxs[0].f = x[0].f;
            int ind = 0;

            for (int i = 1; i < x.Length; i++)
            {
                if (Math.Abs(fxs[0].x - val) > Math.Abs(x[i].x - val))
                {
                    fxs[0].f = x[i].f;
                    fxs[0].x = x[i].x;
                    ind = i;
                }
            }
            #endregion

            #region fsx[1]
            if (ind != 0)
            {
                fxs[1].x = x[0].x;
                fxs[1].f = x[0].f;
            }
            else
            {
                fxs[1].x = x[1].x;
                fxs[1].f = x[1].f;
            }
            for (int i = 0; i < x.Length; i++)
            {
                if (Math.Abs(fxs[1].x - val) > Math.Abs(x[i].x - val) && i != ind)
                {
                    fxs[1].f = x[i].f;
                    fxs[1].x = x[i].x;
                }
            }
            #endregion


            if (fxs[0].x < fxs[1].x)
            {
                double X = fxs[0].x;
                double F = fxs[0].f;

                fxs[0].x = fxs[1].x;
                fxs[0].f = fxs[1].f;

                fxs[1].x = X;
                fxs[1].f = F;

            }
            return fxs;
        }
        static fx[] VectorToFx(Vector x)
        {
            fx[] fxs = new fx[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                fxs[i] = new fx();
                fxs[i].x = i;
                fxs[i].f = x[i];
            }
            return fxs;
        }
        struct fx
        {
            public double f;
            public double x;
        }
        private static double SearchNearestVal(fx[] derivs, double x)
        {
            double value = derivs[0].f;
            double dist = Math.Abs(derivs[0].x - x);
            for (int i = 1; i < derivs.Length; i++)
            {
                if (dist > Math.Abs(derivs[i].x - x))
                {
                    dist = Math.Abs(derivs[i].x - x);
                    value = derivs[i].f;
                }
            }
            return value;
        }
        private static int SearchNearestInd(Vector vec, double x)
        {
            int ind = -1;
            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] < x) ind++;
                else
                {
                    if (ind == vec.Length - 1) return vec.Length - 1;
                    if (ind == -1) return 0;

                    if (x - vec[ind] < vec[ind + 1] - x)
                        return ind;
                    else return ind + 1;
                }
            }
            return ind;
        }
        private static int SearchNearestInd(fx[] vec, double x)
        {
            int ind = -1;
            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i].x < x) ind++;
                else
                {
                    if (ind == vec.Length - 1) return vec.Length - 1;
                    if (ind == -1) return 0;

                    if (x - vec[ind].x < vec[ind + 1].x - x)
                        return ind;
                    else return ind + 1;
                }
            }
            return ind;
        }
        public static Vector Taylor(Vector x, Vector y, Vector X, int accuracy = 3)
        {
            Vector xPreproc = new Vector(accuracy);
            Vector yPreproc = new Vector(xPreproc.Length);
            int sh = x.Length / accuracy;
            int b = 0;
            for (int i = 0; i < accuracy; i++)
            {
                xPreproc[i] = x.GetInterval(b, b + sh - 1).AverageValue();
                yPreproc[i] = y.GetInterval(b, b + sh - 1).AverageValue();
                b += sh;
            }
            //x = xPreproc;
            //y = yPreproc;

            Vector res = new Vector(X.Length);
            if (x.Length < accuracy)
                accuracy = x.Length;

            fx[][] derivs = new fx[accuracy - 1][];
            for (int i = 0; i < derivs.Length; i++)
                derivs[i] = new fx[x.Length - i - 1];


            for (int i = 0; i < derivs.Length; i++)
                for (int j = 0; j < derivs[i].Length; j++)
                    derivs[i][j] = new fx();

            if (derivs.Length > 0)
            {
                for (int i = 1; i < x.Length; i++)
                {
                    derivs[0][i - 1].f = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
                    derivs[0][i - 1].x = x[i - 1] + (x[i] - x[i - 1]) / 2.0;
                }
                for (int i = 1; i < derivs.Length; i++)
                    for (int j = 1; j < derivs[i - 1].Length; j++)
                    {
                        derivs[i][j - 1].x = derivs[i - 1][j - 1].x + (derivs[i - 1][j].x - derivs[i - 1][j - 1].x) / 2.0;
                        derivs[i][j - 1].f = (derivs[i - 1][j].f - derivs[i - 1][j - 1].f) / (derivs[i - 1][j].x - derivs[i - 1][j - 1].x);
                    }
            }
            for (int i = 0; i < X.Length; i++)
            {
                Vector deriv = new Vector(accuracy);
                var index = SearchNearestInd(x, X[i]);
                deriv[0] = y[index];
                double a = x[index];
                Vector A = new Vector(accuracy);
                A[0] = a;
                for (int k = 1; k < A.Length; k++)
                {
                    A[k] = derivs[k - 1][SearchNearestInd(derivs[k - 1], X[i])].x;
                }
                for (int j = 0; j < accuracy - 1; j++)
                {
                    deriv[j + 1] = SearchNearestVal(derivs[j], X[i]);
                }
                //res[i] = matlib.Function.Polynomials.Taylor(X[i], a, deriv);

                if (accuracy < 1) return null;

                double result = 0.0;
                for (int n = 0; n < deriv.Length; n++)
                    result += deriv[n] * Math.Pow(X[i] - A[0], n) / (double)matlib.Factorial(n);

                res[i] = result;
            }
            return res;
        }
    }
}
