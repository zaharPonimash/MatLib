using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Regression
{
    public class LinearRegression
    {
        public double b, k;
        public Vector reg;
        private double gradB, gradK;
        public LinearRegression()
        {
            b = matlib.random.NextDouble();
            k = matlib.random.NextDouble();
        }
        public double Regression(double x)
        {
            return k * x + b;
        }
        public Vector Regression(Vector x, Vector y)
        {
            k = Vector.Cov(x, y) / x.Dispersion();
            b = y.AverageValue() - k * x.AverageValue();
            reg = new Vector(x.Length);
            for (int i = 0; i < reg.Length; i++) { reg[i] = Regression(x[i]);
            }

            return reg;
        }
        double dEdb(Vector x, Vector y)
        {
            double res = 0.0;
            for (int i = 0; i < y.Length; i++)
            {
                res += Regression(x[i]) - y[i];
            }
            return res;
        }
        double dEdk(Vector x, Vector y)
        {
            double res = 0.0;
            for (int i = 0; i < y.Length; i++)
            {
                res += (Regression(x[i]) - y[i]) * x[i];
            }
            return res;
        }
        public void Mds(Vector x, Vector y, int iter = 1000, double n = 0.001)
        {
            for (int i = 0; i < iter; i++)
            {
                var gradB = dEdb(x, y);
                var gradK = dEdk(x, y);
                if (this.gradB == 0 || this.gradK == 0)
                {
                    b -= gradB * n;
                    k -= gradK * n;
                }
                else
                {
                    double w = Math.Sqrt(gradB * gradB + gradK * gradK) / Math.Sqrt(this.gradB * this.gradB + this.gradK * this.gradK);
                    b -= n * (gradB + w * gradB);
                    k -= n * (gradK + w * gradK);
                    //  (gradB)
                }
                this.gradB = gradB;
                this.gradK = gradK;
            }
        }
        public double angle()
        {
            return 180.0 * Math.Atan(k) / Math.PI;
        }
    }
}
