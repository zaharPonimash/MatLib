using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Interpolation
{
    public class Hyperbolic : IFunction
    {

         public int n { get; set; }

        public Vector y { get; set; }

        public double x { get; set; }

        public Hyperbolic(int n)
        {
            this.n = n;
            y = new Vector(n);
        }

        public Vector Calculation(double x)
        {
            this.x = x;
            y = new Vector(y.Length);
            for (int i = 0; i < n; i++) {
                y[i] = x*Math.Sin(2.0 * Math.PI * i);// Math.Sqrt(Math.Abs(x) * i) * Math.Sign(x);//matlib.Function.Func.SigmoidalFunc.HyperbolicTangent(Math.Sqrt(i) * x) * i;
            }

            return y;
        }

    }
}
