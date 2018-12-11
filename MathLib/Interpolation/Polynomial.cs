using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Interpolation
{
    public class Polynomial : IFunction
    {
        public int n { get; set; }

        public Vector y { get; set; }

        public double x { get; set; }

        public Polynomial(int n)
        {
            this.n = n;
            y = new Vector(n);
        }

        public Vector Calculation(double x)
        {
            this.x = x;
            y = new Vector(y.Length);
            for (int i = 0; i < n; i++) {
                y[i] = Math.Pow(x, i);
            }

            return y;
        }
    }
}
