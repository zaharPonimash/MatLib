using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Classifier
{
    class LinearClassifier
    {
        public int n;
        public Vector w;
        public Vector x;
        public LinearClassifier(int n)
        {
            this.n = n;
        }
        public double Calculation(Vector x)
        {
            this.x = x;
            return x * w;
        }

        public Vector Calculation(Vector[] x)
        {
            Vector res = new Vector(x.Length);
            for (int i = 0; i < x.Length; i++)
            {
                res[i] = x[i] * w;
            }
            return w;
        }

        public bool Learn(Vector[] x, Vector y)
        {
            var A = new Matrix(n, n);
            var b = new Vector(n);

                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < x.Length; k++) {
                        b[j] += y[k] * x[k][j];
                    }
                }

                for (int q = 0; q < n; q++)
                    for (int m = 0; m < n; m++)
                        for (int k = 0; k < x.Length; k++) {
                            A[q, m] += x[k][m] * x[k][q];
                        }

                try {
                    w = matlib.SolvingSystems.SolvingSLAY(A, b);

                    return true;

                } catch { return false; }
        }
    }
}
