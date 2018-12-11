using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Interpolation
{
    public class InterpolTrigonom
    {
        public int n;

        public Vector pow = new Vector(21);

        public Vector fi, k, A;
        public  double f;
        const double q = 2.0 * Math.PI;
        public InterpolTrigonom(int n = 10)
        {
            this.n = n;
            k = new Vector(n);
            fi = new Vector(n);
            A = new Vector(n);
            pow[0] = 1;
            for (int i = 1; i < pow.Length; i++) {
                pow[i] = pow[i - 1] * 2;
            }
        }
        int Search(int n)
        {
            int l = -1;
            for (int i = 0; i < pow.Length; i++)
            {
                if (pow[i] >= n)
                {
                    l = (int)pow[i];
                    break;
                }
            }
            return l;
        }
        public double Calculation(double x)
        {
            f = 0.0;
            for (int i = 0; i < n; i++) {
                f += A[i] * Math.Sin(q * k[i] * x + fi[i]);
            }
            return f;
        }
        public Vector Calculation(Vector x)
        {
            var res = new Vector(x.Length);
            for (int j = 0; j < res.Length; j++)
            {
                res[j] = 0;
                for (int i = 0; i < n; i++)
                {
                    res[j] += A[i] * Math.Cos(q * k[i] * x[j] + fi[i]);
                }

            }
            return res;
        }
        public Vector Learn(Vector x, Vector y)
        {
            var y_ = new Vector(y);
            y_.Resize(Search(y.Length));

            var cmplx = matlib.Transofm.FFT.fft(y_);
            var freq = matlib.Convertor.Complex_.ComplexInMagnitude(cmplx);
            freq.Resize(freq.Length / 2);
            //freq = freq.GetInterval(1, freq.Length - 1);
            Vector ind = Vector.GetVectorIndexes(freq.Length, 0);
            var res = matlib.Sort.qsort(ind, freq);
            ind = res[0];

            for (int i = 0; i < n; i++) {
                k[i] = ind[ind.Length - i - 1]; 
            }

            for (int i = 0; i < n; i++)
            {
                A[i] = freq[(int)ind[ind.Length - 1 - i]] * 2.0 / x.Length;
                if ((int)ind[ind.Length - 1 - i] == 0) A[i] /= 2.0;
            }

            for (int i = 0; i < n; i++)
            {
                fi[i] = cmplx[(int)ind[ind.Length - 1 - i]].Phase;// *2.0 / x.Length;
            }
            return k;
        }

    }
}
