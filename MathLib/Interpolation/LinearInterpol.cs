using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MatLib.Classifier;
namespace MatLib.Interpolation
{
    public class LinearInterpol
    {
        LinearClassifier classifier;
        public IFunction func;
        public LinearInterpol(IFunction func = null)
        {
            if(func == null) {
                func = new Polynomial(3);
            }
            this.func = func;
            classifier = new LinearClassifier(func.n);
        }
        public double Interpolation(double x)
        {
            if (classifier.w == null) return double.NaN;

            var inps = func.Calculation(x);

            return classifier.Calculation(inps);
        }
        public Vector Interpolation(Vector x)
        {
            if (classifier.w.IsNaN()) return null;

            Vector y = new Vector(x.Length);
            for (int i = 0; i < x.Length; i++)
            {
                var inps = func.Calculation(x[i]);
                y[i] = classifier.Calculation(inps);
            }
            return y;
        }
        public bool Learn(Vector x, Vector y)
        {
            Vector[] x_ = new Vector[x.Length];
            for (int i = 0; i < x_.Length; i++) {
                x_[i] = func.Calculation(x[i]);
            }

           return classifier.Learn(x_, y);

        }
    }
}
