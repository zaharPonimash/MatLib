using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MatLib;
using MatLib.Interpolation;
namespace Interpolation
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        Vector x, y, X;

        LinearInterpol interp;

        private void buttonBuild_Click(object sender, EventArgs e)
        {
            if (!interp.Learn(x, y)) throw new Exception();
            var Y = interp.Interpolation(X);
            matlib.ZedGraph.DrawZedGraphCurve(zedGraphControl2, Y.elements);
            matlib.ZedGraph.DrawZedGraphCurve(zedGraphControl1, y.elements);
            //MessageBox.Show("" + Y[Y.Length - 1] / Y[Y.Length - 2]);
        }

        private void textBoxX_TextChanged(object sender, EventArgs e)
        {
            try
            {
                var vars = textBoxX.Text.Split(new char[]{' '}, StringSplitOptions.RemoveEmptyEntries );
                x = new Vector(vars.Length);
                for (int i = 0; i < vars.Length; i++)
                    x[i] = double.Parse(vars[i]);
            }
            catch { }
        }

        private void textBoxY_TextChanged(object sender, EventArgs e)
        {
            try
            {
                var vars = textBoxY.Text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                y = new Vector(vars.Length);
                for (int i = 0; i < vars.Length; i++)
                    y[i] = double.Parse(vars[i]);
            }
            catch { }
        }

        private void textBoxTesting_TextChanged(object sender, EventArgs e)
        {
            try
            {
                var vars = textBoxTesting.Text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                X = new Vector(vars.Length);
                for (int i = 0; i < vars.Length; i++)
                    X[i] = double.Parse(vars[i]);
            }
            catch { }
        }
        Vector AverSc(Vector x, double p = 0.97)
        {
            double m = 1.0 - p;
            Vector res = new Vector(x.Length);
            res[0] = x[0];
            for (int i = 1; i < x.Length; i++)
            {
                res[i] = res[i - 1] * p + x[i] * m;
            }
            return res;
        }
        private void Form1_Load(object sender, EventArgs e)
        {
            x = Vector.GetVectorIndexes(25);
            y = new Vector(x.Length);
            X = Vector.GetVectorIndexes(x.Length + 4);
            for (int i = 0; i < y.Length; i++)
            {
                y[i] = Math.Sin(x[i] / 1.5);
            }
            var a = x.AverageValue();
            var d = Math.Sqrt(x.Dispersion());
            x -= a;
           x /= d;

           X -= a;
           X /= d;
           // x /= (double)x.Length;
           // x *= 10.0;
           // x /= (double)x.Length;
           //Hyperbolic(5);//
          // var func = new Polynomial(8);
           // interp = new LinearInterpol(func);
          //  y.Visualize();

            var seq = matlib.Function.Func.Sin(80, 5, 5, 0);
            var seq2 = matlib.Function.Func.Sin(80, 3, 3, 0.0);
            var seq3 = matlib.Function.Func.Sin(80, 4, 7, 0.0);
            seq = seq + seq2 + seq3;

            
            seq = AverSc(seq, 0.91);
            //seq = seq2;// seq.MultiplyElements(seq2) + seq + seq2;

            var meth = new InterpolTrigonom(10);
            for (int i = 0; i < seq.Length; i++) {
                //seq[i] *= Math.Sqrt(i);
            //    seq[i] *= Gausse(i, seq.Length);
            }

         
            seq.VisualizeThread("input");

            //double freq = f.GetInterval(0, f.Length / 2).IndexMax();
          //  var func = new Hyperbolic(40);
          //  LinearInterpol interp = new LinearInterpol(func);
            var inp = seq.GetInterval(0, seq.Length - 20);
            meth.Learn(Vector.GetVectorIndexes(inp.Length) / inp.Length, inp);
            meth.Calculation(Vector.GetVectorIndexes(20 + inp.Length) / inp.Length).Visualize("predict");
          //  InterpolTrigonom trig = new InterpolTrigonom((inp.Length - 1)/2);
         //   trig.Learn(inp);

           // interp.Learn(Vector.GetVectorIndexes(seq.Length), seq);
            //interp.Interpolation(Vector.GetVectorIndexes(seq.Length)).VisualizeThread();
          // trig.Interpolation(getSeq(inp.Length + 10, inp.Length)).VisualizeThread();
          //  var inp = seq.GetInterval(0, seq.Length - 10);

          //  InterpolBase.Taylor(Vector.GetVectorIndexes(inp.Length), inp, Vector.GetVectorIndexes(8+inp.Length), 27).Visualize("predict");
           // phase.Visualize();
            Application.Exit();
        }
        private const double Q = 0.5;
        public static double Gausse(double n, double frameSize)
        {
            var a = (frameSize - 1) / 2;
            var t = (n - a) / (Q * a);
            t = t * t;
            return Math.Exp(-t / 2);
        }
        Vector getSeq(int n, int m)
        {
            Vector x = new Vector(n);
            //double step = 1.0 / Math.PI;
            for (int i = 0; i < n; i++)
            {
                x[i] = -Math.PI + 2.0 * Math.PI / (double)m * i;
            }
            return x;
        }
    }
}
