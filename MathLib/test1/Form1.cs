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
using GeneticAlgorithm;
namespace test1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }
        Genetic genetic;
        private void Form1_Load(object sender, EventArgs e)
        {
            genetic = new Genetic(7);
            double max = genetic.LivingCycle(FitnessFunc);
            for (int i = 0; i < 10000; i++)
            {

                var a = genetic.LivingCycle(FitnessFunc);
                if (max < a) max = a;
            }

                MessageBox.Show("" + Math.Round(max, 4));
        }
        double FitnessFunc(double x)
        {
            return -x * x;//(Math.Sin(x) + 2.0) / (Math.Cos(x) + 2.0);
        }
    }
}
