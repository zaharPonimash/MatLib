using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatLib.Interpolation
{
   public interface IFunction
    {
        int n { get; set; }
        Vector y { get; set; }
        double x { get; set; }
        Vector Calculation(double x);
    }
}
