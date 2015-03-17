using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace GMM_FIELD
{
    class func_vnorm
    {
        public double vnorm(Complex[] vec)
        {
            double tmp = 0; ;
            double vnorm = 0;
            for (int i = 1; i <= 3; i++)
            {
                tmp = Complex.Abs(vec[i - 1]);
                vnorm = vnorm + tmp * tmp;
            }
            return Math.Sqrt(vnorm);
        }
    }
}
