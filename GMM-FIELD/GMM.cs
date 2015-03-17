using System;
using System.Numerics;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace GMM_FIELD
{
    public partial class GMM : Form
    {
        StreamWriter OUT = new StreamWriter("OUT.txt");
        StreamWriter gmm_out = new StreamWriter("gmm_out.txt");
        int lP = 0;//количество частиц
        double[,] cnv, dc;
        int nmp, ni0, nmp0, np, uvmax;
        Complex[] Emn, escati, escatc;
        Complex[,] Nmn3, Mmn3;
        double[] besj, besy, p, pi, tau, cof0, cofsr, bcof, fnr, ga0;
        double[,] smue;
        Complex[, ,] atr;
        int[] nmax, iga0;
        double pih, twopi, d, xt, st, sphi, cphi, cn, dn;
        public GMM()
        {
            InitializeComponent();
        }
        
        private void GMM_Load(object sender, EventArgs e)
        {
      
        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {
            try
            {
                if (textBox1.Text != "")
                {
                    dataGridView1.RowCount = 0;
                    textBox1.BackColor = Color.White;
                    lP = Convert.ToInt32(textBox1.Text);
                    dataGridView1.RowCount = lP;  //количество строк таблицы
                }
                button3.Enabled = true;
            }
            catch (FormatException)
            {
                textBox1.Text = "";
                textBox1.BackColor = Color.Red;
                MessageBox.Show("Необходимо ввести целое положительное число!","Ошибка!");      
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            button2.Enabled = button5.Enabled = button6.Enabled = button3.Enabled=true;
            button1.Enabled = false;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            bool f=true;
            for (int i = 0; i < Convert.ToInt32(textBox1.Text); i++)
                for (int j = 0; j < 6; j++)
                    if (("").CompareTo(Convert.ToString(dataGridView1[j, i].Value))==0)
                        f = false;
            if (f)
            {
                StreamWriter s = new StreamWriter("GMM_out_d.txt");
                for (int i = 0; i < Convert.ToInt32(textBox1.Text); i++)
                {
                    s.WriteLine("Информация о сфере №" + (i + 1) + ":");
                    s.Write("x= {0}; y= {1}; z= {2}; ", dataGridView1[0, i].Value, dataGridView1[1, i].Value, dataGridView1[2, i].Value);
                    s.WriteLine("r= {0}; Re(n)= {1}; Im(n)= {2}", dataGridView1[3, i].Value, dataGridView1[4, i].Value, dataGridView1[5, i].Value);
                    s.WriteLine();
                }
                s.Close();
                MessageBox.Show("Данные успешно сохранены!");
            }
            else
                MessageBox.Show("Введены не все данные!", "Предупреждение!");
        }

        private void button3_Click(object sender, EventArgs e)
        {

            button2.Enabled = button5.Enabled = button6.Enabled = false;
            button1.Enabled = true;
            for (int i = 0; i < Convert.ToInt32(textBox1.Text); i++)
                for (int j = 0; j < 6; j++)
                    dataGridView1[j, i].Value = "";
            dataGridView1.RowCount = 0;
            textBox1.Text = "";
            button3.Enabled = false;
        }

        private void button4_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void button5_Click(object sender, EventArgs e)
        {
            StreamReader gmm = new StreamReader("gmm_in.txt");
            string GM = gmm.ReadToEnd();
            string[] spl = GM.Split('\n');
            string[] obj = spl[0].Split(';');
            np = Convert.ToInt32(obj[0]);
            int nLp = Convert.ToInt32(obj[1]);
            const int NXMAX = 3000, nangmax = 181, MOR = 181, ncmax = 180;
            nmp = np * (np + 2);
            nmp0 = (np + 1) * (np + 4) / 2;
            ni0 = np * (np + 1) * (2 * np + 1) / 3 + np * np;
            int ng0 = (int)(np * (2 * Math.Pow(np, 3) + 10 * Math.Pow(np, 2) + 19 * np + 5) / 6);
            int nrc = 4 * np * (np + 1) * (np + 2) / 3 + np;
            int nij = nLp * (nLp - 1) / 2;
            #region
            atr = new Complex[2, np, nmp];
            cnv = new double[np, np];
            Emn = new Complex[nmp0];
            escati = new Complex[3];
            escatc = new Complex[3];
            Nmn3 = new Complex[3, nmp0];
            Mmn3 = new Complex[2, nmp0];
            besj = new double[2 * np + 2];
            besy = new double[2 * np + 2];
            p = new double[np + 1];
            int u, v, u0;
            int[] uvmax = new int[nLp], ind = new int[nLp];
            nmax = new int[nLp];
            pih = Math.Acos(0);
            double k, gcs = 0, gcv = 0, idpq = 0, pione = 2 * pih;
            double[,] r0 = new double[6, nLp], r00 = new double[3, nLp];
            double[] x = new double[nLp], dang = new double[nangmax], c0i = new double[nLp], c1i = new double[nLp];
            double[] rsr0 = new double[NXMAX], rsi0 = new double[NXMAX], rsx0 = new double[NXMAX], px0 = new double[NXMAX];
            
            double[,] rsr = new double[np, nLp], rsi = new double[np, nLp], rsx = new double[np, nLp], px = new double[np, nLp];
            double[] betar = new double[MOR], thetr = new double[MOR], phair = new double[MOR];
            smue = new double[4, 4];
            double[,] drot = new double[nrc, nij];
            double[, , ,] mue = new double[4, 4, ncmax, nangmax];
            double[] i11 = new double[nangmax], i12 = new double[nangmax], i21 = new double[nangmax], i22 = new double[nangmax];
            double[] pol = new double[nangmax], inat = new double[nangmax];
            double[] cscaxi = new double[nLp], cscayi = new double[nLp], cextxi = new double[nLp], cabsxi = new double[nLp];
            double[] cextyi = new double[nLp], cabsyi = new double[nLp], cscai = new double[nLp], cexti = new double[nLp];
            double[] cabsi = new double[nLp], assymxi = new double[nLp], assymyi = new double[nLp], assymi = new double[nLp];
            double[] cprxi = new double[nLp], cpryi = new double[nLp], cpri = new double[nLp];
            Complex A, B, cmz, A0, B0, Aj, Bj, Aj2, Bj2, A2=new Complex(0,0), B2=new Complex(0,0), ephi, ci, cin;
            Complex[,] atr0 = new Complex[ni0, nij], btr0 = new Complex[ni0, nij], atr1 = new Complex[ni0, nij], btr1 = new Complex[ni0, nij];
            Complex[] at = new Complex[nmp], bt = new Complex[nmp], reff = new Complex[nLp], an = new Complex[np], bn = new Complex[np];
            Complex[,] ek = new Complex[np, nij], p0 = new Complex[nLp, nmp], q0 = new Complex[nLp, nmp];
            Complex[,] aMie = new Complex[nLp, np], bMie = new Complex[nLp, np], ass = new Complex[nLp, nmp], bs = new Complex[nLp, nmp];
            Complex[,] as0 = new Complex[nLp, nmp], bs0 = new Complex[nLp, nmp], asp = new Complex[nLp, nmp], bsp = new Complex[nLp, nmp];
            Complex[,] as1 = new Complex[nLp, nmp], bs1 = new Complex[nLp, nmp], asc = new Complex[nLp, nmp], bsc = new Complex[nLp, nmp];
            Complex[,] ast = new Complex[nLp, nmp], bst = new Complex[nLp, nmp], asv = new Complex[nLp, nmp], bsv = new Complex[nLp, nmp];
            Complex[,] s2x = new Complex[ncmax, nangmax], s4x = new Complex[ncmax, nangmax], s3y = new Complex[ncmax, nangmax], s1y = new Complex[ncmax, nangmax];
            Complex[] atj = new Complex[nmp], btj = new Complex[nmp], py0 = new Complex[NXMAX], py = new Complex[NXMAX], dpy = new Complex[NXMAX];
            bcof = new double[np + 3];
            dc = new double[np * 2 + 1, nmp+1];
            fnr = new double[2 * np + 5];
            pi = new double[nmp0];
            tau = new double[nmp0];
            iga0 = new int[ni0];
            ga0 = new double[ng0];
            cof0 = new double[ni0];
            cofsr = new double[nmp];
            pih = Math.Acos(0);
            twopi = 4 * pih;
            ci = new Complex(0, 1);
            cin = new Complex(0, -1);
            #endregion

            #region
            //файл "gmm01f.in"
            obj = spl[1].Split(';');
            int nbeta = Convert.ToInt32(obj[0]), nthet = Convert.ToInt32(obj[1]), nphai = Convert.ToInt32(obj[2]);

            OUT.WriteLine("nbeta,nthet,nphai: " + nbeta + ", " + nthet + ", " + nphai);

            if (nbeta > MOR || nthet > MOR || nphai > MOR)
            {
                OUT.WriteLine("***  parameter MOR too small  ***");
                OUT.WriteLine("MOR must be >nbeta,nthet,nphai given above");
                OUT.WriteLine("Please change MOR in the parameter line of the");
                OUT.WriteLine("main code, recompile, then try again");
                OUT.Close();
                Close();
            }
            //double idran;
            double nram;
            //if (nbeta * nphai > 1)
            //    idran = 1;
            //else
            //    idran = 0;
            nram = nbeta * nthet * nphai;
            if (nram < 1)
            {
                OUT.WriteLine("please check (nbeta,nthet,nphai) in gmm01f.in");
                OUT.Close();
                Close();
            }
            obj = spl[2].Split(';');
            double betami = Convert.ToDouble(obj[0]), betamx = Convert.ToDouble(obj[1]);
            double thetmi = Convert.ToDouble(obj[2]), thetmx = Convert.ToDouble(obj[3]);
            double phaimi = Convert.ToDouble(obj[4]), phaimx = Convert.ToDouble(obj[5]);
            if (nbeta == 1)
                betamx = betami;
            if (nthet == 1)
                thetmx = thetmi;
            if (nphai == 1)
                phaimx = phaimi;

            OUT.WriteLine("Ranges of Euler angles: " + betami + " " + betamx + " " + thetmi + " " + thetmx + " " + phaimi + " " + phaimx);

            obj = spl[3].Split(';');
            double idMie = Convert.ToDouble(obj[0]);

            OUT.WriteLine("idMie: " + idMie);

            int idd = Convert.ToInt32(obj[1]);
            if(nram==1&&idd==1)
            {
                OUT.WriteLine("Warning: idd=1 while calculating a single");
                OUT.WriteLine("particle-orientation, the results will be an");
                OUT.WriteLine("average over two particle orientations");
                OUT.WriteLine("please check: idd=0 or idd=1?");
            }
            obj = spl[4].Split(';');
            int idc = Convert.ToInt32(obj[0]);
            int iseed = Convert.ToInt32(obj[1]);
            if (idpq == 1)
            {
                nram = nphai + nthet;
                idc = 0;
                idd = 0;
                nbeta = 1;
                betami = 0;
                betamx = betami;
                thetr[nthet] = 0;
                phair[nphai] = 0;
            }

            OUT.WriteLine("idd: " + idd);

            if (idd == 1)
                nram = 2 * nram;

            OUT.WriteLine("# of orientations to be averaged: " + nram);

            if(idc<0)
                OUT.WriteLine("idc,iseed: " + idc+"  "+iseed);
            else
                OUT.WriteLine("idc: " + idc);

            obj = spl[5].Split(';');
            double factor1 = Convert.ToDouble(obj[0]),
                    factor2 = Convert.ToDouble(obj[1]),
                    MXINT = Convert.ToDouble(obj[2]);

            OUT.WriteLine("Numerical factors for convergence: " + factor1 + "  " + factor2);
            OUT.WriteLine("Maximum iterations allowed: " + MXINT);

            obj = spl[6].Split(';');
            int NADD = Convert.ToInt32(obj[0]);

            OUT.WriteLine("Scat. orders added to Wiscombe criterion: " + NADD);

            double eps = 1/Math.Pow(10, 20), small = 1/Math.Pow(10, 10);

            OUT.WriteLine("error tolerance for Mie-expansions: " + eps);
            OUT.WriteLine("Convergence criterion: " + small);

            obj = spl[7].Split(';');
            double fint = Convert.ToDouble(obj[0]);

            if (fint < 0 || fint > 1)
            {
                fint = 0.02;
                OUT.WriteLine("Interaction index: using default 0.02");
            }
            else
                OUT.WriteLine("Interaction index: " + fint);

            obj = spl[8].Split(';');
            double sang = Convert.ToDouble(obj[0]), pang = Convert.ToDouble(obj[1]);

            OUT.WriteLine("scat.-angle-interval in output: " + sang);

            if (sang <= 0)
                sang = 1;
            double nang = 90 / sang + 1;
            double nang2 = 2 * nang - 1;
            if (nang2 > nangmax)
            {
                OUT.WriteLine("sang too small");
                OUT.WriteLine("please increase sang in the input file gmm01f.in");
                OUT.WriteLine("and try again, or");
                OUT.WriteLine("increase nangmax in the parameter line of the");
                OUT.WriteLine("main code, recompile, then try again");
                OUT.Close();
                Close();
            }

            OUT.WriteLine("azimuth-angle-interval in Mueller matrix output: " + pang);

            double npng;
            if (pang < 0.0001)
                npng = 1;
            else
                npng = 360 / pang;
            if (npng > ncmax)
            {
                OUT.WriteLine("pang too small");
                OUT.WriteLine("please increase pang in the input file gmm01f.in");
                OUT.WriteLine("and try again, or increase ncmax in the parameter");
                OUT.WriteLine("line of the main code, recompile, then try again");
                OUT.Close();
                return;
            }
            betami = betami * pih / 90;
            betamx = betamx * pih / 90;
            thetmi = thetmi * pih / 90;
            thetmx = thetmx * pih / 90;
            phaimi = phaimi * pih / 90;
            phaimx = phaimx * pih / 90;
            if (idc > 0)
                orientcd(betami, betamx, thetmi, thetmx, phaimi, phaimx, nbeta, nthet, nphai, betar, thetr, phair);
            else
                orientud(betami, betamx, thetmi, thetmx, phaimi, phaimx, nbeta, nthet, nphai, betar, thetr, phair);
            if(idMie==1)
            {
                OUT.WriteLine("*** Calculating only coherent Mie-scattering ***");
                OUT.WriteLine("*** No interaction included ********************");
            }
            // файл "Ag-Si-2s-405nm.k"
            obj = spl[9].Split(';');
            double w = Convert.ToDouble(obj[0]);
            int nL = 2;    //textBox1.Text
            if (nL > nLp)
            {
                OUT.WriteLine("Parameter nLp too small, must be > " +nL);
                OUT.WriteLine("Change nLp in gmm01f.par, recompile, then try again");
                OUT.Close();
                Close();
            }
            if (nL == 1)
                idMie = 1;

            r0[0,0]=0;
            r0[1,0]=0;
            r0[2,0]=-0.09;
            r0[3,0]=0.006;
            r0[4, 0] = 0.13902;
            r0[5, 0] = 1.38294;
            r0[0,1]=0;
            r0[1,1]=0;
            r0[2, 1] = 0.02;
            r0[3, 1] = 0.06;
            r0[4, 1] = 1.4928;
            r0[5, 1] = 0.5114;
            //for (int i = 0; i < nL; i++)
            //    for (int j = 0; j < 6; j++)
            //        r0[j, i] = Convert.ToDouble(dataGridView1[j, i].Value);
            double x0, y0, z0;
            double gcsr, gcvr;
            for (int i = 1; i <= nL; i++)
            {
                x0 = r0[0, i - 1];
                y0 = r0[1, i - 1];
                z0 = r0[2, i - 1];
                r00[0, i - 1] = x0;
                r00[1, i - 1] = y0;
                r00[2, i - 1] = z0;
                if (r0[5, i - 1] > 0)
                    r0[5, i - 1] = -r0[5, i - 1];
                if (r0[4, i - 1] != 1 || r0[5, i - 1] != 0)
                {
                    gcs = gcs + r0[3, i - 1] * r0[3, i - 1];
                    gcv = gcv + r0[3, i - 1] * r0[3, i - 1] * r0[3, i - 1];
                }
            }
            gcsr = Math.Sqrt(gcs);
            double step_l = (double)1 / 3;
            gcvr = Math.Pow(gcv, step_l);
            double xv, xs;
            k = twopi / w;
            xv = k * gcvr;
            xs = k * gcsr;

            OUT.WriteLine("volume-equiv. xv: " + xv + "   surface-equiv. xs: " + xs);

            double temp1, temp2;
            for (int i = 1; i <= nL; i++)
            {
                x[i - 1] = k * r0[3, i - 1];
                reff[i - 1] = new Complex(r0[4, i - 1], r0[5, i - 1]);
                temp1 = reff[i - 1].Real;
            }
            for (int j = 1; j <= np; j++)
                for (int i = 1; i <= nL; i++)
                {
                    aMie[i - 1, j - 1] = new Complex(0, 0);
                    bMie[i - 1, j - 1] = new Complex(0, 0);
                }
            int nmax0 = 1;
            for (int i = 1; i <= nL; i++)
            {
                if (i != 1 && x[i - 1] == x[i - 2] && reff[i - 1] == reff[i - 2])
                {
                    nmax[i - 1] = nmax[i - 2];
                    uvmax[i - 1] = uvmax[i - 2];
                }
                else
                {
                    OUT.WriteLine("sphere #" + i + "   individual size parameter: "+x[i-1].ToString());
                    abMiexud(x[i - 1], reff[i - 1], np, NXMAX, out nmax[i - 1], an, bn, NADD, rsr0, rsi0, rsx0, px0, eps);         
                    if (nmax[i - 1] > np)
                    {
                        OUT.WriteLine("Parameter np too small, must be > " + nmax[i - 1]);
                        OUT.WriteLine("Please change np in gmm01f.par, recompile, then try again");
                       // OUT.Close();
                        
                    }
                    else
                    {
                        uvmax[i - 1] = nmax[i - 1] * (nmax[i - 1] + 2);

                        OUT.WriteLine("Actual single-sphere expansion truncation: " + nmax[i - 1]);

                        for (int j = 1; j <= nmax[i - 1]; j++)
                        {
                            rsr[j - 1, i - 1] = rsr0[j - 1];
                            rsi[j - 1, i - 1] = rsi0[j - 1];
                            rsx[j - 1, i - 1] = rsx0[j - 1];
                            px[j - 1, i - 1] = px0[j - 1];
                            temp1 = an[j - 1].Real;
                            temp2 = bn[j - 1].Real;
                            if(j == 1 || j == nmax[i-1])
                                OUT.WriteLine(j+"   "+temp1+"   "+an[j-1].Imaginary+"   "+temp2+"   "+bn[j-1].Imaginary);
                        }
                    }
                }
                for (int j = 1; j <= nmax[i - 1]; j++)
                {
                    aMie[i - 1, j - 1] = an[j - 1];
                    bMie[i - 1, j - 1] = bn[j - 1];
                    rsr[j - 1, i - 1] = rsr0[j - 1];
                    rsi[j - 1, i - 1] = rsi0[j - 1];
                    rsx[j - 1, i - 1] = rsx0[j - 1];
                    px[j - 1, i - 1] = px0[j - 1];
                }
                if (nmax[i - 1] > nmax0)
                    nmax0 = nmax[i - 1];
            }
            //StreamWriter gmm01f_w = new StreamWriter("gmm01f.out");
            double cextx = 0, cexty = 0, cabsx = 0, cabsy = 0, cscax = 0;
            double cscay = 0, cprx = 0, cpry = 0, cbakx = 0, cbaky = 0;
            for (int i = 1; i <= nL; i++)
            {
                cextxi[i - 1] = 0;
                cextyi[i - 1] = 0;
                cabsxi[i - 1] = 0;
                cabsyi[i - 1] = 0;
                cscaxi[i - 1] = 0;
                cscayi[i - 1] = 0;
                cprxi[i - 1] = 0;
                cpryi[i - 1] = 0;
            }
            for (int i = 1; i <= nang2; i++)
            {
                i11[i - 1] = 0;
                i21[i - 1] = 0;
                i22[i - 1] = 0;
                i12[i - 1] = 0;
                for (int jc = 1; jc <= npng; jc++)
                    for (int j = 1; j <= 4; j++)
                        for (int m = 1; m <= 4; m++)
                            mue[j - 1, m - 1, jc - 1, i - 1] = 0;
            }
            double iram = 0;
            if (idpq == 1)
                gmm_out.WriteLine("input file: Ag-Si-2s-405nm.k");

            gmm_out.Close();
            
            OUT.WriteLine();
            OUT.WriteLine("original input sphere-positions: ");

            int ii = 1;
            OUT.WriteLine(ii+"   "+ r0[0, 0]+"   "+ r0[1, 0]+"   "+ r0[2, ii-1]);
            ii = nL;
            OUT.WriteLine(ii + "   " + r0[0, ii-1] + "   " + r0[1, ii-1] + "   " + r0[2, ii - 1]);

            int nphaic;
            if (idpq == 1)
            {
                nphaic = nphai + 1;
                phair[nphaic - 1] = 0;
            }
            else
                nphaic = nphai;
            //OUT.Close();
#endregion

            double nthetc;
            double alph;
            double ca, sa, beta, cb, sb, cz, sz;
            int n0;
            double temp, factor;
            int nlarge;
            double xd;
            double cext0, cext1;
            int irc, n1, itrc, indpol, imn;
            for (int ibeta = 1; ibeta <= nbeta; ibeta++)
            {
                for (int iphai = 1; iphai <= nphaic; iphai++)
                {
                    if (idpq == 1 && iphai < nphaic)
                        nthetc = 1;
                    else
                        nthetc = nthet;
                    for (int ithet = 1; ithet <= nthetc; ithet++)
                    {
                        if (idc < 0)
                        {
                            betar[ibeta - 1] = (betamx - betami) * ran1d(iseed);
                            phair[iphai - 1] = (phaimx - phaimi) * ran1d(iseed);
                            thetr[ithet - 1] = (thetmx - thetmi) * ran1d(iseed);
                            OUT.WriteLine();
                            OUT.WriteLine(betar[ibeta - 1] + "   " + phair[iphai - 1] + "   " + thetr[ithet - 1]);
                            OUT.WriteLine();
                        }
                        for (int irot = 1; irot <= 2; irot++)
                        {
                            if (idd == 1 || irot != 2)
                            {
                                #region
                                iram = iram + 1;
                                if (irot == 1)
                                    alph = 0;
                                else
                                    alph = pih;
                                ca = Math.Cos(alph);
                                sa = Math.Sin(alph);
                                beta = betar[ibeta];
                                cb = Math.Cos(beta);
                                sb = Math.Sin(beta);
                                for (int i = 1; i <=nL; i++)
                                {
                                    x0 = r00[0, i - 1];
                                    y0 = r00[1, i - 1];
                                    r0[0, i - 1] = cb * x0 - sb * y0;
                                    r0[1, i - 1] = sb * x0 + cb * y0;
                                }
                                double phai = phair[iphai - 1];
                                double thet = thetr[ithet - 1];
                                if (idpq == 1 && nthetc == 1)
                                    thet = 0;
                                cb = Math.Cos(phai);
                                sb = Math.Sin(phai);
                                cz = Math.Cos(thet);
                                sz = Math.Sin(thet);
                                if (iram == 1 || iram / 50 * 50 == iram)
                                    OUT.WriteLine("iram & nram: "+ iram+"   "+nram);
                                for (int i = 1; i <= nL; i++)
                                {
                                    x0 = r0[0, i - 1];
                                    y0 = r0[1, i - 1];
                                    z0 = r00[2, i - 1];
                                    r0[0, i - 1] = ca * cz * x0 - (ca * sz * sb + sa * cb) * y0 + (ca * sz * cb - sa * sb) * z0;
                                    r0[1, i - 1] = sa * cz * x0 - (sa * sz * sb - ca * cb) * y0 + (sa * sz * cb + ca * sb) * z0;
                                    r0[2, i - 1] = -sz * x0 - cz * sb * y0 + cz * cb * z0;
                                }
                                if (iram == 1 || iram / 50 * 50 == iram)
                                {
                                    ii = 1;
                                    OUT.WriteLine(ii + "   " + r0[0, 0] + "   " + r0[1, 0] + "   " + r0[2, ii - 1]);
                                    ii = nL;
                                    OUT.WriteLine(ii + "   " + r0[0, ii - 1] + "   " + r0[1, ii - 1] + "   " + r0[2, ii - 1]);
                                }  
                                n0 = nmax0 + 2;
                                fnr[0] = 0;
                                for (int n = 1; n <= 2 * n0; n++)
                                    fnr[n] = Math.Sqrt((double)n);
                                bcof[0] = 1;
                                for (int n = 0; n <= n0 - 1; n++)
                                    bcof[n + 1] = fnr[n + n + 2] * fnr[n + n + 1] * bcof[n] / fnr[n + 1] / fnr[n + 1];
                                cofsrd(nmax0);
                                cofd0(nmax0);
                                cofnv0(nmax0);
                                gau0(nmax0);
                                for (int i = 1; i <= nL; i++)
                                {
                                    for (int j = i + 1; j <= nL; j++)
                                    {
                                        int ij = (j - 1) * (j - 2) / 2 + j - i;
                                        x0 = r0[0, i - 1] - r0[0, j - 1];
                                        y0 = r0[1, i - 1] - r0[1, j - 1];
                                        z0 = r0[2, i - 1] - r0[2, j - 1];
                                        carsphd(x0, y0, z0, out d, out sphi, out cphi);
                                        temp = (r0[3, i - 1] + r0[3, j - 1]) / d;
                                        if (temp > fint)
                                        {
                                            ephi = new Complex(cphi, sphi);
                                            nlarge = Math.Max(nmax[i - 1], nmax[j - 1]);
                                            for (int m = 1; m <= nlarge; m++)
                                                ek[m - 1, ij - 1] = Complex.Pow(ephi, m);
                                            xd = k * d;
                                            int nbes = 2 * nlarge + 1;
                                            besseljd(nbes, xd);
                                            besselyd(nbes, xd);
                                            rotcoef(xt, nlarge);
                                            irc = 0;
                                            for (int n = 1; n <= nlarge; n++)
                                            {
                                                n1 = n * (n + 1);
                                                for (u = -n; u <= n; u++)
                                                    for (int m = -n; m <= n; m++)
                                                    {
                                                        imn = n1 + m;
                                                        irc = irc + 1;
                                                        drot[irc - 1, ij - 1] = dc[u +np, imn];
                                                    }
                                            }
                                            itrc = 0;
                                            int nsmall = Math.Min(nmax[i - 1], nmax[j - 1]);
                                            for (int m = -nsmall; m <= nsmall; m++)
                                            {
                                                n1 = Math.Max(1, Math.Abs(m));
                                                for (int n = n1; n <= nlarge; n++)
                                                    for (v = n1; v <= nlarge; v++)
                                                    {
                                                        itrc = itrc + 1;
                                                        /*atr0,btr0,atr1,btr1 (определяются неверно)*/
                                                        cofxuds0(nmax0, m, n, v, besj, besy,out atr0[itrc - 1, ij - 1],out btr0[itrc - 1, ij - 1],out atr1[itrc - 1, ij - 1],out btr1[itrc - 1, ij - 1]);
                                                    }
                                            }
                                        }
                                    }
                                }

                                OUT.WriteLine("atr0[530112, 1]      " + atr0[530111, 0]);
                                OUT.WriteLine("atr0[17668, 1]      " + atr0[17667, 0]);
                                OUT.WriteLine("atr0[300398,1]      " + atr0[300397, 0]);
                                OUT.WriteLine("atr0[ni0,1]      " + atr0[ni0 - 1, 0]);

                                OUT.WriteLine("btr0[132526, 1]      " + btr0[132525, 0]);
                                OUT.WriteLine("btr0[141362, 1]      " + btr0[141361, 0]);
                                OUT.WriteLine("btr0[265056,1]      " + btr0[265055, 0]);
                                OUT.WriteLine("btr0[ni0,1]      " + btr0[ni0 - 1, 0]);

                                OUT.WriteLine("atr1[1, 1]      " + atr1[0, 0]);
                                OUT.WriteLine("atr1[11, 1]      " + atr1[10, 0]);
                                OUT.WriteLine("atr1[121,1]      " + atr1[120, 0]);
                                OUT.WriteLine("atr1[ni0,1]      " + atr1[ni0 - 1, 0]);

                                OUT.WriteLine("btr1[1, 1]      " + btr1[0, 0]);
                                OUT.WriteLine("btr1[11, 1]      " + btr1[10, 0]);
                                OUT.WriteLine("btr1[121,1]      " + btr1[120, 0]);
                                OUT.WriteLine("btr1[ni0,1]      " + btr1[ni0 - 1, 0]);
                         
                                indpol = 0;
                                factor = factor1;

                                OUT.WriteLine("orien.# " + iram+"     Solving for x-pol. inci. state");
                                #endregion
                                #region
                                bool F = true;
                                while (F)
                                {
                                    #region
                                    for (imn = 1; imn <= nmp; imn++)
                                        for (int i = 1; i <= nL; i++)
                                        {
                                            p0[i - 1, imn - 1] = 0;
                                            q0[i - 1, imn - 1] = 0;
                                            ass[i - 1, imn - 1] = 0;
                                            bs[i - 1, imn - 1] = 0;
                                        }
                                    for (int i = 1; i <= nL; i++)
                                    {
                                        cz = Math.Cos(k * r0[2, i - 1]);
                                        sz = Math.Sin(k * r0[2, i - 1]);
                                        cmz = new Complex(cz, sz) * 0.5;
                                        for (int n = 1; n <= nmax[i - 1]; n++)
                                        {
                                            imn = n * n + n + 1;
                                            A = cmz * fnr[2 * n+1];
                                            p0[i - 1, imn - 1] = aMie[i - 1, n - 1] * A;
                                            q0[i - 1, imn - 1] = bMie[i - 1, n - 1] * A;
                                            p0[i - 1, imn - 3] = -p0[i - 1, imn - 1];
                                            q0[i - 1, imn - 3] = q0[i - 1, imn - 1];
                                            if (indpol > 1)
                                            {
                                                p0[i - 1, imn - 1] = p0[i - 1, imn - 1] * cin;
                                                q0[i - 1, imn - 1] = q0[i - 1, imn - 1] * cin;
                                                p0[i - 1, imn - 3] = p0[i - 1, imn - 1];
                                                q0[i - 1, imn - 3] = -q0[i - 1, imn - 1];
                                            }
                                            ass[i - 1, imn - 1] = p0[i - 1, imn - 1];
                                            bs[i - 1, imn - 1] = q0[i - 1, imn - 1];
                                            ass[i - 1, imn - 3] = p0[i - 1, imn - 3];
                                            bs[i - 1, imn - 3] = q0[i - 1, imn - 3];
                                        }
                                    }
                                    bool flag = true;
                                    if (idMie == 1 || nL == 1)
                                        flag = false;
                                    if (flag)
                                    {
                                        for (int i = 1; i <= nL; i++)
                                        {
                                            ind[i - 1] = 0;
                                            c0i[i - 1] = 0;
                                            for (int n = 1; n <= nmax[i - 1]; n++)
                                            {
                                                imn = n * n + n + 1;
                                                c0i[i - 1] = c0i[i - 1] + (p0[i - 1, imn - 1] * Complex.Conjugate(p0[i - 1, imn - 1])).Real;
                                                c0i[i - 1] = c0i[i - 1] + (q0[i - 1, imn - 1] * Complex.Conjugate(q0[i - 1, imn - 1])).Real;
                                                c0i[i - 1] = c0i[i - 1] + (p0[i - 1, imn - 3] * Complex.Conjugate(p0[i - 1, imn - 3])).Real;
                                                c0i[i - 1] = c0i[i - 1] + (q0[i - 1, imn - 3] * Complex.Conjugate(q0[i - 1, imn - 3])).Real;
                                            }
                                        }
                                        double niter = 1;
                                        bool flag2 = true;
                                        if (factor1 < 0.001 && factor2 < 0.001)
                                            flag2 = false;
                                        if (flag2)
                                        {
                                            if (iram == 1)
                                                OUT.WriteLine("Starting iteration solution process");
                                            for (int i = 1; i <= nL; i++)
                                                for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                {
                                                    as0[i - 1, imn - 1] = p0[i - 1, imn - 1];
                                                    bs0[i - 1, imn - 1] = q0[i - 1, imn - 1];
                                                }
                                            while (flag2)
                                            {
                                                trans(nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, as0, bs0, as1, bs1, ind);
                                                for (int i = 1; i <= nL; i++)
                                                {
                                                    if (ind[i - 1] <= 0)
                                                    {
                                                        c1i[i - 1] = 0;
                                                        /*as1,bs1 (определяются неверно)*/
                                                        for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                        {
                                                            int n = (int)Math.Sqrt((double)imn);
                                                            as0[i - 1, imn - 1] = p0[i - 1, imn - 1] - aMie[i - 1, n - 1] * as1[i - 1, imn - 1];
                                                            bs0[i - 1, imn - 1] = q0[i - 1, imn - 1] - bMie[i - 1, n - 1] * bs1[i - 1, imn - 1];
                                                            A = as0[i - 1, imn - 1] - ass[i - 1, imn - 1];
                                                            B = bs0[i - 1, imn - 1] - bs[i - 1, imn - 1];
                                                            c1i[i - 1] = c1i[i - 1] + (A * Complex.Conjugate(A)).Real;
                                                            c1i[i - 1] = c1i[i - 1] + (B * Complex.Conjugate(B)).Real;
                                                            as0[i - 1, imn - 1] = ass[i - 1, imn - 1] + factor * A;
                                                            bs0[i - 1, imn - 1] = bs[i - 1, imn - 1] + factor * B;
                                                            ass[i - 1, imn - 1] = as0[i - 1, imn - 1];
                                                            bs[i - 1, imn - 1] = bs0[i - 1, imn - 1];
                                                        }
                                                    }
                                                }
                                                cext0 = 0;
                                                cext1 = 0;
                                                for (int i = 1; i <= nL; i++)
                                                {
                                                    if (ind[i - 1] <= 0)
                                                    {
                                                        cext0 = cext0 + c0i[i - 1];
                                                        cext1 = cext1 + c1i[i - 1];
                                                        temp = c1i[i - 1] / c0i[i - 1];
                                                        if (temp < small)
                                                            ind[i - 1] = 1;
                                                    }         
                                                }
                                                temp = cext1 / cext0;
                                                if (temp < small)
                                                {
                                                    flag = false;
                                                    flag2 = false;
                                                }
                                                if (flag)
                                                {
                                                    if (iram == 1 || iram == nram)
                                                        OUT.WriteLine("iteration #  " + niter + "   " + temp);
                                                    if (niter > MXINT)
                                                    {
                                                        OUT.WriteLine("*** Maximum iterations exceeded ***");
                                                        OUT.WriteLine("*** Switched to Bi-CGSTAB method***");
                                                        for (int i = 1; i <= nL; i++)
                                                        {
                                                            ind[i - 1] = 0;
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                ass[i - 1, imn - 1] = p0[i - 1, imn - 1];
                                                                bs[i - 1, imn - 1] = q0[i - 1, imn - 1];
                                                            }
                                                        }
                                                        niter = 1;
                                                        flag2 = false;
                                                    }
                                                    niter = niter + 1;
                                                }
                                            }
                                        }
                                        if (flag)
                                        {
                                            if (iram == 1)
                                                OUT.WriteLine("Starting Bi-CGSTAB solution process");
                                            trans(nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, ass, bs, as1, bs1, ind);
                                            for (int i = 1; i <= nL; i++)
                                            {
                                                c1i[i - 1] = 0;
                                                for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                {
                                                    int n = (int)Math.Sqrt((double)imn);
                                                    as1[i - 1, imn - 1] = -aMie[i - 1, n - 1] * as1[i - 1, imn - 1];
                                                    bs1[i - 1, imn - 1] = -bMie[i - 1, n - 1] * bs1[i - 1, imn - 1];
                                                    c1i[i - 1] = c1i[i - 1] + (as1[i - 1, imn - 1] * Complex.Conjugate(as1[i - 1, imn - 1])).Real;
                                                    c1i[i - 1] = c1i[i - 1] + (bs1[i - 1, imn - 1] * Complex.Conjugate(bs1[i - 1, imn - 1])).Real;
                                                }
                                            }
                                            temp = 0;
                                            for (int i = 1; i <= nL; i++)
                                            {
                                                cext0 = c1i[i - 1] / c0i[i - 1];
                                                if (cext0 < small)
                                                    ind[i - 1] = 1;
                                                if (cext0 > temp)
                                                    temp = cext0;
                                            }
                                            if (temp < small)
                                                flag = false;
                                            if (flag)
                                            {


                                                A0 = new Complex(0, 0);
                                                for (int i = 1; i <= nL; i++)
                                                {
                                                    if (ind[i - 1] <= 0)
                                                    {
                                                        for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                        {
                                                            asp[i - 1, imn - 1] = as1[i - 1, imn - 1];
                                                            bsp[i - 1, imn - 1] = bs1[i - 1, imn - 1];
                                                            as0[i - 1, imn - 1] = as1[i - 1, imn - 1];
                                                            bs0[i - 1, imn - 1] = bs1[i - 1, imn - 1];
                                                            A0 = A0 + as1[i - 1, imn - 1] * as1[i - 1, imn - 1];
                                                            A0 = A0 + bs1[i - 1, imn - 1] * bs1[i - 1, imn - 1];
                                                        }
                                                    }
                                                }
                                                bool flag3 = true;
                                                while (flag3)
                                                {
                                                    trans(nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, asp, bsp, asv, bsv, ind);
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                int n = (int)Math.Sqrt((double)imn);
                                                                asv[i - 1, imn - 1] = aMie[i - 1, n - 1] * asv[i - 1, imn - 1] + asp[i - 1, imn - 1];
                                                                bsv[i - 1, imn - 1] = bMie[i - 1, n - 1] * bsv[i - 1, imn - 1] + bsp[i - 1, imn - 1];
                                                            }
                                                        }
                                                    }
                                                    A = new Complex(0, 0);
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                A = A + asv[i - 1, imn - 1] * as1[i - 1, imn - 1];
                                                                A = A + bsv[i - 1, imn - 1] * bs1[i - 1, imn - 1];
                                                            }
                                                        }
                                                    }
                                                    Aj = A0 / A;
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                asc[i - 1, imn - 1] = as0[i - 1, imn - 1] - Aj * asv[i - 1, imn - 1];
                                                                bsc[i - 1, imn - 1] = bs0[i - 1, imn - 1] - Aj * bsv[i - 1, imn - 1];
                                                            }
                                                        }
                                                    }
                                                    trans(nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, asc, bsc, ast, bst, ind);
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                int n = (int)Math.Sqrt((double)imn);
                                                                ast[i - 1, imn - 1] = aMie[i - 1, n - 1] * ast[i - 1, imn - 1] + asc[i - 1, imn - 1];
                                                                bst[i - 1, imn - 1] = bMie[i - 1, n - 1] * bst[i - 1, imn - 1] + bsc[i - 1, imn - 1];
                                                            }
                                                        }
                                                    }
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                A2 = A2 + ast[i - 1, imn - 1] * asc[i - 1, imn - 1];
                                                                A2 = A2 + bst[i - 1, imn - 1] * bsc[i - 1, imn - 1];
                                                                B2 = B2 + ast[i - 1, imn - 1] * ast[i - 1, imn - 1];
                                                                B2 = B2 + bst[i - 1, imn - 1] * bst[i - 1, imn - 1];
                                                            }
                                                        }
                                                    }
                                                    Bj = A2 / B2;
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            c1i[i - 1] = 0;
                                                            for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                            {
                                                                Aj2 = Aj * asp[i - 1, imn - 1] + Bj * asc[i - 1, imn - 1];
                                                                Bj2 = Aj * bsp[i - 1, imn - 1] + Bj * bsc[i - 1, imn - 1];
                                                                c1i[i - 1] = c1i[i - 1] + (Aj2 * Complex.Conjugate(Aj2)).Real;
                                                                c1i[i - 1] = c1i[i - 1] + (Bj2 * Complex.Conjugate(Bj2)).Real;
                                                                ass[i - 1, imn - 1] = ass[i - 1, imn - 1] + Aj2;
                                                                bs[i - 1, imn - 1] = bs[i - 1, imn - 1] + Bj2;
                                                            }
                                                        }
                                                    }
                                                    cext0 = 0;
                                                    cext1 = 0;
                                                    for (int i = 1; i <= nL; i++)
                                                    {
                                                        if (ind[i - 1] <= 0)
                                                        {
                                                            cext0 = cext0 + c0i[i - 1];
                                                            cext1 = cext1 + c1i[i - 1];
                                                        }
                                                    }
                                                    temp = cext1 / cext0;
                                                    if (temp < small)
                                                    {
                                                        flag = false;
                                                        flag3 = false;
                                                    }
                                                    if (niter > MXINT)
                                                    {
                                                        OUT.WriteLine("Caution:");
                                                        OUT.WriteLine("*** Maximum iterations exceeded ***");
                                                        flag = false;
                                                        flag3 = false;
                                                    }
                                                    if (flag)
                                                    {
                                                        if (iram == 1 || iram == nram)
                                                            OUT.WriteLine("iteration #  " + niter + "   " + temp);
                                                        for (int i = 1; i <= nL; i++)
                                                        {
                                                            if (ind[i - 1] <= 0)
                                                            {
                                                                for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                                {
                                                                    as0[i - 1, imn - 1] = asc[i - 1, imn - 1] - Bj * ast[i - 1, imn - 1];
                                                                    bs0[i - 1, imn - 1] = bsc[i - 1, imn - 1] - Bj * bst[i - 1, imn - 1];
                                                                }
                                                            }
                                                        }
                                                        A2 = new Complex(0, 0);
                                                        for (int i = 1; i <= nL; i++)
                                                        {
                                                            ast[i - 1, 0] = new Complex(0, 0);
                                                            if (ind[i - 1] <= 0)
                                                            {
                                                                for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                                {
                                                                    ast[i - 1, 0] = ast[i - 1, 0] + as0[i - 1, imn - 1] * as1[i - 1, imn - 1];
                                                                    ast[i - 1, 0] = ast[i - 1, 0] + bs0[i - 1, imn - 1] * bs1[i - 1, imn - 1];
                                                                }
                                                                A2 = A2 + ast[i - 1, 0];
                                                            }
                                                        }
                                                        B0 = A2 / A0;
                                                        B0 = B0 * Aj / Bj;
                                                        for (int i = 1; i <= nL; i++)
                                                        {
                                                            if (ind[i - 1] <= 0)
                                                            {
                                                                for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                                                {
                                                                    asp[i - 1, imn - 1] = as0[i - 1, imn - 1] + B0 * (asp[i - 1, imn - 1] - Bj * asv[i - 1, imn - 1]);
                                                                    bsp[i - 1, imn - 1] = bs0[i - 1, imn - 1] + B0 * (bsp[i - 1, imn - 1] - Bj * bsv[i - 1, imn - 1]);
                                                                }
                                                            }
                                                        }
                                                        A0 = new Complex(0, 0);
                                                        for (int i = 1; i <= nL; i++)
                                                        {
                                                            if (ind[i - 1] <= 0)
                                                            {
                                                                cext0 = c1i[i - 1] / c0i[i - 1];
                                                                if (cext0 < small)
                                                                    ind[i - 1] = 1;
                                                                else
                                                                    A0 = A0 + ast[i - 1, 0];
                                                            }
                                                        }
                                                        niter = niter + 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    #endregion

                                    /*490*/
                                    for (int i = 1; i <= nL; i++)
                                        ind[i - 1] = 0;
                                    if (indpol == 0)
                                    {
                                        OUT.WriteLine("Calculating near-field. "+indpol);
                                        field(nL, r0, k, nmax, ass, bs);
                                    }
                                    trans(nL, r0, nmax, uvmax, fint, atr1, btr1, ek, drot, ass, bs, as1, bs1, ind);
                                    for (int i = 1; i <= nL; i++)
                                    {
                                        for (imn = 1; imn <= uvmax[i - 1]; imn++)
                                        {
                                            at[imn - 1] = ass[i - 1, imn - 1] + as1[i - 1, imn - 1];
                                            bt[imn - 1] = bs[i - 1, imn - 1] + bs1[i - 1, imn - 1];
                                        }
                                        for (int n = 1; n <= nmax[i - 1]; n++)
                                        {
                                            n1 = n + 1;
                                            int n2 = 2 * n;
                                            double rn = 1 / ((double)(n * n1));
                                            double pp = fnr[n] * fnr[n + 2] / fnr[n2+1] / fnr[n2 + 3] / ((double)n1);
                                            double t = fnr[n - 1] * fnr[n+1] / fnr[n2 - 1] / fnr[n2+1] / ((double)n);
                                            double sc = 0;
                                            temp = 0;
                                            for (int m = -n; m <= n; m++)
                                            {
                                                int iL = n * (n + 1) + m;
                                                sc = sc + (Complex.Conjugate(ass[i - 1, iL - 1]) * at[iL - 1]).Real;
                                                sc = sc + (Complex.Conjugate(bs[i - 1, iL - 1]) * bt[iL - 1]).Real;
                                                double rm = (double)(m) * rn;
                                                A0 = rm * bt[iL - 1];
                                                B0 = rm * at[iL - 1];
                                                if (n != nmax[i - 1])
                                                {
                                                    u = (n + 1) * (n + 2) + m;
                                                    double fnp = fnr[n + m] * fnr[n - m] * pp;
                                                    A0 = A0 + fnp * at[u - 1];
                                                    B0 = B0 + fnp * bt[u - 1];
                                                }
                                                if (n != 1 && Math.Abs(m) <= (n - 1))
                                                {
                                                    u = (n - 1) * n + m;
                                                    double fn = fnr[n + m - 1] * fnr[n - m - 1] * t;
                                                    A0 = A0 + fn * at[u - 1];
                                                    B0 = B0 + fn * bt[u - 1];
                                                }
                                                temp = temp + (Complex.Conjugate(ass[i - 1, iL - 1]) * A0).Real;
                                                temp = temp + (Complex.Conjugate(bs[i - 1, iL - 1]) * B0).Real;
                                            }
                                            if (indpol < 1)
                                            {
                                                cscaxi[i - 1] = cscaxi[i - 1] + sc;
                                                cscax = cscax + sc;
                                                cprxi[i - 1] = cprxi[i - 1] + temp;
                                                cprx = cprx + temp;
                                            }
                                            else
                                            {
                                                cscayi[i - 1] = cscayi[i - 1] + sc;
                                                cscay = cscay + sc;
                                                cpryi[i - 1] = cpryi[i - 1] + temp;
                                                cpry = cpry + temp;
                                            }
                                        }
                                    }
                                    for (int j = 1; j <= nL; j++)
                                    {
                                        cz = Math.Cos(k * r0[2, j - 1]);
                                        sz = Math.Sin(k * r0[2, j - 1]);
                                        cmz = new Complex(cz, -sz);
                                        A = new Complex(0, 0);
                                        B = new Complex(0, 0);
                                        for (int n = 1; n <= nmax[j - 1]; n++)
                                        {
                                            double rn = fnr[2 * n];
                                            int m0 = n * n + n + 1;
                                            u0 = n * n + n - 1;
                                            A = A + rn * (ass[j - 1, m0 - 1] + bs[j - 1, m0 - 1]);
                                            B = B + rn * (ass[j - 1, u0 - 1] - bs[j - 1, u0 - 1]);
                                        }
                                        if (indpol < 1)
                                        {
                                            cextxi[j - 1] = cextxi[j - 1] + ((A - B) * cmz).Real;
                                            cextx = cextx + ((A - B) * cmz).Real;
                                        }
                                        else
                                        {
                                            cextyi[j - 1] = cextyi[j - 1] - ((A + B) * cmz).Imaginary;
                                            cexty = cexty - ((A + B) * cmz).Imaginary;
                                        }
                                    }
                                    for (int j = 1; j <= nL; j++)
                                    {
                                        for (int n = 1; n <= nmax[j - 1]; n++)
                                        {
                                            A = reff[j - 1] * (new Complex(rsr[n - 1, j - 1], -rsi[n - 1, j - 1]));
                                            temp1 = -A.Imaginary;
                                            A = px[n - 1, j - 1] * (reff[j - 1] * rsx[n - 1, j - 1] - (new Complex(rsr[n - 1, j - 1], rsi[n - 1, j - 1])));
                                            temp = Complex.Abs(A) * Complex.Abs(A);

                                            if (temp == 0)
                                                dn = 0;
                                            else
                                                dn = temp1 / temp;
                                            A = (new Complex(r0[4, j - 1], -r0[5, j - 1])) * (new Complex(rsr[n - 1, j - 1], -rsi[n - 1, j - 1]));
                                            temp1 = -A.Imaginary;
                                            A = px[n - 1, j - 1] * (rsx[n - 1, j - 1] - reff[j - 1] * (new Complex(rsr[n - 1, j - 1], rsi[n - 1, j - 1])));
                                            temp = Complex.Abs(A) * Complex.Abs(A);
                                            if (temp == 0)
                                                cn = 0;
                                            else
                                                cn = temp1 / temp;
                                            for (int m = -n; m < n; m++)
                                            {
                                                int i = n * n + n + m;
                                                temp1 = dn * Complex.Abs(ass[j - 1, i - 1]) * Complex.Abs(ass[j - 1, i - 1]) + cn * Complex.Abs(bs[j - 1, i - 1]) * Complex.Abs(bs[j - 1, i - 1]);
                                                if (indpol < 1)
                                                {
                                                    cabsxi[j - 1] = cabsxi[j - 1] + temp1;
                                                    cabsx = cabsx + temp1;
                                                }
                                                else
                                                {
                                                    cabsyi[j - 1] = cabsyi[j - 1] + temp1;
                                                    cabsy = cabsy + temp1;
                                                }
                                            }
                                        }
                                    }
                                    for (int i = 1; i <= nang; i++)
                                    {
                                        int iang = (int)(2 * nang) - i;
                                        dang[i-1] = sang * (double)(i - 1);
                                        dang[iang-1] = 180 - dang[i - 1];
                                        double theta = dang[i - 1] * pione / 180;
                                        xt = Math.Cos(theta);
                                        st = Math.Sin(theta);
                                        tipitaud(nmax0, xt);
                                        for (int jc = 1; jc <= npng; jc++)
                                        {
                                            double azphi = pang * pih * (double)(jc - 1) / 90;
                                            sphi = Math.Sin(azphi);
                                            cphi = Math.Cos(azphi);
                                            for (imn = 1; imn <= nmp; imn++)
                                            {
                                                at[imn - 1] = new Complex(0, 0);
                                                bt[imn - 1] = new Complex(0, 0);
                                                atj[imn - 1] = new Complex(0, 0);
                                                btj[imn - 1] = new Complex(0, 0);
                                            }
                                            for (int j = 1; j <= nL; j++)
                                            {
                                                sb = r0[0, j - 1] * cphi + r0[1, j - 1] * sphi;
                                                sb = sb * st;
                                                cb = r0[2, j - 1] * xt;
                                                cz = k * (sb + cb);
                                                sz = k * (sb - cb);
                                                A = (new Complex(Math.Cos(cz), -Math.Sin(cz)));
                                                B = (new Complex(Math.Cos(sz), -Math.Sin(sz)));
                                                bool flag4 = true;
                                                imn = 1;
                                                while (imn <= uvmax[j - 1] && flag4)
                                                {
                                                    int n = (int)Math.Sqrt((double)imn);
                                                    if (n > nmax[j - 1])
                                                        flag4 = false;
                                                    if (flag4)
                                                    {
                                                        bool flag5 = true;
                                                        if (idMie > 0)
                                                        {
                                                            int m = imn - n * n - n;
                                                            if (Math.Abs(m) == 1)
                                                                flag5 = false;
                                                        }
                                                        if(flag5)
                                                        {
                                                            at[imn - 1] = at[imn - 1] + A * ass[j - 1, imn - 1];
                                                            bt[imn - 1] = bt[imn - 1] + A * bs[j - 1, imn - 1];
                                                            atj[imn - 1] = atj[imn - 1] + B * ass[j - 1, imn - 1];
                                                            btj[imn - 1] = btj[imn - 1] + B * bs[j - 1, imn - 1];
                                                        }
                                                    }
                                                    imn++;
                                                }
                                            }
                                            if (indpol < 1)
                                            {
                                                s2x[jc - 1, i - 1] = new Complex(0, 0);
                                                s4x[jc - 1, i - 1] = new Complex(0, 0);
                                                s2x[jc - 1, iang - 1] = new Complex(0, 0);
                                                s4x[jc - 1, iang - 1] = new Complex(0, 0);
                                            }
                                            else
                                            {
                                                s3y[jc - 1, i - 1] = new Complex(0, 0);
                                                s1y[jc - 1, i - 1] = new Complex(0, 0);
                                                s3y[jc - 1, iang - 1] = new Complex(0, 0);
                                                s1y[jc - 1, iang - 1] = new Complex(0, 0);
                                            }
                                            A = new Complex(0, 0);
                                            B = new Complex(0, 0);
                                            Aj = new Complex(0, 0);
                                            Bj = new Complex(0, 0);
                                            for (int j = 1; j <= nmax0; j++)
                                            {
                                                imn = (j - 1) * (j + 2) / 2 + 1;
                                                u = j * j + j;
                                                A = A + at[u - 1] * tau[imn - 1];
                                                B = B + bt[u - 1] * tau[imn - 1];
                                                if (i != iang)
                                                {
                                                    double t = Math.Pow((-1), (j + 1));
                                                    Aj = Aj + atj[u - 1] * tau[imn - 1] * t;
                                                    Bj = Bj + btj[u - 1] * tau[imn - 1] * t;
                                                }
                                            }
                                            if (indpol < 1)
                                            {
                                                s2x[jc - 1, i - 1] = s2x[jc - 1, i - 1] + A * cphi;
                                                s4x[jc - 1, i - 1] = s4x[jc - 1, i - 1] + B * (new Complex(0, -1)) * cphi;
                                                if (i != iang)
                                                {
                                                    s2x[jc - 1, iang - 1] = s2x[jc - 1, iang - 1] + Aj * cphi;
                                                    s4x[jc - 1, iang - 1] = s4x[jc - 1, iang - 1] + Bj * (new Complex(0, -1)) * cphi;
                                                }
                                            }
                                            else
                                            {
                                                s3y[jc - 1, i - 1] = s3y[jc - 1, i - 1] - A * cphi;
                                                s1y[jc - 1, i - 1] = s1y[jc - 1, i - 1] + B * (new Complex(0, 1)) * cphi;
                                                if (i != iang)
                                                {
                                                    s3y[jc - 1, iang - 1] = s3y[jc - 1, iang - 1] - Aj * cphi;
                                                    s1y[jc - 1, iang - 1] = s1y[jc - 1, iang - 1] + Bj * (new Complex(0, 1)) * cphi;
                                                }
                                            }
                                            double rm = 1;
                                            for (int m = 1; m <= nmax0; m++)
                                            {
                                                A = new Complex(0, 0);
                                                B = new Complex(0, 0);
                                                A2 = new Complex(0, 0);
                                                B2 = new Complex(0, 0);
                                                Aj = new Complex(0, 0);
                                                Bj = new Complex(0, 0);
                                                Aj2 = new Complex(0, 0);
                                                Bj2 = new Complex(0, 0);
                                                rm = -rm;
                                                for (int j = m; j <= nmax0; j++)
                                                {
                                                    imn = (j - 1) * (j + 2) / 2 + m + 1;
                                                    u = j * j + j + m;
                                                    v = u - 2 * m;
                                                    A0 = at[u - 1] * tau[imn - 1] + bt[u - 1] * pi[imn - 1];
                                                    B0 = rm * (at[v - 1] * tau[imn - 1] - bt[v - 1] * pi[imn - 1]);
                                                    A = A + A0 + B0;
                                                    A2 = A2 + A0 - B0;
                                                    A0 = at[u - 1] * pi[imn - 1] + bt[u - 1] * tau[imn - 1];
                                                    B0 = rm * (at[v - 1] * pi[imn - 1] - bt[v - 1] * tau[imn - 1]);
                                                    B = B + A0 - B0;
                                                    B2 = B2 + A0 + B0;
                                                    if (i != iang)
                                                    {
                                                        double t = Math.Pow((-1), (j + m + 1));
                                                        double pp = -t;
                                                        A0 = atj[u - 1] * tau[imn - 1] * t + btj[u - 1] * pi[imn - 1] * pp;
                                                        B0 = rm * (atj[v - 1] * tau[imn - 1] * t - btj[v - 1] * pi[imn - 1] * pp);
                                                        Aj = Aj + A0 + B0;
                                                        Aj2 = Aj2 + A0 - B0;
                                                        A0 = atj[u - 1] * pi[imn - 1] * pp + btj[u - 1] * tau[imn - 1] * t;
                                                        B0 = rm * (atj[v - 1] * pi[imn - 1] * pp - btj[v - 1] * tau[imn - 1] * t);
                                                        Bj = Bj + A0 - B0;
                                                        Bj2 = Bj2 + A0 + B0;
                                                    }
                                                }
                                                temp = (double)(m - 1) * azphi;
                                                sb = Math.Sin(temp);
                                                cb = Math.Cos(temp);
                                                if (indpol < 1)
                                                {
                                                    s2x[jc - 1, i - 1] = s2x[jc - 1, i - 1] + A * cb + A2 * (new Complex(0, 1)) * sb;
                                                    s4x[jc - 1, i - 1] = s4x[jc - 1, i - 1] + B * (new Complex(0, -1)) * cb + B2 * sb;
                                                    if (i != iang)
                                                    {
                                                        s2x[jc - 1, iang - 1] = s2x[jc - 1, iang - 1] + Aj * cb;
                                                        s2x[jc - 1, iang - 1] = s2x[jc - 1, iang - 1] + Aj2 * (new Complex(0, 1)) * sb;
                                                        s4x[jc - 1, iang - 1] = s4x[jc - 1, iang - 1] + Bj * (new Complex(0, -1)) * cb;
                                                        s4x[jc - 1, iang - 1] = s4x[jc - 1, iang - 1] + Bj2 * sb;
                                                    }
                                                }
                                                else
                                                {
                                                    s3y[jc - 1, i - 1] = s3y[jc - 1, i - 1] - A * cb;
                                                    s3y[jc - 1, i - 1] = s3y[jc - 1, i - 1] + A2 * (new Complex(0, -1)) * sb;
                                                    s1y[jc - 1, i - 1] = s1y[jc - 1, i - 1] + B * (new Complex(0, 1)) * cb;
                                                    s1y[jc - 1, i - 1] = s1y[jc - 1, i - 1] - B2 * sb;
                                                    if (i != iang)
                                                    {
                                                        s3y[jc - 1, iang - 1] = s3y[jc - 1, iang - 1] - Aj * cb;
                                                        s3y[jc - 1, iang - 1] = s3y[jc - 1, iang - 1] + Aj2 * (new Complex(0, -1)) * sb;
                                                        s1y[jc - 1, iang - 1] = s1y[jc - 1, iang - 1] + Bj * (new Complex(0, 1)) * cb;
                                                        s1y[jc - 1, iang - 1] = s1y[jc - 1, iang - 1] - Bj * sb;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    indpol = indpol + 2;
                                    factor = factor2;
                                    if (indpol >= 3)
                                        F = false;
                                    else
                                        OUT.WriteLine("orien.# " + iram + "     Solving for y-pol. inci. state");
                                }
                                #endregion
                                for (int i = 1; i <= nang2; i++)
                                {
                                    i22[i - 1] = i22[i - 1] + Complex.Abs(s2x[0, i - 1]) * Complex.Abs(s2x[0, i - 1]);
                                    i21[i - 1] = i21[i - 1] + Complex.Abs(s4x[0, i - 1]) * Complex.Abs(s4x[0, i - 1]);
                                    i11[i - 1] = i11[i - 1] + Complex.Abs(s1y[0, i - 1]) * Complex.Abs(s1y[0, i - 1]);
                                    i12[i - 1] = i12[i - 1] + Complex.Abs(s3y[0, i - 1]) * Complex.Abs(s3y[0, i - 1]);
                                    for (int jc = 1; jc <= npng; jc++)
                                    {
                                        mueller(s1y[jc - 1, i - 1], s2x[jc - 1, i - 1], s3y[jc - 1, i - 1], s4x[jc - 1, i - 1], smue);                                 
                                        for (int j = 1; j <= 4; j++)
                                            for (int m = 1; m <= 4; m++)
                                                mue[j - 1, m - 1, jc - 1, i - 1] = mue[j - 1, m - 1, jc - 1, i - 1] + smue[j - 1, m - 1];
                                    }
                                }
                                cbakx = cbakx + Complex.Abs(s2x[0, (int)nang2 - 1]) * Complex.Abs(s2x[0, (int)nang2 - 1]);
                                cbaky = cbaky + Complex.Abs(s1y[0, (int)nang2 - 1]) * Complex.Abs(s1y[0, (int)nang2 - 1]);
                                cz = 4 / (gcs * k * k);
                                if (idpq == 1)
                                {
                                    temp1 = s2x[0, 0].Real * cz;
                                    temp2 = s1y[0, 0].Real * cz;
                                    OUT.WriteLine(iram.ToString() + " " + thet / pih * 90 + " " + phai / pih * 90 + " " + temp1 + " " + s2x[0, 0].Imaginary * cz + " " + temp2 + " " + s1y[0, 0].Imaginary * cz);
                                    if(idc>0)
                                        OUT.WriteLine(iram.ToString() + " " + dang[0] + " " + Math.Cos(thet) + " " + phai*90/pih);
                                    else
                                        OUT.WriteLine(iram.ToString() + " " + dang[0] + " " + thet*90/pih + " " + phai*90/pih);
                                }
                                if (idMie != 1)
                                {
                                    //if(iram<10)
                                    //{
                                    //    write(cnr1,'(i1)') iram
                                    //    tailn='00'//cnr1
                                    //}
                                    //else
                                    //{
                                    //    if(iram<100)
                                    //    {
                                    //        write(cnr2,'(i2)') iram
                                    //        tailn='0'//cnr2
                                    //    }
                                    //    else
                                    //    {
                                    //        write(cnr3,'(i3)') iram
                                    //        tailn=cnr3
                                    //    }
                                    //}
                                    if (nram == 1)
                                    {
                                        if (iram == 1)
                                        {
                                            OUT.WriteLine();
                                            OUT.WriteLine("Output file for scattering amplitude matrix: gmm01f.Aout");
                                            OUT.WriteLine();
                                        }
                                        StreamWriter gmm01fAout_w = new StreamWriter("gmm01f.Aout");
                                        gmm01fAout_w.WriteLine("gmm01f.Aout      (Scattering amplitude matrix)");
                                        gmm01fAout_w.WriteLine("The results are for the x-z plane of phi=0 only");
                                        gmm01fAout_w.WriteLine("wavelength: " + w + "      input filename: Ag-Si-2s-405nm.k");
                                        gmm01fAout_w.WriteLine("sphere#,x,y,z,radius,complex refractive index:");
                                        for (int i = 1; i <= nL; i++)
                                            gmm01fAout_w.WriteLine(i.ToString() + " " + r0[0, i-1].ToString() + " " + r0[1, i-1].ToString() + " " + r0[2, i-1].ToString() + " " + r0[3, i-1].ToString() + " " + r0[4, i-1].ToString() + " " + r0[5, i-1].ToString());
                                        gmm01fAout_w.WriteLine("scattering angle, s2x(complex), s3y(complex)");
                                        gmm01fAout_w.WriteLine("                  s4x(complex), s1y(complex)");
                                        for (int i = 1; i <= nang2; i++)
                                        {
                                            gmm01fAout_w.WriteLine(dang[i - 1].ToString() + " " + s2x[0, i - 1].Real.ToString() + " " + s2x[0, i - 1].Imaginary.ToString() + " " + s3y[0, i - 1].Real.ToString() + " " + s3y[0, i - 1].Imaginary.ToString());
                                            gmm01fAout_w.WriteLine(s4x[0, i - 1].Real.ToString() + " " + s4x[0, i - 1].Imaginary.ToString() + " " + s1y[0, i - 1].Real.ToString() + " " + s1y[0, i - 1].Imaginary.ToString());
                                        }
                                        gmm01fAout_w.Close();
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (iram != nram)
            {
                OUT.WriteLine("Note: iram not equal to nram!");
                OUT.WriteLine("iram: "+iram+"   nram: "+nram);
            }

            if (idpq == 1)
                gmm_out.Close();
            cz = nram;
            for (int i = 1; i <= nang2; i++)
            {
                i11[i - 1] = i11[i - 1] / cz;
                i21[i - 1] = i21[i - 1] / cz;
                i22[i - 1] = i22[i - 1] / cz;
                i12[i - 1] = i12[i - 1] / cz;
                inat[i - 1] = i11[i - 1] + i22[i - 1] + i12[i - 1] + i21[i - 1];
                pol[i - 1] = (i11[i - 1] - i22[i - 1]) / inat[i - 1];
                for (int jc = 1; jc <= npng; jc++)
                    for (int j = 1; j <= 4; j++)
                        for (int m = 1; m <= 4; m++)
                            mue[j - 1, m - 1, jc - 1, i - 1] = mue[j - 1, m - 1, jc - 1, i - 1] / cz;
            }
            cz = cz * k * k;
            cscax = 2 * twopi * cscax / cz;
            cscay = 2 * twopi * cscay / cz;
            double csca = 0.5 * (cscax + cscay);
            cextx = twopi * cextx / cz;
            cexty = twopi * cexty / cz;
            double cext = 0.5 * (cextx + cexty);
            cabsx = 2 * twopi * cabsx / cz;
            cabsy = 2 * twopi * cabsy / cz;
            double cabs = 0.5 * (cabsx + cabsy);
            cprx = 2 * twopi * cprx / cz;
            cpry = 2 * twopi * cpry / cz;
            double cpr = 0.5 * (cprx + cpry);
            double assym = (cprx + cpry) / (cscax + cscay);
            double assym0 = 0.5 * (cprx / cscax + cpry / cscay);
            cbakx = 2 * twopi * cbakx / cz;
            cbaky = 2 * twopi * cbaky / cz;
            double cbak = 0.5 * (cbakx + cbaky);

            OUT.WriteLine(assym + "    " + assym0);
            OUT.WriteLine("Cext,Cabs,Csca,Cbak,Cpr,<cos(theta)>");
            OUT.WriteLine(cext + "   " + cabs + "   " + csca+"   "+cbak+"   "+(cext - cpr)+"   "+assym);

            cscax = 0;
            cscay = 0;
            cextx = 0;
            cexty = 0;
            cabsx = 0;
            cabsy = 0;
            cprx = 0;
            cpry = 0;
            for (int i = 1; i <= nL; i++)
            {
                cscax = cscax + cscaxi[i - 1];
                cscay = cscay + cscayi[i - 1];
                cextx = cextx + cextxi[i - 1];
                cabsx = cabsx + cabsxi[i - 1];
                cexty = cexty + cextyi[i - 1];
                cabsy = cabsy + cabsyi[i - 1];
                cprx = cprx + cprxi[i - 1];
                cpry = cpry + cpryi[i - 1];
            }
            double assymx = cprx / cscax;
            double assymy = cpry / cscay;
            assym0 = 0.5 * (assymx + assymy);
            cscax = 2 * twopi * cscax / cz;
            cscay = 2 * twopi * cscay / cz;
            csca = 0.5 * (cscax + cscay);
            cextx = twopi * cextx / cz;
            cexty = twopi * cexty / cz;
            cext = 0.5 * (cextx + cexty);
            cabsx = 2 * twopi * cabsx / cz;
            cabsy = 2 * twopi * cabsy / cz;
            cabs = 0.5 * (cabsx + cabsy);
            cprx = 2 * twopi * cprx / cz;
            cpry = 2 * twopi * cpry / cz;
            cpr = 0.5 * (cprx + cpry);
            assym = cpr / csca;

            OUT.WriteLine(cext + "   " + cabs + "   " + csca + "   " + cbak + "   " + (cext - cpr) + "   " + assym);

            for (int i = 1; i <= nL; i++)
            {
                cabsxi[i - 1] = 4 * pione * cabsxi[i - 1] / cz;
                cabsyi[i - 1] = 4 * pione * cabsyi[i - 1] / cz;
                cextxi[i - 1] = 2 * pione * cextxi[i - 1] / cz;
                cextyi[i - 1] = 2 * pione * cextyi[i - 1] / cz;
                cscaxi[i - 1] = 4 * pione * cscaxi[i - 1] / cz;
                cscayi[i - 1] = 4 * pione * cscayi[i - 1] / cz;
                cprxi[i - 1] = 4 * pione * cprxi[i - 1] / cz;
                cpryi[i - 1] = 4 * pione * cpryi[i - 1] / cz;
                cscai[i - 1] = 0.5 * (cscaxi[i - 1] + cscayi[i - 1]);
                cexti[i - 1] = 0.5 * (cextxi[i - 1] + cextyi[i - 1]);
                cabsi[i - 1] = 0.5 * (cabsxi[i - 1] + cabsyi[i - 1]);
                cpri[i - 1] = 0.5 * (cprxi[i - 1] + cpryi[i - 1]);
                cpri[i - 1] = cscai[i - 1] + cabsi[i - 1] - cpri[i - 1];
                assymxi[i - 1] = cprxi[i - 1] / cscaxi[i - 1];
                assymyi[i - 1] = cpryi[i - 1] / cscayi[i - 1];
                assymi[i - 1] = 0.5 * (cprxi[i - 1] + cpryi[i - 1]) / csca;
                cprxi[i - 1] = cscaxi[i - 1] + cabsxi[i - 1] - cprxi[i - 1];
                cpryi[i - 1] = cscayi[i - 1] + cabsyi[i - 1] - cpryi[i - 1];

                OUT.WriteLine(i + "   " + cexti[i - 1] + "   " + cabsi[i - 1] + "   " + cscai[i - 1] + "   " + cpri[i - 1] + "   " + assymi[i - 1]);

            }

            OUT.WriteLine("efficiencies for radiation pressure");
            OUT.WriteLine(assym+"   "+assym0);
            OUT.WriteLine(assym + "   " + (cext - cpr) + "   " + assymx + "   " + (cextx - cprx) + "   " + assymy + "   " + (cexty - cpry));

            for (int i = 1; i <= nL; i++)
                OUT.WriteLine(i + "   " + assymi[i - 1] + "   " + cpri[i - 1] + "   " + assymxi[i - 1] + "   " + cprxi[i - 1] + "   " + assymyi[i - 1] + "   " + cpryi[i - 1]);

            betami = betami * 90 / pih;
            betamx = betamx * 90 / pih;
            thetmi = thetmi * 90 / pih;
            thetmx = thetmx * 90 / pih;
            phaimi = phaimi * 90 / pih;
            phaimx = phaimx * 90 / pih;
            StreamWriter crgmm01f_w = new StreamWriter("crgmm01f.out");
            crgmm01f_w.WriteLine("crgmm01f.out                 Total and individual-particle cross sections");
            crgmm01f_w.WriteLine("input sphere-aggregate filename: Ag-Si-2s-405nm.k");
            crgmm01f_w.WriteLine(nbeta.ToString() + " " + nthet.ToString() + " " + nphai.ToString());
            crgmm01f_w.WriteLine("Ranges of Euler angles: {0} {1} {2} {3} {4} {5}", betami.ToString(), betamx.ToString(), thetmi.ToString(), thetmx.ToString(), phaimi.ToString(), phaimx.ToString());
            crgmm01f_w.WriteLine("# of orientations averaged: {0}", nram.ToString());
            crgmm01f_w.WriteLine("Cext       Cabs      Csca     Cpr    <cos(theta)>");
            crgmm01f_w.WriteLine("total {0} {1} {2} {3} {4}", cext.ToString(), cabs.ToString(), csca.ToString(), (cext - cpr).ToString(), assym.ToString());
            for (int i = 1; i <= nL; i++)
                crgmm01f_w.WriteLine("{0} {1} {2} {3} {4} {5}", i.ToString(), cexti[i - 1].ToString(), cabsi[i - 1].ToString(), cscai[i - 1].ToString(), cpri[i - 1].ToString(), assymi[i - 1].ToString());
            crgmm01f_w.Close();
            cz = pione * gcvr * gcvr;
            assym = cpr / csca;
            assymx = cprx / cscax;
            assymy = cpry / cscay;
            double cabsxv = cabsx / cz;
            double cabsyv = cabsy / cz;
            double cextxv = cextx / cz;
            double cextyv = cexty / cz;
            double cscaxv = cscax / cz;
            double cscayv = cscay / cz;
            double cprxv = cprx / cz;
            cprxv = cextxv - cprxv;
            double cpryv = cpry / cz;
            cpryv = cextyv - cpryv;
            double cscav = 0.5 * (cscaxv + cscayv);
            double cextv = 0.5 * (cextxv + cextyv);
            double cabsv = 0.5 * (cabsxv + cabsyv);
            double cprv = 0.5 * (cprxv + cpryv);
            double cbakxv = cbakx / cz;
            double cbakyv = cbaky / cz;
            double cbakv = 0.5 * (cbakxv + cbakyv);
            temp = gcvr * gcvr / gcs;
            double cabsxs = cabsxv * temp;
            double cabsys = cabsyv * temp;
            double cextxs = cextxv * temp;
            double cextys = cextyv * temp;
            double cscaxs = cscaxv * temp;
            double cscays = cscayv * temp;
            double cprxs = cprxv * temp;
            double cprys = cpryv * temp;
            double cscas = cscav * temp;
            double cexts = cextv * temp;
            double cabss = cabsv * temp;
            double cprs = cprv * temp;
            double cbakxs = cbakxv * temp;
            double cbakys = cbakyv * temp;
            double cbaks = cbakv * temp;
            // 222      format(1x,a1,6e13.5)
            // 221      format(6x,a5,8x,a5,8x,a5,8x,a5,8x,a4,5x,a12)

            OUT.WriteLine("Qextv,Qabsv,Qscav,Qbakv,Qprv,<cos(theta)>");
            OUT.WriteLine("t" + "   " + cextv + "   " + cabsv + "   " + cscav + "   " + cbakv + "   " + cprv + "   " + assym);
            OUT.WriteLine("x" + "   " + cextxv + "   " + cabsxv + "   " + cscaxv + "   " + cbakxv + "   " + cprxv + "   " + assymx);
            OUT.WriteLine("y" + "   " + cextyv + "   " + cabsyv + "   " + cscayv + "   " + cbakyv + "   " + cpryv + "   " + assymy);
            OUT.WriteLine("Qexts,Qabss,Qscas,Qbaks,Qprs,<cos(theta)>");
            OUT.WriteLine("t" + "   " + cexts + "   " + cabss + "   " + cscas + "   " + cbaks + "   " + cprs + "   " + assym);
            OUT.WriteLine("x" + "   " + cextxs + "   " + cabsxs + "   " + cscaxs + "   " + cbakxs + "   " + cprxs + "   " + assymx);
            OUT.WriteLine("y" + "   " + cextys + "   " + cabsys + "   " + cscays + "   " + cbakys + "   " + cprys + "   " + assymy);

            temp = -(cabs + csca - cext) / cext;

            OUT.WriteLine("Accuracy of this numerical solution: "+temp);
            OUT.Close();

            StreamWriter gmm01f_w = new StreamWriter("gmm01f.out");
            gmm01f_w.WriteLine("gmm01f.out      --- input file: Ag-Si-2s-405nm.k xv: {0}   xs: {1}", xv, xs);
            gmm01f_w.WriteLine("Ranges of Euler angles: " + betami + "   " + betamx + "   " + thetmi + "   " + thetmx + "   " + phaimi + "   " + phaimx);
            gmm01f_w.WriteLine("nbeta,nthet,nphai: " + "   " + nbeta + "   " + nthet + "   " + nphai + "   # of orientations averaged: "+nram);
            gmm01f_w.WriteLine("Cext,Cabs,Csca,Cbak,Cpr,<cos(theta)>");
            gmm01f_w.WriteLine("t" + "   " + cext + "   " + cabs + "   " + csca + "   " + cbak + "   " + (cext - cpr) + "   " + assym);
            gmm01f_w.WriteLine("x" + "   " + cextx + "   " + cabsx + "   " + cscax + "   " + cbakx + "   " + (cextx - cprx) + "   " + assymx);
            gmm01f_w.WriteLine("y" + "   " + cexty + "   " + cabsy + "   " + cscay + "   " + cbaky + "   " + (cexty - cpry) + "   " + assymy);
            gmm01f_w.WriteLine("Qextv,Qabsv,Qscav,Qbakv,Qprv,<cos(theta)>");
            gmm01f_w.WriteLine("t" + "   " + cextv + "   " + cabsv + "   " + cscav + "   " + cbakv + "   " + cprv + "   " + assym);
            gmm01f_w.WriteLine("x" + "   " + cextxv + "   " + cabsxv + "   " + cscaxv + "   " + cbakxv + "   " + cprxv + "   " + assymx);
            gmm01f_w.WriteLine("y" + "   " + cextyv + "   " + cabsyv + "   " + cscayv + "   " + cbakyv + "   " + cpryv + "   " + assymy);
            gmm01f_w.WriteLine("Qexts,Qabss,Qscas,Qbaks,Qprs,<cos(theta)>");
            gmm01f_w.WriteLine("t" + "   " + cexts + "   " + cabss + "   " + cscas + "   " + cbaks + "   " + cprs + "   " + assym);
            gmm01f_w.WriteLine("x" + "   " + cextxs + "   " + cabsxs + "   " + cscaxs + "   " + cbakxs + "   " + cprxs + "   " + assymx);
            gmm01f_w.WriteLine("y" + "   " + cextys + "   " + cabsys + "   " + cscays + "   " + cbakys + "   " + cprys + "   " + assymy);
            gmm01f_w.WriteLine("s.a.,i11+i22,pol.,i11,i21,i12,i22");
            for (int i = 1; i <= nang2; i++)
                gmm01f_w.WriteLine(dang[i - 1] + "   " + inat[i - 1] + "   " + pol[i - 1] + "   " + i11[i - 1] + "   " + i21[i - 1] + "   " + i12[i - 1] + "   " + i22[i - 1]);
            gmm01f_w.Close();
            StreamWriter mueller_out = new StreamWriter("mueller.out");
            mueller_out.WriteLine("mueller.out     (Mueller matrix)");
            mueller_out.WriteLine("Input filename:    Ag-Si-2s-405nm.k");
            mueller_out.WriteLine("nbeta,nthet,nphai: " + nbeta + "   " + nthet + "   " + nphai);
            mueller_out.WriteLine("Ranges of Euler angles: " + betami + "   " + betamx + "   " + thetmi + "   " + thetmx + "   " + phaimi + "   " + phaimx);
            mueller_out.WriteLine("# of orientations averaged: " + nram);
            for (int jc = 1; jc <= npng; jc++)
            {
                double t = pang * (double)(jc - 1);
                mueller_out.WriteLine("phi (in degrees): " + t);
                for(int i=1;i<=nang2;i++)
                {
                    mueller_out.WriteLine(dang[i - 1] + "   " + mue[0, 0, jc - 1, i - 1] + "   " + mue[0, 1, jc - 1, i - 1] + "   " + mue[0, 2, jc - 1, i - 1] + "   " + mue[0, 3, jc - 1, i - 1]);
                    mueller_out.WriteLine(mue[1, 0, jc - 1, i - 1] + "   " + mue[1, 1, jc - 1, i - 1] + "   " + mue[1, 2, jc - 1, i - 1] + "   " + mue[1, 3, jc - 1, i - 1]);
                    mueller_out.WriteLine(mue[2, 0, jc - 1, i - 1] + "   " + mue[2, 1, jc - 1, i - 1] + "   " + mue[2, 2, jc - 1, i - 1] + "   " + mue[2, 3, jc - 1, i - 1]);
                    mueller_out.WriteLine(mue[3, 0, jc - 1, i - 1] + "   " + mue[3, 1, jc - 1, i - 1] + "   " + mue[3, 2, jc - 1, i - 1] + "   " + mue[3, 3, jc - 1, i - 1]);
                }
            }
            mueller_out.WriteLine("phi (in degrees): 360");
            for (int i = 1; i <= nang2; i++)
            {
                mueller_out.WriteLine(dang[i - 1] + "   " + mue[0, 0, 0, i - 1] + "   " + mue[0, 1, 0, i - 1] + "   " + mue[0, 2, 0, i - 1] + "   " + mue[0, 3, 0, i - 1]);
                mueller_out.WriteLine(mue[1, 0, 0, i - 1] + "   " + mue[1, 1, 0, i - 1] + "   " + mue[1, 2, 0, i - 1] + "   " + mue[1, 3, 0, i - 1]);
                mueller_out.WriteLine(mue[2, 0, 0, i - 1] + "   " + mue[2, 1, 0, i - 1] + "   " + mue[2, 2, 0, i - 1] + "   " + mue[2, 3, 0, i - 1]);
                mueller_out.WriteLine(mue[3, 0, 0, i - 1] + "   " + mue[3, 1, 0, i - 1] + "   " + mue[3, 2, 0, i - 1] + "   " + mue[3, 3, 0, i - 1]);
            }
            mueller_out.Close();
            OUT.Close();
        }

        private void button6_Click(object sender, EventArgs e)
        {
            button2.Enabled = button5.Enabled = button6.Enabled = button3.Enabled = false;
            button1.Enabled = true;
        }

        private void carsphd(double x, double y, double z,out double r,out double sphi,out double cphi)
        {
            r = Math.Sqrt(x * x + y * y + z * z);
            if (r == 0)
            {
                xt = 1;
                sphi = 0;
                cphi = 1;
                return;
            }
            xt = z / r;
            if (x==0 && y==0)
            {
                sphi = 0;
                cphi = 1;
                return;
            }
            sphi = Math.Sqrt(x * x + y * y);
            cphi = x / sphi;
            sphi = y / sphi;
        }
        private void cofd0(int nmax)
        {
            int i = 0;
            int inm;
            double c = 0, c0 = 0, c1 = 0;
            double sm = -0.5 * Math.Pow((-1), nmax);
            for (int m = -nmax; m <= nmax; m++)
            {
                int ns = Math.Max(1, Math.Abs(m));
                sm = -sm;
                for (int n = ns; n <= nmax; n++)
                {
                    inm = n * (n + 1) - m;
                    for (int v = ns; v <= nmax; v++)
                    {
                        i = i + 1;
                        int ivm = v * (v + 1) + m;
                        c = cofsr[inm - 1] + cofsr[ivm - 1];
                        c = sm * Math.Exp(c);
                        c0 = fnr[2 * n+1] * fnr[2 * v+1];
                        c1 = fnr[n] * fnr[v] * fnr[n+1] * fnr[v+1];
                        c0 = c0 / c1;
                        cof0[i - 1] = c * c0;
                    }
                }
            }
        }
        private void cofnv0(int nmax)
        {
            double c1 = 0;
            for (int n = 1; n <= nmax; n++)
            {
                for (int v = n; v <= nmax; v++)
                {
                    c1 = lnfacd(2 * n) + lnfacd(2 * v);
                    c1 = c1 - lnfacd(2 * n + 2 * v);
                    c1 = c1 + 2 * lnfacd(n + v);
                    c1 = c1 - lnfacd(n) - lnfacd(v);
                    cnv[n - 1, v - 1] = c1;
                }
            }
        }
        private void cofsrd(int nmax)
        {
            double c;
            int i = 0;
            for (int n = 1; n <= nmax; n++)
            {
                for (int m = -n; m <= n; m++)
                {
                    i = i + 1;
                    c = lnfacd(n - m) - lnfacd(n + m);
                    cofsr[i - 1] = 0.5 * c;
                }
            }
        }
        private double lnfacd(double z)
        {
            double[] c0 = new double[11] { 0.16427423239836267*Math.Pow(10,5), -0.48589401600331902*Math.Pow(10,5),
                                            0.55557391003815523*Math.Pow(10,5), -0.30964901015912058*Math.Pow(10,5),
                                            0.87287202992571788*Math.Pow(10,4), -0.11714474574532352*Math.Pow(10,4),
                                            0.63103078123601037*Math.Pow(10,2), -0.93060589791758878,
                                            0.13919002438227877*Math.Pow(10,-2),-0.45006835613027859*Math.Pow(10,-8),
                                            0.13069587914063262*Math.Pow(10,-9)};
            double a = 1;
            double cp = 2.5066282746310005;
            double b = z + 10.5;
            b = (z + 0.5) * Math.Log(b, Math.E) - b;
            for (int i = 1; i <= 11; i++)
            {
                z = z + 1;
                a = a + c0[i - 1] / z;
            }
            return b + Math.Log(cp * a, Math.E);
        }
        private void field(int nL, double[,] r0, double k, int[] nmax, Complex[,] ass, Complex[,] bs)
        {
            int[] ngrd = new int[3];
            double[] grdmin = new double[3];
            double[] grdstp = new double[3];
            double x = 0, y = 0, z = 0;
            Complex cplxi = new Complex(0, 1);
            Complex cplx0 = new Complex(0, 0);
            Complex eiphi;

            Complex[] etot = new Complex[3];
            Complex cimn;
            normlz(nmax[nmax.Length-1]);
            double r = 0;
            int imn = 0;
            int jj = 0;

            StreamReader grid_in = new StreamReader("grid_in.txt");
            string Grid = grid_in.ReadToEnd();
            string[] spl = Grid.Split('\n');
            string[] obj = spl[0].Split(';');
            ngrd[0] = Convert.ToInt32(obj[0]);
            ngrd[1] = Convert.ToInt32(obj[1]);
            ngrd[2] = Convert.ToInt32(obj[2]);
            obj = spl[1].Split(';');
            grdmin[0] = Convert.ToDouble(obj[0]);
            grdmin[1] = Convert.ToDouble(obj[1]);
            grdmin[2] = Convert.ToDouble(obj[2]);
            grid_in.Close();
            for (int i = 1; i <= 3; i++)
                grdstp[i - 1] = -grdmin[i - 1] * 2 / (ngrd[i - 1] - 1);
            StreamWriter grid_w = new StreamWriter("grid.out");
            grid_w.WriteLine("number of grid points (nx, ny, nz): " + ngrd[0].ToString() + ", " + ngrd[1].ToString() + ", " + ngrd[2].ToString());
            grid_w.WriteLine("grid corner (x0, y0, z0): " + grdmin[0].ToString() + ", " + grdmin[1].ToString() + ", " + grdmin[2].ToString());
            grid_w.Write("grid step (dx, dy, dz): " + grdstp[0].ToString() + ", " + grdstp[1].ToString() + ", " + grdstp[2].ToString());
            grid_w.Close();
            StreamWriter field_w = new StreamWriter("field.dat");
            z = grdmin[2];
            for (int iz = 1; iz <= ngrd[2]; iz++)
            {
                OUT.WriteLine(iz + "        /       " + ngrd[2]);
                y = grdmin[1];
                for (int iy = 1; iy <= ngrd[1]; iy++)
                {
                    x = grdmin[0];
                    for (int ix = 1; ix <= ngrd[0]; ix++)
                    {
                        etot[0] = new Complex(Math.Cos(k * z), Math.Sin(k * z));
                        etot[1] = cplx0;
                        etot[2] = cplx0;
                        for (int i = 1; i <= nL; i++)
                        {
                            carsphd(x - r0[0, i - 1], y - r0[1, i - 1], z - r0[2, i - 1],out r,out sphi,out cphi);
                            eiphi = new Complex(cphi, sphi);
                            if (r < r0[3, i - 1])
                            {
                                etot[0] = cplx0;
                                etot[1] = cplx0;
                                etot[2] = cplx0;
                            }
                            else
                            {
                                vswf(nmax[i - 1], k * r, sphi, cphi, Nmn3, Mmn3);
                                escati[0] = cplx0;
                                escati[1] = cplx0;
                                escati[2] = cplx0;
                                imn = 1;
                                jj = 0;
                                for (int n = 1; n <= nmax[i - 1]; n++)
                                {
                                    imn = imn + n;
                                    for (int m = -n; m <= -1; m++)
                                    {
                                        jj++;
                                        cimn = cplxi * ass[i - 1, jj - 1] * Emn[imn - 1] * Math.Pow((-1), m) * Complex.Pow(eiphi, (2 * m));
                                        escati[0] = escati[0] + cimn * Nmn3[0, imn - 1];
                                        escati[1] = escati[1] - cimn * Nmn3[1, imn - 1];
                                        escati[2] = escati[2] + cimn * Nmn3[2, imn - 1];
                                        cimn = cplxi * bs[i - 1, jj - 1] * Emn[imn - 1] * Math.Pow((-1), m) * Complex.Pow(eiphi, (2 * m));
                                        escati[1] = escati[1] + cimn * Mmn3[0, imn - 1];
                                        escati[2] = escati[2] - cimn * Mmn3[1, imn - 1];
                                        imn--;
                                    }
                                    for (int m = 0; m <= n; m++)
                                    {
                                        jj++;
                                        cimn = cplxi * Emn[imn - 1] * ass[i - 1, jj - 1];
                                        escati[0] = escati[0] + cimn * Nmn3[0, imn - 1];
                                        escati[1] = escati[1] - cimn * Nmn3[1, imn - 1];
                                        escati[2] = escati[2] + cimn * Nmn3[2, imn - 1];
                                        cimn = cplxi * Emn[imn - 1] * bs[i - 1, jj - 1];
                                        escati[1] = escati[1] + cimn * Mmn3[0, imn - 1];
                                        escati[2] = escati[2] - cimn * Mmn3[1, imn - 1];
                                        imn++;
                                    }
                                }
                                sphcrtv(sphi, cphi, escati, escatc);
                                etot[0] = etot[0] + escatc[0];
                                etot[1] = etot[1] + escatc[1];
                                etot[2] = etot[2] + escatc[2];
                            }
                        }
                        field_w.Write(x.ToString() + " " + y.ToString() + " " + z.ToString());
                        field_w.Write("  ( " + etot[0].Real.ToString() + "," + etot[0].Imaginary.ToString() + " )  ");
                        field_w.Write("  ( " + etot[1].Real.ToString() + "," + etot[1].Imaginary.ToString() + " )  ");
                        field_w.Write("  ( " + etot[2].Real.ToString() + "," + etot[2].Imaginary.ToString() + " )  ");
                        field_w.WriteLine("  " + vnorm(etot));
                        x = x + grdstp[0];
                    }
                    y = y + grdstp[1];
                }
                z = z + grdstp[2];
            }
            field_w.Close();
        } 
        private void vswf(int nmax, double kr, double sphi, double cphi, Complex[,] Nmn3, Complex[,] Mmn3)
        {
            pitaud(nmax, xt);
            legdre(nmax, xt);
            besseljd(nmax, kr);
            besselyd(nmax, kr);
            Complex eiphi = new Complex(cphi, sphi);
            Complex cplxi = new Complex(0, 1);
            Complex hankln, psinpr, eimphi;
            double pimn, taumn, pmn;
            int imn = 0;
            for (int n = 1; n <= nmax; n++)
            {
                hankln = new Complex(besj[n], besy[n]);
                psinpr = new Complex(kr * besj[n - 1] - n * besj[n], kr * besy[n - 1] - n * besy[n]);
                for (int m = 0; m <= n; m++)
                {
                    eimphi = Complex.Pow(eiphi, m);
                    imn = imn + 1;
                    pimn = pi[imn - 1];
                    taumn = tau[imn - 1];
                    if (m == 0)
                        pmn = p[n];
                    else
                        pmn = pimn * Math.Sqrt(1 - xt * xt) / m;
                    Mmn3[0, imn - 1] = -taumn * hankln * eimphi;
                    Mmn3[1, imn - 1] = cplxi * pimn * hankln * eimphi;
                    Nmn3[0, imn - 1] = n * (n + 1) * pmn * hankln * eimphi / kr;
                    Nmn3[1, imn - 1] = cplxi * pimn * psinpr * eimphi / kr;
                    Nmn3[2, imn - 1] = taumn * psinpr * eimphi / kr;
                }
            }
        } 
        private void sphcrtv(double sphi, double cphi, Complex[] vsph, Complex[] vcrt)
        {
            double st = 0;
            st = Math.Sqrt(1 - xt * xt);
            vcrt[0] = vsph[0] * cphi * st - vsph[1] * sphi + vsph[2] * cphi * xt;
            vcrt[1] = vsph[0] * sphi * st + vsph[1] * cphi + vsph[2] * sphi * xt;
            vcrt[2] = vsph[0] * xt - vsph[2] * st;
        }
        private void pitaud(int nmax, double x)
        {
            int n, i1, m, i2;
            double t, tx;
            int nt = (nmax + 1) * (nmax + 4) / 2;
            if (nt > nmp0 || Math.Abs(x) > 1)
            {
                MessageBox.Show("dimension or argument wrong in sub. tipitaud! argument: " + x.ToString());
                Close();
            }
            double sx = Math.Sqrt(1 - x * x);
            pi[0] = 0;
            pi[1] = Math.Sqrt(0.75);
            pi[2] = 0;
            double t125 = Math.Sqrt(1.25);
            pi[3] = t125 * x;
            pi[4] = t125 * sx;
            int imn = 5;
            for (int i = 3; i <= nmax + 1; i++)
            {
                n = i;
                imn = imn + 1;
                pi[imn - 1] = 0;
                for (int j = 2; j <= n; j++)
                {
                    m = j - 1;
                    imn = imn + 1;
                    i1 = (n - 2) * (n + 1) / 2 + m + 1;
                    if (m == (n - 1))
                        pi[imn - 1] = fnr[n - 1] * fnr[2 * n + 1] / fnr[n + 1] * x * pi[i1 - 1];
                    else
                    {
                        i2 = (n - 3) * n / 2 + m + 1;
                        t = fnr[n] * fnr[2 * n - 3];
                        t = fnr[n - 2] * fnr[n - m - 1] * fnr[n + m - 1] / t;
                        pi[imn - 1] = fnr[2 * n - 1] * x * pi[i1 - 1] - t * pi[i2 - 1];
                        t = fnr[n + 1] * fnr[n - m] * fnr[n + m];
                        t = fnr[n - 1] * fnr[2 * n + 1] / t;
                        pi[imn - 1] = t * pi[imn - 1];
                    }
                }
                imn = imn + 1;
                i1 = (n - 2) * (n + 1) / 2 + n;
                t = fnr[n - 1] * fnr[n + 1];
                t = Math.Sqrt(0.5) * fnr[n] * fnr[2 * n + 1] / t;
                pi[imn - 1] = t * sx * pi[i1 - 1];
            }
            tx = x * sx;
            tau[0] = -Math.Sqrt(1.5) * sx;
            tau[1] = pi[1] * x;
            tau[2] = -Math.Sqrt(7.5) * tx;
            tau[3] = t125 * (2 * x * x - 1);
            tau[4] = t125 * tx;
            imn = 5;
            for (int i = 3; i <= nmax; i++)
            {
                n = i;
                for (int j = 1; j <= n + 1; j++)
                {
                    m = j - 1;
                    imn = imn + 1;
                    if (m == 0)
                    {
                        i1 = (n - 2) * (n + 1) / 2 + 1;
                        i2 = (n - 3) * n / 2 + 1;
                        t = fnr[2 * n - 3];
                        t = fnr[n - 2] * fnr[n] / t;
                        tau[imn - 1] = fnr[2 * n - 1] * x * tau[i1 - 1] - t * tau[i2 - 1];
                        t = fnr[n - 1] * fnr[n + 1];
                        t = fnr[2 * n + 1] / t;
                        tau[imn - 1] = t * tau[imn - 1];
                    }
                    else
                    {
                        i1 = n * (n + 3) / 2 + m + 1;
                        t = fnr[n] * fnr[2 * n + 3];
                        t = fnr[n + 2] * fnr[2 * n + 1] * fnr[n - m + 1] * fnr[n + m + 1] / t;
                        tau[imn - 1] = t * pi[i1 - 1] - (double)(n + 1) * x * pi[imn - 1];
                        tau[imn - 1] = tau[imn - 1] / (double)(m);
                    }
                }
            }
        }
        private void legdre(int nmax, double xt)
        {
            int npp, nm, n2p;
            double en;
            p[0] = 1;
            p[1] = xt;
            npp = 2;
            nm = 0;
            n2p = 3;
            for (int n = 1; n <= nmax - 1; n++)
            {
                p[npp] = (double)(n2p * xt * p[n] - n * p[nm]) / (double)(npp);
                npp = npp + 1;
                nm = nm + 1;
                n2p = n2p + 2;
            }
            for (int n = 1; n <= nmax - 1; n++)
            {
                en = Math.Sqrt((double)(2 * n + 1) / (double)(n * (n + 1)));
                p[n] = p[n] * en;
            }
        }
        private void normlz(int nmax)
        {
            Complex cplxi = new Complex(0, 1);
            double en;
            int imn = 0;
            Complex cin = new Complex(1, 0);
            for (int n = 1; n <= nmax; n++)
            {
                cin = cin * cplxi;
                en = Math.Sqrt((2 * n + 1) / (double)(n * (n + 1)));
                for (int m = 0; m <= n; m++)
                {
                    imn = imn + 1;
                    Emn[imn - 1] = new Complex(1, 0) * cin;
                }
            }
        }
        private double vnorm(Complex[] vec)
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
        private void gau0(int nmax)
        {
            int na = 0;
            int qmax = 0;
            uvmax = nmax * (nmax + 2);
            int i = 0;
            for (int m = -nmax; m <= nmax; m++)
            {
                int ns = Math.Max(1, Math.Abs(m));
                for (int n = ns; n <= nmax; n++)
                {
                    for (int v = ns; v <= nmax; v++)
                    {
                        gxurcd0(-m, n, v,out qmax, na);
                        i = i + 1;
                        iga0[i - 1] = na;
                        na = na + qmax + 1;
                    }
                }
            }
        }
        private void gxurcd0(int m, int n, int v,out int qmax, int na)
        {
            
            if (Math.Abs(m) > n || Math.Abs(m) > v)
            {
                MessageBox.Show("warning: |m|>n or v in gxurcd0");
                qmax = -1;
                return;
            }
            double c1 = 0, c2 = 0, c3 = 0;
            int nq;
            qmax = Math.Min(n, v);
            nq = qmax + 1;
            if (n <= v)
                c1 = cnv[n-1, v-1];
            else
                c1 = cnv[v-1, n-1];
            c1 = c1 - lnfacd((double)(n - m)) - lnfacd((double)(v + m));
            ga0[na] = Math.Exp(c1);
            if (qmax < 1)
                return;
            int p = n + v;
            for (int i = 2; i <= nq; i++)
            {
                p = p - 2;
                if (m == 0)
                {
                    c1=fb(n,v,p+1);
                    c2=fb(n,v,p+2);
                    ga0[na + i - 1] = c2 * ga0[na + i - 2] / c1;
                }
                else
                {
                    c1=fb(n,v,p+1);
                    c2=(double)(4*m*m)+fb(n,v,p+2)+fb(n,v,p+3);
                    if (i == 2)
                        ga0[na + i - 1] = c2 * ga0[na + i - 2] / c1;
                    else
                    {
                        c3 = -fb(n, v, p + 4);
                        ga0[na + i - 1] = (c2 * ga0[na + i - 2] + c3 * ga0[na + i - 3]) / c1;
                    }
                }
            }
        }
        public double fb(int n, int v, int p)
        {
            return (double)(p - (n + v + 1)) * (double)(p + (n + v + 1)) * (double)(p - (n - v)) * (double)(p + (n - v)) / ((double)(2 * p + 1) * (double)(2 * p - 1));
        }
        private void besseljd(int NC, double X)
        {
            int NX, N = 0;
            double PN, CN, X2;
            for (int K = 1; K <= NC; K++)
                besj[K] = 0;
            if (Math.Abs(X) < Math.Pow(1, -12))
            {
                besj[0] = 1;
                besj[1] = 1 / 3 * X;
                return;
            }
            NX = (int)(1.1 * X) + 10;
            NX = NC + NX;
            PN = X / (double)(2 * NX + 3);
            for (int K = 1; K <= NX - 1; K++)
            {
                N = NX - K + 1;
                CN = (double)(N);
                PN = X / ((double)(2 * N + 1) - PN * X);
                if (N <= NC)
                    besj[N] = PN;
            }
            if (Math.Abs(X) - 0.1 < 0)
            {
                X2 = X * X;
                besj[0] = 1 - X2 / 72;
                besj[0] = 1 - X2 / 42 * besj[0];
                besj[0] = 1 - X2 / 20 * besj[0];
                besj[0] = 1 - X2 / 6 * besj[0];
                besj[1] = 1 / 45360 - X2 / 3991680;
                besj[1] = 1 / 840 - X2 * besj[1];
                besj[1] = 1 / 30 - X2 * besj[1];
                besj[1] = X * (1 / 3 - X2 * besj[1]);
                for (N = 2; N <= NC; N++)
                    besj[N] = besj[N] * besj[N - 1];
            }
            else
            {
                besj[0] = Math.Sin(X) / X;
                besj[1] = (Math.Sin(X) / X - Math.Cos(X)) / X;
                for (N = 2; N <= NC; N++)
                    besj[N] = besj[N] * besj[N - 1];
            }
        }
        private void besselyd(int n, double x)
        {
            double besyn, x2;
            if (x == 0)
            {
                MessageBox.Show("bad argument in sub. besselyd");
                Close();
            }
            if (Math.Abs(x) - 0.1 < 0)
            {
                x2 = x * x;
                besyn = 1 - x2 / 72;
                besyn = 1 - x2 / 42 * besyn;
                besyn = 1 - x2 / 20 * besyn;
                besyn = 1 - x2 / 6 * besyn;
                besy[0] = 1 - x2 / 56;
                besy[0] = 1 - x2 / 30 * besy[0];
                besy[0] = 1 - x2 / 12 * besy[0];
                besy[0] = 1 - 0.5 * x2 * besy[0];
                besy[0] = -besy[0] / x;
                besy[1] = besy[0] / x - besyn;
                for (int i = 2; i <= n; i++)
                    besy[i] = (double)(2 * i - 1) / x * besy[i - 1] - besy[i - 2];
            }
            else
            {
                besyn = Math.Sin(x) / x;
                besy[0] = -Math.Cos(x) / x;
                besy[1] = besy[0] / x - besyn;
                for (int i = 2; i <= n; i++)
                    besy[i] = (double)(2 * i - 1) / x * besy[i - 1] - besy[i - 2];
            }
        }
        private void cofxuds0(int nmax, int m, int n, int v, double[] sja, double[] sya, out Complex a, out Complex B, out Complex Aj, out Complex Bj)
        {
            //if (Math.Abs(m) > n || Math.Abs(m) > v)
            //{
            //    MessageBox.Show("|m|>n or v in subroutine cofxuds0.f");
            //    Close();
            //}
            int qmax = 0;
            double cp = 0;
            double sj, sy,c4;
            int id;
            a = 0;
            B = 0;
            Aj = 0;
            Bj = 0;
            id= gid0(nmax, m, n, v);
            double c = cof0[id - 1];
            int ig = iga0[id - 1];
            int nv2 = v * (v + 1) + n * (n + 1);
            Complex signz = Complex.Pow(new Complex(0, 1), (n + v));  //на fortran результат другой
            int pp = n + v + 2;
            qmax = Math.Min(n, v);
            for (int i = 1; i <= qmax + 1; i++)
            {
                pp = pp - 2;
                cp = (double)(nv2 - pp * (pp + 1)) * ga0[ig + i-1];
                sj = sja[pp];
                sy = sya[pp];
                a = a + (new Complex(sj, sy)) * signz * cp;
                Aj = Aj + sj * signz * cp;
                signz = -signz;
            }
            a = a * c;
            Aj = Aj * c;
            if (m != 0)
            {
                signz = Complex.Pow(new Complex(0, 1), (n + v + 1));
                pp = n + v;
                for (int i = 1; i <= qmax; i++)
                {
                    pp = pp - 2;
                    signz = -signz;
                    if (i == 1)
                    {
                        cp = (double)(2 * pp + 3) * fa(m, pp + 3);
                        cp = cp * ga0[ig] / (double)((pp + 3) * (pp + 2));
                        sj = sja[pp + 1];
                        sy = sya[pp + 1];
                        B = B + (new Complex(sj, sy)) * signz * cp;
                        Bj = Bj + sj * signz * cp;
                    }
                    else
                    {
                        if (i == qmax)
                        {
                            if (pp == 0)
                            {
                                c4 = fa(m, pp + 2);
                                cp = -(double)((pp + 1) * (pp + 2)) * fa(n, v, pp + 2) * ga0[ig + i - 1];
                                cp = cp + (double)((pp + 2) * (pp + 1)) * fa(n, v, pp + 1) * ga0[ig + i];
                                cp = cp * (double)(2 * pp + 3) / c4;
                                sj = sja[pp + 1];
                                sy = sya[pp + 1];
                                B = B + (new Complex(sj, sy)) * signz * cp;
                                Bj = Bj + sj * signz * cp;
                            }
                            else
                            {
                                nv2 = pp * (pp + 1);
                                cp = (double)(2 * pp + 3) * fa(m, pp + 1);
                                cp = -cp * ga0[ig + i] / (double)(nv2);
                                sj = sja[pp + 1];
                                sy = sya[pp + 1];
                                B = B + (new Complex(sj, sy)) * signz * cp;
                                Bj = Bj + sj * signz * cp;
                            }
                        }
                        else
                        {
                            c4 = fa(m, pp + 2);
                            cp = -(double)((pp + 1) * (pp + 2)) * fa(n, v, pp + 2) * ga0[ig + i - 1];
                            cp = cp + (double)((pp + 2) * (pp + 1)) * fa(n, v, pp + 1) * ga0[ig + i];
                            cp = cp * (double)(2 * pp + 3) / c4;
                            sj = sja[pp + 1];
                            sy = sya[pp + 1];
                            B = B + (new Complex(sj, sy)) * signz * cp;
                            Bj = Bj + sj * signz * cp;
                        }
                    }    
                }
                B = B * c;
                Bj = Bj * c;
            }
        }
        public double fa(int m,int p)
        {
            return (double)(-2 * m * p * (p - 1));
        }
        public double fa(int n, int v, int p)
        {
            return (double)(p * p - (n + v + 1) * (n + v + 1)) * (double)(p * p - (n - v) * (n - v)) / (double)(4 * p * p - 1);
        }
        private int gid0(int nmax, int m, int n, int iv)
        {
            int id = 0;
            int nt = nmax * (nmax + 1) * (2 * nmax + 1) / 3 + nmax * nmax;
            int ns = Math.Max(1, Math.Abs(m));
            int nc0 = nmax - Math.Abs(m);
            id = nc0 * (nc0 + 1) * (2 * nc0 + 1) / 6;
            if (m < 0)
                id = id + (n - ns) * (nc0 + 1) + iv - ns + 1;
            if (m == 0)
                id = id + (n - ns) * nmax + iv;
            if (m > 0)
            {
                id = id + (nmax - n) * (nc0 + 1) + nmax - iv;
                id = nt - id;
            }
            return id;
        }
        private void mueller(Complex s1, Complex s2, Complex s3, Complex s4, double[,] s)
        {
            double s1s, s2s, s3s, s4s;
            Complex s2s3c, s1s4c, s2s4c, s1s3c, s1s2c, s3s4c, s2s1c, s4s3c, s2cs4, s3cs1;
            s1s = Math.Pow(Complex.Abs(s1), 2);
            s2s = Math.Pow(Complex.Abs(s2), 2);
            s3s = Math.Pow(Complex.Abs(s3), 2);
            s4s = Math.Pow(Complex.Abs(s4), 2);
            s2s3c = s2 * Complex.Conjugate(s3);
            s1s4c = s1 * Complex.Conjugate(s4);
            s2s4c = s2 * Complex.Conjugate(s4);
            s1s3c = s1 * Complex.Conjugate(s3);
            s1s2c = s1 * Complex.Conjugate(s2);
            s3s4c = s3 * Complex.Conjugate(s4);
            s2s1c = s2 * Complex.Conjugate(s1);
            s4s3c = s4 * Complex.Conjugate(s3);
            s2cs4 = Complex.Conjugate(s2) * s4;
            s3cs1 = Complex.Conjugate(s3) * s1;
            s[0, 0] = 0.5 * (s1s + s2s + s3s + s4s);
            s[0, 1] = 0.5 * (s2s - s1s + s4s - s3s);
            s[0, 2] = (s2s3c + s1s4c).Real;
            s[0, 3] = (s2s3c - s1s4c).Imaginary;
            s[1, 0] = 0.5 * (s2s - s1s - s4s + s3s);
            s[1, 1] = 0.5 * (s2s + s1s - s4s - s3s);
            s[1, 2] = (s2s3c - s1s4c).Real;
            s[1, 3] = (s2s3c + s1s4c).Imaginary;
            s[2, 0] = (s2s4c + s1s3c).Real;
            s[2, 1] = (s2s4c - s1s3c).Real;
            s[2, 2] = (s1s2c + s3s4c).Real;
            s[2, 3] = (s2s1c + s4s3c).Imaginary;
            s[3, 0] = (s2cs4 + s3cs1).Imaginary;
            s[3, 1] = (s2cs4 - s3cs1).Imaginary;
            s[3, 2] = (s1s2c - s3s4c).Imaginary;
            s[3, 3] = (s1s2c - s3s4c).Real;
        }
        private void orientcd(double BETAMI, double BETAMX, double THETMI, double THETMX, double PHIMIN, double PHIMAX, double NBETA, double NTHETA, double NPHI, double[] BETA, double[] THETA, double[] PHI)
        {
            double delta;
            BETA[0] = BETAMI;
            if (NBETA > 1)
            {
                delta = (BETAMX - BETAMI) / (double)(NBETA - 1);
                for (int j = 2; j <= NBETA; j++)
                    BETA[j - 1] = BETA[0] + delta * (double)(j - 1);
            }
            if (NPHI == 1 && NTHETA > 1)
            {
                delta = (Math.Cos(THETMX) - Math.Cos(THETMI)) / (double)(NTHETA - 1);
                THETA[0] = THETMI;
            }
            else
            {
                delta = (Math.Cos(THETMX) - Math.Cos(THETMI)) / (double)(NTHETA);
                THETA[0] = Math.Acos(Math.Cos(THETMI) + 0.5 * delta);
            }
            if (NTHETA > 1)
            {
                for (int j = 2; j <= NTHETA; j++)
                    THETA[j - 1] = Math.Acos(Math.Cos(THETA[0]) + delta * (double)(j - 1));
            }
            PHI[0] = PHIMIN;
            if (NPHI > 1)
            {
                delta = (PHIMAX - PHIMIN) / (double)(NPHI - 1);
                for (int j = 2; j <= NPHI; j++)
                    PHI[j - 1] = PHI[0] + delta * (double)(j - 1);
            }
        }
        private void orientud(double BETAMI, double BETAMX, double THETMI, double THETMX, double PHIMIN, double PHIMAX, double NBETA, double NTHETA, double NPHI, double[] BETA, double[] THETA, double[] PHI)
        {
            double delta;
            BETA[0] = BETAMI;
            if (NBETA > 1)
            {
                delta = (BETAMX - BETAMI) / (double)(NBETA - 1);
                for (int j = 2; j <= NBETA; j++)
                    BETA[j - 1] = BETA[0] + delta * (double)(j - 1);
            }
            THETA[0] = THETMI;
            if (NTHETA > 1)
            {
                delta = (THETMX - THETMI) / (double)(NTHETA - 1);
                for (int j = 2; j <= NTHETA; j++)
                    THETA[j - 1] = THETA[0] + delta * (double)(j - 1);
            }
            PHI[0] = PHIMIN;
            if (NPHI > 1)
            {
                delta = (PHIMAX - PHIMIN) / (double)(NPHI - 1);
                for (int j = 2; j <= NPHI; j++)
                    PHI[j - 1] = PHI[0] + delta * (double)(j - 1);
            }
        }
        private double ran1d(int idum)
        {
            int NTAB = 32;
            int IM = 2147483647;
            int IR = 2836;
            int IA = 16807;
            int IQ = 127773;
            double AM = 1 / (double)IM, EPS = Math.Pow(1.2, -7), RNMX = 1 - EPS;
            int NDIV = 1 + (IM - 1) / NTAB;
            int iy = 0;
            int k, j;
            int[] iv = new int[NTAB];
            if (idum <= 0 || iy == 0)
            {
                idum = Math.Max(-idum, 1);
                for (j = NTAB + 8; j >= 1; j--)
                {
                    k = idum / IQ;
                    idum = IA * (idum - k * IQ) - IR * k;
                    if (idum < 0)
                        idum = idum + IM;
                    if (j <= NTAB)
                        iv[j - 1] = idum;
                }
                iy = iv[0];
            }
            k = idum / IQ;
            idum = IA * (idum - k * IQ) - IR * k;
            if (idum < 0)
                idum = idum + IM;
            j = 1 + iy / NDIV;
            iy = iv[j - 1];
            iv[j - 1] = idum;
            return Math.Min(AM * iy, RNMX);
        }
        private void rotcoef(double cbe, int nmax)
        {
            double[] dk0 = new double[4 * np + 1];
            double[] dk01 = new double[4 * np + 1];
            double sbe = Math.Sqrt((1 + cbe) * (1 - cbe));
            double cbe2 = 0.5 * (1 + cbe);
            double sbe2 = 0.5 * (1 - cbe);
            int inn = 1;
            dk0[0+2 * np] = 1;
            double sben = 1;
            dc[0 + np, 0] = 1;
            dk01[0+2 * np] = 0;
            for (int n = 1; n <= nmax; n++)
            {
                int nn1 = n * (n + 1);
                inn = -inn;
                sben = sben * sbe / 2;
                dk0[n + 2 * np] = (double)(inn) * sben * bcof[n];
                dk0[-n + 2 * np] = (double)(inn) * dk0[n + 2 * np];
                dk01[n + 2 * np] = 0;
                dk01[-n + 2 * np] = 0;
                dc[0 + np, nn1 + n] = dk0[n + 2 * np];
                dc[0 + np, nn1 - n] = dk0[-n + 2 * np];
                for (int k = -n + 1; k <= n - 1; k++)
                {
                    int kn = nn1 + k;
                    double dkt = dk01[k + 2 * np];
                    dk01[k + 2 * np] = dk0[k + 2 * np];
                    dk0[k + 2 * np] = (cbe * (double)(n + n - 1) * dk01[k + 2 * np] - fnr[n - k - 1] * fnr[n + k - 1] * dkt) / (fnr[n + k] * fnr[n - k]);
                    dc[0 + np, kn] = dk0[k + 2 * np];
                }
                int im = 1;
                for (int m = 1; m <= n; m++)
                {
                    im = -im;
                    double fmn = 1 / fnr[n - m + 1] / fnr[n + m];
                    int m1 = m - 1;
                    double dkm0 = 0;
                    for (int k = -n; k <= n; k++)
                    {
                        int kn = nn1 + k;
                        double dkm1 = dkm0;
                        double dkn1;
                        dkm0 = dc[m1 + np, kn];
                        if (k == n)
                            dkn1 = 0;
                        else
                            dkn1 = dc[m1 + np, kn + 1];
                        dc[m + np, kn] = (fnr[n + k] * fnr[n - k + 1] * cbe2 * dkm1 - fnr[n - k] * fnr[n + k + 1] * sbe2 * dkn1 - (double)(k) * sbe * dc[m1 + np, kn]) * fmn;
                        dc[-m + np, nn1 - k] = (double)(Math.Pow((-1), k) * im) * dc[m + np, kn];
                    }
                }
            }
        }
        private void rtr(Complex[,] Anpt, int nodrj, int nodri, Complex ekt, double drot)
        {
            Complex[] ek = new Complex[2 * np + 1];
            Complex[,] amt = new Complex[2, 2 * np + 1];
            Complex[,] ant = new Complex[2, 2 * np];
            Complex a, b;
            ek[0+np] = 1;
            int mmax = Math.Max(nodrj, nodri);
            for (int m = 1; m <= mmax; m++)
            {
                //ek[m+np]=ekt[m+np];
                ek[-m+np] = Complex.Conjugate(ek[m+np]);
            }
            int irc = 0;
            for (int n = 1; n <= nodrj; n++)
            {
                int n1 = (np + n) * (np + n + 1);
                for (int m = -n; m <= n; m++)
                {
                    amt[0, m+np] = 0;
                    amt[1, m+np] = 0;
                }
                for (int k = -n; k <= n; k++)
                {
                    int kn = n1 + np + k;
                    a = ek[k+np] * Anpt[0, kn];
                    b = ek[k+np] * Anpt[1, kn];
                    for (int m = -n; m <= n; m++)
                    {
                        irc = irc + 1;
                        //amt[0,m+np]=amt[0,m+np]+a*drot[irc-1];
                        //amt[1,m+np]=amt[1,m+np]+b*drot[irc-1];
                    }
                }
                for (int m = -n; m <= n; m++)
                {
                    int imn = n1 + m;
                    Anpt[0, imn - 1] = amt[0, m+np];
                    Anpt[1, imn - 1] = amt[1, m+np];
                }
            }
            mmax = Math.Min(nodrj, nodri);
            for (int m = -mmax; m <= mmax; m++)
            {
                int n1 = Math.Max(1, Math.Abs(m));
                for (int n = n1; n <= nodrj; n++)
                {
                    int imn = np + n * (n + 1) + m;
                    for (int ip = 1; ip <= 2; ip++)
                        ant[ip - 1, n - 1] = Anpt[ip - 1, imn - 1];
                }
                for (int n = n1; n <= nodri; n++)
                {
                    int imn = n * (n + 1) + m;
                    a = 0;
                    b = 0;
                    for (int l = n1; l <= nodrj; l++)
                    {
                        int ml = l * (l + 1) + m;
                        a = a + atr[0, n - 1, ml - 1] * ant[0, l - 1];
                        //1            +atr(2,n,ml)*ant(2,l)
                        b = b + atr[0, n - 1, ml - 1] * ant[1, l - 1];
                        //1            +atr(2,n,ml)*ant(1,l)
                    }
                    Anpt[0, imn - 1] = a;
                    Anpt[1, imn - 1] = b;
                }
            }
            int inn = 1;
            irc = 0;
            for (int n = 1; n <= nodri; n++)
            {
                inn = -inn;
                int n1 = n * (n + 1);
                for (int m = -n; m <= n; m++)
                {
                    amt[0, m+np] = 0;
                    amt[1, m + np] = 0;
                }
                int sik = -inn;
                for (int k = -n; k <= n; k++)
                {
                    sik = -sik;
                    int kn = n1 + np + k;
                    a = sik * Anpt[0, kn + np];
                    b = sik * Anpt[1, kn + np];
                    for (int m = -n; m <= n; m++)
                    {
                        irc = irc + 1;
                        //amt[0, m + np] = amt[0, m + np] + a * drot[irc - 1];
                        //amt[1, m + np] = amt[1, m + np] + b * drot[irc - 1];
                    }
                }
                sik = -inn;
                for (int m = -n; m <= n; m++)
                {
                    sik = -sik;
                    int imn = n1 + m;
                    Anpt[0, imn - 1] = amt[0, m + np] * ek[-m + np] * sik;
                    Anpt[1, imn - 1] = amt[1, m + np] * ek[-m + np] * sik;
                }
            }
        }
        private void tipitaud(int nmax, double x)
        {
            int i1, n, m, i2;
            double t;
            double nt = (nmax + 1) * (nmax + 4) / 2;
            if (nt > nmp0 || Math.Abs(x) > 1)
            {
                MessageBox.Show("dimension or argument wrong in sub. tipitaud! argument: " + x);
                Close();
            }
            double sx = Math.Sqrt(1 - x * x);
            pi[0] = 0;
            pi[1] = Math.Sqrt(0.75);
            pi[2] = 0;
            double t125 = Math.Sqrt(1.25);
            pi[3] = t125 * x;
            pi[4] = t125 * sx;
            int imn = 5;
            for (int i = 3; i <= nmax + 1; i++)
            {
                n = i;
                imn = imn + 1;
                pi[imn - 1] = 0;
                for (int j = 2; j <= n; j++)
                {
                    m = j - 1;
                    imn = imn + 1;
                    i1 = (n - 2) * (n + 1) / 2 + m + 1;
                    if (m == n - 1)
                        pi[imn - 1] = fnr[n - 1] * fnr[2 * n + 1] / fnr[n + 1] * x * pi[i1 - 1];
                    else
                    {
                        i2 = (n - 3) * n / 2 + m + 1;
                        t = fnr[n] * fnr[2 * n - 3];
                        t = fnr[n - 2] * fnr[n - m - 1] * fnr[n + m - 1] / t;
                        pi[imn - 1] = fnr[2 * n - 1] * x * pi[i1 - 1] - t * pi[i2 - 1];
                        t = fnr[n + 1] * fnr[n - m] * fnr[n + m];
                        t = fnr[n - 1] * fnr[2 * n + 1] / t;
                        pi[imn - 1] = t * pi[imn - 1];
                    }
                }
                imn = imn + 1;
                i1 = (n - 2) * (n + 1) / 2 + n;
                t = fnr[n - 1] * fnr[n + 1];
                t = Math.Sqrt(0.5) * fnr[n] * fnr[2 * n + 1] / t;
                pi[imn - 1] = t * sx * pi[i1 - 1];
            }
            double tx = x * sx;
            tau[0] = -Math.Sqrt(1.5) * sx;
            tau[1] = pi[1] * x;
            tau[2] = -Math.Sqrt(7.5) * tx;
            tau[3] = t125 * (2 * x * x - 1);
            tau[4] = t125 * tx;
            imn = 5;
            for (int i = 3; i <= nmax; i++)
            {
                n = i;
                for (int j = 1; j <= n + 1; j++)
                {
                    m = j - 1;
                    imn = imn + 1;
                    if (m == 0)
                    {
                        i1 = (n - 2) * (n + 1) / 2 + 1;
                        i2 = (n - 3) * n / 2 + 1;
                        t = fnr[2 * n - 3];
                        t = fnr[n - 2] * fnr[n] / t;
                        tau[imn - 1] = fnr[2 * n - 1] * x * tau[i1 - 1] - t * tau[i2 - 1];
                        t = fnr[n - 1] * fnr[n + 1];
                        t = fnr[2 * n + 1] / t;
                        tau[imn - 1] = t * tau[imn - 1];
                    }
                    else
                    {
                        i1 = n * (n + 3) / 2 + m + 1;
                        t = fnr[n] * fnr[2 * n + 3];
                        t = fnr[n + 2] * fnr[2 * n + 1] * fnr[n - m + 1] * fnr[n + m + 1] / t;
                        tau[imn - 1] = t * pi[i1 - 1] - (double)(n + 1) * x * pi[imn - 1];
                        tau[imn - 1] = tau[imn - 1] / (double)m;
                    }
                }
            }
        }
        private void trv(Complex[,] Anpt, int nodrj, int nodri)
        {
            int imn, n1, ml;
            Complex a, b;
            Complex[,] ant = new Complex[2, 2 * np];
            int mmax = Math.Min(nodrj, nodri);
            for (int m = -mmax; m <= mmax; m++)
            {
                n1 = Math.Max(1, Math.Abs(m));
                for (int n = n1; n <= nodrj; n++)
                {
                    imn = n * (n + 1) + m;
                    for (int ip = 1; ip <= 2; ip++)
                    {
                        ant[ip - 1, n - 1] = Anpt[ip - 1, imn - 1];
                    }
                }
                for (int n = n1; n <= nodri; n++)
                {
                    imn = n * (n + 1) + m;
                    a = 0;
                    b = 0;
                    for (int l = n1; l <= nodrj; l++)
                    {
                        ml = l * (l + 1) + m;
                        a = a + atr[0, n - 1, ml - 1] * ant[0, l - 1] + atr[1, n - 1, ml - 1] * ant[1, l - 1];
                        b = b + atr[0, n - 1, ml - 1] * ant[1, l - 1] + atr[1, n - 1, ml - 1] * ant[0, l - 1];
                    }
                    Anpt[0, imn - 1] = a;
                    Anpt[1, imn - 1] = b;
                }
            }
        }
        private void abMiexud(double X, Complex REFREL, int NP, int NMAX, out int NM, Complex[] AN, Complex[] BN, int NADD, double[] RSR, double[] RSI, double[] RSX, double[] px, double EPS)
        {
            double[] AR = new double[np], AI = new double[np], BR = new double[np], BI = new double[np];
            double XM = 0, YM = 0, SNX = 0, CNM1X = 0, CNX = 0, CNM2X = 0;
            double ALN, BEN, PZD, SNM1X;
            int M = 0;
            int N;
            double CN;
            if (EPS > 1 || EPS < 0)
                EPS = Math.Pow(1, -20);
            if (NADD != 0)
                EPS = 0;
            double CTC = EPS;
            XM = REFREL.Real;
            YM = REFREL.Imaginary;
            double XMX = X * XM;
            double YMX = X * YM;
            double RP2 = XMX * XMX + YMX * YMX;
            int NSTOP = (int)(X + 4 * Math.Pow(X, 0.3333));
            NSTOP = NSTOP + 2;
            NM = NSTOP + NADD;
            double XN = Math.Sqrt(XM * XM + YM * YM) * X;
            int NX = (int)(1.1 * XN + 10);
            if ((NX - NM) < 10)
                NX = NM + 10;
            OUT.WriteLine("Wiscombe criterion: " + NSTOP);
            OUT.WriteLine("NADD: " + NADD);
            OUT.WriteLine("NX: " + NX);
            if (NX > NMAX)
            {
                OUT.WriteLine("Parameter NXMAX too small");
                OUT.WriteLine("NXMAX must be greater than " + NX);
                OUT.WriteLine("Please correct NXMAX in main code, recompile, then try again");
                Close();
            }
            else
            {
                if (NM > NP)
                {
                    OUT.WriteLine("Parameter np too small");
                    OUT.WriteLine("np must be greater than " + NM);
                    OUT.WriteLine("Please correct np in gmm01f.par, recompile the code, then try again");
                    //NM = 0;
                    Close();
                }
                else
                {
                    double PNX = X / (double)(2 * NX + 3);
                    double PNR = XMX / (double)(2 * NX + 3);
                    double PNI = YMX / (double)(2 * NX + 3);
                    int k = 1;
                    bool flag1 = true;
                    while (k <= NX && flag1)
                    {
                        N = NX - k + 1;
                        CN = (double)(N);
                        ALN = (2 * CN + 1) * XMX / RP2 - PNR;
                        BEN = (2 * CN + 1) * YMX / RP2 + PNI;
                        RSR[N - 1] = -CN * XMX / RP2 + ALN;
                        RSI[N - 1] = CN * YMX / RP2 - BEN;
                        PZD = ALN * ALN + BEN * BEN;
                        PNR = ALN / PZD;
                        PNI = BEN / PZD;
                        RSX[N - 1] = (CN + 1) / X - PNX;
                        if (N != 1)
                        {
                            PNX = X / (2 * CN + 1 - PNX * X);
                            px[N - 1] = PNX;
                        }
                        else
                            flag1 = false;
                        k++;
                    }
                    SNM1X = Math.Sin(X);
                    CNM1X = Math.Cos(X);
                    if ((X - 0.1) < 0)
                        SNX = Math.Pow(X,2) / 3 - Math.Pow(X,4) / 30 + Math.Pow(X, 6) / 840 - Math.Pow(X, 8) / 45360;
                    else
                        SNX = SNM1X / X - CNM1X;
                    CNX = CNM1X / X + SNM1X;
                    N = 1;
                    flag1 = true;
                    while (N <= NX && flag1)
                    {
                        px[N - 1] = SNX;
                        double C = (double)(N);
                        double DCNX = CNM1X - C * CNX / X;
                        double DSNX = RSX[N - 1] * SNX;
                        double ANNR = RSR[N - 1] * SNX - XM * DSNX;
                        double ANNI = RSI[N - 1] * SNX - YM * DSNX;
                        double TA1 = RSR[N - 1] * SNX - RSI[N - 1] * CNX;
                        double TA2 = RSI[N - 1] * SNX + RSR[N - 1] * CNX;
                        double ANDR = TA1 - XM * DSNX + YM * DCNX;
                        double ANDI = TA2 - XM * DCNX - YM * DSNX;
                        double AND = ANDR * ANDR + ANDI * ANDI;
                        double BNNR = (XM * RSR[N - 1] - YM * RSI[N - 1]) * SNX - DSNX;
                        double BNNI = (XM * RSI[N - 1] + YM * RSR[N - 1]) * SNX;
                        double TB1 = RSR[N - 1] * SNX - RSI[N - 1] * CNX;
                        double TB2 = RSR[N - 1] * CNX + RSI[N - 1] * SNX;
                        double BNDR = XM * TB1 - YM * TB2 - DSNX;
                        double BNDI = XM * TB2 + YM * TB1 - DCNX;
                        double BND = BNDR * BNDR + BNDI * BNDI;
                        
                        AR[N - 1] = (ANNR * ANDR + ANNI * ANDI) / AND;
                        AI[N - 1] = (ANNI * ANDR - ANNR * ANDI) / AND;
                        BR[N - 1] = (BNNR * BNDR + BNNI * BNDI) / BND;
                        BI[N - 1] = (BNNI * BNDR - BNNR * BNDI) / BND;
                        double TI = AR[N - 1] * AR[N - 1] + AI[N - 1] * AI[N - 1] + BR[N - 1] * BR[N - 1] + BI[N - 1] * BI[N - 1];
                        TI = TI / (AR[0] * AR[0] + AI[0] * AI[0] + BR[0] * BR[0] + BI[0] * BI[0]);
                        if (TI - CTC < 0)
                        {
                            OUT.WriteLine("*** NOTE THAT THE FIELD-EXPANSION TRANCATION");
                            OUT.WriteLine("*** IS DETERMINED BY eps GIVEN IN THE INPUT");
                            OUT.WriteLine("*** FILE gmm01f.in");
                            OUT.WriteLine("*** IN CASE YOU NEED A HIGHER ORDER, eps MUST");
                            OUT.WriteLine(" *** BE SMALLER THAN THE CURRENT VALUE  " + EPS);
                            flag1 = false;
                        }
                        else
                        {
                            if (NM - N <= 0)
                                flag1 = false;
                            else
                            {
                                if (N - NX >= 0)
                                    flag1 = false;
                                else
                                {
                                    M = N + 1;
                                    SNX = px[M - 1] * SNX;
                                    CNM2X = CNM1X;
                                    CNM1X = CNX;
                                    CNX = (2 * C + 1) * CNM1X / X - CNM2X;
                                }
                            }
                        }
                        N++;
                    }
                    NM = N - 1;
                    for (int i = 1; i <= NM; i++)
                    {
                        AN[i - 1] = new Complex(AR[i - 1], -AI[i - 1]);
                        BN[i - 1] = new Complex(BR[i - 1], -BI[i - 1]);
                    }
                }
            }
        }
        private void trans(int nL, double[,] r0, int[] nmax, int[] uvmax, double fint, Complex[,] atr0, Complex[,] btr0, Complex[,] ek, double[,] drot, Complex[,] ass, Complex[,] bs, Complex[,] as1, Complex[,] bs1, int[] ind)
        {
            int ij = 0;
            Complex[,] at1 = new Complex[2, nmp];
            for (int i = 1; i <= nL; i++)
            {
                if (ind[i - 1] <= 0)
                {
                    for (int imn = 1; imn <= uvmax[i - 1]; imn++)
                    {
                        as1[i - 1, imn - 1] = new Complex(0, 0);
                        bs1[i - 1, imn - 1] = new Complex(0, 0);
                    }
                    for (int j = 1; j <= nL; j++)
                    {
                        if (j != i)
                        {
                            double x0 = r0[0, i - 1] - r0[0, j - 1];
                            double y0 = r0[1, i - 1] - r0[1, j - 1];
                            double z0 = r0[2, i - 1] - r0[2, j - 1];
                            d = Math.Sqrt(x0 * x0 + y0 * y0 + z0 * z0);
                            double temp = (r0[3, i - 1] + r0[3, j - 1]) / d;
                            if (temp > fint)
                            {
                                if (i < j)
                                    ij = (j - 1) * (j - 2) / 2 + j - i;
                                else
                                    ij = (i - 1) * (i - 2) / 2 + i - j;
                                int nlarge = Math.Max(nmax[i - 1], nmax[j - 1]);
                                int itrc = 0;
                                int nsmall = Math.Min(nmax[i - 1], nmax[j - 1]);

                                for (int m = -nsmall; m <= nsmall; m++)
                                {
                                    int n1 = Math.Max(1, Math.Abs(m));
                                    for (int n = n1; n <= nlarge; n++)
                                    {
                                        for (int v = n1; v <= nlarge; v++)
                                        {
                                            itrc = itrc + 1;
                                            int iuv = v * (v + 1) + m;
                                            atr[0, n - 1, iuv - 1] = atr0[itrc - 1, ij - 1];
                                            atr[1, n - 1, iuv - 1] = btr0[itrc - 1, ij - 1];
                                            if (x0 == 0 && y0 == 0)
                                            {
                                                if (z0 < 0 || j < i)
                                                {
                                                    double sic = (double)(Math.Pow((-1), (n + v)));
                                                    atr[0, n - 1, iuv - 1] = sic * atr[0, n - 1, iuv - 1];
                                                    atr[1, n - 1, iuv - 1] = -sic * atr[1, n - 1, iuv - 1];
                                                }
                                            }
                                        }
                                    }
                                }
                                for (int iuv = 1; iuv <= uvmax[j - 1]; iuv++)
                                {
                                    at1[0, iuv - 1] = ass[j - 1, iuv - 1];
                                    at1[1, iuv - 1] = bs[j - 1, iuv - 1];
                                }
                                if (x0 == 0 && y0 == 0)
                                    trv(at1, nmax[j - 1], nmax[i - 1]);
                                else
                                    rtr(at1, nmax[j - 1], nmax[i - 1], ek[0, ij - 1], drot[0, ij - 1]);
                                for (int imn = 1; imn <= uvmax[j - 1]; imn++)
                                {
                                    as1[i - 1, imn - 1] = as1[i - 1, imn - 1] + at1[0, imn - 1];
                                    bs1[i - 1, imn - 1] = bs1[i - 1, imn - 1] + at1[1, imn - 1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
