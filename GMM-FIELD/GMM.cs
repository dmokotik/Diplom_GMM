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
        int lP = 0;//количество частиц
        double[,] cnv, dc;
        int nmp, ni0, nmp0, np;
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
            int u=0, v, u0;
            int[] uvmax = new int[nLp], ind = new int[nLp];
            nmax = new int[nLp];
            pih = Math.Acos(0);
            double k, gcs = 0, gcv = 0, idpq = 0, pione = 2 * pih;
            double[,] r0 = new double[6, nLp], r00 = new double[3, nLp];
            double[] x = new double[nLp], dang = new double[nangmax], c0i = new double[nLp], c1i = new double[nLp];
            double[] rsr0 = new double[NXMAX], rsi0 = new double[NXMAX], rsx0 = new double[NXMAX], px0 = new double[NXMAX];
            double[] w1 = new double[np], w2 = new double[np], w3 = new double[np], w4 = new double[np];
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
            Complex A, B, cmz, A0=0, B0=0, Aj, Bj, Aj2, Bj2, A2=new Complex(0,0), B2=new Complex(0,0), ephi, ci, cin;
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
                new orientcd(betami, betamx, thetmi, thetmx, phaimi, phaimx, nbeta, nthet, nphai, betar, thetr, phair);
            else
                new orientud(betami, betamx, thetmi, thetmx, phaimi, phaimx, nbeta, nthet, nphai, betar, thetr, phair);
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
            r0[2,0]=-0.04;
            r0[3,0]=0.04;
            r0[4, 0] = 1.15402;
            r0[5, 0] = 1.88224;
            r0[0,1]=0;
            r0[1,1]=0;
            r0[2, 1] = 0.36;
            r0[3, 1] = -0.16;
            r0[4, 1] = 1.54928;
            r0[5, 1] = 0.21614;
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
                    new abMiexud(OUT, np, x[i - 1], reff[i - 1], np, NXMAX, out nmax[i - 1], an, bn, NADD, rsr0, rsi0, rsx0, px0,w1,w2,w3,w4, eps);
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
                OUT.WriteLine("input file: Ag-Si-2s-405nm.k");
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
            ran1 ran = new ran1();
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
                            betar[ibeta - 1] = (betamx - betami) * ran.ran1d(iseed);
                            phair[iphai - 1] = (phaimx - phaimi) * ran.ran1d(iseed);
                            thetr[ithet - 1] = (thetmx - thetmi) * ran.ran1d(iseed);
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
                                new cofsrd(cofsr, nmax0);
                                new cofd0(fnr, cof0, cofsr, nmax0);
                                new cofnv0(cnv, nmax0);
                                new gau0(cnv,ga0,iga0, nmax0);
                                for (int i = 1; i <= nL; i++)
                                {
                                    for (int j = i + 1; j <= nL; j++)
                                    {
                                        int ij = (j - 1) * (j - 2) / 2 + j - i;
                                        x0 = r0[0, i - 1] - r0[0, j - 1];
                                        y0 = r0[1, i - 1] - r0[1, j - 1];
                                        z0 = r0[2, i - 1] - r0[2, j - 1];
                                        new carsphd(out xt,x0, y0, z0, out d, out sphi, out cphi);
                                        
                                        temp = (r0[3, i - 1] + r0[3, j - 1]) / d;
                                        if (temp > fint)
                                        {
                                            ephi = new Complex(cphi, sphi);
                                            nlarge = Math.Max(nmax[i - 1], nmax[j - 1]);
                                            for (int m = 1; m <= nlarge; m++)
                                                ek[m - 1, ij - 1] = Complex.Pow(ephi, m);
                                            xd = k * d;
                                            int nbes = 2 * nlarge + 1;
                                            new besseljd(besj, nbes, xd);
                                            new besselyd(besy, nbes, xd);
                                            new rotcoef(np, dc,bcof, fnr, xt, nlarge);
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
                                                        new cofxuds0(cof0, iga0, ga0,nmax0, m, n, v, besj, besy,out atr0[itrc - 1, ij - 1],out btr0[itrc - 1, ij - 1],out atr1[itrc - 1, ij - 1],out btr1[itrc - 1, ij - 1]);
                                                    }
                                            }
                                        }
                                    }
                                }
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
                                                new trans(np, atr, d, nmp, nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, as0, bs0, as1, bs1, ind);
                                                for (int i = 1; i <= nL; i++)
                                                {
                                                    if (ind[i - 1] <= 0)
                                                    {
                                                        c1i[i - 1] = 0;
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
                                            new trans(np, atr, d, nmp, nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, ass, bs, as1, bs1, ind);
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
                                                    new trans(np, atr, d, nmp, nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, asp, bsp, asv, bsv, ind);
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
                                                    new trans(np, atr, d, nmp, nL, r0, nmax, uvmax, fint, atr0, btr0, ek, drot, asc, bsc, ast, bst, ind);
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
                                    for (int i = 1; i <= nL; i++)
                                        ind[i - 1] = 0;
                                    if (indpol == 0)
                                    {
                                        OUT.WriteLine("Calculating near-field. "+indpol);
                                        new field(Nmn3, Mmn3,escati, escatc, nmp0, fnr, besj, besy, pi, tau, p, OUT, Emn, xt, sphi, cphi, nL, r0, k, nmax, ass, bs);
                                    }
                                    new trans(np, atr, d, nmp, nL, r0, nmax, uvmax, fint, atr1, btr1, ek, drot, ass, bs, as1, bs1, ind);
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
                                                    double fnp = fnr[n + m+1] * fnr[n - m+1] * pp;
                                                    A0 = A0 + fnp * at[u - 1];
                                                    B0 = B0 + fnp * bt[u - 1];
                                                }
                                                if (n != 1 && Math.Abs(m) <= (n - 1))
                                                {
                                                    u = (n - 1) * n + m;
                                                    double fn = fnr[n + m] * fnr[n - m] * t;
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
                                            double rn = fnr[2 * n+1];
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
                                            for (int m = -n; m <= n; m++)
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
                                        new tipitaud(pi, tau, fnr, nmax0, xt, nmp0);
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
                                                Aj2 = 0;
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
                                        new mueller(s1y[jc - 1, i - 1], s2x[jc - 1, i - 1], s3y[jc - 1, i - 1], s4x[jc - 1, i - 1], smue);                                 
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
                                    OUT.WriteLine(iram + " " + thet / pih * 90 + " " + phai / pih * 90 + " " + temp1 + " " + s2x[0, 0].Imaginary * cz + " " + temp2 + " " + s1y[0, 0].Imaginary * cz);
                                    if(idc>0)
                                        OUT.WriteLine(iram + " " + dang[0] + " " + Math.Cos(thet) + " " + phai*90/pih);
                                    else
                                        OUT.WriteLine(iram + " " + dang[0] + " " + thet*90/pih + " " + phai*90/pih);
                                }
                                if (idMie != 1)
                                {
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
                                            gmm01fAout_w.WriteLine("{0,5:d} {1,10:f4} {2,10:f4} {3,10:f4} {4,10:f4} {5,10:f4} {6,10:f4}", i, r0[0, i - 1], r0[1, i - 1], r0[2, i - 1], r0[3, i - 1], r0[4, i - 1], r0[5, i - 1]);
                                        gmm01fAout_w.WriteLine("scattering angle, s2x(complex), s3y(complex)");
                                        gmm01fAout_w.WriteLine("                  s4x(complex), s1y(complex)");
                                        for (int i = 1; i <= nang2; i++)
                                        {
                                            gmm01fAout_w.WriteLine("{0,8:f2} {1,14:f6} {2,14:f6} {3,14:f6} {4,14:f6}", dang[i - 1], s2x[0, i - 1].Real, s2x[0, i - 1].Imaginary, s3y[0, i - 1].Real, s3y[0, i - 1].Imaginary);
                                            gmm01fAout_w.WriteLine("         {0,14:f6} {1,14:f6} {2,14:f6} {3,14:f6}", s4x[0, i - 1].Real, s4x[0, i - 1].Imaginary, s1y[0, i - 1].Real, s1y[0, i - 1].Imaginary);
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
                OUT.Close();
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
            OUT.WriteLine("     {0,14:f6} {1,14:f6}",assym, assym0);
            OUT.WriteLine("      Cext         Cabs         Csca         Cbak         Cpr        <cos(theta)>");
            OUT.WriteLine("  {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cext, cabs, csca, cbak, (cext - cpr), assym);
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
            
            OUT.WriteLine("  {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cext, cabs, csca, cbak, (cext - cpr), assym);
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
                OUT.WriteLine("{0,5:d} {1,14:f6} {2,14:f6} {3,14:f6} {4,14:f6} {5,14:f6}",i, cexti[i - 1], cabsi[i - 1], cscai[i - 1], cpri[i - 1], assymi[i - 1]);
            }
            OUT.WriteLine("efficiencies for radiation pressure");
            OUT.WriteLine("     {0,12:f6} {1,12:f6}",assym,assym0);
            OUT.WriteLine("     {0,12:f6} {1,12:f6} {2,12:f6} {3,12:f6} {4,12:f6} {5,12:f6}", assym, (cext - cpr), assymx, (cextx - cprx), assymy, (cexty - cpry));
            for (int i = 1; i <= nL; i++)
                OUT.WriteLine("{0,5:d} {1,12:f6} {2,12:f6} {3,12:f6} {4,12:f6} {5,12:f6} {6,12:f6}", i, assymi[i - 1], cpri[i - 1], assymxi[i - 1], cprxi[i - 1], assymyi[i - 1], cpryi[i - 1]);
            betami = betami * 90 / pih;
            betamx = betamx * 90 / pih;
            thetmi = thetmi * 90 / pih;
            thetmx = thetmx * 90 / pih;
            phaimi = phaimi * 90 / pih;
            phaimx = phaimx * 90 / pih;
            StreamWriter crgmm01f_w = new StreamWriter("crgmm01f.out");
            crgmm01f_w.WriteLine("crgmm01f.out                 (Total and individual-particle cross sections)");
            crgmm01f_w.WriteLine("input sphere-aggregate filename: Ag-Si-2s-405nm.k");
            crgmm01f_w.WriteLine("nbeta,nthet,nphai:    "+nbeta + " " + nthet + " " + nphai);
            crgmm01f_w.WriteLine("Ranges of Euler angles: {0} {1} {2} {3} {4} {5}", betami, betamx, thetmi, thetmx, phaimi, phaimx);
            crgmm01f_w.WriteLine("# of orientations averaged: {0}", nram.ToString());
            crgmm01f_w.WriteLine("            Cext           Cabs           Csca           Cpr        <cos(theta)>");
            crgmm01f_w.WriteLine("total {0,14:f9} {1,14:f9} {2,14:f9} {3,14:f9} {4,14:f9}", cext, cabs, csca, (cext - cpr), assym);
            for (int i = 1; i <= nL; i++)
                crgmm01f_w.WriteLine("{0,5:d} {1,14:f9} {2,14:f9} {3,14:f9} {4,14:f9} {5,14:f9}", i, cexti[i - 1], cabsi[i - 1], cscai[i - 1], cpri[i - 1], assymi[i - 1]);
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
            OUT.WriteLine("      Qextv        Qabsv        Qscav        Qbakv        Qprv     <cos(theta)>");
            OUT.WriteLine("t {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextv, cabsv, cscav, cbakv, cprv, assym);
            OUT.WriteLine("x {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextxv, cabsxv, cscaxv, cbakxv, cprxv, assymx);
            OUT.WriteLine("y {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextyv, cabsyv, cscayv, cbakyv, cpryv, assymy);
            OUT.WriteLine("      Qexts        Qabss        Qscas        Qbaks        Qprs     <cos(theta)>");
            OUT.WriteLine("t {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cexts, cabss, cscas, cbaks, cprs, assym);
            OUT.WriteLine("x {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextxs, cabsxs, cscaxs, cbakxs, cprxs, assymx);
            OUT.WriteLine("y {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextys, cabsys, cscays, cbakys, cprys, assymy);
            temp = -(cabs + csca - cext) / cext;

            OUT.WriteLine("Accuracy of this numerical solution: "+temp);
            OUT.Close();

            StreamWriter gmm01f_w = new StreamWriter("gmm01f.out");
            gmm01f_w.WriteLine("gmm01f.out        --- input file: Ag-Si-2s-405nm.k     xv: {0}      xs: {1}", xv, xs);
            gmm01f_w.WriteLine("Ranges of Euler angles: " + betami + "   " + betamx + "   " + thetmi + "   " + thetmx + "   " + phaimi + "   " + phaimx);
            gmm01f_w.WriteLine("nbeta,nthet,nphai:    " + nbeta + "    " + nthet + "    " + nphai + "      # of orientations averaged: "+nram);
            gmm01f_w.WriteLine("       Cext         Cabs         Csca         Cbak         Cpr     <cos(theta)>");
            gmm01f_w.WriteLine(" t {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cext, cabs, csca, cbak, (cext - cpr), assym);
            gmm01f_w.WriteLine(" x {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextx, cabsx, cscax, cbakx, (cextx - cprx), assymx);
            gmm01f_w.WriteLine(" y {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cexty, cabsy, cscay, cbaky, (cexty - cpry), assymy);
            gmm01f_w.WriteLine("      Qextv        Qabsv        Qscav        Qbakv        Qprv     <cos(theta)>");
            gmm01f_w.WriteLine(" t {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextv, cabsv, cscav, cbakv, cprv, assym);
            gmm01f_w.WriteLine(" x {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextxv, cabsxv, cscaxv, cbakxv, cprxv, assymx);
            gmm01f_w.WriteLine(" y {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextyv, cabsyv, cscayv, cbakyv, cpryv, assymy);
            gmm01f_w.WriteLine("      Qexts        Qabss        Qscas        Qbaks        Qprs     <cos(theta)>");
            gmm01f_w.WriteLine(" t {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cexts, cabss, cscas, cbaks, cprs, assym);
            gmm01f_w.WriteLine(" x {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextxs, cabsxs, cscaxs, cbakxs, cprxs, assymx);
            gmm01f_w.WriteLine(" y {0,12:f5} {1,12:f5} {2,12:f5} {3,12:f5} {4,12:f5} {5,12:f5}", cextys, cabsys, cscays, cbakys, cprys, assymy);
            gmm01f_w.WriteLine("    s.a.    i11+i22     pol.       i11          i21          i12          i22");
            for (int i = 1; i <= nang2; i++)
                gmm01f_w.WriteLine("{0,7:f1} {1,12:f5} {2,8:f4} {3,12:f5} {4,12:f5} {5,12:f5} {6,12:f5}", dang[i - 1], inat[i - 1], pol[i - 1], i11[i - 1], i21[i - 1], i12[i - 1], i22[i - 1]);
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
                    mueller_out.WriteLine("{0,7:f1} {1,16:e7} {2,16:e7} {3,16:e7} {4,16:e7}", dang[i - 1], mue[0, 0, jc - 1, i - 1], mue[0, 1, jc - 1, i - 1], mue[0, 2, jc - 1, i - 1], mue[0, 3, jc - 1, i - 1]);
                    mueller_out.WriteLine("        {0,16:e7} {1,16:e7} {2,16:e7} {3,16:e7}", mue[1, 0, jc - 1, i - 1], mue[1, 1, jc - 1, i - 1], mue[1, 2, jc - 1, i - 1], mue[1, 3, jc - 1, i - 1]);
                    mueller_out.WriteLine("        {0,16:e7} {1,16:e7} {2,16:e7} {3,16:e7}", mue[2, 0, jc - 1, i - 1], mue[2, 1, jc - 1, i - 1], mue[2, 2, jc - 1, i - 1], mue[2, 3, jc - 1, i - 1]);
                    mueller_out.WriteLine("        {0,16:e7} {1,16:e7} {2,16:e7} {3,16:e7}", mue[3, 0, jc - 1, i - 1], mue[3, 1, jc - 1, i - 1], mue[3, 2, jc - 1, i - 1], mue[3, 3, jc - 1, i - 1]);
                }
            }
            mueller_out.WriteLine("phi (in degrees): 360");
            for (int i = 1; i <= nang2; i++)
            {
                mueller_out.WriteLine("{0,7:f1} {1,16:e7} {2,16:e7} {3,16:e7} {4,16:e7}", dang[i - 1], mue[0, 0, 0, i - 1], mue[0, 1, 0, i - 1], mue[0, 2, 0, i - 1], mue[0, 3, 0, i - 1]);
                mueller_out.WriteLine("        {0,16:e7} {1,16:e7} {2,16:e7} {3,16:e7}", mue[1, 0, 0, i - 1], mue[1, 1, 0, i - 1], mue[1, 2, 0, i - 1], mue[1, 3, 0, i - 1]);
                mueller_out.WriteLine("        {0,16:e7} {1,16:e7} {2,16:e7} {3,16:e7}", mue[2, 0, 0, i - 1], mue[2, 1, 0, i - 1], mue[2, 2, 0, i - 1], mue[2, 3, 0, i - 1]);
                mueller_out.WriteLine("        {0,16:e7} {1,16:e7} {2,16:e7} {3,16:e7}", mue[3, 0, 0, i - 1], mue[3, 1, 0, i - 1], mue[3, 2, 0, i - 1], mue[3, 3, 0, i - 1]);
            }
            mueller_out.Close();
            OUT.Close();
        }
        private void button6_Click(object sender, EventArgs e)
        {
            button2.Enabled = button5.Enabled = button6.Enabled = button3.Enabled = false;
            button1.Enabled = true;
        }
    }
}