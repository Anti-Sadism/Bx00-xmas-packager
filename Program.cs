using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;

/*
 * Student ID: 000185550/00711016
 * KSU Undergraduate Research
 * 
 * 
 * 
 * Wu! It works, I told you it would work! -RC
 * 
 * 
 * 
 * 
 * Some notes, so I don't have to re-google these:
 * 
 * 
 *      The range of DCT for 4x4 blocks 
 *      of pixel intensities between (0,255) 
 *      is (-4080,4080) (255*4*4)
 *      
 * 
 *      The range of DCT for 4x4x4 blocks 
 *      of pixel intensities between (0,255) 
 *      is (-16320,16320) (255*4*4*4)
 *      
 * 
 *      Look at the size of the transition matrix:
 *      
 *      16,646,400 matrix entities for DCT 4x4 === (HUGE)
 *      266,342,400 matrix entities DCT 4x4 === (TOO BIG)
 *      
 *      So we use markov_t = 10
 *      and truncate coefficients at 10
 *      to govern the matrix's size 
 *      so that it doesn't grow too large
 *      our transition matrix is 10x10
 *      
 *      With [1] and markov_t =10, each group gives us 
 *      
 *      10x10 intra frame features
 *      10x10 inter frame features
 *      3 abs central moment features
 *      1 kurtosis feature
 *      1 skewness feature
 *      
 *      For a 205-dimensional feature vector
 * 
 *      [1] performs YouTube Steganalysis on a frame by frame (Group by Group) basis NOT a video by video basis.
 *          This means that the cluster representatives are not entire videos, but instead, groups (4-frames).
 *          
 *      Explained by authors:
 *      
 *                   * 4.2 Experiment Setup
 *                  As shown in [18], unsupervised classification is more
 *                  practical, even the performance is not as good as supervised
 *                  classification. In our experiments, K-means clustering,
 *                  which is a very basic and widely used clustering algorithm
 *                  with cosine distance was considered. The performance is
 *                  measured with True Positive Rate (TPR) and True Negative
 *                  Rate (TNR). TPR represents the probability of that stego is
 *                  classified into stego. And TNR means the probability of that
 *                  cover is classified into cover. One thing should be
 *                  remembered is that our proposed scheme is based on a
 *                  frame-by-frame or group-by-group method.
 * 
 * 
 *      [1] uses k-means clustering with 2 clusters, "Steg" and "Clean". 
 *      It is up to you to decide the initial clusters and k.
 *      
 * 
 *      Try to decide initial centers as the most extreme example of high activity steganography  group for "Steg", 
 *      and the smoothest low-motion video group for "Clean"
 *      
 * 
 *      The steganography algorithms that we use are different than those used in [1] but steganalysis is the same.
 *      
 *      The TPR (Total Positive Rate) detection is based on how many groups are correctly classified in the "Steg" cluster per total frames of a video
 *      The TNR (Total Negative Rate) detection is based on how many groups are correctly classified in the "Clean" cluster per total frames of a video
 *      
 *      A Markov() represents a single "Group", mentioned in the paper
 *      A Markov()[] represents a collection of "Group" objects.
 *      
 *      A Markov()[] is used in this case to represent a collection of groups in a cluster.
 
        See Main() for a full and detailed explanation of Steganalysis.
 
 
 */

namespace H264Stego
{
    class Program
    {
        static long filesizespecified;
        static bool loop = true;
        static long pos;
        static double[] acc;

        static long pos_bit;
        static long pos_rep;
        static int tmod = 0;
        static byte[] finalFile;
        static int markov_t = 10; // as in [1]

        public class Frequency_4x4
        {
            public double[,] frequencies;

            public Frequency_4x4()
            {
                frequencies = new double[4, 4];
            }
        }

        public class Markov
        {
            public double[,,] markov;

            public Markov(int t1, int t2)
            {
                markov = new double[t1, t2, 3];
            }

            public Markov()
            {
                markov = new double[markov_t, markov_t, 3];
            }
        
        }

        public class FeatureVector
        {
            public double[] vec;
            public int id = 0;
            public FeatureVector(int size, int id)
            {
                vec = new double[size];
                this.id = id;
            }

        }

        public class Frequency_4x4x4
        {
            public double[,,] frequencies;

            public Frequency_4x4x4()
            {
                frequencies = new double[4, 4, 4];
            }
        }

        /// <summary>
        /// Generates a DCT3D "group" of 4 adjacent video frames' transformed into 4 adjacent coefficient planes, see [1]
        /// </summary>
        /// <param name="b0">First Frame</param>
        /// <param name="b1">Second Frame</param>
        /// <param name="b2">Fourth Frame</param>
        /// <param name="b3">Third Frame</param>
        /// <returns></returns>
        static double[,,] M_Group(Bitmap b0, Bitmap b1, Bitmap b2, Bitmap b3)
        {
            double[, ,] coefficients = new double[b0.Width, b0.Height, 4];
            Frequency_4x4x4 f = new Frequency_4x4x4();
            Bitmap src0 = new Bitmap(4, 4);
            Bitmap src1 = new Bitmap(4, 4);
            Bitmap src2 = new Bitmap(4, 4);
            Bitmap src3 = new Bitmap(4, 4);

            for (int i = 0; i < b0.Height; i+=4)
                for (int j = 0; j < b0.Width; j+=4)
                {
                
                    src0 = b0.Clone(new Rectangle(j, i, 4, 4), b0.PixelFormat);
                    src1 = b1.Clone(new Rectangle(j, i, 4, 4), b1.PixelFormat);
                    src2 = b2.Clone(new Rectangle(j, i, 4, 4), b2.PixelFormat);
                    src3 = b3.Clone(new Rectangle(j, i, 4, 4), b3.PixelFormat);

                    f = DCT3D_4x4x4(src0, src1, src2, src3);

                    for (int z = 0; z < 4; z++)
                        for (int y = 0; y < src0.Height; y++)
                            for (int x = 0; x < src0.Width; x++)
                                coefficients[j+x,i+y,z] = f.frequencies[x,y,z];

                    src0.Dispose();
                    src1.Dispose();
                    src2.Dispose();
                    src3.Dispose();
                }
            
            return coefficients;
        }

                /// <summary>
                /// Alternative DCT, faster than 
                /// </summary>
                /// <param name="b0"></param>
                /// <param name="b1"></param>
                /// <param name="b2"></param>
                /// <param name="b3"></param>
                /// <returns></returns>
        static Frequency_4x4x4 Fast_DCT3D_4x4x4(Bitmap b0, Bitmap b1, Bitmap b2, Bitmap b3)
        {

            Bitmap[] b = {b0, b1, b2, b3};

            Frequency_4x4x4 f = new Frequency_4x4x4();
            int x = 0;
            int y = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int z = 0;

            double summand = 0;

            for (z = 0; z < 4; z++)
            {
                for (y = 0; y < 4; y++)
                {
                    for (x = 0; x < 4; x++)
                    {
                        for (k = 0; k < 4; k++)
                        {
                            for (i = 0; i < 4; i++)
                            {
                                for (j = 0; j < 4; j++)
                                {
                                    summand += (b[z].GetPixel(j, i).ToArgb() * Math.Cos(Math.PI * z * ((2 * k) + 1) / 8) * Math.Cos(Math.PI * y * ((2 * i) + 1) / 8) * Math.Cos(Math.PI * x * ((2 * j) + 1) / 8));
                                }
                            }
                        }
                        //summand = 0.5 * summand;
                        if (x == 0) summand *= (0.5);
                        if (y == 0) summand *= (0.5);
                        if (z == 0) summand *= (0.5);
                        if (x > 0) summand *= (1.0 / Math.Sqrt(2));
                        if (y > 0) summand *= (1.0 / Math.Sqrt(2));
                        if (z > 0) summand *= (1.0 / Math.Sqrt(2));
                        f.frequencies[x, y, z] = summand;
                        summand = 0;
                    }
                }
            }
            return f;
        }
        
        static Frequency_4x4x4 DCT3D_4x4x4(Bitmap b0, Bitmap b1, Bitmap b2, Bitmap b3)
        {

            Bitmap[] b = {b0, b1, b2, b3};

            Frequency_4x4x4 f = new Frequency_4x4x4();
            int x = 0;
            int y = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int z = 0;

            double summand = 0;

            for (z = 0; z < 4; z++)
            {
                for (y = 0; y < 4; y++)
                {
                    for (x = 0; x < 4; x++)
                    {
                        for (k = 0; k < 4; k++)
                        {
                            for (i = 0; i < 4; i++)
                            {
                                for (j = 0; j < 4; j++)
                                {
                                    summand += (b[z].GetPixel(j, i).ToArgb() * Math.Cos(Math.PI * z * ((2 * k) + 1) / 8) * Math.Cos(Math.PI * y * ((2 * i) + 1) / 8) * Math.Cos(Math.PI * x * ((2 * j) + 1) / 8));
                                }
                            }
                        }
                        //summand = 0.5 * summand;
                        if (x == 0) summand *= (0.5);
                        if (y == 0) summand *= (0.5);
                        if (z == 0) summand *= (0.5);
                        if (x > 0) summand *= (1.0 / Math.Sqrt(2));
                        if (y > 0) summand *= (1.0 / Math.Sqrt(2));
                        if (z > 0) summand *= (1.0 / Math.Sqrt(2));
                        f.frequencies[x, y, z] = summand;
                        summand = 0;
                    }
                }
            }
            return f;
        }

        static Frequency_4x4 DCT2D_4x4(Bitmap b)
        {

            Frequency_4x4 f = new Frequency_4x4();
            int x = 0;
            int y = 0;
            int i = 0;
            int j = 0;

            double summand = 0;

            for (y = 0; y < 4; y++)
            {
                for (x = 0; x < 4; x++)
                {
                    for (i = 0; i < 4; i++)
                    {
                        for (j = 0; j < 4; j++)
                        {
                            summand += (b.GetPixel(j, i).ToArgb() * Math.Cos(Math.PI * y * ((2 * i) + 1) / 8) * Math.Cos(Math.PI * x * ((2 * j) + 1) / 8));
                        }
                    }
                    //summand = 0.5 * summand;
                    if (x == 0) summand *= (0.5);
                    if (y == 0) summand *= (0.5);
                    if (x > 0) summand *= (1.0 / Math.Sqrt(2));
                    if (y > 0) summand *= (1.0 / Math.Sqrt(2));
                    f.frequencies[x, y] = summand;
                    summand = 0;
                }
            }
            return f;
        }

        static Bitmap iDCT2D_4x4(Frequency_4x4 f)
        {
            Bitmap b = new Bitmap(4, 4);
            int x = 0;
            int y = 0;
            int i = 0;
            int j = 0;

            double summand = 0;
            double temp = 0;
            for (y = 0; y < 4; y++)
            {
                for (x = 0; x < 4; x++)
                {
                    for (i = 0; i < 4; i++)
                    {
                        for (j = 0; j < 4; j++)
                        {
                            temp = f.frequencies[j, i] * Math.Cos(Math.PI * i * ((2 * y) + 1) / 8) * Math.Cos(Math.PI * j * ((2 * x) + 1) / 8);
                            if (i == 0) temp *= (0.5);
                            if (j == 0) temp *= (0.5);
                            if (i > 0) temp *= (1.0 / Math.Sqrt(2));
                            if (j > 0) temp *= (1.0 / Math.Sqrt(2));
                            //summand *= 0.5;
                            summand += temp;
                        }
                    }
                    //summand = 0.5 * summand;

                    b.SetPixel(x, y, Color.FromArgb((int)Math.Round(summand)));
                    summand = 0;
                }
            }
            return b;
        }


        static Bitmap[] iDCT3D_4x4x4(Frequency_4x4x4 f)
        {
            Bitmap[] b = new Bitmap[4];
            b[0] = new Bitmap(4, 4);
            b[1] = new Bitmap(4, 4);
            b[2] = new Bitmap(4, 4);
            b[3] = new Bitmap(4, 4);

            int x = 0;
            int y = 0;
            int z = 0;
            int i = 0;
            int j = 0;
            int k = 0;


            double summand = 0;
            double temp = 0;
            for (z = 0; z < 4; z++)
            {
                for (y = 0; y < 4; y++)
                {
                    for (x = 0; x < 4; x++)
                    {
                        for (k = 0; k < 4; k++)
                        {
                            for (i = 0; i < 4; i++)
                            {
                                for (j = 0; j < 4; j++)
                                {
                                    temp = f.frequencies[j, i, k] * Math.Cos(Math.PI * k * ((2 * z) + 1) / 8) * Math.Cos(Math.PI * i * ((2 * y) + 1) / 8) * Math.Cos(Math.PI * j * ((2 * x) + 1) / 8);
                                    if (i == 0) temp *= (0.5);
                                    if (j == 0) temp *= (0.5);
                                    if (i > 0) temp *= (1.0 / Math.Sqrt(2));
                                    if (j > 0) temp *= (1.0 / Math.Sqrt(2));
                                    //summand *= 0.5;
                                    summand += temp;
                                }
                            }
                        }
                        //summand = 0.5 * summand;

                        b[z].SetPixel(x, y, Color.FromArgb((int)Math.Round(summand)));
                        summand = 0;
                    }
                }
            }
            return b;
        }


        /// <summary>
        /// Embed the bit within the 8 dimensional vector
        /// </summary>
        /// <param name="f"></param>
        /// <returns></returns>
        static Frequency_4x4 vector_Embed8D(Frequency_4x4 f, int T, byte bit)
        {
            double distance = Math.Pow(f.frequencies[0, 0], 2) + Math.Pow(f.frequencies[1, 0], 2) + Math.Pow(f.frequencies[0, 1], 2) + Math.Pow(f.frequencies[0, 2], 2) + Math.Pow(f.frequencies[1, 1], 2) + Math.Pow(f.frequencies[2, 0], 2) + Math.Pow(f.frequencies[3, 0], 2) + Math.Pow(f.frequencies[2, 1], 2);
            distance = Math.Sqrt(distance);
            double distance_prime = Math.Round(distance / (double)T);
            if (bit > 0)
                distance_prime += 0.25;
            else
                distance_prime -= 0.25;
            distance_prime = distance_prime * (double)T;
            double coeff = distance_prime / distance;
            f.frequencies[0, 0] *= coeff;
            f.frequencies[1, 0] *= coeff;
            f.frequencies[0, 1] *= coeff;
            f.frequencies[0, 2] *= coeff;
            f.frequencies[1, 1] *= coeff;
            f.frequencies[2, 0] *= coeff;
            f.frequencies[3, 0] *= coeff;
            f.frequencies[2, 1] *= coeff;

            return f;
        }

        /// <summary>
        /// Retrieves the bit from within the 8 dimensional vector
        /// </summary>
        /// <param name="f"></param>
        /// <returns></returns>
        static byte vector_Retrieve8D(Frequency_4x4 f, int T)
        {
            byte bit = 0;
            double distance = Math.Pow(f.frequencies[0, 0], 2) + Math.Pow(f.frequencies[1, 0], 2) + Math.Pow(f.frequencies[0, 1], 2) + Math.Pow(f.frequencies[0, 2], 2) + Math.Pow(f.frequencies[1, 1], 2) + Math.Pow(f.frequencies[2, 0], 2) + Math.Pow(f.frequencies[3, 0], 2) + Math.Pow(f.frequencies[2, 1], 2);
            distance = Math.Sqrt(distance);
            double distance_prime = distance / (double)T;
            distance_prime = distance_prime - Math.Round(distance_prime);
            if (distance_prime >= 0)
                bit = 1;
            else
                bit = 0;

            return bit;
        }

        /// <summary>
        /// Calculates the alternating sign matrix for a spatial domain block
        /// </summary>
        /// <param name="block">The block of pixel values</param>
        /// <returns></returns>
        static int[,] matrix_SGN(double[,] block)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);

            int[,] sgn = new int[N, M];
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    sgn[j, i] = ((((i) % 2) == ((j) % 2)) ? 1 : -1);
            return sgn;
        }

        /// <summary>
        /// Calculates the arithmetic difference \alpha of the block
        /// </summary>
        /// <param name="block"></param>
        static int arithmetic_Difference(double[,] block, int[,] sgn)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);

            int alpha = 0;
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    alpha += (sgn[j, i] * (int)block[j, i]);
            return alpha;
        }

        /// <summary>
        /// Calculates the arithmetic threshold of the block
        /// </summary>
        /// <param name="block">The Block</param>
        /// <param name="alpha">The Arithmetic Difference</param>
        /// <param name="gamma">Gamma Value</param>
        /// <returns></returns>
        static double[,] arithmetic_Threshold(double[,] block, int alpha, double gamma, int T)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                {
                    if ((((i + 1) % 2) == ((j + 1) % 2)) && alpha > T) block[j, i] += gamma;
                    if ((((i + 1) % 2) != ((j + 1) % 2)) && alpha < -T) block[j, i] += gamma;
                }
            return block;
        }

        /// <summary>
        /// Embeds a bit within the current block. Returns a block to be stored.
        /// </summary>
        /// <param name="bit"></param>
        /// <param name="alpha"></param>
        /// <param name="block"></param>
        /// <param name="gamma_prime"></param>
        /// <param name="T"></param>
        /// <returns></returns>
        static double[,] arithmetic_Embed(byte bit, int alpha, double[,] block, double gamma_prime, int T)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);

            if (bit > 0)
            {
                for (int i = 0; i < M; i++)
                    for (int j = 0; j < N; j++)
                    {
                        if (((((i + 1) % 2) == ((j + 1) % 2)) && alpha >= 0 && alpha <= T)) block[j, i] += gamma_prime;
                        if (((((i + 1) % 2) != ((j + 1) % 2)) && alpha < 0 && alpha >= -T)) block[j, i] += gamma_prime;
                    }
            }
            return block;
        }

        /// <summary>
        /// Retrieves a bit from the current block
        /// </summary>
        /// <param name="block"></param>
        /// <param name="T"></param>
        /// <param name="G"></param>
        /// <returns></returns>
        static byte arithmetic_Retrieve(double[,] block, int T, int G)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);

            int alpha = arithmetic_Difference(block, matrix_SGN(block));

            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                {
                    if (alpha >= -T && alpha <= T) return 0;
                    if ((alpha > T && alpha <= (2 * T) + G) || (alpha < -T && alpha >= -((2 * T) + G))) return 1;
                }
            Console.Write("?"); //salt
            return 1;
        }

        static double[,] Array2DFromBitmap(Bitmap b)
        {
            double[,] a = new double[b.Width, b.Height];
            for (int i = 0; i < b.Height; i++)
                for (int j = 0; j < b.Width; j++)
                    a[j, i] = (b.GetPixel(j, i).ToArgb());
            return a;
        }

        static Bitmap BitmapFromArray2D(double[,] a)
        {
            Bitmap b = new Bitmap(a.GetLength(0), a.GetLength(1));
            for (int i = 0; i < b.Height; i++)
                for (int j = 0; j < b.Width; j++)
                    b.SetPixel(j, i, Color.FromArgb((int)a[j, i]));
            return b;
        }

        static byte[] Array1DFromBitmap(Bitmap b, int size)
        {
            byte[] a = new byte[size];
            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            int stride = data.Stride;
            unsafe
            {
                byte* ptr = (byte*)data.Scan0;
                for (int i = 0; i < b.Height; i++)
                    for (int j = 0; j < b.Width; j++)
                    {
                        //a[(b.Width * i) + j] = (byte)(b.GetPixel(j, i).ToArgb());
                        // layer.GetBitmap().SetPixel(x, y, m_colour);
                        a[(b.Width * i) + j] = ptr[(j * 3) + i * stride];
                        //ptr[(j * 3) + i * stride + 1] = m_colour.G;
                        //ptr[(j * 3) + i * stride + 2] = m_colour.R;
                    }
            }
            b.UnlockBits(data);
            return a;
        }

        static byte[] Array1DFromBitmapSafe(Bitmap b, int size)
        {
            byte[] a = new byte[size];
            //BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            //int stride = data.Stride;
            //unsafe
            //{
            //byte* ptr = (byte*)data.Scan0;    
            for (int i = 0; i < b.Height; i++)
                for (int j = 0; j < b.Width; j++)
                {
                    a[(b.Width * i) + j] = (byte)(b.GetPixel(j, i).ToArgb());
                    // layer.GetBitmap().SetPixel(x, y, m_colour);
                    //a[(b.Width * i) + j] = ptr[(j * 3) + i * stride];
                    //ptr[(j * 3) + i * stride + 1] = m_colour.G;
                    //ptr[(j * 3) + i * stride + 2] = m_colour.R;
                }
            //}
            //b.UnlockBits(data);
            return a;
        }

        static Bitmap BitmapFromArray1DSafe(byte[] a, int width, int height)
        {
            Bitmap b = new Bitmap(width, height);
            //BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            //int stride = data.Stride;
            //unsafe
            //{
            //byte* ptr = (byte*)data.Scan0;
            for (int i = 0; i < a.Length; i++)
                if (i / width < 720)
                    //for (int j = 0; j < b.Width; j++)
                    //{
                    //a[(b.Width * i) + j] = (byte)(b.GetPixel(j, i).ToArgb());
                    // layer.GetBitmap().SetPixel(x, y, m_colour);
                    //a[(b.Width * i % width) + i / width] = ptr[(j * 3) + i * stride];
                    //ptr[(i % width * 3) + ((i / width) * stride)] = a[i];
                    //ptr[(j * 3) + i * stride + 1] = m_colour.G;
                    //ptr[(j * 3) + i * stride + 2] = m_colour.R;
                    b.SetPixel(i % width, i / width, Color.FromArgb((int)a[i]));
            //}
            //
            //}
            //b.UnlockBits(data);
            return b;
        }

        /// <summary>
        /// Gives a running sum of an array
        /// </summary>
        /// <param name="input"></param>
        /// <param name="frames"></param>
        /// <returns></returns>
        static double[] AccArray(byte[] input, double[] prev)
        {
            byte[] result = new byte[input.Length*8];

            for (int i = 0; i < input.Length; i++)
                for (int k = 0; k < 8; k++)
                {
                    result[(8*i) + k] = (byte)(((input[i] >> k) % 2));
                    prev[(8 * i) + k] += result[(8 * i) + k];
                }
                    // each slot of result contains either a zero or one
            return prev;            
        }

        /// <summary>
        /// Gives a running average of an array
        /// </summary>
        /// <param name="input"></param>
        /// <param name="frames"></param>
        /// <returns></returns>
        static double[] AvgArray(double[] input, double frames)
        {
            double[] result = new double[input.Length];

            for (int i = 0; i < input.Length; i++)
                    result[i] = input[i] / frames;
                
            // each slot of result contains either a zero or one
            return result;
        }

        /// <summary>
        /// Calculates the P value
        /// </summary>
        /// <param name="input"></param>
        /// <param name="frames"></param>
        /// <returns></returns>
        static double P(double[] input)
        {
            double result = 0;

            for (int i = 0; i < input.Length; i++)
                    result += input[i];

            result /= input.Length;
            return Math.Pow(((2 * result) - 1), 2);
        }





        static Bitmap BitmapFromArray1D(byte[] a, int width, int height)
        {
            Bitmap b = new Bitmap(width, height);
            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            int stride = data.Stride;
            unsafe
            {
                byte* ptr = (byte*)data.Scan0;
                for (int i = 0; i < a.Length; i++)
                    if (i / width < 720)
                        //for (int j = 0; j < b.Width; j++)
                        //{
                        //a[(b.Width * i) + j] = (byte)(b.GetPixel(j, i).ToArgb());
                        // layer.GetBitmap().SetPixel(x, y, m_colour);
                        //a[(b.Width * i % width) + i / width] = ptr[(j * 3) + i * stride];
                        ptr[(i % width * 3) + ((i / width) * stride)] = a[i];
                //ptr[(j * 3) + i * stride + 1] = m_colour.G;
                //ptr[(j * 3) + i * stride + 2] = m_colour.R;
                //}
                //b.SetPixel(i%width, i/width, Color.FromArgb((int)a[i]));
            }
            b.UnlockBits(data);
            return b;
        }

        /// <summary>
        /// Embeds a byte[] within the spatial domain of the video and returns the bitmap.
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        static Bitmap spatial_Embed(Bitmap b, int block_size, int T, int G, byte[] bytes, int rep, int trep)
        {
            double[,] block = new double[block_size, block_size];
            double gamma = Gamma(block, T, G);
            double gamma_prime = Gamma_Prime(block, T, G);
            int[,] sgn = matrix_SGN(block);

            long temppos = pos;
            long temppos_rep = pos_rep;
            long temppos_bit = pos_bit;

            if (loop) pos = 0;

            if (loop) pos_bit = 0;

            if (loop) pos_rep = 0;

            Bitmap src = new Bitmap(block_size, block_size);
            int x = 0;
            int y = 0;

            for (; ; pos++)
            {

                if (pos == bytes.Length && loop == true) return b;
                if (pos == bytes.Length) pos = 0;




                for (pos_bit = 0; pos_bit < 8; pos_bit++)
                {
                    //int y = (block_size * (((8 *  i) + (  k))) / b.Width );
                    //int x = (block_size * (((8 *  i) + (  k))) % b.Width );




                    for (pos_rep = 0; pos_rep < rep; pos_rep++)
                    {
                        src = b.Clone(new Rectangle(x, y, block_size, block_size), b.PixelFormat);
                        block = Array2DFromBitmap(src);
                        byte bit = (byte)((bytes[pos] >> (byte)pos_bit) % 2);
                        int alpha = arithmetic_Difference(block, sgn);
                        block = arithmetic_Threshold(block, alpha, gamma, T);

                        alpha = arithmetic_Difference(block, sgn);//See below

                        if (alpha >= -T && alpha <= T) //Thank you Alavianamehr! for leaving important details out after plagiari- "paraphrasing" Xian
                            block = arithmetic_Embed(bit, alpha, block, gamma_prime, T);
                        if (alpha > T && alpha < -T)
                            pos_bit--;

                        src = BitmapFromArray2D(block);

                        //g.DrawImage(src, x, y, new Rectangle(0, 0, block_size, block_size), GraphicsUnit.Pixel);
                        for (int m = 0; m < block_size; m++)
                            for (int n = 0; n < block_size; n++)
                                b.SetPixel(x + n, y + m, src.GetPixel(n, m));
                        //src.Dispose();

                        x += block_size;
                        if (x >= b.Width) { y += block_size; x = 0; }
                
                        
                        if (y > b.Height+block_size)
                        {
                            tmod++;
                            if (!loop)
                            {
                                if (tmod % (trep + 1) != 0)
                                {
                                    //reset the values
                                    pos = temppos;
                                    pos_bit = temppos_bit;
                                    pos_rep = temppos_rep;
                                    tmod = 0;
                                }
                                else tmod = 0;
                            }
                            //otherwise, allow the parameters to advance

                            return b;
                        }
                    }
                    tmod++;



                }
            }
            //return b;
        }

        /// <summary>
        /// Embeds a byte[] within the transform domain of the video and returns the bitmap.
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        static Bitmap transform_Embed(Bitmap b, int T, byte[] bytes, int rep)
        {

            int block_size = 4;
            double[,] block = new double[block_size, block_size];

            Bitmap src = new Bitmap(block_size, block_size);
            Frequency_4x4 f = new Frequency_4x4();
            int x = 0;
            int y = 0;

            if (loop) pos = 0;

            if (loop) pos_bit = 0;

            if (loop) pos_rep = 0;

            for (; ; pos++)
            {
                if (pos == bytes.Length && loop == true) return b;
                if (pos == bytes.Length) pos = 0;

                for (pos_bit = 0; pos_bit < 8; pos_bit++)
                {
                    for (pos_rep = 0; pos_rep < rep; pos_rep++)
                    {
                        src = b.Clone(new Rectangle(x, y, block_size, block_size), b.PixelFormat);
                        f = DCT2D_4x4(src);
                        byte bit = (byte)((bytes[pos] >> (byte)(pos_bit)) % 2);
                        f = vector_Embed8D(f, T, bit);

                        src = iDCT2D_4x4(f);
                        for (int m = 0; m < block_size; m++)
                            for (int n = 0; n < block_size; n++)
                                b.SetPixel(x + n, y + m, src.GetPixel(n, m));

                        //int foo = src.GetPixel(0, 0).ToArgb();
                        //foo = 0;

                        x += block_size;
                        if (x >= b.Width) { y += block_size; x = 0; }
                        if (y >= b.Height) return b;

                    }
                }
            }
            //return b;
        }

        /// <summary>
        /// Retrieves a byte[] within the spatial domain of the video and returns the byte[].
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        static byte[] spatial_Retrieve(Bitmap b, int block_size, int T, int G, byte[] bytes, int rep, int trep)
        {
            double[,] block = new double[block_size, block_size];
            double gamma = Gamma(block, T, G);
            double gamma_prime = Gamma_Prime(block, T, G);

            Bitmap src = new Bitmap(block_size, block_size);
            int x = 0;
            int y = 0;
            byte bit = 0;
            byte avg = 0;

            long temppos = pos;
            long temppos_rep = pos_rep;
            long temppos_bit = pos_bit;


            if (loop) pos = 0;


            if (loop) pos_bit = 0;

            if (loop) pos_rep = 0;

            byte[] newbytes = new byte[bytes.Length];

            for (; ; pos++)
            {
                if (pos == bytes.Length && loop == true) return newbytes;
                if (pos == bytes.Length && loop == false) pos = 0;


                for (pos_bit = 0; pos_bit < 8; pos_bit++)
                {
                    for (pos_rep = 0; pos_rep < rep; pos_rep++)
                    {
                        src = b.Clone(new Rectangle(x, y, block_size, block_size), b.PixelFormat);
                        block = Array2DFromBitmap(src);

                        int alpha = arithmetic_Difference(block, matrix_SGN(block));
                        //if (alpha <= ((2*T) + G) && alpha >= -((2*T)+ G)) //see my frustration above at Alavianmehr et. als' "original" work
                        //{
                        avg += arithmetic_Retrieve(block, T, G);


                        //}
                        //if (alpha > ((2 * T) + G) || alpha < -((2 * T) + G)) k--;

                        x += block_size;
                        if (x >= b.Width) { y += block_size; x = 0; }
                        if (y+block_size > b.Height)
                        {
                            tmod++;
                            if (!loop)
                            {
                                if (tmod % trep != 0)
                                {
                                    //reset the values
                                    pos = temppos;
                                    pos_bit = temppos_bit;
                                    pos_rep = temppos_rep;
                                    tmod = 0;
                                }
                                else
                                    tmod = 0;
                            }
                            //otherwise, allow the parameters to advance
                            return newbytes;
                        }
                    }
                    bit = (byte)Math.Round((double)avg / (double)rep);
                    newbytes[pos] = (byte)(newbytes[pos] >> 1);
                    newbytes[pos] += (byte)(128 * bit);
                    avg = 0;
                }
            }
            //return newbytes;
        }



        /// <summary>
        /// Retrieves a byte[] within the transform domain of the video and returns the byte[].
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        static byte[] transform_Retrieve(Bitmap b, byte[] bytes, int T, int rep)
        {
            double[,] block = new double[4, 4];

            if (loop) pos = 0;

            if (loop) pos_bit = 0;

            if (loop) pos_rep = 0;

            Bitmap src = new Bitmap(4, 4);
            int x = 0;
            int y = 0;
            byte bit = 0;
            byte avg = 0;
            Frequency_4x4 f = new Frequency_4x4();
            byte[] newbytes = new byte[bytes.Length];

            for (; ; pos++)
            {
                if (pos == bytes.Length && loop == true) return newbytes;
                if (pos == bytes.Length) pos = 0;

                for (pos_bit = 0; pos_bit < 8; pos_bit++)
                {
                    for (pos_rep = 0; pos_rep < rep; pos_rep++)
                    {
                        src = b.Clone(new Rectangle(x, y, 4, 4), b.PixelFormat);
                        f = DCT2D_4x4(src);
                        //block = Array2DFromBitmap(src);

                        //int alpha = arithmetic_Difference(block, matrix_SGN(block));
                        //if (alpha <= ((2*T) + G) && alpha >= -((2*T)+ G)) //see my frustration above at Alavianmehr et. als' "original" work
                        //{
                        avg += vector_Retrieve8D(f, T);


                        //}
                        //if (alpha > ((2 * T) + G) || alpha < -((2 * T) + G)) k--;

                        x += 4;
                        if (x >= b.Width) { y += 4; x = 0; }
                        if (y >= b.Height)
                            return newbytes;
                    }
                    bit = (byte)Math.Round((double)avg / (double)rep);
                    newbytes[pos] = (byte)(newbytes[pos] >> 1);
                    newbytes[pos] += (byte)(128 * bit);
                    avg = 0;
                }
            }
            //return newbytes;
        }


        static double MSE(Bitmap b1, Bitmap b2)
        {
            double mse = 0;
            for (int i = 0; i < b1.Height; i++)
                for (int j = 0; j < b1.Width; j++)
                    mse += Math.Pow(b1.GetPixel(j, i).ToArgb() - b2.GetPixel(j, i).ToArgb(), 2);
            mse /= (b1.Width * b1.Height);
            return mse;
        }

        static double PSNR(double MSE)
        {
            return 10 * Math.Log10(Math.Pow(255, 2) / MSE);
        }

        static Bitmap GetBMPFromFrame(string filename, int frame, int VideoWidth, int VideoHeight)
        {
            Bitmap b = new Bitmap(VideoWidth, VideoHeight);

            return b;
        }

        /// <summary>
        /// Calculates Gamma Value
        /// </summary>
        /// <param name="block">The block</param>
        /// <returns></returns>
        static double Gamma(double[,] block, int T, int G)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);

            return (double)Math.Ceiling((2 * ((2 * G) + T)) / (double)(N * M));
        }

        /// <summary>
        /// Calculates Gamma_Prime Value
        /// </summary>
        /// <param name="block"></param>
        /// <param name="T"></param>
        /// <param name="G"></param>
        /// <returns></returns>
        static double Gamma_Prime(double[,] block, int T, int G)
        {
            int N = block.GetLength(0);
            int M = block.GetLength(1);


            return (double)Math.Ceiling((2 * (T + G)) / (double)(N * M));
        }

        /// <summary>
        /// Gets a byte[] from a file.
        /// </summary>
        /// <param name="fullFilePath"></param>
        /// <param name="offset"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static byte[] GetBytesFromFile(string fullFilePath, int offset, int length)
        {
            FileStream fs = null;
            try
            {
                fs = File.OpenRead(fullFilePath);
                byte[] bytes = new byte[length];

                fs.Seek(offset, SeekOrigin.Begin);
                fs.Read(bytes, 0, length);
                return bytes;
            }
            finally
            {
                if (fs != null)
                {
                    fs.Close();
                    fs.Dispose();
                }
            }

        }

        /// <summary>
        /// Saves a byte[] to a file.
        /// </summary>
        /// <param name="fullFilePath"></param>
        /// <param name="offset"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static byte[] SaveBytesToFile(string fullFilePath, int offset, int length, byte[] bytes)
        {
            FileStream fs = null;
            try
            {
                fs = File.OpenWrite(fullFilePath);
                //bytes = new byte[length];
                fs.Seek(offset, SeekOrigin.Begin);
                fs.Write(bytes, 0, length);

                return bytes;
            }
            finally
            {
                if (fs != null)
                {
                    fs.Close();
                    fs.Dispose();
                }
            }

        }



        /// <summary>
        /// Gets a byte[] from a file.
        /// </summary>
        /// <param name="fullFilePath"></param>
        /// <param name="offset"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static byte[] GetBytesFromFile(string fullFilePath, int offset)
        {
            FileStream fs = null;
            try
            {
                fs = File.OpenRead(fullFilePath);
                byte[] bytes = new byte[fs.Length];

                fs.Seek(offset, SeekOrigin.Begin);
                fs.Read(bytes, 0, (int)bytes.Length);
                return bytes;
            }
            finally
            {
                if (fs != null)
                {
                    fs.Close();
                    fs.Dispose();
                }
            }

        }

        static double BER(byte[] a, byte[] b)
        {
            int i = 0;
            int k = 0;
            double ber = 0;

           

            for (i = 0; i < a.Length; i++)
                for (k = 0; k < 8; k++)
                    if (((a[i] >> k) % 2) != ((b[i] >> k) % 2)) 
                        ber++;


            return 100.0 * ber / a.Length / 8;

        }

        /// <summary>
        /// Takes the BER between a chunk of recovered data and the host data at a certain position and length
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        static double BER(byte[] src, byte[] rec, long position, int width)
        {
            int i = 0;
            int k = 0;
            double ber = 0;

            for (i = 0; i < width; i++)
                for (k = 0; k < 8; k++)
                    if (((src[i+position] >> k) % 2) != ((rec[i+position] >> k) % 2)) 
                        ber++;

            return 100.0 * ber / width / 8;

        }

        /// <summary>
        /// Expectation operator
        /// </summary>
        /// <param name="f"></param>
        /// <returns></returns>
        static double Expectation(double[] f)
        {
            
            return m(f);

            /* wow this is embarrasing,  so much code for something entirely simple
            double sum = 0;
            int acc = 1;

            if (f.Length < 500) //O(n^2)
            {
                for (int i = 0; i < f.Length - 1; i++)
                {
                    for (int j = (i + 1); j < f.Length; j++)
                        if (f[i] == f[j]) acc++;
                    f[i] = acc * f[i];
                    acc = 1; // :-)
                }
                for (int i = 0; i < f.Length; i++)
                    sum += f[i];

                sum /= f.Length;
            }

            if (f.Length >= 500)  //O(nlogn)
            {
                double tweedledee, tweedledum = 0;
                
                tweedledee = f[0];
                Array.Sort(f); //O(nlogn)

                for (int i = 1; i < f.Length; i++) //O(n)
                {
                    tweedledum = f[i];
                    acc++;
                    if (tweedledum != tweedledee)
                    {
                        tweedledee = f[i];
                        acc = 0;
                        sum += (acc * tweedledum);
                    }
                }
                    
            }

            return sum;
             */
        }

        /// <summary>
        /// Skewness
        /// </summary>
        /// <param name="mean"></param>
        /// <param name="std"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        static double Skewness(double mean, double std, double[] f)
        {
            for (int i = 0; i < f.Length; i++)
                f[i] = Math.Pow(((f[i] - mean)/std), 3);
            return Expectation(f);
        }

        /// <summary>
        /// Kurtosis
        /// </summary>
        /// <param name="mean"></param>
        /// <param name="std"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        static double Kurtosis(double mean, double std, double[] f)
        {
            for (int i = 0; i < f.Length; i++)
                f[i] = Math.Pow(((f[i] - mean) / std), 4);
            return Expectation(f);
        }

        /// <summary>
        /// Simply calculates the mean
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        static double m(double[] a)
        {
            double sum = 0; for (int i = 0; i < a.Length; i++)
                sum += a[i];
            sum /= a.Length; return sum;
        }

        /// <summary>
        /// Simply calculates the mean
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        static double m(double[,,] a)
        {
            double sum = 0;
                                    for (int i = 0; i < a.GetLength(0); i++)
                                for (int j = 0; j < a.GetLength(1); j++)
                            for (int k = 0; k < a.GetLength(2); k++)
                        sum += a[i, j, k];

            sum /= (a.GetLength(0)*a.GetLength(1)*a.GetLength(2)); return sum;
        }

        /// <summary>
        /// Simply calculates the standard deviation sigma
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        static double std(double[, ,] a, double mean)
        {
            double sum = 0;
            for (int i = 0; i < a.GetLength(0); i++)
                for (int j = 0; j < a.GetLength(1); j++)
                    for (int k = 0; k < a.GetLength(2); k++)
                        sum += Math.Pow((a[i, j, k] - mean), 2);

            sum /= (a.GetLength(0) * a.GetLength(1) * a.GetLength(2));
            sum = Math.Sqrt(sum);

            return sum;
        }

        /// <summary>
        /// Flattens out a double[,,] to a double[]
        /// </summary>
        /// <returns></returns>
        static double[] flatten(double[,,] d)
        {

            double[] foo = new double[d.GetLength(0) * d.GetLength(1) * d.GetLength(2)];

            for (int z = 0; z < d.GetLength(2); z++)
                for (int y = 0; y < d.GetLength(1); y++)
                    for (int x = 0; x < d.GetLength(0); x++)
                        foo[x + (d.GetLength(1)*y) + (z*d.GetLength(2))] = d[x, y, z];

            return foo;
        }


        /// <summary>
        /// Finds the centroid feature of a feature cluster
        /// </summary>
        /// <param name="c">The cluster of feature vectors</param>
        /// <returns></returns>
        static FeatureVector centroid(List<FeatureVector> c)
        {
            FeatureVector fv = new FeatureVector(c[0].vec.Length, -1);
            //double[] elements = new double[1];
            /*
            for (int z = 0; z < c[0].markov.GetLength(2); z++)
                for (int y = 0; y < c[0].markov.GetLength(1); y++)
                    for (int x = 0; x < c[0].markov.GetLength(0); x++)
                    {
                            elements = new double[c.Count];
                            for (int i = 0; i < c.Count(); i++) elements[i] = c[i].markov[x, y, z];
                            mv.markov[x, y, z] = m(elements);
                    }
             * */
            for (int x = 0; x < c[0].vec.Length; x++)
                fv.vec[x] = m(c[x].vec);


            return fv;
        }

        /// <summary>
        /// Finds the Absolute Central Moment of a group
        /// </summary>
        /// <param name="p"></param>
        /// <param name="height">The height of the video frame</param>
        /// <param name="width">The width of the video frame</param>
        /// <param name="f"></param>
        /// <param name="mean"></param>
        /// <returns></returns>
        static double AbsCentralMoment(int p, int height, int width, double[] f, double mean)
        {
            double sum = 0;
            for (int i = 0; i < f.Length; i++)
                sum += Math.Pow((Math.Abs(f[i] - mean)), p);
            sum /= (height * width);
            return sum;
        }


        

        /// <summary>
        /// Reduces the feature dimension of M_Group markov feature transition matrix to (markov_t X markov_t, 2) by reducing the range of coefficients to [0, markov_t]. If any values in the frequency block are above the threshold, then set them to the threshold value.
        /// </summary>
        /// <param name="T"></param>
        /// <param name="f"></param>
        /// <returns></returns>
        static double[,,] Truncate(double[,,] f)
        {
            int T = markov_t;
            for (int i = 0; i < f.GetLength(2); i++)
                for (int j = 0; j < f.GetLength(1); j++)
                    for (int k = 0; k < f.GetLength(0); k++)
                        if (f[k, j, i] > T || f[k, j, i] < -T) f[k, j, i] = (T * (f[k, j, i] / Math.Abs(f[k, j, i])));
            return f;
        }

        /// <summary>
        /// Calculates the cosine distance between two feature vectors.
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        static double cosine(FeatureVector a, FeatureVector b)
        {
            double numerator = 0;
            double denominator_a = 0;
            double denominator_b = 0;

            for (int i = 0; i < a.vec.Length; i++)
                numerator += a.vec[i] * b.vec[i];

            for (int i = 0; i < a.vec.Length; i++)
                denominator_a += Math.Pow(a.vec[i], 2);

            for (int i = 0; i < a.vec.Length; i++)
                denominator_b += Math.Pow(b.vec[i], 2);

            try
            {
                return (1 - (numerator / (Math.Sqrt(denominator_a) * Math.Sqrt(denominator_b))));
            }
            catch
            {
                Console.WriteLine("(^^^) Division by 0: rawr");
            }
            return 0;
        }





        /// <summary>
        /// Generates the 205 dimensional feature vector [Markov, abs(p), skew, kurt], p = {1,2,3} essential for classification
        /// *Applies weights (10x) to Kurtosis and Skewness
        /// </summary>
        /// <param name="mv">Markov Object</param>
        /// <param name="abs1">Absolute Central Moment, p = 1</param>
        /// <param name="abs2">Absolute Central Moment, p = 2</param>
        /// <param name="abs3">Absolute Central Moment, p = 3</param>
        /// <param name="skew">Skewness</param>
        /// <param name="kurt">Kurtosis</param>
        /// <returns></returns>
        static FeatureVector featureVector1(Markov mv, double abs1, double abs2, double abs3, double skew, double kurt, int id)
        {
            double[] vec = new double[2 * (markov_t * markov_t) + 5];
            for (int y = 0; y < markov_t; y++) //intra
                for (int x = 0; x < markov_t; x++)
                    vec[((markov_t) * y) + x] = mv.markov[x, y, 0];
            for (int y = 0; y < markov_t; y++) //inter
                for (int x = 0; x < markov_t; x++)
                    vec[((markov_t) * markov_t) + ((markov_t) * y) + x] = mv.markov[x, y, 1];

            vec[2 * (markov_t * markov_t) + 0] = abs1;
            vec[2 * (markov_t * markov_t) + 1] = abs2;
            vec[2 * (markov_t * markov_t) + 2] = abs3;
            vec[2 * (markov_t * markov_t) + 3] = 10 * skew;
            vec[2 * (markov_t * markov_t) + 4] = 10 * kurt;

            FeatureVector f = new FeatureVector((markov_t * markov_t) + 5, id);

            f.vec = vec; // ur pushing the envelope on ambiguity here chief
            return f;
        }




        /// <summary>
        /// Random initialization k-means with cosine distance
        /// </summary>
        static List<FeatureVector>[] k_means_cosine(FeatureVector[] o, int k, int terminateIn)
        {
            //k means
            FeatureVector[] means = new FeatureVector[k];
            FeatureVector mean = new FeatureVector(markov_t, markov_t);
            List<FeatureVector>[] clusters = new List<FeatureVector>[k];
            double min = 0;

            bool notDone = true;

            //use random initialization
            o = randomShuffle(o);

            Console.WriteLine("\n ..Clustering with " + k + " randomly selected cluster centers");
            Console.WriteLine("First cluster center selected is: Feature Vector " + means[0].id);
            Console.WriteLine("Second cluster center selected is: Feature Vector " + means[1].id);

            //select k random observations and add them to the cluster centers
            for (int i = 0; i < k; i++)
            {
                clusters[i].Add(o[i]);
                means[i] = o[i];
            }
            
            //k means Algorithm
            for (int q = 0; q < terminateIn; q++)
            {
                Console.WriteLine("Iteration: " + q);
                notDone = false;

                //put each element in o[] to the centroid that its closer to
                for (int i = 0; i < o.Length; i++)
                {
                    min = cosine(o[i], means[0]);

                    for (int kk = 0; kk < k; kk++)
                        min = Math.Min(cosine(o[i], means[kk]), min);

                    for (int kk = 0; kk < k; kk++)
                        if (cosine(o[i], means[kk]) == min)
                        {
                            //check if o[i] is not already in the proper cluster
                            if (!clusters[kk].Contains(o[i]))
                            {
                                //o[i] has found a new cluster, move on
                                Console.WriteLine("Feature Vector " + i + " has moved to cluster " + kk);
                                clusters[kk].Add(o[i]);

                                //remove o[i] from other clusters
                                for (int kkk = 0; kkk < kk; kkk++)
                                    if (kkk != kk)
                                        clusters[kkk].Remove(o[i]);
                            
                                //calculate the new centroid after inserting the new member

                                //make sure everything checks out
                                for (int ui = 0; ui < mean.vec.Length; ui++)
                                    if (!(mean.vec[ui] == means[kk].vec[ui])) notDone = true;
                                means[kk] = mean;

                                /*mean = centroid(clusters[kk]);
                                if (mean != means[kk]) notDone = true;
                                means[kk] = mean;
                                */

                                


                                continue;
                            }
                        }
                }

                for (int s = 0; s < k; s++)
                {
                    //Output..
                    Console.WriteLine("= Elements in cluster " + s + " are:");
                    for (int ss = 0; ss < clusters[k].Count(); ss++)
                        Console.WriteLine(ss + ". Feature Vector " + clusters[k][ss].id);
                    Console.WriteLine("= Centroid for cluster " + s + " is: \n <");

                    for (int ss = 0; ss < means[k].vec.Length; ss++)
                        Console.Write(means[k].vec[ss] + ", ");
                    Console.Write("> \n");
                }

                if (!notDone)
                {
                    Console.WriteLine("K means algorithm has been completed after " + q + " steps.");
                    return clusters;
                }
            }

            Console.WriteLine("Did not converge in specified steps: Terminating.");

            return clusters;
        }



        /// <summary>
        /// k means with cosine distance
        /// 
        /// Use cosine distance (1- cosine similarity)
        /// </summary>
        static List<FeatureVector>[] k_means_cosine(FeatureVector[] o, FeatureVector[] initialMeans, int terminateIn)
        {
            //k means
            FeatureVector[] means = initialMeans;
            int k = initialMeans.Length;

            FeatureVector mean = new FeatureVector(markov_t, markov_t);
            List<FeatureVector>[] clusters = new List<FeatureVector>[k];
            double min = 0;

            bool notDone = true;

           
            Console.WriteLine("\n ..Clustering with " + k + " specified cluster centers");


            //k means Algorithm
            for (int q = 0; q < terminateIn; q++)
            {
                Console.WriteLine("Iteration: " + q);
                notDone = false;

                //put each element in o[] to the centroid that its closer to
                for (int i = 0; i < o[0].vec.Length; i++)
                {
                    min = cosine(o[i], means[0]);

                    for (int kk = 0; kk < k; kk++)
                        min = Math.Min(cosine(o[i], means[kk]), min);

                    for (int kk = 0; kk < k; kk++)
                        if (cosine(o[i], means[kk]) == min)
                        {
                            //check if o[i] is not already in the proper cluster
                            if (!clusters[kk].Contains(o[i]))
                            {
                                //o[i] has found a new cluster, move on
                                Console.WriteLine("Feature Vector " + i + " has moved to cluster " + kk);
                                clusters[kk].Add(o[i]);

                                //remove o[i] from other clusters
                                for (int kkk = 0; kkk < kk; kkk++)
                                    if (kkk != kk)
                                        clusters[kkk].Remove(o[i]);
                            
                                //calculate the new centroid after inserting the new member
                                mean = centroid(clusters[kk]);
                                //make sure everything checks out
                                for (int ui = 0; ui < mean.vec.Length; ui++)
                                    if (!(mean.vec[ui] == means[kk].vec[ui])) notDone = true;
                                means[kk] = mean;

                                continue;
                            }
                        }
                }

                for (int s = 0; s < k; s++)
                {
                    //Output..
                    Console.WriteLine("= Elements in cluster " + s + " are:");
                    for (int ss = 0; ss < clusters[k].Count(); ss++)
                        Console.WriteLine(ss + ". Feature Vector " + clusters[k][ss].id);
                    Console.WriteLine("= Centroid for cluster " + s + " is: \n <");

                    for (int ss = 0; ss < means[k].vec.Length; ss++)
                        Console.Write(means[k].vec[ss] + ", ");
                    Console.Write("> \n");
                }

                if (!notDone)
                {
                    Console.WriteLine("K means algorithm has been completed after " + q + " steps.");
                    return clusters;
                }
            }

            Console.WriteLine("Did not converge in specified steps: Terminating.");

            return clusters;

        }

        /// <summary>
        /// Pseudo-randomly shuffles a feature cluster with divide and conquer
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        static FeatureVector[] randomShuffle(FeatureVector[] a)
        {
            int[] u = new int[a.Length];
            int n = 0;
            FeatureVector o = new FeatureVector((markov_t*markov_t)+5, -1);
            for (int i = 0; i < u.Length; i++)
                u[i] = i; //populate
            randomSelect(u, ref n, ref u);

            for (int i = 0; i < u.Length; i++)
            {
                o = a[i];
                a[i] = a[u[i]]; // rawr
                a[u[i]] = o;
            }

            return a;
        }


        private static void randomSelect(int[] a, ref int n, ref int[] u)
        {
           
            if (a.Length == 1)
            {
                u[n] = a[0];
                n++;
                return;
            }

            Random r = new Random();
            int pivot = r.Next(1, a.Length - 1);
            u[n] = a[pivot];
            n++;

            int[] c = new int[pivot];
            int[] d = new int[a.Length - pivot - 1];
            Array.Copy(a, c, pivot);
            Array.Copy(a, pivot + 1, d, 0, a.Length - pivot - 1);
     
            if (c.Length > 0)
                randomSelect(c, ref n, ref u);
            if (d.Length > 0)
                randomSelect(d, ref n, ref u);
            
            return; 
        }

        
        /// <summary>
        /// 
        /// Returns a state transition tensor (multidimensional array) with
        /// the markov features in the horizontal, vertical, and main diagonal directions
        /// as well as markov features in the time direction
        /// 
        /// 
        /// </summary>
        /// <param name="f"> The output of MGroup(), which is a Rectangular Prism Y(width, height, depth) : Y(FrameWidth, FrameHeight, 4) 
        /// representing 3d DCT coefficients for the current group</param>
        /// <returns>A transition tensor (double[,,) where the x and y positions are average intra-frame markov features, and z position is inter-frame markov feature</returns>
        static double[,,] transition_Tensor(double[,,] f)
        {
            int t1 = markov_t;
            int t2 = markov_t;

            //For a 4x4x4 cube of 3D DCT frequencies f[x,y,z]:

            //
            //
            //
            //      EXPLANATION
            //
            //
            //
            //      Position [a,b,c] of the tensor represents:
            //      a is the row position of the transition tensor
            //      b is the column position of the transition tensor
            //      c is the z position in the transition tensor (inter frame)
            //
            //      Where kron(*) casts the result of a binary statement * to a bit {0,1}
            //
            ///     
            ///     Implementation:
            /// 
            /// [1] STEGANALYSIS OF YOUTUBE COMPRESSED VIDEO USING HIGH-ORDER STATISTICS
            /// IN 3D DCT DOMAIN
            /// Hong Zhao1, 2, Hongxia Wang1, Hafiz
            ///
            /// Info:
            /// 
            /// [2] JPEG Image Steganalysis Utilizing both Intrablock
            /// and Interblock Correlations
            /// 

            double[,,] markov = new double[t1, t2, 2];
            double horizontal = 0;
            double vertical = 0;
            double diagonal = 0;



                    for (int m = 0; m < t2; m++)
                    {
                        for (int n = 0; n < t1; n++)
                        {
                            for (int i = 0; i < f.GetLength(1); i++)
                            {
                                for (int j = 0; j < f.GetLength(0); j++)
                                {

                                    if (j < f.GetLength(0) - 1)
                                        horizontal += ((int)f[j, i, 3] == m ? (((int)f[j + 1, i, 3] == n && (int)f[j, i, 3] == m) ? 1.0 : 0.0) : 0.0);

                                    if (i < f.GetLength(1) - 1)
                                        vertical += ((int)f[j, i, 3] == m ? (((int)f[j, i + 1, 3] == n && (int)f[j, i, 3] == m) ? 1.0 : 0.0) : 0.0);

                                    if (j < f.GetLength(0) - 1 && i < f.GetLength(1) - 1)
                                        diagonal += ((int)f[j, i, 3] == m ? (((int)f[j + 1, i + 1, 3] == n && (int)f[j, i, 3] == m) ? 1.0 : 0.0) : 0.0);
                                    
                                    //inter frame markov feature
                                    markov[m, n, 1] += ((int)f[j, i, 3] == m ? (((int)f[j, i, 2] == n && (int)f[j, i, 3] == m) ? 1.0 : 0.0) : 0.0);   
                                
                                }
                            }
                            // intra frame markov feature
                                    markov[m, n, 0] = (horizontal + vertical + diagonal / 3.0);

                                    horizontal = 0;
                                    vertical = 0;
                                    diagonal = 0;
                                    
                        }
                    }

            return markov;
        }

        /*
        /// <summary>
        /// Helper 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        static double[] _f(Frequency_4x4x4 f, int i)
        {
            double[] d = { f.frequencies[0, 0, i], f.frequencies[0, 1, i], f.frequencies[0, 2, i], f.frequencies[0, 3, i], f.frequencies[1, 0, i], f.frequencies[1, 1, i], f.frequencies[1, 2, i], f.frequencies[1, 3, i], f.frequencies[2, 0, i], f.frequencies[2, 1, i], f.frequencies[2, 2, i], f.frequencies[2, 3, i], f.frequencies[3, 0, i], f.frequencies[3, 1, i], f.frequencies[3, 2, i], f.frequencies[3, 3, i] };
            return d; 
        }



        /// <summary>
        /// Helper 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="num"></param>
        /// <returns></returns>
        static Frequency_4x4x4 f_(Frequency_4x4x4 f,  double[] f0, int i)
        {
            Frequency_4x4x4 f1 = new Frequency_4x4x4();
            Array.Copy( f.frequencies, f1.frequencies, f1.frequencies.Length);

            f1.frequencies[0, 0, i] = f.frequencies[0, 0, i];
            f1.frequencies[0, 1, i] = f.frequencies[0, 1, i];
            f1.frequencies[0, 2, i] = f.frequencies[0, 2, i];
            f1.frequencies[0, 3, i] = f.frequencies[0, 3, i];
            f1.frequencies[1, 0, i] = f.frequencies[1, 0, i];
            f1.frequencies[1, 1, i] = f.frequencies[1, 1, i];
            f1.frequencies[1, 2, i] = f.frequencies[1, 2, i];
            f1.frequencies[1, 3, i] = f.frequencies[1, 3, i];
            f1.frequencies[2, 0, i] = f.frequencies[2, 0, i];
            f1.frequencies[2, 1, i] = f.frequencies[2, 1, i];
            f1.frequencies[2, 2, i] = f.frequencies[2, 2, i];
            f1.frequencies[2, 3, i] = f.frequencies[2, 3, i];
            f1.frequencies[3, 0, i] = f.frequencies[3, 0, i];
            f1.frequencies[3, 1, i] = f.frequencies[3, 1, i];
            f1.frequencies[3, 2, i] = f.frequencies[3, 2, i];
            f1.frequencies[3, 3, i] = f.frequencies[3, 3, i];

            return f1;
        }
         * 
         * 
        */


        static void Main(string[] args)
        {
            byte[] videoData = new byte[0];
            byte[] bytes = new byte[0];
            Bitmap b = new Bitmap(1, 1);
            Bitmap b1 = new Bitmap(1, 1);
            Bitmap b2 = new Bitmap(1, 1);
            Bitmap b3 = new Bitmap(1, 1);
            
            string savefilename = "file";
            byte[] saveFile = new byte[1];
            Frequency_4x4x4 steg = new Frequency_4x4x4();

            Console.WriteLine("XMAS PACKAGER -- (0x00)");
            Console.WriteLine("------------");
            Console.WriteLine("Digital Data Hiding and Forensic Steganalysis for Online Videos");
            Console.WriteLine("What is the video width?");
            int videoWidth = Int32.Parse(Console.ReadLine());
            Console.WriteLine("What is the video height?");
            int videoHeight = Int32.Parse(Console.ReadLine());
            Console.WriteLine("What is the number of frames of the video?");
            int frameLength = Int32.Parse(Console.ReadLine());
            Console.WriteLine("What is the video's filename?");
            string filename = Console.ReadLine();
            Console.WriteLine("What is the stegano file's filename?");
            string sfilename = Console.ReadLine();
            Console.WriteLine("Are you hiding data, reconstructing it, or performing steganalysis? <H/R/S>");
            string hiding = Console.ReadLine();

            
            
            if (hiding == "S")
            {

                Console.WriteLine("\n > == Blind Steganalysis using the 3D DCT's high frequencies and K-means Clustering ==");
                Console.WriteLine("\n 1. Specify Cluster Centroids (Recommended) \n 2. CPU-Generated Random Cluster Centroids \n ");
                string cent = Console.ReadLine();
                if (cent == "2")
                {
                    Console.WriteLine("The CPU will generate k random cluster centroids. This may or may not give the best output. \nWhat is k?");
                    int k = Int32.Parse(Console.ReadLine());
                    
                }
                if (cent == "1")
                {
                    Console.WriteLine("To initialize the k means algorithm, you need to select two representative videos so that the program has an understanding of the difference between the two classes. "
                    + "\n\n The first video should have data hidden inside of all frames (Stegano), while the second should have no data hidden inside of any frame. (Clean) \nWhat is k?");
                    int k = Int32.Parse(Console.ReadLine());
                    Console.WriteLine("Specify the filename of the stegano video");
                    string steganovid = Console.ReadLine();
                    Console.WriteLine("Specify the filename of the clean video");
                    string cleanvid = Console.ReadLine();
                    Console.WriteLine("Choose a parameter T for (Typically 0-4)");
                    int T = Int32.Parse(Console.ReadLine());

                    Bitmap[] sData = new Bitmap[4];
                    Bitmap[] cData = new Bitmap[4];
                    
                    //
                    // An explanation of how this all works with apologies for false rates
                    //
                    // With k = 2, we harvest the first four frames of the videos
                    // and group them together into two separate classes (or clusters)
                    // To create two ideal classes:
                    //
                    //
                    //
                    // Ideal Dirty Video
                    // Ideal Clean Video
                    //
                    //
                    //
                    // We then check the similarity/distance between incoming group and existing classes in k-means.
                    // We group by high similarity / low distance
                    // Cosine Distance = (1 - Cosine similarity)
                    //
                    //
                    //
                    // A group contains 4 frames.
                    // Each class is initially represented by a group of first 4 frames from a separate video.
                    //
                    // Really and truly, the "Dirty" and "Clean" are meaningless
                    // Watch, for example, we use k = 2
                    // and set each cluster representative to be a "Ideally Clean" (strong correlations) and "Ideally Dirty" (weak correlations) video
                    // We imagine a Clean video has strong spatial and temporal correlations
                    // and a Dirty video has weak correlations of both.
                    //
                    //
                    // But our attacker's crafty data hiding algorithm might be different from our own
                    // and it might preserve either temporal and spatial correlations or both.
                    //
                    // This would mean that not all dirty videos may exhibit weak correlations.
                    //
                    // We could use larger values of k
                    // To create classes where unique steganography only modifies certain correlations
                    //
                    // All this algorithm is doing is measuring if an input video group has strong and weak correlations
                    // which it captures using intra-frame (spatial) markov features, and inter-frame (temporal) markov features
                    // Which only capture the unique probabilities of transitioning from DCT coefficient to DCT coefficient.
                    //
                    //
                    // The algorithm classifies due to weak spatial and temporal correlations in DCT domain.
                    // It checks to see how similar DCT frequencies are to each other.
                    //
                    // BUT
                    //
                    // One could set k to anything
                    // And have
                    //
                    // Strong Temporal Correlations with Weak Spatial Correlations,
                    // Strong Spatial Correlations with Weak Temporal Correlations, 
                    // ...
                    // etc
                    //
                    // (any number of classes)
                    //
                    // and other features are exhibited from YouTube videos, such as the Absolute Central Moment, Kurtosis, Skewness, etc
                    // which introduce the possibility for even more classes
                    //
                    //
                    //
                    // Because Steganography algorithms are very different, some create weak spatial correlations, some, weak temporal correlations.
                    // There are a seemingly limitless number of ways to hide data in videos which all have unique footprints.
                    //
                    // The steganalysis scheme then
                    // but with contiguous spacing (scaled by T) to compensate 
                    // for temporal error correcting codes.
                    //
                    // Why?
                    //
                    // Data hiding in online video is unreliable without "temporal error correction",
                    // This means we repeat the message in multiple neighboring frames to improve the reliability of data hiding.
                    //
                    // YouTube compression can be seen as a form of steganographic attack called a Collusion Attack
                    // A Collusion Attack on a video weakens a hidden message by altering the frames around it and creating weak temporal correlation.
                    // 
                    // This only means that many of our frames in time series become too different from each other, 
                    // which would convince us that there is hidden data in the video using this detection scheme.
                    //
                    // In reality YouTube messes up our hidden data by introducing errors in it and around it by dropping frames and encoding them
                    // with motion vector data, which is why it is necessary to use temporal error correction.
                    //
                    // But with temporal error correction, if we try to perform steganalysis on neighboring frames, 
                    // we may get false negatives thanks to error correcting codes.
                    //
                    // With the strong temporal correlation, or, the very close similarity between frames, 
                    // temporal error correcting codes fool us into thinking that nothing is there.
                    // Alternatively, we may get true positives thanks to noisy collusion attack
                    // but this is too difficult to determine in advance what we will get with YouTube.
                    // All it will likely do is interfere with our detection scheme
                    //
                    //
                    // Temporal error correction lowers the bitrate in the temporal direction if it repeats the message in time
                    // This creates strong temporal correlations between DCT blocks.
                    // With the model in [1] this would give us a false negative.
                    // This algorithm classifies using WEAK temporal correlations between DCT blocks
                    // But the temporal corrected message has STRONG temporal correction.
                    //
                    //
                    // We just want to get an honest detection reading that isn't affected by error correction or collusion.
                    //
                    //
                    // Likewise, fast moving, choppy, and noisy videos 
                    // give many false positives for this reason because they have weak temporal correlations.
                    // and the authors mention this.
                    //
                    // YouTube does not collude videos which have no information as long as they are smooth and slow-moving,
                    // and they usually remain this way after compression
                    //
                    // Slow moving videos with hidden data that isn't temporally corrected are colluded, 
                    // but it should not negatively affect the detection performance.
                    //
                    // It is the error correction which is annoying
                    // because it reduces bitrate
                    // a very widely held fact of Forensic Steganalysis is that bitrate reduction can harm Steganalysis
                    // because hiding less information lowers detection rates :(
                    //
                    //
                    Console.WriteLine("A. Initializing Steganalytic Tool (8 Steps)");

                    Console.WriteLine("1. Loading Data");
                    for (int i = 0; i < 4; i++)
                        sData[i] = BitmapFromArray1DSafe(GetBytesFromFile(steganovid, (((videoWidth * videoHeight) * i * T) + (((videoWidth * videoHeight) * i * T) / 2)), videoWidth * videoHeight), videoWidth, videoHeight);;

                    for (int i = 0; i < 4; i++)
                        cData[i] = BitmapFromArray1DSafe(GetBytesFromFile(steganovid, (((videoWidth * videoHeight) * i * T) + (((videoWidth * videoHeight) * i * T) / 2)), videoWidth * videoHeight), videoWidth, videoHeight); ;

                    Console.WriteLine("2. Performing 3D DCT Transform, Acquiring Groups Y(k1,k2,k3,k) for Initial Centroids, calculating Low-Order Statistics");
                    double[,,] sDataGroup = M_Group(sData[0], sData[1], sData[2], sData[3]);
                    double[,,] cDataGroup = M_Group(cData[0], cData[1], cData[2], cData[3]);

                    Console.WriteLine("2. Assembling transform histograms"); ;
                    double smean = m(sDataGroup); Console.WriteLine("Mean of stegano video is" + smean);
                    double cmean = m(cDataGroup); Console.WriteLine("Mean of clean video is " + cmean);


                    double sstd = std(sDataGroup, smean); Console.WriteLine("Sigma of stegano video is " + sstd);
                    double cstd = std(cDataGroup, cmean); Console.WriteLine("Sigma of clean video is " + cstd);


                    Console.WriteLine("\n3. Calculating Kurtosis for Initial Centroids");
                    double sKurtosis = Kurtosis(smean, sstd, flatten(sDataGroup));
                    double cKurtosis = Kurtosis(cmean, cstd, flatten(cDataGroup));
                    Console.WriteLine("Kurtosis for Clean Videos: " + cKurtosis);
                    Console.WriteLine("Kurtosis for Dirty Videos: " + sKurtosis);

                    Console.WriteLine("\n4. Calculating Skewness for Initial Centroids");
                    double sSkewness = Skewness(smean, sstd, flatten(sDataGroup));
                    double cSkewness = Skewness(cmean, cstd, flatten(cDataGroup));
                    Console.WriteLine("Skewness for Clean Videos: " + cSkewness);
                    Console.WriteLine("Skewness for Dirty Videos: " + sSkewness);

                    Console.WriteLine("\n5. Calculating 1st order Absolute Central Moment for Initial Centroids");
                    double s1stabs = AbsCentralMoment(1, videoHeight, videoWidth, flatten(sDataGroup), smean);
                    double c1stabs = AbsCentralMoment(1, videoHeight, videoWidth, flatten(cDataGroup), cmean);
                    Console.WriteLine("1st order Absolute Moment for Clean Videos: " + c1stabs);
                    Console.WriteLine("1st order Absolute Moment for Dirty Videos: " + s1stabs);

                    Console.WriteLine("\n6. Calculating 2nd order Absolute Central Moment for Initial Centroids");
                    double s2ndabs = AbsCentralMoment(2, videoHeight, videoWidth, flatten(sDataGroup), smean);
                    double c2ndabs = AbsCentralMoment(2, videoHeight, videoWidth, flatten(cDataGroup), cmean);
                    Console.WriteLine("2nd order Absolute Moment for Clean Videos: " + c2ndabs);
                    Console.WriteLine("2nd order Absolute Moment for Dirty Videos: " + s2ndabs);

                    Console.WriteLine("\n7. Calculating 3rd order Absolute Central Moment for Initial Centroids");
                    double s3rdabs = AbsCentralMoment(3, videoHeight, videoWidth, flatten(sDataGroup), smean);
                    double c3rdabs = AbsCentralMoment(3, videoHeight, videoWidth, flatten(cDataGroup), cmean);
                    Console.WriteLine("3rd order Absolute Moment for Clean Videos: " + c3rdabs);
                    Console.WriteLine("3rd order Absolute Moment for Dirty Videos: " + s3rdabs);

                    Console.WriteLine("\n8. Calculating Markov Features for Initial Centroids");
                    Markov smv = new Markov();
                    Markov cmv = new Markov();

                    //Reduce feature dimension by T, which is 10
                    sDataGroup = Truncate(sDataGroup);
                    cDataGroup = Truncate(cDataGroup);

                    smv.markov = transition_Tensor(sDataGroup); 
                    cmv.markov = transition_Tensor(cDataGroup);

                    Console.WriteLine("\n8. Assembling essential feature Vectors for initializing 2-means");
                    FeatureVector fv1 = featureVector1(smv, s1stabs, s2ndabs, s3rdabs, sSkewness, sKurtosis, -1);
                    FeatureVector fv2 = featureVector1(cmv, c1stabs, c2ndabs, c3rdabs, cSkewness, cKurtosis, -1);

                    Console.WriteLine("Done Initializing! \n\n");

                    Console.WriteLine("B. Gathering Low-order and High-order statistics from target video frames.");
                    


                    Console.ReadLine();
                }

                 
                for (int i = 0; i < frameLength - 3; i++)
                {


                    Console.Write("\n Scanning the video with " + (100 * (((double)i / (double)frameLength))) + "% Completed \n");

                    videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                    b = BitmapFromArray1D(videoData, videoWidth, videoHeight);

                    videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * (i+1)) + (((videoWidth * videoHeight) * (i+1)) / 2)), videoWidth * videoHeight);
                    b1 = BitmapFromArray1D(videoData, videoWidth, videoHeight);

                    videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * (i+2)) + (((videoWidth * videoHeight) * (i+2)) / 2)), videoWidth * videoHeight);
                    b2 = BitmapFromArray1D(videoData, videoWidth, videoHeight);

                    videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * (i+3)) + (((videoWidth * videoHeight) * (i+3)) / 2)), videoWidth * videoHeight);
                    b3 = BitmapFromArray1D(videoData, videoWidth, videoHeight);

                    steg = DCT3D_4x4x4(b, b1, b2, b3); 

                    videoData = Array1DFromBitmap(b, videoWidth * videoHeight);
                    SaveBytesToFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight, videoData);

                }
                  
            }

            if (hiding == "R")
            {
                Console.WriteLine("Enter the filename to save the reconstructed file to");
                savefilename = Console.ReadLine();
                Console.WriteLine("What is the filesize of the file in bytes? If unknown, type 42");
                filesizespecified = long.Parse(Console.ReadLine());
            }

            Console.WriteLine("What algorithm will you use? <Choose Number> \n \n \n 1. (checkers): Spatial Domain PVD - Very Robust, uses computationally efficient Pixel Value Differencing and Histogram Shifting algorithms. Easily detectable visually however. Payload capacity is large. Reccomended if video is being uploaded to YouTube or highly compressed on disk. Not recommended for covert operation and is weak to known steganalysis attacks.");
            Console.WriteLine("\n \n \n 2. (smiling face): Transform Domain DCT2D 8-D Vector Quantization - Fragile, with borderline robustness, uses computationally intensive Discrete Cosine Transform and eight dimensional vector quantization. Invisible with distortion < 3db (6% visual distortion calculated with PSNR). Payload capacity is large but small once error correction used. Highly recommended if video is stored on hard disk and will not be re-encoded. Recommended for covert operation as embedded data is invisible.");
            string response = Console.ReadLine();
            
            while (true)
            {
                if (response == "2")
                {
                    Console.WriteLine("What is your thresholding parameter T?");
                    int T = Int32.Parse(Console.ReadLine());
                    Console.WriteLine("What is the error correcting codeword length you'll use? (Set this to 900 or greater for YouTube.)");
                    int rep = Int32.Parse(Console.ReadLine());
                    bytes = GetBytesFromFile(sfilename, 0);
                    acc = new double[8 * bytes.Length];

                    finalFile = new byte[bytes.Length];

                    if (bytes.Length * 8 * rep > videoWidth * videoHeight / 4 / 4)
                    {
                        Console.WriteLine(" * File is much larger than embedding capacity of frame with chosen error code, would you like to embed it across all frames? <Y/n>");
                        string allframes = Console.ReadLine();

                        if (allframes != "Y")

                        {
                            Console.WriteLine("File too large, aborting...");
                        }
                        else
                        {
                            loop = false;
                        }
                    }
                    Frequency_4x4 f = new Frequency_4x4();



                    for (int i = 0; i < frameLength; i++)
                    {
                        if (hiding == "H")
                        {
                            
                            Console.Write("\n Processing.." + (100 * (((double)i / (double)frameLength))) + "% Completed \n");
                            
                            videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                            b = BitmapFromArray1DSafe(videoData, videoWidth, videoHeight);
                            b = transform_Embed(b, T, bytes, rep);

                            videoData = Array1DFromBitmapSafe(b, videoWidth * videoHeight);
                            SaveBytesToFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight, videoData);

                        }
                        if (hiding == "R")
                        {

                                if (bytes.Length * rep >= (videoWidth * videoHeight / 8 / 16))
                                    loop = false;

                                filesizespecified = bytes.Length;
                             
                            Console.Write("\n Processing.." + (100 * (((double)i / (double)frameLength))) + "% Completed \n");
                            
                            videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                            b = BitmapFromArray1DSafe(videoData, videoWidth, videoHeight);
                            saveFile = transform_Retrieve(b, bytes, T, rep);
                            
                            if (!loop)
                                Console.WriteLine("BER for chunk " + i + " (Bit Error Rate) is " + BER(bytes, saveFile, (int)(pos), (int)(videoWidth * videoHeight / 16 / rep)));
                            if (loop)
                                Console.WriteLine("BER (Bit Error Rate) is " + BER(bytes, saveFile));

                            SaveBytesToFile(savefilename, 0, bytes.Length, saveFile);
                            if (loop)
                                acc = AccArray(saveFile, acc);

                        }

                    }
                    acc = AvgArray(acc, frameLength);
                    Console.WriteLine("P value is " + P(acc));
                    Console.WriteLine("Done!");
                    Console.ReadLine();
                }
                if (response == "1")
                {
                    Console.WriteLine("What is your thresholding parameter T?");
                    int T = Int32.Parse(Console.ReadLine());
                    Console.WriteLine("What is your thresholding parameter G?");
                    int G = Int32.Parse(Console.ReadLine());
                    Console.WriteLine("What is the block size");
                    int block_size = Int32.Parse(Console.ReadLine());
                    Console.WriteLine("What is the spatial error correcting codeword length you'll use? (Set this to 50 or greater for YouTube.)");
                    int rep = Int32.Parse(Console.ReadLine());
                    Console.WriteLine("What is the temporal error correcting codeword length you'll use? (Set this to 50 or greater for YouTube.)");
                    int trep = Int32.Parse(Console.ReadLine());
                    bytes = GetBytesFromFile(sfilename, 0);
                    finalFile = new byte[bytes.Length];
                    acc = new double[8 * bytes.Length];

                    long loc = 0;
                    if (bytes.Length * 8 * rep > videoWidth * videoHeight / block_size / block_size)
                    {
                        Console.WriteLine(" * File is much larger than embedding capacity of frame with chosen error code, would you like to embed it across all frames? <Y/n>");
                        string allframes = Console.ReadLine();

                        if (allframes != "Y")
                        {
                            Console.WriteLine("File too large, aborting...");
                        
                        }

                        else
                        {
                            loop = false;
                        }

                    }

                    for (int i = 0; i < frameLength; i++)
                    {

                        if (hiding == "H")
                        {
                            
                            Console.Write("\n Processing.." + (100 * (((double)i / (double)frameLength))) + "% Completed \n");
                            
                            videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                            b = BitmapFromArray1D(videoData, videoWidth, videoHeight);
                            b = spatial_Embed(b, block_size, T, G, bytes, rep, trep);

                            videoData = Array1DFromBitmap(b, videoWidth * videoHeight);
                            SaveBytesToFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight, videoData);

                        }
                        if (hiding == "R")
                        {

                                if (bytes.Length * rep >= (videoWidth * videoHeight / 8 / block_size / block_size))
                                    loop = false;

                                filesizespecified = bytes.Length;
                            
                            
                            Console.Write("\n Processing.." + (100 * (((double)i / (double)frameLength))) + "% Completed.");
                            loc = pos;     
                            videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                            b = BitmapFromArray1D(videoData, videoWidth, videoHeight);
                            saveFile = spatial_Retrieve(b, block_size, T, G, bytes, rep, trep);

                            if (!loop)
                                if ((i+1) % (trep) == 0)
                                {
                                    Console.WriteLine("BER for chunk " + i + " (Bit Error Rate) is " + BER(bytes, saveFile, (int)(loc), (int)(pos - loc)));
                                    Array.Copy(saveFile, pos, finalFile, pos, pos - loc); //this is the cumulative file
                                    Console.WriteLine((100 * (((double)i / (double)frameLength))) + " % of the file ( " + pos + " ) bytes : total filesize of ( " + bytes.Length + " ) bytes has been saved to " + savefilename);
                                    SaveBytesToFile(savefilename, 0, bytes.Length, finalFile); //save the whole file
                                }
                            
                            if (loop)
                                Console.WriteLine("BER (Bit Error Rate) is " + BER(bytes, saveFile));
                            if (loop)
                                SaveBytesToFile(savefilename, 0, bytes.Length, saveFile); //save this
                            if (loop)
                                acc = AccArray(saveFile, acc);
                          
                        }

                    }
                    acc = AvgArray(acc, frameLength);
                    Console.WriteLine("P value is " + P(acc));
                    Console.WriteLine("SUCCESS! Completed 100%. Please be careful with this video.");
                    Console.ReadLine();
                }
                Console.WriteLine("Try a different response.");
                Console.ReadLine();
            }
        }
    }
}

/*
            Cited

        [1] : Hong Zhao; Hongxia Wang; Malik, H., 
 * 
 * 
 * "Steganalysis of Youtube Compressed Video Using High-order Statistics in 3D DCT Domain," 
 *      in Intelligent Information Hiding and Multimedia Signal Processing (IIH-MSP), 
 *      2012 Eighth International Conference on , vol., no., pp.191-194, 18-20 July 2012
        doi: 10.1109/IIH-MSP.2012.52
 
 * 
 * 
 * 
 * keywords: {
 * Markov processes;
 * discrete cosine transforms;
 * image classification;
 * pattern clustering;
 * social networking (online);
 * spatiotemporal phenomena;
 * statistical analysis;
 * steganography;
 * video coding;
 * 3D DCT domain;
 * 3D discrete cosine transform domain;
 * Markov features;
 * YouTube compressed video steganalysis;
 * absolute central moments;
 * cover-videos;
 * data hiding;
 * hidden message detection;
 * high-order statistics classification;
 * kurtosis;
 * skewness;
 * spatial-temporal correlation;
 * stego-videos;
 * unsupervised k-means clustering;
 * video frames;
 * Correlation;
 * Discrete cosine transforms;
 * Feature extraction;
 * Image coding;
 * Markov processes;
 * Video compression;
 * YouTube;
 * 3D-DCT;
 * Video Steganalysis},

 
 */
