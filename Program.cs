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
 * Wu! It works, I told you it would work! -RC
 */

namespace H264Stego
{
    class Program
    {
        static long filesizespecified;
        static bool loop = true;
        static long pos;
        static long pos_bit;
        static long pos_rep;

        public class Frequency_4x4
        {
            public double[,] frequencies;
            
            public Frequency_4x4()
            {
                frequencies = new double[4, 4];
            }
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
                            summand += (b.GetPixel(j, i).ToArgb() * Math.Cos(Math.PI*y*((2*i)+1)/8) * Math.Cos(Math.PI*x*((2*j)+1)/8));
                        }
                    }
                    summand = 0.5 * summand;
                    if (x == 0) summand *= (1.0 / Math.Sqrt(2));
                    if (y == 0) summand *= (1.0 / Math.Sqrt(2));
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
                            temp = f.frequencies[j, i] * Math.Cos(Math.PI*i*((2*y)+1)/8)*Math.Cos(Math.PI*j*((2*x)+1)/8);
                            if (i == 0) temp *= (1.0 / Math.Sqrt(2));
                            if (j == 0) temp *= (1.0 / Math.Sqrt(2));
                            //summand *= 0.5;
                            summand += temp;
                        }
                    }
                    summand = 0.5 * summand;

                    b.SetPixel(x, y, Color.FromArgb((int)Math.Round(summand)));
                    summand = 0;
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
            double distance_prime = Math.Round(distance/(double)T);
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
            double distance_prime = distance / T;
            distance_prime = distance_prime - Math.Round(distance_prime);
            if (bit >= 0)
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
                    if ((((i+1) % 2) == ((j+1) % 2)) && alpha > T) block[j, i] += gamma;
                    if ((((i+1) % 2) != ((j+1) % 2)) && alpha < -T) block[j, i] += gamma;                
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
            Console.WriteLine(" >> ERROR.. assuming salt");
            return 1;
        }

        static double[,] Array2DFromBitmap(Bitmap b)
        {
            double[,] a = new double[b.Width, b.Height];
            for (int i = 0; i < b.Height; i++)
                for (int j = 0; j < b.Width; j++)
                    a[j,i] = (b.GetPixel(j,i).ToArgb());
            return a;
        }

        static Bitmap BitmapFromArray2D(double[,] a)
        {
            Bitmap b = new Bitmap(a.GetLength(0),a.GetLength(1));
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
        static Bitmap spatial_Embed(Bitmap b, int block_size, int T, int G, byte[] bytes, int rep)
        {
            double[,] block = new double[block_size, block_size];
            double gamma = Gamma(block, T, G);
            double gamma_prime = Gamma_Prime(block, T, G);
            int[,] sgn = matrix_SGN(block);

            if (loop) pos = 0;

            if (loop) pos_bit = 0;

            if (loop) pos_rep = 0;

            Bitmap src = new Bitmap(block_size, block_size);
            int x = 0;
            int y = 0;

            for (; ; pos++)
            {

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
                        if (y >= b.Height) return b;

                    }

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

                        int foo = src.GetPixel(0, 0).ToArgb();
                        foo = 0;

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
        static byte[] spatial_Retrieve(Bitmap b, int block_size, int T, int G, byte[] bytes, int rep)
        {
            double[,] block = new double[block_size, block_size];
            double gamma = Gamma(block, T, G);
            double gamma_prime = Gamma_Prime(block, T, G);

            Bitmap src = new Bitmap(block_size, block_size);
            int x = 0;
            int y = 0;
            byte bit = 0;
            byte avg = 0;

            if (loop) pos = 0;

            if (loop) pos_bit = 0;

            if (loop) pos_rep = 0;

            byte[] newbytes = new byte[bytes.Length];

            for (; ; pos++)
            {
                if (pos == bytes.Length) pos = 0;

                for (pos_bit = 0; pos_bit < 8; pos_bit++)
                {
                    for (pos_rep = 0; pos_rep < rep; pos_rep++)
                    {
                        src = b.Clone(new Rectangle(x, y, block_size, block_size), b.PixelFormat);
                        block = Array2DFromBitmap(src);

                        int alpha = arithmetic_Difference(block, matrix_SGN(block));
                        //if (alpha <= ((2*T) + G) && alpha >= -((2*T)+ G)) //see my frustration above at Alavianmehr et. als' "original" work
                        //{
                        avg+= arithmetic_Retrieve(block, T, G);


                        //}
                        //if (alpha > ((2 * T) + G) || alpha < -((2 * T) + G)) k--;

                        x += block_size;
                        if (x >= b.Width) { y += block_size; x = 0; }
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
                    if (((a[i] >> k) % 2) != ((b[i] >> k) % 2)) ber++;

                        return 100.0*ber/a.Length/8;

        }

        static void Main(string[] args)
        {
            byte[] videoData = new byte[0];
            byte[] bytes = new byte[0];
            Bitmap b = new Bitmap(1,1);
            string savefilename = "file";
            byte[] saveFile = new byte[1];
            Console.WriteLine("XMAS PACKAGER-- (0x00)");
                Console.WriteLine("------------");
                Console.WriteLine("Silent Steganography on YouTube while preserving Peak Video Signal : Noise");
                Console.WriteLine("What is the video width?");
                int videoWidth = Int32.Parse(Console.ReadLine());
                Console.WriteLine("What is the video height?");
                int videoHeight = Int32.Parse(Console.ReadLine());
                Console.WriteLine("What are the number of frames of the video?");
                int frameLength = Int32.Parse(Console.ReadLine());
                Console.WriteLine("What is the video's filename?");
                string filename = Console.ReadLine();
                Console.WriteLine("What is the stegano file's filename?");
                string sfilename = Console.ReadLine();
                Console.WriteLine("Are you hiding data or reconstructing it? <H/R>");
                string hiding = Console.ReadLine();
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
                        if (bytes.Length * 8 * rep > videoWidth * videoHeight / 4 / 4)
                        {
                            Console.WriteLine(" * File is much larger than embedding capacity of frame with chosen error code, would you like to embed it across all frames? <Y/n>");
                            string allframes = Console.ReadLine();
                            if (allframes != "Y")
                            {
                                Console.WriteLine("File too large, aborting...");
                                loop = false;
                            }
                        }
                        Frequency_4x4 f = new Frequency_4x4();
                        


                        for (int i = 0; i < frameLength; i++)
                        {
                            if (hiding == "H")
                            {
                                Console.Write("\n Processing.. " + i / frameLength + "% completed");

                                videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                                b = BitmapFromArray1DSafe(videoData, videoWidth, videoHeight);
                                b = transform_Embed(b, T, bytes, rep);

                                videoData = Array1DFromBitmapSafe(b, videoWidth * videoHeight);
                                SaveBytesToFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight, videoData);

                            }
                            if (hiding == "R")
                            {
                                if (filesizespecified == 42)
                                {
                                    if (bytes.Length * rep >= (videoWidth * videoHeight / 8 / 16))
                                        loop = false;

                                    filesizespecified = bytes.Length;
                                }

                                Console.Write("\n Processing.." + i / frameLength + "% Done \n");

                                videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                                b = BitmapFromArray1D(videoData, videoWidth, videoHeight);
                                saveFile = transform_Retrieve(b, bytes, T, rep);
                                Console.WriteLine("BER (Bit Error Rate) is " + BER(bytes, saveFile));
                                SaveBytesToFile(savefilename, 0, bytes.Length, saveFile);
                            }
                        }
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
                        Console.WriteLine("What is the error correcting codeword length you'll use? (Set this to 50 or greater for YouTube.)");
                        int rep = Int32.Parse(Console.ReadLine());
                        bytes = GetBytesFromFile(sfilename, 0);
                        if (bytes.Length * 8 * rep > videoWidth * videoHeight / block_size / block_size)
                        {
                            Console.WriteLine(" * File is much larger than embedding capacity of frame with chosen error code, would you like to embed it across all frames? <Y/n>");
                            string allframes = Console.ReadLine();
                            if (allframes != "Y")
                                Console.WriteLine("File too large, aborting...");
                        }

                        for (int i = 0; i < frameLength; i++)
                        {

                            if (hiding == "H")
                            {
                                Console.Write("\n Processing.. " + i/frameLength + "% completed");

                                videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                                b = BitmapFromArray1D(videoData, videoWidth, videoHeight);
                                b = spatial_Embed(b, block_size, T, G, bytes, rep);

                                videoData = Array1DFromBitmap(b, videoWidth * videoHeight);
                                SaveBytesToFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight, videoData);
                       
                            }
                            if (hiding == "R")
                            {
                                if (filesizespecified == 42)
                                {
                                    if (bytes.Length * rep >= (videoWidth * videoHeight /8 /block_size /block_size))
                                    loop = false;

                                    filesizespecified = bytes.Length;
                                }

                                Console.Write("\n Processing.." + i / frameLength + "% Done");

                                videoData = GetBytesFromFile(filename, (((videoWidth * videoHeight) * i) + (((videoWidth * videoHeight) * i) / 2)), videoWidth * videoHeight);
                                b = BitmapFromArray1D(videoData, videoWidth, videoHeight);
                                saveFile = spatial_Retrieve(b, 2, T, G, bytes, rep);
                                Console.WriteLine("BER (Bit Error Rate) is " + BER(bytes, saveFile));
                                SaveBytesToFile(savefilename, 0, bytes.Length, saveFile);
                            }

                        }
                        Console.WriteLine("SUCCESS! Completed 100%. Please be careful with this video.");
                        Console.ReadLine();
                    }
                    Console.WriteLine("Try a different response.");
                    Console.ReadLine();
                }
        }
    }
}
