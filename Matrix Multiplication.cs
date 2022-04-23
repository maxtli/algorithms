using System;
using System.IO;
using System.Linq;
using System.Diagnostics;
class App
{
    //item 0 is an unused flag, item 1 is the size of the matrix, item 2 is the file name of a file with one matrix entry on each line
    static void Main(string[] cmdInputs)
    {
        if (cmdInputs.Length < 3)
            return;
        int n = int.Parse(cmdInputs[1]);
        int[,] m1 = new int[n, n];
        int[,] m2 = new int[n, n];
        int r = 0;
        int c = 0;
        string[] lines = File.ReadAllLines(cmdInputs[2]);

        int i = 0;
        while (r < n)
        {
            while (c < n)
            {
                m1[r, c] = int.Parse(lines[i]);
                i++;
                c++;
            }
            c = 0;
            r++;
        }
        r = 0;
        c = 0;
        while (r < n)
        {
            while (c < n)
            {
                m2[r, c] = int.Parse(lines[i]);
                i++;
                c++;
            }
            c = 0;
            r++;
        }
        //throw new Exception(String.Join(", ", lines) + "\n" + String.Join(", ", m1.Cast<int>()) + "\n" + String.Join(", ", m2.Cast<int>()));

        Strassen sts = new Strassen();
        int[,] result = sts.strass(n, m1, m2);
        for (int j = 0; j < n; j++)
            Console.WriteLine(result[j, j]);
        Console.WriteLine();
    }

    static void checkThreshold()
    {
        Strassen sts = new Strassen();
        Stopwatch sw = new Stopwatch();
        int trials = 200;
        for (int n = 2; n < 150; n++)
        {
            int multMs = 0;
            int strassMs = 0;
            sts.Threshold = n - 1;
            for (int i = 0; i < trials; i++)
            {
                int[,] m1 = sts.randomMatrix(n);
                int[,] m2 = sts.randomMatrix(n);
                sw.Restart();
                int[,] result = sts.multiplyMatrices(n, m1, m2);
                sw.Stop();
                multMs += (int)sw.ElapsedTicks;

                sw.Restart();
                int[,] stsresult = sts.strass(n, m1, m2);
                sw.Stop();
                strassMs += (int)sw.ElapsedTicks;

                sts.checkMatrices(n, stsresult, result);
            }
            Console.WriteLine(n + ", " + multMs + ", " + strassMs);
        }
    }

    static void findTriangles()
    {
        Strassen sts = new Strassen();
        double[] probs = { 0.01, 0.02, 0.03, 0.04, 0.05 };
        int trials = 100;
        int s = 1024;
        foreach (double p in probs)
        {
            Console.WriteLine("p=" + p);
            for (int i = 0; i < trials; i++)
            {
                int triangles = 0;
                int[,] rg = sts.randomGraph(s, p);
                int[,] result = sts.strass(s, sts.strass(s, rg, rg), rg);
                for (int j = 0; j < s; j++)
                    triangles += result[j, j];
                Console.WriteLine(triangles/6);
            }
        }
    }
}

class Strassen {
    Random random;
    public Strassen()
    {
        random = new Random(909089732);

        //empirically found
        this.Threshold = 68;
    }
    public int Threshold { get; set; }
    public int[,] randomMatrix(int n)
    {
        int[,] result = new int[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i, j] = random.Next(-50, 50);
        return result;
    }

    public int[,] randomGraph(int n, double p)
    {
        int[,] result = new int[n, n];
        for(int i = 0; i < n; i++)
            for (int j = i+1; j < n; j++)
                if(random.NextDouble() < p)
                {
                    result[i, j] = 1;
                    result[j, i] = 1;
                }
        return result;
    }

    public int[,] multiplyMatrices(int n, int[,] m1, int[,] m2, int m1r = 0, int m1c = 0, int m2r = 0, int m2c = 0)
    {
        int[,] result = new int[n, n];
        int imax = Math.Min(n, m1.GetLength(0) - m1r);
        int jmax = Math.Min(n, m2.GetLength(1) - m2c);
        int kmax = Math.Min(Math.Min(n, m1.GetLength(1) - m1c), m2.GetLength(0) - m2r);
        for (int i = 0; i < imax; i++)
            for (int j = 0; j < jmax; j++)
                for (int k = 0; k < kmax; k++)
                    result[i, j] += m1[m1r + i, m1c + k] * m2[m2r + k, m2c + j];
        return result;
    }

    int[,] addMatrices(int n, int[,] m1, int[,] m2, bool subtract = false, int m1r = 0, int m1c = 0, int m2r = 0, int m2c = 0)
    {
        int[,] result = new int[n, n];
        for (int i = 0; i < Math.Min(n, m1.GetLength(0) - m1r); i++)
            for (int j = 0; j < Math.Min(n, m1.GetLength(1) - m1c); j++)
                result[i, j] += m1[m1r + i, m1c + j];
        for (int i = 0; i < Math.Min(n, m2.GetLength(0) - m2r); i++)
            for (int j = 0; j < Math.Min(n, m2.GetLength(1) - m2c); j++)
                result[i, j] += (subtract ? -1 : 1) * m2[m2r + i, m2c + j];
        return result;
    }

    public int[,] strass(int n, int[,] m1, int[,] m2, int m1r = 0, int m1c = 0, int m2r = 0, int m2c = 0)
    {
        if (n <= this.Threshold)
            return multiplyMatrices(n, m1, m2, m1r, m1c, m2r, m2c);

        int remainder = n % 2;
        int newSize = (n + remainder) / 2;

        int[,] p1 = strass(newSize, m1,
            addMatrices(newSize, m2, m2, true, m2r, m2c + newSize, m2r + newSize, m2c + newSize),
            m1r, m1c, 0, 0);
        int[,] p2 = strass(newSize,
            addMatrices(newSize, m1, m1, false, m1r, m1c, m1r, m1c + newSize),
            m2, 0, 0, m2r + newSize, m2c + newSize);
        int[,] p3 = strass(newSize,
            addMatrices(newSize, m1, m1, false, m1r + newSize, m1c, m1r + newSize, m1c + newSize),
            m2, 0, 0, m2r, m2c);
        int[,] p4 = strass(newSize, m1,
            addMatrices(newSize, m2, m2, true, m2r + newSize, m2c, m2r, m2c), 
            m1r + newSize, m1c + newSize, 0, 0);
        int[,] p5 = strass(newSize,
            addMatrices(newSize, m1, m1, false, m1r, m1c, m1r + newSize, m1c + newSize),
            addMatrices(newSize, m2, m2, false, m2r, m2c, m2r + newSize, m2c + newSize),
            0, 0, 0, 0);
        int[,] p6 = strass(newSize,
            addMatrices(newSize, m1, m1, true, m1r, m1c + newSize, m1r + newSize, m1c + newSize),
            addMatrices(newSize, m2, m2, false, m2r + newSize, m2c, m2r + newSize, m2c + newSize),
            0, 0, 0, 0);
        int[,] p7 = strass(newSize,
            addMatrices(newSize, m1, m1, true, m1r + newSize, m1c, m1r, m1c),
            addMatrices(newSize, m2, m2, false, m2r, m2c, m2r, m2c + newSize),
            0, 0, 0, 0);

        int[,] result = new int[n, n];

        for (int i = 0; i < newSize; i++)
        {
            for (int j = 0; j < newSize; j++)
                result[i, j] = p4[i, j] - p2[i, j] + p5[i, j] + p6[i, j];
            for (int j = 0; j < newSize - remainder; j++)
                result[i, newSize + j] = p1[i, j] + p2[i, j];
        }
        for (int i = 0; i < newSize - remainder; i++)
        {
            for (int j = 0; j < newSize; j++)
                result[i + newSize, j] = p3[i, j] + p4[i, j];
            for (int j = 0; j < newSize - remainder; j++)
                result[i + newSize, j + newSize] = p1[i, j] - p3[i, j] + p5[i, j] + p7[i, j];
        }
        return result;
    }

    public void printMatrix(int n, int[,] matrix)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                Console.Write(matrix[i, j] + ", ");
            Console.WriteLine();
        }
    }
    public void checkMatrices(int n, int[,] correct, int[,] check)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (correct[i, j] != check[i, j])
                    Console.Write("Oh no");
    }
}