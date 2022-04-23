using System;
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Collections.Generic;
class App
{
    static long longRandom(Random rand)
    {
        byte[] buf = new byte[8];
        rand.NextBytes(buf);
        long longRand = BitConverter.ToInt64(buf, 0);

        return (Math.Abs(longRand) % 1000000000000 + 1);
    }

    static void Main(string[] cmdInputs)
    {
        Random rand = new Random();
        List<long> l;
        NumPar np = new NumPar(rand);
        if (cmdInputs.Length < 3)
        {
            for (int i = 0; i < 50; i++)
            {
                l = new List<long>();
                for (int j = 0; j < 100; j++)
                    l.Add(longRandom(rand));

                Stopwatch sw = Stopwatch.StartNew();

                string s1 = np.runKKarp(l).ToString();
                string t1 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                string s2 = np.runRegAlgo(l).ToString();
                string t2 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                string s3 = np.runRegAlgo(l, true).ToString();
                string t3 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                string s4 = np.runRegAlgo(l, true, true).ToString();
                string t4 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                string s5 = np.runPrePartitionedAlgo(l).ToString();
                string t5 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                string s6 = np.runPrePartitionedAlgo(l, true).ToString();
                string t6 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                string s7 = np.runPrePartitionedAlgo(l, true, true).ToString();
                string t7 = sw.ElapsedMilliseconds.ToString();
                sw.Restart();

                Console.WriteLine(String.Join("; ", new string[] {s1, t1, s2, t2, s3, t3, s4, t4, s5, t5, s6, t6, s7, t7}));
            }
        }
        else
        {
            //command line inputs: 0 is an unused flag, 1 is an enumeration for the algorithm to use, and 2 is the file name of a file containing a list with one number on each line to run the algorithm
            l = File.ReadAllLines(cmdInputs[2]).ToList()
                .ConvertAll<long>(line => long.Parse(line));
            Console.WriteLine(
                (cmdInputs[1] == "0" ? np.runKKarp(l) :
                cmdInputs[1] == "1" ? np.runRegAlgo(l) :
                cmdInputs[1] == "2" ? np.runRegAlgo(l, true) :
                cmdInputs[1] == "3" ? np.runRegAlgo(l, true, true) :
                cmdInputs[1] == "11" ? np.runPrePartitionedAlgo(l) :
                cmdInputs[1] == "12" ? np.runPrePartitionedAlgo(l, true) :
                cmdInputs[1] == "13" ? np.runPrePartitionedAlgo(l, true, true) :
                0).ToString()
            );
        }
    }
}


class NumPar
{
    int maxiter = 25000;
    Random rand;
    public NumPar(Random r)
    {
        rand = r;
    }
    public int MaxIter
    {
        set { maxiter = value; }
    }
    public long runKKarp(List<long> l)
    {
        Heap heap = new Heap(l);
        while (heap.Length > 1)
            heap.insert(heap.extractMax() - heap.extractMax());
        return heap.extractMax();
    }

    //no prepartition
    public long runRegAlgo(List<long> l, bool hillClimb = false, bool simAnn = false)
    {
        int i = 0;
        int[] signs = randomSigns(l.Count);

        //for an allocation, the residue is the absolute value of the dot product between the original list and list of signs
        long minres = Math.Abs(l.Sum(num => num * (signs[i++])));

        //prevres keeps track of the most recently visited allocation
        long prevres = minres;
        for (int j = 0; j < maxiter; j++)
        {
            i = 0;

            //repeated random method
            if (!hillClimb)
            {
                //re-draw random signs
                signs = randomSigns(l.Count);

                //if the residue is below the minimum, update it
                minres = Math.Min(minres, Math.Abs(l.Sum(num => num * (signs[i++]))));
                continue;
            }

            //hill climb + simulated annealing
            //draw distint nb1, nb2 in [0, count-1]
            int nb1 = rand.Next(0, l.Count);
            int nb2 = rand.Next(0, l.Count - 1);
            if (nb2 >= nb1)
                nb2 += 1;

            //change nb1 with probability 1, change nb2 eith probability 1/2
            bool changeNb2 = rand.NextDouble() < .5;

            //residue under new allocation
            long newres = Math.Abs(l.Sum(num => num * ((i == nb1
                || (changeNb2 && i == nb2) ? -1 : 1) * signs[i++])));

            //if the new residue beats the most recently visited residue (or with some probility), "visit" the new allocation by changing the signs in the original list
            if (newres < prevres || (
                simAnn && rand.NextDouble() < Math.Exp(
                (prevres - newres) / cooling(j))
                ))
            {
                prevres = newres;
                signs[nb1] *= 1;
                signs[nb2] *= changeNb2 ? -1 : 1;
            }

            //update the minimum
            minres = Math.Min(minres, newres);
        }
        return minres;
    }
    public long runPrePartitionedAlgo(List<long> l, bool hillClimb = false, bool simAnn = false)
    {
        //prepartition the list with a random partition
        int[] partition = randomSigns(l.Count, true);

        //plist stores the sum of each group of numbers in the prepartition
        List<long> pList = new List<long>(new long[l.Count]);
        for (int i = 0; i < l.Count; i++)
            pList[partition[i]] += l[i];

        //kkarp -- best way to divide the prepartition without splitting any sets
        long minres = runKKarp(pList);
        long prevres = minres;

        for(int j=0; j< maxiter; j++)
        {
            //repeated random method
            if (!hillClimb)
            {
                //make a new random partition and update the minimum
                partition = randomSigns(l.Count, true);
                pList = new List<long>(new long[l.Count]);
                for (int i = 0; i < l.Count; i++)
                    pList[partition[i]] += l[i];
                minres = Math.Min(minres, runKKarp(pList));
                continue;
            }

            //hill climb + simulated annealing
            //draw distinct nb1, nb2 in [0, count-1]
            int nb1 = rand.Next(0, l.Count);
            int nb2 = rand.Next(0, l.Count - 1);
            if (nb2 >= partition[nb1])
                nb2 += 1;

            //update the partition with the random change by removing a random element nb1 from its original set in the prepartition and adding it to a different set nb2
            pList[partition[nb1]] -= l[nb1];
            pList[nb2] += l[nb1];
            long newres = runKKarp(pList);

            //if the new residue beats the previously visited residue (or with some probability), "visit" the new allocation
            if (newres < prevres || (
                simAnn && rand.NextDouble() < Math.Exp(
                (prevres - newres) / cooling(j))
                ))
                prevres = newres;
            //otherwise, don't visit the new allocation, and reset the prepartition back to its previous value
            else
            {
                pList[partition[nb1]] += l[nb1];
                pList[nb2] -= l[nb1];
            }
            
            //update the minimum
            minres = Math.Min(minres, newres);
        }
        return minres;
    }

    //partition = true: outputs random partition rather than random sign
    int[] randomSigns(int length, bool partition = false)
    {
        int[] signs = new int[length];
        for (int i = 0; i < length; i++)
            signs[i] = partition ? rand.Next(0,100) : rand.Next(0, 2) * 2 - 1;
        return signs;
    }
    double cooling(int iter)
    {
        return Math.Pow(10, 10) * Math.Pow(0.8, Math.Floor((double) iter / 300));
    }
}

class Heap
{
    List<long> heap;

    public Heap(List<long> l)
    {
        heap = new List<long>(l);
        for (int i = (int) Math.Floor((double) l.Count / 2); i>= 0; i--)
            maxHeapify(i);

    }
    public int Length
    {
        get { return heap.Count; }
    }
    int parent(int i)
    {
        return (int) Math.Floor((double) (i+1) / 2) - 1;
    }

    int left(int i)
    {
        return 2 * (i + 1) - 1;
    }
    
    int right(int i)
    {
        return 2 * (i + 1);
    }

    void maxHeapify(int n)
    {
        int largest = n;
        int l = left(n);
        int r = right(n);

        //bottom layer node automatically a heap
        if (heap.Count <= left(n))
            return;

        //get largest of node and its children
        if (heap[l] > heap[n])
            largest = l;
        if (heap.Count > right(n) && heap[r] > heap[largest])
            largest = r;
        if (largest != n)
        {
            //swap node with largest child
            long lg = heap[largest];
            heap[largest] = heap[n];
            heap[n] = lg;

            //reheapify child heap
            maxHeapify(largest);
        }
    }

    public void printHeap()
    {
        Console.WriteLine(String.Join(", ", heap));
    }

    public long extractMax()
    {
        if (heap.Count == 0)
            throw new Exception("empty heap");

        //extract first element
        long max = heap[0];

        //put last element first
        heap[0] = heap[heap.Count -1];
        heap.RemoveAt(heap.Count - 1);

        //re-heapify
        maxHeapify(0);
        return max;
    }

    public void insert(long v)
    {
        //add element at end of heap
        int n = heap.Count;
        heap.Add(v);

        //re-heapify all parent nodes of this element as necessary
        int np = parent(n);
        while (n > 0 && heap[np] < heap[n])
        {
            long lg = heap[n];
            heap[n] = heap[np];
            heap[np] = lg;
            n = np;
            np = parent(n);
        }
    }

}
