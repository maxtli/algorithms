using System.Diagnostics;

class App
{
    static void Main(string[] cmdInputs)
    {
        MST mst = new MST();
        if (cmdInputs.Length == 0)
        {
            foreach((int n, int trials) in new (int, int)[] {
                (128, 100), (256, 100), (512, 100),
                (1024, 100), (2048, 50), (4096, 50),
                (8192, 25), (16384, 25), (32768, 10),
                (65536, 10), (131072, 5), (262144, 5)
            })
            {
                foreach (int dimension in new int[] { 0, 2, 3, 4 })
                    mst.runTrials(n, dimension, trials);
            }
        }
        else {
            //support for command line arguments
            int[] newArgs = Array.ConvertAll(cmdInputs, int.Parse);
            mst.runTrials(newArgs[1], newArgs[3], newArgs[2]);
        }
        Console.ReadLine();
    }
}

class MST
{
    private readonly Random r = new();
    public void runTrials(int n, int dimension, int trials)
    {
        Console.WriteLine(n + ", " + dimension + ", " + trials);
        Stopwatch sw = new Stopwatch();

        //list containing tree weights for each trial
        List<double> totalWeights = new List<double>();

        //perform multiple trials
        for (int trial = 0; trial < trials; trial++)
        {
            //show progress on large n cases
            if (n > 8192)
            {
                Console.WriteLine("Trial " + trial);

                //generate edges
                sw.Restart();
                List<(int, int, double weight)> edges = generateEdges(n, dimension);
                sw.Stop();
                if (n > 8192)
                    Console.WriteLine("Edges Produced (" + edges.Count + ") in " + sw.Elapsed.TotalSeconds + " seconds");

                //find MST
                sw.Restart();
                totalWeights.Add(findMST(n, edges, dimension > 0));
                sw.Stop();

                if (n > 8192)
                    Console.WriteLine("MST produced in " + sw.Elapsed.TotalSeconds + " seconds");
            } else
            {
                //generate edges
                List<(int, int, double weight)> edges = generateEdges(n, dimension);

                //find MST
                totalWeights.Add(findMST(n, edges, dimension > 0));
            }
        }
        Console.WriteLine();

        //show individual trials if sufficiently small number
        if (trials < 32)
        {
            Console.WriteLine("Trial tree weights");
            Console.WriteLine(String.Join("\n ", totalWeights.ToArray()));
            Console.WriteLine();
        }

        Console.WriteLine("Average tree weight");
        Console.WriteLine(totalWeights.Average());
        Console.WriteLine();
    }
    public List<(int, int, double)> generateEdges(int n, int dimension)
    {
        //list of relevant edges to search through
        List<(int x, int y, double weight)> edges = new List<(int, int, double)>();

        //empirical edge bounds: MST will not contain any edge greater than the edge bound 99.99% of the time
        double edgeBound = dimension;
        if (dimension == 0)
        {
            edgeBound = 3.3607 * Math.Pow(n, -.686);

            //generate all edges among n vertices with random weights
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    //exclude edge if greater than edge bound
                    if (r.NextDouble() is double nxt && nxt < edgeBound)
                        edges.Add((i, j, nxt));

            return edges;
        }
        switch (dimension)
        {
            case 2: edgeBound = 2.0207 * Math.Pow(n, -.384 ); break;
            case 3: edgeBound = 1.9626 * Math.Pow(n, -.295); break;
            case 4: edgeBound = 1.9917 * Math.Pow(n, -.246); break;
        }

        //square edge bounds (monotonic transformation)
        edgeBound = edgeBound * edgeBound;

        //array is faster than list for vertices because we know there are exactly n
        double[,] vertices = new double[n,dimension];

        for (int i = 0; i < n; i++)
        {
            //generate vertex i
            for (int coord = 0; coord < dimension; coord++)
                vertices[i,coord] = r.NextDouble();

            //generate edge connecting vertex i with all existing vertices
            for (int j = 0; j < i; j++)
            {
                double d = 0;
                for (int e = 0; e < dimension; e++)
                    d += (vertices[i,e] - vertices[j,e]) * (vertices[i,e] - vertices[j,e]);

                //compare squared edges to squared edge bound
                if (d < edgeBound)
                    edges.Add((i, j, d));
            }
        }
        return edges;
    }
    public double findMST(int n, List<(int, int, double weight)> edges, bool squared)
    {
        //sort edges
        edges.Sort((x, y) => x.weight.CompareTo(y.weight));

        //union find data structure -- each vertex included in any iteration of the MST is identified with the set of vertices it is already connected with
        Dictionary<int, int> unionfind = new Dictionary<int, int>();

        //find the root of the tree of a vertex
        int findSet(int key)
        {
            if (!unionfind.ContainsKey(key))
                unionfind[key] = key;
            //path compaction
            else if (unionfind[key] != key)
                unionfind[key] = findSet(unionfind[key]);
            return unionfind[key];
        }
        double totalWeight = 0;
        int edgeCount = 0;
        int resetCounter = n / 2;

        //loop through edges starting from least weight
        foreach ((int x, int y, double weight) in edges)
        {
            int rootx = findSet(x);
            int rooty = findSet(y);

            //edge between x and y does not create a cycle
            if (findSet(x) != findSet(y))
            {
                //merge sets of x and y
                unionfind[rootx] = rooty;

                //account for squared edges (monotonic transformation)
                totalWeight += squared ? Math.Sqrt(weight) : weight;
                edgeCount++;

                //when there are n-1 vertices, all vertices must be connected
                if (edgeCount == n - 1)
                    break;
            }
        }

        //there were not enough edges to complete the MST
        if (edgeCount < n - 1)
            throw new Exception("couldn't find MST due to bad edge bound");
        return totalWeight;
    }
}