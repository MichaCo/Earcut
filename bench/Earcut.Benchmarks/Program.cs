using BenchmarkDotNet.Running;

namespace Earcut.Benchmarks;

public class Program
{
    public static void Main(string[] args)
    {
        BenchmarkRunner.Run<TriangulationBenchmarks>();
    }
}
