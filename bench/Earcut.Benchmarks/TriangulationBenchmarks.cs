using BenchmarkDotNet.Attributes;
using ModernEarcut;

namespace Earcut.Benchmarks;

[MemoryDiagnoser]
public class TriangulationBenchmarks
{
    private double[] _square = null!;
    private double[] _complexPolygon = null!;

    [GlobalSetup]
    public void Setup()
    {
        _square = [0, 0, 10, 0, 10, 10, 0, 10];
        
        // Create a 100-vertex polygon
        var vertices = new List<double>();
        for (int i = 0; i < 100; i++)
        {
            double angle = 2 * Math.PI * i / 100;
            vertices.Add(Math.Cos(angle) * 100);
            vertices.Add(Math.Sin(angle) * 100);
        }
        _complexPolygon = vertices.ToArray();
    }

    [Benchmark]
    public int[] TriangulateSquare()
    {
        return Earcut.Triangulate(_square);
    }

    [Benchmark]
    public int[] TriangulateComplexPolygon()
    {
        return Earcut.Triangulate(_complexPolygon);
    }
}
