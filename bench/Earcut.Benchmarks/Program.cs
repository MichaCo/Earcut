// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using BenchmarkDotNet.Running;

namespace EarcutDotNet.Benchmarks;

public class Program
{
    public static void Main(string[] args)
    {
        var m = new TriangulationBenchmarks();
        m.Setup();
        m.TriangulateSquare();
        m.TriangulateWater();
        m.TriangulateWaterHuge();
        m.TriangulateDude();
        m.TriangulateComplexPolygon();

        BenchmarkSwitcher.FromAssembly(typeof(TriangulationBenchmarks).Assembly).Run(args);
    }
}
