// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using BenchmarkDotNet.Running;

namespace EarcutDotNet.Benchmarks;

public class Program
{
    public static void Main()
    {
        var x = new TriangulationBenchmarks();
        x.Setup();

        Console.WriteLine("start");

        for (var i = 0; i < 100; i++)
        {
            var m = x.TriangulateComplexPolygon();
        }

        Console.WriteLine("end");

        BenchmarkRunner.Run<TriangulationBenchmarks>();
    }
}
