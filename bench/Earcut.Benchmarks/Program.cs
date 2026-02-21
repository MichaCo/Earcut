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
        BenchmarkRunner.Run<TriangulationBenchmarks>();
    }
}
