// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Text.Json;
using BenchmarkDotNet.Attributes;

namespace EarcutDotNet.Benchmarks;

[MemoryDiagnoser]
[ShortRunJob]
public class TriangulationBenchmarks
{
    private double[] _square = null!;
    private double[] _complexPolygon = null!;
    private double[] _dudeVertices = null!;
    private int[] _dudeHoles = null!;
    private double[] _waterVertices = null!;
    private int[] _waterHoles = null!;
    private double[] _waterHugeVertices = null!;
    private int[] _waterHugeHoles = null!;

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

        // Load dude fixture
        LoadFixture("dude", out _dudeVertices, out _dudeHoles);

        // Load water fixture
        LoadFixture("water", out _waterVertices, out _waterHoles);

        // Load water-huge fixture
        LoadFixture("water-huge", out _waterHugeVertices, out _waterHugeHoles);
    }

    private static void LoadFixture(string name, out double[] vertices, out int[] holes)
    {
        string path = Path.Combine(
            AppDomain.CurrentDomain.BaseDirectory, "fixtures", $"{name}.json");

        var json = File.ReadAllText(path);
        var coords = JsonSerializer.Deserialize<double[][][]>(json)
            ?? throw new InvalidDataException($"Fixture '{name}' deserialized to null.");

        (vertices, holes, _) = Earcut.Flatten(coords);
    }

    [Benchmark]
    public IReadOnlyList<int> TriangulateSquare()
    {
        var r = Earcut.Triangulate(_square);
        if (r.Length == 0)
        {
            throw new Exception();
        }

        return r;
    }

    [Benchmark]
    public IReadOnlyList<int> TriangulateComplexPolygon()
    {
        var r =  Earcut.Triangulate(_complexPolygon);
        if (r.Length == 0)
        {
            throw new Exception();
        }

        return r;
    }

    [Benchmark]
    public IReadOnlyList<int> TriangulateDude()
    {
        var r = Earcut.Triangulate(_dudeVertices, _dudeHoles);
        if (r.Length == 0)
        {
            throw new Exception();
        }

        return r;
    }

    [Benchmark]
    public IReadOnlyList<int> TriangulateWater()
    {
        var r = Earcut.Triangulate(_waterVertices, _waterHoles);
        if (r.Length == 0)
        {
            throw new Exception();
        }

        return r;
    }

    [Benchmark]
    public IReadOnlyList<int> TriangulateWaterHuge()
    {
        var r = Earcut.Triangulate(_waterHugeVertices, _waterHugeHoles);
        if (r.Length == 0)
        {
            throw new Exception();
        }

        return r;
    }
}
