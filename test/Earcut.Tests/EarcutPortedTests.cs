// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Text.Json;
using Xunit;

namespace EarcutDotNet.Tests;

public class EarcutPortedTests
{
    private static readonly JsonSerializerOptions s_cachedJsonOptions = new()
    {
        PropertyNameCaseInsensitive = true
    };

    private class ExpectedResults
    {
        public Dictionary<string, int> Triangles { get; set; } = new();
        public Dictionary<string, double> Errors { get; set; } = new();

        [System.Text.Json.Serialization.JsonPropertyName("errors-with-rotation")]
        public Dictionary<string, double> ErrorsWithRotation { get; set; } = new();
    }

    private static ExpectedResults LoadExpectedResults()
    {
        var expectedPath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "expected.json");
        var json = File.ReadAllText(expectedPath);
        return JsonSerializer.Deserialize<ExpectedResults>(json, s_cachedJsonOptions) ?? new ExpectedResults();
    }

    public static TheoryData<string, int, int, double> GetFixtureTestData()
    {
        var expected = LoadExpectedResults();
        var rotations = new[] { 0, 90, 180, 270 };
        var theoryData = new TheoryData<string, int, int, double>();

        foreach (var kvp in expected.Triangles)
        {
            string fixtureName = kvp.Key;
            int expectedTriangles = kvp.Value;

            foreach (int rotation in rotations)
            {
                double expectedDeviation = 0.0;

                // Determine expected deviation based on rotation
                if (rotation != 0 && expected.ErrorsWithRotation.TryGetValue(fixtureName, out var value))
                {
                    expectedDeviation = value;
                }
                else if (expected.Errors.TryGetValue(fixtureName, out var value1))
                {
                    expectedDeviation = value1;
                }

                theoryData.Add(fixtureName, rotation, expectedTriangles, expectedDeviation);
            }
        }

        return theoryData;
    }

    [Theory]
    [MemberData(nameof(GetFixtureTestData))]
    public void FixtureTest(string fixtureName, int rotation, int expectedTriangles, double expectedDeviation)
    {
        // Load fixture
        var fixturePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "fixtures", $"{fixtureName}.json");
        Assert.True(File.Exists(fixturePath), $"Fixture file not found: {fixtureName}.json");

        var coordsJson = File.ReadAllText(fixturePath);
        var coords = JsonSerializer.Deserialize<double[][][]>(coordsJson);
        Assert.NotNull(coords);

        // Apply rotation if needed
        if (rotation != 0)
        {
            double theta = rotation * Math.PI / 180;
            int xx = (int)Math.Round(Math.Cos(theta));
            int xy = (int)Math.Round(-Math.Sin(theta));
            int yx = (int)Math.Round(Math.Sin(theta));
            int yy = (int)Math.Round(Math.Cos(theta));

            foreach (var ring in coords)
            {
                foreach (var coord in ring)
                {
                    double x = coord[0];
                    double y = coord[1];
                    coord[0] = xx * x + xy * y;
                    coord[1] = yx * x + yy * y;
                }
            }
        }

        // Flatten and triangulate
        var (Vertices, Holes, Dimensions) = Earcut.Flatten(coords);
        var indices = Earcut.Triangulate(Vertices, Holes, Dimensions);
        var err = Earcut.Deviation(Vertices, Holes, Dimensions, indices);

        int numTriangles = indices.Length / 3;

        // Validate triangle count (only for rotation = 0)
        if (rotation == 0)
        {
            Assert.Equal(expectedTriangles, numTriangles);
        }

        // Validate deviation
        if (expectedTriangles > 0)
        {
            Assert.True(err <= expectedDeviation,
                $"Fixture {fixtureName} rotation {rotation}: deviation {err} > expected {expectedDeviation}");
        }
    }

    [Fact]
    public void Indices2D()
    {
        double[] data = [10, 0, 0, 50, 60, 60, 70, 10];
        var indices = Earcut.Triangulate(data);
        int[] expected = [1, 0, 3, 1, 3, 2];
        Assert.Equal(expected, indices);
    }

    [Fact]
    public void Indices3D()
    {
        double[] data = [10, 0, 0, 0, 50, 0, 60, 60, 0, 70, 10, 0];
        var indices = Earcut.Triangulate(data, [], 3);
        int[] expected = [1, 0, 3, 1, 3, 2];
        Assert.Equal(expected, indices);
    }

    [Fact]
    public void Empty()
    {
        double[] data = [];
        Assert.Empty(Earcut.Triangulate(data));
    }

    [Fact]
    public void InfiniteLoop()
    {
        double[] data = [1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 4, 1, 5, 1, 3, 2, 4, 2, 4, 1];
        int[] holes = [5];
        // Should not hang
        Earcut.Triangulate(data, holes, 2);
    }

    [Fact]
    public void Refine_ImprovesABadQuadDiagonal()
    {
        double[] vertices = [0, 0, 3, 0, 10, 1, 0, 2];
        int[] triangles = [2, 3, 0, 2, 0, 1];
        double beforePerimeter = TrianglePerimeter(triangles, vertices);
        Earcut.Refine(triangles, vertices);
        double afterPerimeter = TrianglePerimeter(triangles, vertices);

        Assert.Equal(new int[] { 2, 3, 1, 3, 0, 1 }, triangles);
        Assert.True(afterPerimeter < beforePerimeter * 0.7);
        Assert.Equal(0.0, Earcut.Deviation(vertices, [], 2, triangles));
    }

    [Fact]
    public void Refine_LeavesAGoodQuadDiagonalAlone()
    {
        double[] vertices = [0, 0, 5, 0, 4, 1, 0, 4];
        int[] triangles = [2, 3, 0, 2, 0, 1];
        Earcut.Refine(triangles, vertices);

        Assert.Equal(new int[] { 2, 3, 0, 2, 0, 1 }, triangles);
        Assert.Equal(0.0, Earcut.Deviation(vertices, [], 2, triangles));
    }

    [Fact]
    public void Refine_PreservesAConcavePolygon()
    {
        double[] vertices = [0, 0, 4, 0, 4, 1, 1, 1, 1, 4, 0, 4];
        int[] triangles = Earcut.Triangulate(vertices);
        int length = triangles.Length;
        double beforePerimeter = TrianglePerimeter(triangles, vertices);
        Earcut.Refine(triangles, vertices);
        double afterPerimeter = TrianglePerimeter(triangles, vertices);

        Assert.Equal(length, triangles.Length);
        Assert.True(afterPerimeter < beforePerimeter * 0.9);
        Assert.Equal(0.0, Earcut.Deviation(vertices, [], 2, triangles));
    }

    [Fact]
    public void BlockIndexCollinear()
    {
        const int N = 30;
        var outer = new List<double[]>();
        for (int x = 0; x <= N; x++) outer.Add([x, 0]);
        for (int y = 1; y <= N; y++) outer.Add([N, y]);
        for (int x = N - 1; x >= 0; x--) outer.Add([x, N]);
        for (int y = N - 1; y >= 1; y--) outer.Add([0, y]);

        static double[][] Rect(int x0, int y0, int w, int h) =>
        [
            [x0, y0], [x0, y0 + h], [x0 + w, y0 + h], [x0 + w, y0]
        ];

        var rings = new List<double[][]> { outer.ToArray(), Rect(5, 5, 2, 4), Rect(2, 23, 1, 1) };

        foreach (int rotation in new[] { 0, 90, 180, 270 })
        {
            double theta = rotation * Math.PI / 180;
            int xx = (int)Math.Round(Math.Cos(theta));
            int xy = (int)Math.Round(-Math.Sin(theta));
            int yx = (int)Math.Round(Math.Sin(theta));
            int yy = (int)Math.Round(Math.Cos(theta));

            var rotated = rings.Select(ring =>
                ring.Select(coord => new double[] { xx * coord[0] + xy * coord[1], yx * coord[0] + yy * coord[1] })
                    .ToArray()).ToArray();

            var (vertices, holes, dimensions) = Earcut.Flatten(rotated);
            int[] indices = Earcut.Triangulate(vertices, holes, dimensions);
            double err = Earcut.Deviation(vertices, holes, dimensions, indices);
            Assert.True(err < 1e-9, $"rotation {rotation}: deviation {err} (hole dropped?)");
        }
    }

    private static double TrianglePerimeter(int[] triangles, double[] vertices, int dim = 2)
    {
        double perimeter = 0;
        for (int i = 0; i < triangles.Length; i += 3)
        {
            double ax = vertices[triangles[i] * dim],     ay = vertices[triangles[i] * dim + 1];
            double bx = vertices[triangles[i + 1] * dim], by = vertices[triangles[i + 1] * dim + 1];
            double cx = vertices[triangles[i + 2] * dim], cy = vertices[triangles[i + 2] * dim + 1];
            perimeter += Math.Sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by)) +
                         Math.Sqrt((bx - cx) * (bx - cx) + (by - cy) * (by - cy)) +
                         Math.Sqrt((cx - ax) * (cx - ax) + (cy - ay) * (cy - ay));
        }
        return perimeter;
    }
}
