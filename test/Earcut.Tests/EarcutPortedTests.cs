// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Text.Json;
using Xunit;

namespace Earcut.Tests;

public class EarcutPortedTests
{
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
        var options = new JsonSerializerOptions
        {
            PropertyNameCaseInsensitive = true
        };
        return JsonSerializer.Deserialize<ExpectedResults>(json, options) ?? new ExpectedResults();
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
                if (rotation != 0 && expected.ErrorsWithRotation.ContainsKey(fixtureName))
                {
                    expectedDeviation = expected.ErrorsWithRotation[fixtureName];
                }
                else if (expected.Errors.ContainsKey(fixtureName))
                {
                    expectedDeviation = expected.Errors[fixtureName];
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
        var data = Earcut.Flatten(coords);
        var indices = Earcut.Triangulate(data.Vertices, data.Holes, data.Dimensions);
        var err = Earcut.Deviation(data.Vertices, data.Holes, data.Dimensions, indices.ToArray());

        int numTriangles = indices.Count / 3;

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
        int[] expected = [1, 0, 3, 3, 2, 1];
        Assert.Equal(expected, indices);
    }

    [Fact]
    public void Indices3D()
    {
        double[] data = [10, 0, 0, 0, 50, 0, 60, 60, 0, 70, 10, 0];
        var indices = Earcut.Triangulate(data, [], 3);
        int[] expected = [1, 0, 3, 3, 2, 1];
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
}
