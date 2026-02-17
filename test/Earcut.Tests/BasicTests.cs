// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using Xunit;

namespace Earcut.Tests;

public class BasicTests
{
    [Fact]
    public void Triangulate_NUll_ReturnsEmpty()
    {
        var result = Earcut.Triangulate(null);
        Assert.Empty(result);
    }

    [Fact]
    public void Triangulate_EmptyArray_ReturnsEmpty()
    {
        var result = Earcut.Triangulate([], [1]);
        Assert.Empty(result);
    }

    [Fact]
    public void Triangulate_Triangle_ReturnsThreeIndices()
    {
        double[] coords = [0, 0, 1, 0, 0.5, 1];
        var result = Earcut.Triangulate(coords);
        Assert.Equal(3, result.Count);
    }

    [Fact]
    public void Triangulate_Square_ReturnsTwoTriangles()
    {
        double[] coords = [0, 0, 1, 0, 1, 1, 0, 1];
        var result = Earcut.Triangulate(coords);
        Assert.Equal(6, result.Count);
    }

    [Fact]
    public void Flatten_SpanOverload_MatchesTupleOverload()
    {
        double[][][] data =
        [
            [[0, 0], [10, 0], [10, 10], [0, 10]],   // outer ring
            [[2, 2], [2,  8], [ 8,  8], [8,  2]],   // hole
        ];

        // Reference result from the tuple-returning overload.
        var (refVertices, refHoles, refDim) = Earcut.Flatten(data);

        // Span overload: caller provides pre-allocated buffers.
        double[] vertices = new double[refVertices.Length];
        int[]    holes    = new int[refHoles.Length];

        int dim = Earcut.Flatten(data, vertices, holes);

        Assert.Equal(refDim, dim);
        Assert.Equal(refVertices, vertices);
        Assert.Equal(refHoles, holes);
    }
}
