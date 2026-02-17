// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using Xunit;

namespace Earcut.Tests;

public class BasicTests
{
    [Fact]
    public void Triangulate_EmptyArray_ReturnsEmpty()
    {
        var result = Earcut.Triangulate([]);
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
}
