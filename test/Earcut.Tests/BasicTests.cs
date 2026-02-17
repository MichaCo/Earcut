using Xunit;
using ModernEarcut;

namespace Earcut.Tests;

public class BasicTests
{
    [Fact]
    public void Triangulate_EmptyArray_ReturnsEmpty()
    {
        var result = ModernEarcut.Earcut.Triangulate([]);
        Assert.Empty(result);
    }

    [Fact]
    public void Triangulate_Triangle_ReturnsThreeIndices()
    {
        double[] coords = [0, 0, 1, 0, 0.5, 1];
        var result = ModernEarcut.Earcut.Triangulate(coords);
        Assert.Equal(3, result.Length);
    }

    [Fact]
    public void Triangulate_Square_ReturnsTwoTriangles()
    {
        double[] coords = [0, 0, 1, 0, 1, 1, 0, 1];
        var result = ModernEarcut.Earcut.Triangulate(coords);
        Assert.Equal(6, result.Length);
    }
}
