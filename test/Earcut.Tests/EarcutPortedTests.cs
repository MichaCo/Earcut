using System.Text.Json;
using Xunit;
using ModernEarcut;

namespace Earcut.Tests;

public class EarcutPortedTests
{
    private static readonly string FixturesPath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "fixtures");
    private static readonly string ExpectedPath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "expected.json");

    [Fact]
    public void Indices2D()
    {
        double[] data = [10, 0, 0, 50, 60, 60, 70, 10];
        var indices = ModernEarcut.Earcut.Triangulate(data);
        int[] expected = [1, 0, 3, 3, 2, 1];
        Assert.Equal(expected, indices);
    }

    [Fact]
    public void Indices3D()
    {
        double[] data = [10, 0, 0, 0, 50, 0, 60, 60, 0, 70, 10, 0];
        var indices = ModernEarcut.Earcut.Triangulate(data, [], 3);
        int[] expected = [1, 0, 3, 3, 2, 1];
        Assert.Equal(expected, indices);
    }

    [Fact]
    public void Empty()
    {
        double[] data = [];
        Assert.Empty(ModernEarcut.Earcut.Triangulate(data));
    }

    [Fact]
    public void InfiniteLoop()
    {
        double[] data = [1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 4, 1, 5, 1, 3, 2, 4, 2, 4, 1];
        int[] holes = [5];
        // Should not hang
        ModernEarcut.Earcut.Triangulate(data, holes, 2);
    }

    [Fact]
    public void Building()
    {
        var fixturePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "fixtures", "building.json");
        if (!File.Exists(fixturePath)) return;
        
        var coordsJson = File.ReadAllText(fixturePath);
        var coords = JsonSerializer.Deserialize<double[][][]>(coordsJson);
        if (coords == null) return;

        var data = ModernEarcut.Earcut.Flatten(coords);
        var indices = ModernEarcut.Earcut.Triangulate(data.vertices, data.holes, data.dimensions);
        
        Assert.Equal(13, indices.Length / 3);
    }

    [Fact]
    public void Dude()
    {
        var fixturePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "fixtures", "dude.json");
        if (!File.Exists(fixturePath)) return;
        
        var coordsJson = File.ReadAllText(fixturePath);
        var coords = JsonSerializer.Deserialize<double[][][]>(coordsJson);
        if (coords == null) return;

        var data = ModernEarcut.Earcut.Flatten(coords);
        var indices = ModernEarcut.Earcut.Triangulate(data.vertices, data.holes, data.dimensions);
        var err = ModernEarcut.Earcut.Deviation(data.vertices, data.holes, data.dimensions, indices);
        
        Assert.Equal(106, indices.Length / 3);
        Assert.True(err <= 2e-15);
    }

    [Fact]
    public void Water()
    {
        var fixturePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "fixtures", "water.json");
        if (!File.Exists(fixturePath)) return;
        
        var coordsJson = File.ReadAllText(fixturePath);
        var coords = JsonSerializer.Deserialize<double[][][]>(coordsJson);
        if (coords == null) return;

        var data = ModernEarcut.Earcut.Flatten(coords);
        var indices = ModernEarcut.Earcut.Triangulate(data.vertices, data.holes, data.dimensions);
        var err = ModernEarcut.Earcut.Deviation(data.vertices, data.holes, data.dimensions, indices);
        
        Assert.Equal(2482, indices.Length / 3);
        Assert.True(err <= 0.0008);
    }

    [Fact]
    public void BadHole()
    {
        var fixturePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "fixtures", "bad-hole.json");
        if (!File.Exists(fixturePath)) return;
        
        var coordsJson = File.ReadAllText(fixturePath);
        var coords = JsonSerializer.Deserialize<double[][][]>(coordsJson);
        if (coords == null) return;

        var data = ModernEarcut.Earcut.Flatten(coords);
        var indices = ModernEarcut.Earcut.Triangulate(data.vertices, data.holes, data.dimensions);
        var err = ModernEarcut.Earcut.Deviation(data.vertices, data.holes, data.dimensions, indices);
        
        Assert.Equal(42, indices.Length / 3);
        Assert.True(err <= 0.019);
    }
}
