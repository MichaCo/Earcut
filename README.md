# Earcut.NET

A C# port of the [mapbox/earcut](https://github.com/mapbox/earcut) polygon triangulation library.

Fast and robust ear-clipping polygon triangulation library for .NET 8+.

## Features

- Fast polygon triangulation using ear-clipping algorithm
- Supports holes
- Zero-order curve spatial indexing for performance
- Modern C# implementation with .NET 8+ optimizations

## Usage

```csharp
using ModernEarcut;

// Simple polygon
double[] coords = [0, 0, 1, 0, 0.5, 1];
int[] triangles = Earcut.Triangulate(coords);

// Polygon with holes
double[] coords = [0, 0, 10, 0, 10, 10, 0, 10,  // outer ring
                   2, 2, 2, 8, 8, 8, 8, 2];      // hole
int[] holeIndices = [4];  // hole starts at vertex 4
int[] triangles = Earcut.Triangulate(coords, holeIndices);
```

## Building

```bash
dotnet build
```

## Testing

```bash
dotnet test
```

## Benchmarking

```bash
dotnet run -c Release --project bench/Earcut.Benchmarks
```

## License

ISC License - Same as the original earcut library