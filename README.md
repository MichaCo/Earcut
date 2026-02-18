# Earcut.NET

[![Build](https://github.com/MichaCo/Earcut/actions/workflows/ci-cd.yml/badge.svg?branch=main)](https://github.com/MichaCo/Earcut/actions/workflows/ci-cd.yml)
[![codecov](https://codecov.io/gh/MichaCo/Earcut/graph/badge.svg)](https://codecov.io/gh/MichaCo/Earcut)
[![NuGet](https://img.shields.io/nuget/v/Earcut.NET.svg)](https://www.nuget.org/packages/Earcut.NET)

A C# port of the [mapbox/earcut](https://github.com/mapbox/earcut) polygon triangulation library.

> **Note:** This is a 100% AI-driven port, created entirely by AI agents to demonstrate modern .NET development practices.

Fast and robust ear-clipping polygon triangulation library for .NET 8+.

## Features

- Fast polygon triangulation using ear-clipping algorithm
- Supports holes
- Zero-order curve spatial indexing for performance
- Modern C# implementation with .NET 8+ optimizations
- Uses ReadOnlySpan for zero-copy performance

## Installation

```
dotnet add package Earcut.NET
```

## Usage

```csharp
using EarcutDotNet;

// Simple triangle polygon (3 vertices)
double[] coords = [0, 0, 1, 0, 0.5, 1];
var triangles = Earcut.Triangulate(coords);
// triangles: [1, 0, 2] — indices into coords (each vertex is x,y pair)

// Square polygon (4 vertices) — produces 2 triangles
double[] square = [0, 0, 1, 0, 1, 1, 0, 1];
var squareTriangles = Earcut.Triangulate(square);
// squareTriangles: [3, 0, 1, 3, 1, 2]

// Polygon with a hole
double[] polyWithHole = [
    0, 0, 10, 0, 10, 10, 0, 10,   // outer ring (4 vertices)
    2, 2, 2, 8,  8,  8, 8,  2,    // hole ring  (4 vertices, starting at vertex index 4)
];
int[] holeIndices = [4]; // hole starts at vertex index 4
var result = Earcut.Triangulate(polyWithHole, holeIndices);

// Optional: verify triangulation quality (deviation should be close to 0)
double deviation = Earcut.Deviation(polyWithHole, holeIndices, 2, result.ToArray());
```

### 3D coordinates

Pass `dim: 3` when your data contains x, y, z per vertex (only x and y are used for triangulation):

```csharp
using EarcutDotNet;

double[] coords3d = [0, 0, 1,  1, 0, 0,  0.5, 1, 0]; // x,y,z per vertex
var triangles = Earcut.Triangulate(coords3d, dim: 3);
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

MIT License