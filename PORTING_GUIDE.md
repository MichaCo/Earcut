# Porting Guide: mapbox/earcut JS → Earcut.NET (C#)

This document describes how constructs in the upstream [mapbox/earcut](https://github.com/mapbox/earcut) JavaScript library map to the C# implementation in this repository. It is maintained alongside the port and updated automatically when upstream changes are ported.

## Upstream Reference

- **Upstream repository:** https://github.com/mapbox/earcut
- **Upstream file:** `src/earcut.js`
- **Last synced commit:** see [`.last-sync-commit`](.last-sync-commit)

---

## Language Concept Mappings

### Module / Exports

| JavaScript | C# |
|---|---|
| `export default function earcut(...)` | `public static IReadOnlyList<int> Triangulate(...)` on `public static class Earcut` |
| `export function deviation(...)` | `public static double Deviation(...)` on `Earcut` |
| `export function flatten(...)` | `public static (double[] Vertices, int[] Holes, int Dimensions) Flatten(...)` on `Earcut` |
| Internal module-private functions | `private static` methods |

### Types

| JavaScript | C# |
|---|---|
| `number` (float64) | `double` |
| `number` (integer) | `int` or `long` |
| `Array` / `number[]` | `ReadOnlySpan<double>` / `ReadOnlySpan<int>` (zero-copy) or `double[]` / `int[]` |
| `null` / `undefined` | `null` (nullable reference type: `Node?`) |
| `boolean` | `bool` |
| `{vertices, holes, dimensions}` object | C# value tuple `(double[] Vertices, int[] Holes, int Dimensions)` |
| `Infinity` | `double.PositiveInfinity` |
| `-Infinity` | `double.NegativeInfinity` |

### Linked List Node Object

The JS library uses plain objects for linked list nodes. In C# these are represented as a compact `private struct Node` (pure data, no linkage) stored inside a `private struct NodePool`.

`Node` holds only the vertex data fields:

| JavaScript property | C# field |
|---|---|
| `i` | `int I` |
| `x` | `double X` |
| `y` | `double Y` |
| `z` | `int Z` |
| `steiner` | `bool Steiner` |

The linkage fields (`prev`, `next`, `prevZ`, `nextZ`) are stored separately in `NodePool` as parallel `int[]` arrays (indices into the node array), with `NodePool.Nil = -1` as the sentinel for null:

| JavaScript property | C# accessor |
|---|---|
| `node.prev` | `pool.Prev(idx)` / `pool.SetPrev(idx, value)` |
| `node.next` | `pool.Next(idx)` / `pool.SetNext(idx, value)` |
| `node.prevZ` | `pool.PrevZ(idx)` / `pool.SetPrevZ(idx, value)` |
| `node.nextZ` | `pool.NextZ(idx)` / `pool.SetNextZ(idx, value)` |

`pool.Data(idx)` returns a `ref Node` for zero-copy field access. `pool.Allocate(i, x, y)` replaces `new Node(i, x, y)` and returns an `int` index.

### Control Flow

| JavaScript | C# |
|---|---|
| `for (let i = 0; ...)` | `for (int i = 0; ...)` |
| `for (const item of collection)` | `foreach (var item in collection)` |
| `do { ... } while (condition)` | `do { ... } while (condition);` |
| `while (condition)` | `while (condition)` |
| `if (!x)` (falsy check) | `if (x is null)` or `if (x == 0)` (explicit) |
| `const x = ...` | `var x = ...` or typed declaration |
| `let x` | `var x` or typed declaration |

### Arithmetic and Bitwise Operations

| JavaScript | C# |
|---|---|
| `` x \| 0 `` (truncate to int, also converts NaN/Infinity to 0) | `(int)x` cast — note: throws `OverflowException` for out-of-range values; the library's internal use is safe as inputs are validated |
| `` i / dim \| 0 `` | `i / dim` (integer division is already truncating in C#) |
| `` x \| (x << 8) `` | `` x \| (x << 8) `` (same syntax, operates on `int`) |
| `Math.min(a, b)` | `Math.Min(a, b)` |
| `Math.max(a, b)` | `Math.Max(a, b)` |
| `Math.abs(x)` | `Math.Abs(x)` |

### Functions

| JavaScript | C# |
|---|---|
| `function foo(a, b) { ... }` | `private static ReturnType Foo(ParamType a, ParamType b) { ... }` |
| Arrow functions `(a) => expr` | Lambda `(a) => expr` or local function |
| Default parameter `dim = 2` | Default parameter `int dim = 2` |
| Destructuring `const {vertices, holes} = flatten(...)` | Deconstruction `var (vertices, holes, dim) = Flatten(...)` |

### Collections

| JavaScript | C# |
|---|---|
| `const arr = []` | `var list = new List<int>()` |
| `arr.push(a, b, c)` | `list.Add(a); list.Add(b); list.Add(c);` |
| `arr.length` | `list.Count` |
| `arr.sort(compareFn)` | `list.Sort(Comparison<T>)` |
| `arr[i]` | `list[i]` (list) or `span[i]` (span) |

### Null / Undefined Checks

| JavaScript | C# |
|---|---|
| `if (!x)` | `if (x == NodePool.Nil)` (for node indices) or `if (x is null)` (for other refs) |
| `if (x)` | `if (x != NodePool.Nil)` (for node indices) or `if (x is not null)` (for other refs) |
| `x && x.y` | `x?.Y` (null-conditional, for non-node refs) |
| `` x \|\| default `` | `x ?? default` (null-coalescing, for non-node refs) |

### String Handling

Not applicable — the library is purely numeric.

---

## Performance Adaptations

The C# port applies several .NET-specific performance improvements beyond a mechanical translation:

| Concern | JS approach | C# approach |
|---|---|---|
| Zero-copy input | Arrays passed by reference | `ReadOnlySpan<double>` / `ReadOnlySpan<int>` |
| Memory allocation | Dynamic arrays | Pre-allocated `List<int>` with capacity hint |
| Node storage | Plain objects (heap, GC-traced) | `NodePool` struct-of-arrays: one `Node[]` + four `int[]` link arrays; nodes are `int` indices (no GC references) |
| Cache locality | Random heap access per node | Contiguous parallel arrays, sequential traversal |
| Unsafe stack frames | N/A | `[module: SkipLocalsInit]` to skip zero-init of locals |
| JIT warmup | JIT always runs | `TieredPGO=true`, `TieredCompilationQuickJit=false` |
| Math intrinsics | `Math.min/max` | `Math.Min/Max` (maps to CPU intrinsics via JIT) |

---

## API Surface Mapping

### `earcut(data, holeIndices?, dim?)` → `Earcut.Triangulate`

```js
// JS
import earcut from 'earcut';
const triangles = earcut([0,0, 1,0, 0.5,1]);
```

```csharp
// C#
using Earcut;
int[] triangles = Earcut.Triangulate([0,0, 1,0, 0.5,1]).ToArray();
```

### `flatten(data)` → `Earcut.Flatten`

```js
// JS
import { flatten } from 'earcut';
const { vertices, holes, dimensions } = flatten([[[0,0],[1,0],[0.5,1]]]);
```

```csharp
// C#
var (vertices, holes, dimensions) = Earcut.Flatten([[[0,0],[1,0],[0.5,1]]]);
```

### `deviation(data, holeIndices, dim, triangles)` → `Earcut.Deviation`

```js
// JS
import { deviation } from 'earcut';
const err = deviation(data, holeIndices, 2, triangles);
```

```csharp
// C#
double err = Earcut.Deviation(data, holeIndices, 2, triangles);
```

---

## Porting Checklist (for future upstream changes)

When a new upstream commit is detected and a sync PR is auto-generated, follow this process:

1. Review the upstream diff (`src/earcut.js`)
2. Map each JS change to its C# equivalent using the tables above
3. Update `src/Earcut/Earcut.cs` with the ported changes
4. If new test fixtures were added upstream, add them to `test/Earcut.Tests/fixtures/` and update `test/Earcut.Tests/expected.json`
5. Run `dotnet test` — all tests must pass
6. Update the `.last-sync-commit` file to the new upstream commit SHA
7. Update this guide if new language patterns were introduced
