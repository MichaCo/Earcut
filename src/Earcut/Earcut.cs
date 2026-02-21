// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Runtime.CompilerServices;

[module: SkipLocalsInit]
namespace EarcutDotNet;

/// <summary>
/// Fast polygon triangulation using the ear-clipping technique.
/// Port of the mapbox/earcut JavaScript library.
/// </summary>
public static class Earcut
{
    // ───────────────────────────── public API ──────────────────────────────

    /// <summary>
    /// Triangulates a polygon with optional holes.
    /// </summary>
    /// <param name="data">
    /// Flat array of vertex coordinates: [x0, y0, x1, y1, …].
    /// </param>
    /// <param name="holeIndices">
    /// Hole starting indices expressed in <em>vertex</em> counts, not coordinate counts.
    /// </param>
    /// <param name="dim">Number of coordinates per vertex (default 2).</param>
    /// <returns>Flat array of triangle vertex indices (length is a multiple of 3).</returns>
    public static int[] Triangulate(
        ReadOnlySpan<double> data,
        ReadOnlySpan<int> holeIndices = default,
        int dim = 2)
    {
        if (data.IsEmpty)
        {
            return [];
        }

        bool hasHoles = !holeIndices.IsEmpty;
        int outerLen = hasHoles ? holeIndices[0] * dim : data.Length;

        // Pre-allocate: a simple polygon with V vertices produces at most V-2 triangles.
        int vertexCount = data.Length / dim;
        var pool = new NodePool(Math.Max(8, vertexCount * 2 + 16));

        int outerNode = BuildLinkedList(ref pool, data, 0, outerLen, dim, clockwise: true);

        if (outerNode == NodePool.Nil || pool.Next(outerNode) == pool.Prev(outerNode))
        {
            return [];
        }

        var triangles = new TriangleList(Math.Max(0, (vertexCount - 2) * 3));

        if (hasHoles)
        {
            outerNode = EliminateHoles(ref pool, data, holeIndices, outerNode, dim);
        }

        double minX = 0, minY = 0, invSize = 0;

        // For non-trivial shapes use a z-order curve hash for faster look-ups.
        if (data.Length > 80 * dim)
        {
            minX = data[0];
            minY = data[1];
            double maxX = minX;
            double maxY = minY;

            // Scalar loop with Math.Min/Max intrinsics
            for (int i = dim; i < outerLen; i += dim)
            {
                double x = data[i];
                double y = data[i + 1];
                minX = Math.Min(minX, x);
                minY = Math.Min(minY, y);
                maxX = Math.Max(maxX, x);
                maxY = Math.Max(maxY, y);
            }

            invSize = Math.Max(maxX - minX, maxY - minY);
            invSize = invSize != 0.0 ? 32767.0 / invSize : 0.0;
        }

        EarcutLinked(ref pool, outerNode, triangles, dim, minX, minY, invSize, pass: 0);

        return triangles.ToArray();
    }

    /// <summary>
    /// Returns the percentage difference between the polygon area and its
    /// triangulation area — useful for verifying triangulation correctness.
    /// </summary>
    public static double Deviation(
        ReadOnlySpan<double> data,
        ReadOnlySpan<int> holeIndices,
        int dim,
        ReadOnlySpan<int> triangles)
    {
        bool hasHoles = !holeIndices.IsEmpty;
        int outerLen = hasHoles ? holeIndices[0] * dim : data.Length;

        double polygonArea = Math.Abs(SignedArea(data, 0, outerLen, dim));

        if (hasHoles)
        {
            for (int i = 0; i < holeIndices.Length; i++)
            {
                int start = holeIndices[i] * dim;
                int end = i < holeIndices.Length - 1
                    ? holeIndices[i + 1] * dim
                    : data.Length;
                polygonArea -= Math.Abs(SignedArea(data, start, end, dim));
            }
        }

        double trianglesArea = 0.0;

        for (int i = 0; i < triangles.Length; i += 3)
        {
            int a = triangles[i] * dim;
            int b = triangles[i + 1] * dim;
            int c = triangles[i + 2] * dim;
            trianglesArea += Math.Abs(
                (data[a] - data[c]) * (data[b + 1] - data[a + 1]) -
                (data[a] - data[b]) * (data[c + 1] - data[a + 1]));
        }

        return polygonArea == 0.0 && trianglesArea == 0.0
            ? 0.0
            : Math.Abs((trianglesArea - polygonArea) / polygonArea);
    }

    /// <summary>
    /// Converts a multi-dimensional coordinate array (e.g. GeoJSON rings)
    /// into the flat form that <see cref="Triangulate"/> expects.
    /// </summary>
    public static (double[] Vertices, int[] Holes, int Dimensions) Flatten(
        double[][][] data)
    {
        int dimensions = data[0][0].Length;

        // Pre-compute total capacity.
        int totalCoords = 0;
        foreach (double[][] ring in data)
        {
            totalCoords += ring.Length * dimensions;
        }

        // Allocate arrays directly - no intermediate List<T>
        var vertices = new double[totalCoords];
        var holes = new int[data.Length - 1];

        int vertexIndex = 0;
        int holeIndex = 0;
        int holeCount = 0;
        int prevLen = 0;

        foreach (double[][] ring in data)
        {
            foreach (double[] point in ring)
            {
                point.AsSpan(0, dimensions).CopyTo(vertices.AsSpan(vertexIndex));
                vertexIndex += dimensions;
            }

            if (prevLen > 0)
            {
                holeIndex += prevLen;
                holes[holeCount++] = holeIndex;
            }

            prevLen = ring.Length;
        }

        return (vertices, holes, dimensions);
    }

    // ──────────────────────── linked-list helpers ───────────────────────────

    /// <summary>
    /// Creates a circular doubly-linked list from polygon points in the
    /// specified winding order.
    /// </summary>
    private static int BuildLinkedList(
        ref NodePool pool,
        ReadOnlySpan<double> data,
        int start,
        int end,
        int dim,
        bool clockwise)
    {
        int last = NodePool.Nil;

        if (clockwise == (SignedArea(data, start, end, dim) > 0))
        {
            for (int i = start; i < end; i += dim)
            {
                last = InsertNode(ref pool, i / dim, data[i], data[i + 1], last);
            }
        }
        else
        {
            for (int i = end - dim; i >= start; i -= dim)
            {
                last = InsertNode(ref pool, i / dim, data[i], data[i + 1], last);
            }
        }

        if (last != NodePool.Nil && CoordsEqual(in pool.Data(last), in pool.Data(pool.Next(last))))
        {
            RemoveNode(ref pool, last);
            last = pool.Next(last);
        }

        return last;
    }

    /// <summary>Eliminates collinear or duplicate points.</summary>
    private static int FilterPoints(ref NodePool pool, int start, int end = NodePool.Nil)
    {
        if (start == NodePool.Nil)
        {
            return NodePool.Nil;
        }

        if (end == NodePool.Nil)
        {
            end = start;
        }

        int p = start;
        bool again;

        do
        {
            again = false;

            if (!pool.Data(p).Steiner &&
                (CoordsEqual(in pool.Data(p), in pool.Data(pool.Next(p))) ||
                 Area(in pool.Data(pool.Prev(p)), in pool.Data(p), in pool.Data(pool.Next(p))) == 0.0))
            {
                RemoveNode(ref pool, p);
                p = end = pool.Prev(p);
                if (p == pool.Next(p))
                {
                    break;
                }

                again = true;
            }
            else
            {
                p = pool.Next(p);
            }
        }
        while (again || p != end);

        return end;
    }

    // ───────────────────────── ear-clipping core ───────────────────────────

    /// <summary>Main ear-slicing loop.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void EarcutLinked(
        ref NodePool pool,
        int ear,
        TriangleList triangles,
        int dim,
        double minX,
        double minY,
        double invSize,
        int pass)
    {
        if (ear == NodePool.Nil)
        {
            return;
        }

        if (pass == 0 && invSize != 0.0)
        {
            IndexCurve(ref pool, ear, minX, minY, invSize);
        }

        int stop = ear;

        while (pool.Prev(ear) != pool.Next(ear))
        {
            int prev = pool.Prev(ear);
            int next = pool.Next(ear);

            if (invSize != 0.0
                ? IsEarHashed(ref pool, ear, minX, minY, invSize)
                : IsEar(ref pool, ear))
            {
                // Emit triangle.
                triangles.Add(pool.Data(prev).I, pool.Data(ear).I, pool.Data(next).I);

                RemoveNode(ref pool, ear);

                // Skip the next vertex — produces fewer sliver triangles.
                ear = pool.Next(next);
                stop = pool.Next(next);
                continue;
            }

            ear = next;

            if (ear == stop)
            {
                switch (pass)
                {
                    case 0:
                        EarcutLinked(
                            ref pool, FilterPoints(ref pool, ear), triangles, dim,
                            minX, minY, invSize, pass: 1);
                        break;

                    case 1:
                        ear = CureLocalIntersections(ref pool, FilterPoints(ref pool, ear), triangles);
                        EarcutLinked(
                            ref pool, ear, triangles, dim,
                            minX, minY, invSize, pass: 2);
                        break;

                    case 2:
                        SplitEarcut(ref pool, ear, triangles, dim, minX, minY, invSize);
                        break;
                }

                break;
            }
        }
    }

    /// <summary>Checks whether a polygon node forms a valid ear.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEar(ref NodePool pool, int ear)
    {
        int a = pool.Prev(ear);
        int b = ear;
        int c = pool.Next(ear);

        if (Area(in pool.Data(a), in pool.Data(b), in pool.Data(c)) >= 0.0)
        {
            return false;
        }

        double ax = pool.Data(a).X, ay = pool.Data(a).Y;
        double bx = pool.Data(b).X, by = pool.Data(b).Y;
        double cx = pool.Data(c).X, cy = pool.Data(c).Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int p = pool.Next(c);

        while (p != a)
        {
            double px = pool.Data(p).X, py = pool.Data(p).Y;
            if (px >= x0 && px <= x1 &&
                py >= y0 && py <= y1 &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, px, py) &&
                Area(in pool.Data(pool.Prev(p)), in pool.Data(p), in pool.Data(pool.Next(p))) >= 0.0)
            {
                return false;
            }

            p = pool.Next(p);
        }

        return true;
    }

    /// <summary>Ear check optimised with z-order spatial indexing.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEarHashed(
        ref NodePool pool,
        int ear,
        double minX,
        double minY,
        double invSize)
    {
        int a = pool.Prev(ear);
        int b = ear;
        int c = pool.Next(ear);

        if (Area(in pool.Data(a), in pool.Data(b), in pool.Data(c)) >= 0.0)
        {
            return false;
        }

        double ax = pool.Data(a).X, ay = pool.Data(a).Y;
        double bx = pool.Data(b).X, by = pool.Data(b).Y;
        double cx = pool.Data(c).X, cy = pool.Data(c).Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int minZ = ZOrder(x0, y0, minX, minY, invSize);
        int maxZ = ZOrder(x1, y1, minX, minY, invSize);

        int p = pool.PrevZ(ear);
        int n = pool.NextZ(ear);

        while (p != NodePool.Nil && pool.Data(p).Z >= minZ &&
               n != NodePool.Nil && pool.Data(n).Z <= maxZ)
        {
            double px = pool.Data(p).X, py = pool.Data(p).Y;
            if (px >= x0 && px <= x1 && py >= y0 && py <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, px, py) &&
                Area(in pool.Data(pool.Prev(p)), in pool.Data(p), in pool.Data(pool.Next(p))) >= 0.0)
            {
                return false;
            }

            p = pool.PrevZ(p);

            double nx = pool.Data(n).X, ny = pool.Data(n).Y;
            if (nx >= x0 && nx <= x1 && ny >= y0 && ny <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, nx, ny) &&
                Area(in pool.Data(pool.Prev(n)), in pool.Data(n), in pool.Data(pool.Next(n))) >= 0.0)
            {
                return false;
            }

            n = pool.NextZ(n);
        }

        while (p != NodePool.Nil && pool.Data(p).Z >= minZ)
        {
            double px = pool.Data(p).X, py = pool.Data(p).Y;
            if (px >= x0 && px <= x1 && py >= y0 && py <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, px, py) &&
                Area(in pool.Data(pool.Prev(p)), in pool.Data(p), in pool.Data(pool.Next(p))) >= 0.0)
            {
                return false;
            }

            p = pool.PrevZ(p);
        }

        while (n != NodePool.Nil && pool.Data(n).Z <= maxZ)
        {
            double nx = pool.Data(n).X, ny = pool.Data(n).Y;
            if (nx >= x0 && nx <= x1 && ny >= y0 && ny <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, nx, ny) &&
                Area(in pool.Data(pool.Prev(n)), in pool.Data(n), in pool.Data(pool.Next(n))) >= 0.0)
            {
                return false;
            }

            n = pool.NextZ(n);
        }

        return true;
    }

    // ───────────────────────── intersection cures ──────────────────────────

    /// <summary>
    /// Walks the polygon and fixes small local self-intersections by
    /// emitting a triangle at each crossing.
    /// </summary>
    private static int CureLocalIntersections(ref NodePool pool, int start, TriangleList triangles)
    {
        int p = start;

        do
        {
            int a = pool.Prev(p);
            int pNext = pool.Next(p);
            int b = pool.Next(pNext);

            if (!CoordsEqual(in pool.Data(a), in pool.Data(b)) &&
                Intersects(ref pool, a, p, pNext, b) &&
                LocallyInside(ref pool, a, b) &&
                LocallyInside(ref pool, b, a))
            {
                triangles.Add(pool.Data(a).I, pool.Data(p).I, pool.Data(b).I);

                RemoveNode(ref pool, p);
                RemoveNode(ref pool, pNext);

                p = start = b;
            }

            p = pool.Next(p);
        }
        while (p != start);

        return FilterPoints(ref pool, p);
    }

    /// <summary>
    /// Last-resort: find a valid diagonal, split the polygon in two, and
    /// triangulate each half independently.
    /// </summary>
    private static void SplitEarcut(
        ref NodePool pool,
        int start,
        TriangleList triangles,
        int dim,
        double minX,
        double minY,
        double invSize)
    {
        int a = start;

        do
        {
            int b = pool.Next(pool.Next(a));

            while (b != pool.Prev(a))
            {
                if (pool.Data(a).I != pool.Data(b).I && IsValidDiagonal(ref pool, a, b))
                {
                    int c = SplitPolygon(ref pool, a, b);

                    a = FilterPoints(ref pool, a, pool.Next(a));
                    c = FilterPoints(ref pool, c, pool.Next(c));

                    EarcutLinked(ref pool, a, triangles, dim, minX, minY, invSize, pass: 0);
                    EarcutLinked(ref pool, c, triangles, dim, minX, minY, invSize, pass: 0);
                    return;
                }

                b = pool.Next(b);
            }

            a = pool.Next(a);
        }
        while (a != start);
    }

    // ────────────────────────── hole elimination ───────────────────────────

    /// <summary>
    /// Links every hole into the outer loop, producing a single-ring polygon
    /// without holes.
    /// </summary>
    private static int EliminateHoles(
        ref NodePool pool,
        ReadOnlySpan<double> data,
        ReadOnlySpan<int> holeIndices,
        int outerNode,
        int dim)
    {
        var queue = new int[holeIndices.Length];
        int queueCount = 0;

        for (int i = 0; i < holeIndices.Length; i++)
        {
            int start = holeIndices[i] * dim;
            int end = i < holeIndices.Length - 1
                ? holeIndices[i + 1] * dim
                : data.Length;

            int list = BuildLinkedList(ref pool, data, start, end, dim, clockwise: false);

            if (list != NodePool.Nil)
            {
                if (list == pool.Next(list))
                {
                    pool.Data(list).Steiner = true;
                }

                queue[queueCount++] = GetLeftmost(ref pool, list);
            }
        }

        // capture a copy for the sort lambda (arrays are shared, so reads are safe)
        var poolForSort = pool;
        queue.AsSpan(0, queueCount).Sort((a, b) =>
        {
            double ax = poolForSort.Data(a).X, ay = poolForSort.Data(a).Y;
            double bx = poolForSort.Data(b).X, by = poolForSort.Data(b).Y;
            int cmp = ax.CompareTo(bx);
            if (cmp != 0)
            {
                return cmp;
            }

            cmp = ay.CompareTo(by);
            if (cmp != 0)
            {
                return cmp;
            }

            int aNext = poolForSort.Next(a);
            int bNext = poolForSort.Next(b);
            double aSlope = (poolForSort.Data(aNext).Y - ay) / (poolForSort.Data(aNext).X - ax);
            double bSlope = (poolForSort.Data(bNext).Y - by) / (poolForSort.Data(bNext).X - bx);
            return aSlope.CompareTo(bSlope);
        });

        for (int i = 0; i < queueCount; i++)
        {
            outerNode = EliminateHole(ref pool, queue[i], outerNode);
        }

        return outerNode;
    }

    private static int EliminateHole(ref NodePool pool, int hole, int outerNode)
    {
        int bridge = FindHoleBridge(ref pool, hole, outerNode);

        if (bridge == NodePool.Nil)
        {
            return outerNode;
        }

        int bridgeReverse = SplitPolygon(ref pool, bridge, hole);
        FilterPoints(ref pool, bridgeReverse, pool.Next(bridgeReverse));
        return FilterPoints(ref pool, bridge, pool.Next(bridge));
    }

    /// <summary>
    /// David Eberly's algorithm for finding a bridge between a hole and the
    /// outer polygon.
    /// </summary>
    private static int FindHoleBridge(ref NodePool pool, int hole, int outerNode)
    {
        int p = outerNode;
        double hx = pool.Data(hole).X;
        double hy = pool.Data(hole).Y;
        double qx = double.NegativeInfinity;
        int m = NodePool.Nil;

        if (CoordsEqual(in pool.Data(hole), in pool.Data(p)))
        {
            return p;
        }

        do
        {
            int pNext = pool.Next(p);

            if (CoordsEqual(in pool.Data(hole), in pool.Data(pNext)))
            {
                return pNext;
            }

            double py = pool.Data(p).Y;
            double pNextY = pool.Data(pNext).Y;
            double px = pool.Data(p).X;
            double pNextX = pool.Data(pNext).X;

            if (hy <= py && hy >= pNextY && pNextY != py)
            {
                double x = px + (hy - py) * (pNextX - px) / (pNextY - py);

                if (x <= hx && x > qx)
                {
                    qx = x;
                    m = px < pNextX ? p : pNext;
                    if (x == hx)
                    {
                        return m;
                    }
                }
            }

            p = pNext;
        }
        while (p != outerNode);

        if (m == NodePool.Nil)
        {
            return NodePool.Nil;
        }

        int stop = m;
        double mx = pool.Data(m).X;
        double my = pool.Data(m).Y;
        double tanMin = double.PositiveInfinity;

        p = m;

        do
        {
            double px = pool.Data(p).X, py = pool.Data(p).Y;
            if (hx >= px && px >= mx && hx != px &&
                PointInTriangle(
                    hy < my ? hx : qx, hy,
                    mx, my,
                    hy < my ? qx : hx, hy,
                    px, py))
            {
                double tan = Math.Abs(hy - py) / (hx - px);

                if (LocallyInside(ref pool, p, hole) &&
                    (tan < tanMin ||
                     (tan == tanMin &&
                      (px > pool.Data(m).X || (px == pool.Data(m).X && SectorContainsSector(ref pool, m, p))))))
                {
                    m = p;
                    tanMin = tan;
                }
            }

            p = pool.Next(p);
        }
        while (p != stop);

        return m;
    }

    // ───────────────────── z-order spatial indexing ─────────────────────────

    /// <summary>Interlinks polygon nodes in z-order.</summary>
    private static void IndexCurve(
        ref NodePool pool,
        int start,
        double minX,
        double minY,
        double invSize)
    {
        int p = start;

        do
        {
            ref Node pData = ref pool.Data(p);
            if (pData.Z == 0)
            {
                pData.Z = ZOrder(pData.X, pData.Y, minX, minY, invSize);
            }

            pool.SetPrevZ(p, pool.Prev(p));
            pool.SetNextZ(p, pool.Next(p));
            p = pool.Next(p);
        }
        while (p != start);

        int prevZOfStart = pool.PrevZ(p);
        pool.SetNextZ(prevZOfStart, NodePool.Nil);
        pool.SetPrevZ(p, NodePool.Nil);

        SortLinked(ref pool, p);
    }

    /// <summary>
    /// Simon Tatham's linked-list merge sort.
    /// <see href="http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html"/>
    /// </summary>
    private static int SortLinked(ref NodePool pool, int list)
    {
        int numMerges;
        int inSize = 1;

        do
        {
            int p = list;
            list = NodePool.Nil;
            int tail = NodePool.Nil;
            numMerges = 0;

            while (p != NodePool.Nil)
            {
                numMerges++;
                int q = p;
                int pSize = 0;

                for (int i = 0; i < inSize; i++)
                {
                    pSize++;
                    q = pool.NextZ(q);
                    if (q == NodePool.Nil)
                    {
                        break;
                    }
                }

                int qSize = inSize;

                while (pSize > 0 || (qSize > 0 && q != NodePool.Nil))
                {
                    int e;

                    if (pSize != 0 &&
                        (qSize == 0 || q == NodePool.Nil || pool.Data(p).Z <= pool.Data(q).Z))
                    {
                        e = p;
                        p = pool.NextZ(p);
                        pSize--;
                    }
                    else
                    {
                        e = q;
                        q = pool.NextZ(q);
                        qSize--;
                    }

                    if (tail != NodePool.Nil)
                    {
                        pool.SetNextZ(tail, e);
                    }
                    else
                    {
                        list = e;
                    }

                    pool.SetPrevZ(e, tail);
                    tail = e;
                }

                p = q;
            }

            pool.SetNextZ(tail, NodePool.Nil);
            inSize *= 2;
        }
        while (numMerges > 1);

        return list;
    }

    /// <summary>
    /// Computes the z-order (Morton code) of a point.
    /// Coords are mapped to a 15-bit unsigned integer range before
    /// bit-interleaving.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ZOrder(
        double x, double y,
        double minX, double minY,
        double invSize)
    {
        uint ix = (uint)((x - minX) * invSize);
        uint iy = (uint)((y - minY) * invSize);

        ix = (ix | (ix << 8)) & 0x00FF00FFu;
        ix = (ix | (ix << 4)) & 0x0F0F0F0Fu;
        ix = (ix | (ix << 2)) & 0x33333333u;
        ix = (ix | (ix << 1)) & 0x55555555u;

        iy = (iy | (iy << 8)) & 0x00FF00FFu;
        iy = (iy | (iy << 4)) & 0x0F0F0F0Fu;
        iy = (iy | (iy << 2)) & 0x33333333u;
        iy = (iy | (iy << 1)) & 0x55555555u;

        return (int)(ix | (iy << 1));
    }

    // ──────────────────────── geometry predicates ──────────────────────────

    /// <summary>Signed area of triangle (p → q → r).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Area(in Node p, in Node q, in Node r) =>
        (q.Y - p.Y) * (r.X - q.X) - (q.X - p.X) * (r.Y - q.Y);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool CoordsEqual(in Node a, in Node b) => a.X == b.X && a.Y == b.Y;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool PointInTriangle(
        double ax, double ay,
        double bx, double by,
        double cx, double cy,
        double px, double py) =>
        (cx - px) * (ay - py) >= (ax - px) * (cy - py) &&
        (ax - px) * (by - py) >= (bx - px) * (ay - py) &&
        (bx - px) * (cy - py) >= (cx - px) * (by - py);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool PointInTriangleExceptFirst(
        double ax, double ay,
        double bx, double by,
        double cx, double cy,
        double px, double py) =>
        !(ax == px && ay == py) &&
        PointInTriangle(ax, ay, bx, by, cx, cy, px, py);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool SectorContainsSector(ref NodePool pool, int m, int p) =>
        Area(in pool.Data(pool.Prev(m)), in pool.Data(m), in pool.Data(pool.Prev(p))) < 0.0 &&
        Area(in pool.Data(pool.Next(p)), in pool.Data(m), in pool.Data(pool.Next(m))) < 0.0;

    private static bool IsValidDiagonal(ref NodePool pool, int a, int b) =>
        pool.Data(pool.Next(a)).I != pool.Data(b).I &&
        pool.Data(pool.Prev(a)).I != pool.Data(b).I &&
        !IntersectsPolygon(ref pool, a, b) &&
        (LocallyInside(ref pool, a, b) && LocallyInside(ref pool, b, a) && MiddleInside(ref pool, a, b) &&
         (Area(in pool.Data(pool.Prev(a)), in pool.Data(a), in pool.Data(pool.Prev(b))) != 0.0 ||
          Area(in pool.Data(a), in pool.Data(pool.Prev(b)), in pool.Data(b)) != 0.0) ||
         CoordsEqual(in pool.Data(a), in pool.Data(b)) &&
         Area(in pool.Data(pool.Prev(a)), in pool.Data(a), in pool.Data(pool.Next(a))) > 0.0 &&
         Area(in pool.Data(pool.Prev(b)), in pool.Data(b), in pool.Data(pool.Next(b))) > 0.0);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool Intersects(ref NodePool pool, int p1, int q1, int p2, int q2)
    {
        int o1 = Math.Sign(Area(in pool.Data(p1), in pool.Data(q1), in pool.Data(p2)));
        int o2 = Math.Sign(Area(in pool.Data(p1), in pool.Data(q1), in pool.Data(q2)));
        int o3 = Math.Sign(Area(in pool.Data(p2), in pool.Data(q2), in pool.Data(p1)));
        int o4 = Math.Sign(Area(in pool.Data(p2), in pool.Data(q2), in pool.Data(q1)));

        return (o1 != o2 && o3 != o4) ||
               (o1 == 0 && OnSegment(ref pool, p1, p2, q1)) ||
               (o2 == 0 && OnSegment(ref pool, p1, q2, q1)) ||
               (o3 == 0 && OnSegment(ref pool, p2, p1, q2)) ||
               (o4 == 0 && OnSegment(ref pool, p2, q1, q2));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool OnSegment(ref NodePool pool, int p, int q, int r) =>
        pool.Data(q).X <= Math.Max(pool.Data(p).X, pool.Data(r).X) &&
        pool.Data(q).X >= Math.Min(pool.Data(p).X, pool.Data(r).X) &&
        pool.Data(q).Y <= Math.Max(pool.Data(p).Y, pool.Data(r).Y) &&
        pool.Data(q).Y >= Math.Min(pool.Data(p).Y, pool.Data(r).Y);

    private static bool IntersectsPolygon(ref NodePool pool, int a, int b)
    {
        int p = a;

        do
        {
            int pNext = pool.Next(p);
            if (pool.Data(p).I != pool.Data(a).I && pool.Data(pNext).I != pool.Data(a).I &&
                pool.Data(p).I != pool.Data(b).I && pool.Data(pNext).I != pool.Data(b).I &&
                Intersects(ref pool, p, pNext, a, b))
            {
                return true;
            }

            p = pNext;
        }
        while (p != a);

        return false;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool LocallyInside(ref NodePool pool, int a, int b) =>
        Area(in pool.Data(pool.Prev(a)), in pool.Data(a), in pool.Data(pool.Next(a))) < 0.0
            ? Area(in pool.Data(a), in pool.Data(b), in pool.Data(pool.Next(a))) >= 0.0 &&
              Area(in pool.Data(a), in pool.Data(pool.Prev(a)), in pool.Data(b)) >= 0.0
            : Area(in pool.Data(a), in pool.Data(b), in pool.Data(pool.Prev(a))) < 0.0 ||
              Area(in pool.Data(a), in pool.Data(pool.Next(a)), in pool.Data(b)) < 0.0;

    private static bool MiddleInside(ref NodePool pool, int a, int b)
    {
        int p = a;
        bool inside = false;
        double px = (pool.Data(a).X + pool.Data(b).X) * 0.5;
        double py = (pool.Data(a).Y + pool.Data(b).Y) * 0.5;

        do
        {
            int pNext = pool.Next(p);
            double pY = pool.Data(p).Y;
            double pNextY = pool.Data(pNext).Y;
            double pNextX = pool.Data(pNext).X;
            double pX = pool.Data(p).X;
            if (((pY > py) != (pNextY > py)) &&
                pNextY != pY &&
                px < (pNextX - pX) * (py - pY) / (pNextY - pY) + pX)
            {
                inside = !inside;
            }

            p = pNext;
        }
        while (p != a);

        return inside;
    }

    private static int GetLeftmost(ref NodePool pool, int start)
    {
        int p = start;
        int leftmost = start;

        do
        {
            double px = pool.Data(p).X, py = pool.Data(p).Y;
            double lx = pool.Data(leftmost).X, ly = pool.Data(leftmost).Y;
            if (px < lx || (px == lx && py < ly))
            {
                leftmost = p;
            }

            p = pool.Next(p);
        }
        while (p != start);

        return leftmost;
    }

    // ──────────────────── linked-list manipulation ─────────────────────────

    private static int SplitPolygon(ref NodePool pool, int a, int b)
    {
        // Read values before Allocate which may grow the pool
        int aI = pool.Data(a).I;
        double aX = pool.Data(a).X;
        double aY = pool.Data(a).Y;
        int bI = pool.Data(b).I;
        double bX = pool.Data(b).X;
        double bY = pool.Data(b).Y;

        int a2 = pool.Allocate(aI, aX, aY);
        int b2 = pool.Allocate(bI, bX, bY);
        int an = pool.Next(a);
        int bp = pool.Prev(b);

        pool.SetNext(a, b);
        pool.SetPrev(b, a);

        pool.SetNext(a2, an);
        pool.SetPrev(an, a2);

        pool.SetNext(b2, a2);
        pool.SetPrev(a2, b2);

        pool.SetNext(bp, b2);
        pool.SetPrev(b2, bp);

        return b2;
    }

    private static int InsertNode(ref NodePool pool, int i, double x, double y, int last)
    {
        int p = pool.Allocate(i, x, y);

        if (last == NodePool.Nil)
        {
            pool.SetPrev(p, p);
            pool.SetNext(p, p);
        }
        else
        {
            int lastNext = pool.Next(last);
            pool.SetNext(p, lastNext);
            pool.SetPrev(p, last);
            pool.SetPrev(lastNext, p);
            pool.SetNext(last, p);
        }

        return p;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void RemoveNode(ref NodePool pool, int p)
    {
        int pNext = pool.Next(p);
        int pPrev = pool.Prev(p);

        pool.SetPrev(pNext, pPrev);
        pool.SetNext(pPrev, pNext);

        int pPrevZ = pool.PrevZ(p);
        int pNextZ = pool.NextZ(p);

        if (pPrevZ != NodePool.Nil) pool.SetNextZ(pPrevZ, pNextZ);
        if (pNextZ != NodePool.Nil) pool.SetPrevZ(pNextZ, pPrevZ);
    }

    private static double SignedArea(
        ReadOnlySpan<double> data,
        int start,
        int end,
        int dim)
    {
        double sum = 0.0;

        for (int i = start, j = end - dim; i < end; i += dim)
        {
            sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
            j = i;
        }

        return sum;
    }

    // ─────────────────────────── tiny helpers ──────────────────────────────

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Min3(double a, double b, double c) =>
        a < b
            ? (a < c ? a : c)
            : (b < c ? b : c);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Max3(double a, double b, double c) =>
        a > b
            ? (a > c ? a : c)
            : (b > c ? b : c);

    // ─────────────────────── triangle accumulator ──────────────────────────

    /// <summary>
    /// Growable <see langword="int[]"/>-backed buffer for accumulating triangle indices.
    /// </summary>
    private sealed class TriangleList(int capacity)
    {
        private int[] _buffer = capacity > 0 ? new int[capacity] : [];
        private int _count;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Add(int a, int b, int c)
        {
            if (_count + 3 > _buffer.Length)
            {
                Array.Resize(ref _buffer, Math.Max(_buffer.Length * 2, _count + 3));
            }

            _buffer[_count] = a;
            _buffer[_count + 1] = b;
            _buffer[_count + 2] = c;
            _count += 3;
        }

        public int[] ToArray() => _count == _buffer.Length ? _buffer : _buffer[.._count];
    }

    // ──────────────────────────── node pool ────────────────────────────────

    private struct NodePool
    {
        public const int Nil = -1;

        private Node[] _nodes;
        private int[] _prev;
        private int[] _next;
        private int[] _prevZ;
        private int[] _nextZ;
        private int _count;

        public NodePool(int capacity)
        {
            _nodes = new Node[capacity];
            _prev = new int[capacity];
            _next = new int[capacity];
            _prevZ = new int[capacity];
            _nextZ = new int[capacity];
            _count = 0;
        }

        public ref Node Data(int idx) => ref _nodes[idx];
        public int Prev(int idx) => _prev[idx];
        public int Next(int idx) => _next[idx];
        public int PrevZ(int idx) => _prevZ[idx];
        public int NextZ(int idx) => _nextZ[idx];
        public void SetPrev(int idx, int value) => _prev[idx] = value;
        public void SetNext(int idx, int value) => _next[idx] = value;
        public void SetPrevZ(int idx, int value) => _prevZ[idx] = value;
        public void SetNextZ(int idx, int value) => _nextZ[idx] = value;

        public int Allocate(int i, double x, double y)
        {
            if (_count == _nodes.Length)
            {
                Grow();
            }
            int idx = _count++;
            _nodes[idx] = new Node { I = i, X = x, Y = y, Z = 0, Steiner = false };
            _prev[idx] = Nil;
            _next[idx] = Nil;
            _prevZ[idx] = Nil;
            _nextZ[idx] = Nil;
            return idx;
        }

        private void Grow()
        {
            int newCapacity = _nodes.Length * 2;
            Array.Resize(ref _nodes, newCapacity);
            Array.Resize(ref _prev, newCapacity);
            Array.Resize(ref _next, newCapacity);
            Array.Resize(ref _prevZ, newCapacity);
            Array.Resize(ref _nextZ, newCapacity);
        }
    }

    // ──────────────────────────── node type ────────────────────────────────

    private struct Node
    {
        public int I;        // vertex index in the coordinate array
        public double X;     // X coordinate
        public double Y;     // Y coordinate
        public int Z;        // z-order curve value
        public bool Steiner; // is this a Steiner point?
    }
}
