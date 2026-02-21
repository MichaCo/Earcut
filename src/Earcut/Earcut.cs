// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Buffers;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

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

        // Every hole-bridge split adds 2 cloned nodes (SplitPolygon); extra headroom.
        int nodeCapacity = vertexCount + holeIndices.Length * 2 + 2;
        var pool = new NodePool(nodeCapacity);
        try
        {
            int outerNode = BuildLinkedList(data, 0, outerLen, dim, clockwise: true, ref pool);

            if (outerNode == -1 || pool[outerNode].Next == pool[outerNode].Prev)
            {
                return [];
            }

            var triangles = new TriangleList(Math.Max(0, (vertexCount - 2) * 3));

            if (hasHoles)
            {
                outerNode = EliminateHoles(data, holeIndices, outerNode, dim, ref pool);
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

            EarcutLinked(outerNode, triangles, dim, minX, minY, invSize, pass: 0, ref pool);

            return triangles.ToArray();
        }
        finally
        {
            pool.Dispose();
        }
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
        ReadOnlySpan<double> data,
        int start,
        int end,
        int dim,
        bool clockwise,
        ref NodePool pool)
    {
        int last = -1;

        if (clockwise == (SignedArea(data, start, end, dim) > 0))
        {
            for (int i = start; i < end; i += dim)
            {
                last = InsertNode(i / dim, data[i], data[i + 1], last, ref pool);
            }
        }
        else
        {
            for (int i = end - dim; i >= start; i -= dim)
            {
                last = InsertNode(i / dim, data[i], data[i + 1], last, ref pool);
            }
        }

        if (last != -1)
        {
            int lastNext = pool[last].Next;
            if (pool[last].X == pool[lastNext].X && pool[last].Y == pool[lastNext].Y)
            {
                RemoveNode(last, ref pool);
                last = lastNext;
            }
        }

        return last;
    }

    /// <summary>Eliminates collinear or duplicate points.</summary>
    private static int FilterPoints(int start, ref NodePool pool, int end = -1)
    {
        if (start == -1)
        {
            return -1;
        }

        if (end == -1)
        {
            end = start;
        }

        int p = start;
        bool again;

        do
        {
            again = false;

            ref NodeData pNode = ref pool[p];
            if (!pNode.Steiner &&
                (pool[pNode.Next].X == pNode.X && pool[pNode.Next].Y == pNode.Y ||
                 Area(pNode.Prev, p, pNode.Next, ref pool) == 0.0))
            {
                int pPrev = pNode.Prev;
                RemoveNode(p, ref pool);
                p = end = pPrev;
                if (p == pool[p].Next)
                {
                    break;
                }

                again = true;
            }
            else
            {
                p = pNode.Next;
            }
        }
        while (again || p != end);

        return end;
    }

    // ───────────────────────── ear-clipping core ───────────────────────────

    /// <summary>Main ear-slicing loop.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void EarcutLinked(
        int ear,
        TriangleList triangles,
        int dim,
        double minX,
        double minY,
        double invSize,
        int pass,
        ref NodePool pool)
    {
        if (ear == -1)
        {
            return;
        }

        if (pass == 0 && invSize != 0.0)
        {
            IndexCurve(ear, minX, minY, invSize, ref pool);
        }

        int stop = ear;

        while (pool[ear].Prev != pool[ear].Next)
        {
            int prev = pool[ear].Prev;
            int next = pool[ear].Next;

            if (invSize != 0.0
                ? IsEarHashed(ear, minX, minY, invSize, ref pool)
                : IsEar(ear, ref pool))
            {
                // Emit triangle.
                triangles.Add(pool[prev].I, pool[ear].I, pool[next].I);

                RemoveNode(ear, ref pool);

                // Skip the next vertex — produces fewer sliver triangles.
                ear = pool[next].Next;
                stop = pool[next].Next;
                continue;
            }

            ear = next;

            if (ear == stop)
            {
                switch (pass)
                {
                    case 0:
                        EarcutLinked(
                            FilterPoints(ear, ref pool), triangles, dim,
                            minX, minY, invSize, pass: 1, ref pool);
                        break;

                    case 1:
                        ear = CureLocalIntersections(FilterPoints(ear, ref pool), triangles, ref pool);
                        EarcutLinked(
                            ear, triangles, dim,
                            minX, minY, invSize, pass: 2, ref pool);
                        break;

                    case 2:
                        SplitEarcut(ear, triangles, dim, minX, minY, invSize, ref pool);
                        break;
                }

                break;
            }
        }
    }

    /// <summary>Checks whether a polygon node forms a valid ear.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEar(int ear, ref NodePool pool)
    {
        int a = pool[ear].Prev;
        int b = ear;
        int c = pool[ear].Next;

        if (Area(a, b, c, ref pool) >= 0.0)
        {
            return false;
        }

        double ax = pool[a].X, ay = pool[a].Y;
        double bx = pool[b].X, by = pool[b].Y;
        double cx = pool[c].X, cy = pool[c].Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int p = pool[c].Next;

        while (p != a)
        {
            if (pool[p].X >= x0 && pool[p].X <= x1 &&
                pool[p].Y >= y0 && pool[p].Y <= y1 &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, pool[p].X, pool[p].Y) &&
                Area(pool[p].Prev, p, pool[p].Next, ref pool) >= 0.0)
            {
                return false;
            }

            p = pool[p].Next;
        }

        return true;
    }

    /// <summary>Ear check optimised with z-order spatial indexing.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEarHashed(
        int ear,
        double minX,
        double minY,
        double invSize,
        ref NodePool pool)
    {
        int a = pool[ear].Prev;
        int b = ear;
        int c = pool[ear].Next;

        if (Area(a, b, c, ref pool) >= 0.0)
        {
            return false;
        }

        double ax = pool[a].X, ay = pool[a].Y;
        double bx = pool[b].X, by = pool[b].Y;
        double cx = pool[c].X, cy = pool[c].Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int minZ = ZOrder(x0, y0, minX, minY, invSize);
        int maxZ = ZOrder(x1, y1, minX, minY, invSize);

        int p = pool[ear].PrevZ;
        int n = pool[ear].NextZ;

        while (p != -1 && pool[p].Z >= minZ &&
               n != -1 && pool[n].Z <= maxZ)
        {
            if (pool[p].X >= x0 && pool[p].X <= x1 && pool[p].Y >= y0 && pool[p].Y <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, pool[p].X, pool[p].Y) &&
                Area(pool[p].Prev, p, pool[p].Next, ref pool) >= 0.0)
            {
                return false;
            }

            p = pool[p].PrevZ;

            if (pool[n].X >= x0 && pool[n].X <= x1 && pool[n].Y >= y0 && pool[n].Y <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, pool[n].X, pool[n].Y) &&
                Area(pool[n].Prev, n, pool[n].Next, ref pool) >= 0.0)
            {
                return false;
            }

            n = pool[n].NextZ;
        }

        while (p != -1 && pool[p].Z >= minZ)
        {
            if (pool[p].X >= x0 && pool[p].X <= x1 && pool[p].Y >= y0 && pool[p].Y <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, pool[p].X, pool[p].Y) &&
                Area(pool[p].Prev, p, pool[p].Next, ref pool) >= 0.0)
            {
                return false;
            }

            p = pool[p].PrevZ;
        }

        while (n != -1 && pool[n].Z <= maxZ)
        {
            if (pool[n].X >= x0 && pool[n].X <= x1 && pool[n].Y >= y0 && pool[n].Y <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, pool[n].X, pool[n].Y) &&
                Area(pool[n].Prev, n, pool[n].Next, ref pool) >= 0.0)
            {
                return false;
            }

            n = pool[n].NextZ;
        }

        return true;
    }

    // ───────────────────────── intersection cures ──────────────────────────

    /// <summary>
    /// Walks the polygon and fixes small local self-intersections by
    /// emitting a triangle at each crossing.
    /// </summary>
    private static int CureLocalIntersections(int start, TriangleList triangles, ref NodePool pool)
    {
        int p = start;

        do
        {
            int a = pool[p].Prev;
            int pNext = pool[p].Next;
            int b = pool[pNext].Next;

            if (!(pool[a].X == pool[b].X && pool[a].Y == pool[b].Y) &&
                Intersects(a, p, pNext, b, ref pool) &&
                LocallyInside(a, b, ref pool) &&
                LocallyInside(b, a, ref pool))
            {
                triangles.Add(pool[a].I, pool[p].I, pool[b].I);

                RemoveNode(p, ref pool);
                RemoveNode(pNext, ref pool);

                p = start = b;
            }

            p = pool[p].Next;
        }
        while (p != start);

        return FilterPoints(p, ref pool);
    }

    /// <summary>
    /// Last-resort: find a valid diagonal, split the polygon in two, and
    /// triangulate each half independently.
    /// </summary>
    private static void SplitEarcut(
        int start,
        TriangleList triangles,
        int dim,
        double minX,
        double minY,
        double invSize,
        ref NodePool pool)
    {
        int a = start;

        do
        {
            int b = pool[pool[a].Next].Next;

            while (b != pool[a].Prev)
            {
                if (pool[a].I != pool[b].I && IsValidDiagonal(a, b, ref pool))
                {
                    int c = SplitPolygon(a, b, ref pool);

                    a = FilterPoints(a, ref pool, pool[a].Next);
                    c = FilterPoints(c, ref pool, pool[c].Next);

                    EarcutLinked(a, triangles, dim, minX, minY, invSize, pass: 0, ref pool);
                    EarcutLinked(c, triangles, dim, minX, minY, invSize, pass: 0, ref pool);
                    return;
                }

                b = pool[b].Next;
            }

            a = pool[a].Next;
        }
        while (a != start);
    }

    // ────────────────────────── hole elimination ───────────────────────────

    /// <summary>
    /// Links every hole into the outer loop, producing a single-ring polygon
    /// without holes.
    /// </summary>
    private static int EliminateHoles(
        ReadOnlySpan<double> data,
        ReadOnlySpan<int> holeIndices,
        int outerNode,
        int dim,
        ref NodePool pool)
    {
        int queueCount = holeIndices.Length;
        int[]? rentedQueue = queueCount > StackAllocHoleThreshold ? ArrayPool<int>.Shared.Rent(queueCount) : null;
        Span<int> queue = rentedQueue != null
            ? rentedQueue.AsSpan(0, queueCount)
            : stackalloc int[queueCount];

        int actualQueueCount = 0;

        for (int i = 0; i < holeIndices.Length; i++)
        {
            int start = holeIndices[i] * dim;
            int end = i < holeIndices.Length - 1
                ? holeIndices[i + 1] * dim
                : data.Length;

            int list = BuildLinkedList(data, start, end, dim, clockwise: false, ref pool);

            if (list != -1)
            {
                if (list == pool[list].Next)
                {
                    pool[list].Steiner = true;
                }

                queue[actualQueueCount++] = GetLeftmost(list, ref pool);
            }
        }

        // Sort holes by leftmost point X (then Y, then slope) so we process
        // left-to-right. Capture pool by value so the lambda can read node
        // coordinates without a ref parameter. The copy shares the same
        // backing array as the original and is never written to during the sort.
        NodePool poolSnap = pool;
        queue[..actualQueueCount].Sort((a, b) =>
        {
            int cmp = poolSnap[a].X.CompareTo(poolSnap[b].X);
            if (cmp != 0)
            {
                return cmp;
            }

            cmp = poolSnap[a].Y.CompareTo(poolSnap[b].Y);
            if (cmp != 0)
            {
                return cmp;
            }

            double aSlope = (poolSnap[poolSnap[a].Next].Y - poolSnap[a].Y) /
                            (poolSnap[poolSnap[a].Next].X - poolSnap[a].X);
            double bSlope = (poolSnap[poolSnap[b].Next].Y - poolSnap[b].Y) /
                            (poolSnap[poolSnap[b].Next].X - poolSnap[b].X);
            return aSlope.CompareTo(bSlope);
        });

        for (int i = 0; i < actualQueueCount; i++)
        {
            outerNode = EliminateHole(queue[i], outerNode, ref pool);
        }

        if (rentedQueue != null)
        {
            ArrayPool<int>.Shared.Return(rentedQueue);
        }

        return outerNode;
    }

    private static int EliminateHole(int hole, int outerNode, ref NodePool pool)
    {
        int bridge = FindHoleBridge(hole, outerNode, ref pool);

        if (bridge == -1)
        {
            return outerNode;
        }

        int bridgeReverse = SplitPolygon(bridge, hole, ref pool);
        FilterPoints(bridgeReverse, ref pool, pool[bridgeReverse].Next);
        return FilterPoints(bridge, ref pool, pool[bridge].Next);
    }

    /// <summary>
    /// David Eberly's algorithm for finding a bridge between a hole and the
    /// outer polygon.
    /// </summary>
    private static int FindHoleBridge(int hole, int outerNode, ref NodePool pool)
    {
        int p = outerNode;
        double hx = pool[hole].X;
        double hy = pool[hole].Y;
        double qx = double.NegativeInfinity;
        int m = -1;

        if (pool[hole].X == pool[p].X && pool[hole].Y == pool[p].Y)
        {
            return p;
        }

        do
        {
            int pNext = pool[p].Next;

            if (pool[hole].X == pool[pNext].X && pool[hole].Y == pool[pNext].Y)
            {
                return pNext;
            }

            if (hy <= pool[p].Y && hy >= pool[pNext].Y && pool[pNext].Y != pool[p].Y)
            {
                double x = pool[p].X + (hy - pool[p].Y) * (pool[pNext].X - pool[p].X) / (pool[pNext].Y - pool[p].Y);

                if (x <= hx && x > qx)
                {
                    qx = x;
                    m = pool[p].X < pool[pNext].X ? p : pNext;
                    if (x == hx)
                    {
                        return m;
                    }
                }
            }

            p = pool[p].Next;
        }
        while (p != outerNode);

        if (m == -1)
        {
            return -1;
        }

        int stop = m;
        double mx = pool[m].X;
        double my = pool[m].Y;
        double tanMin = double.PositiveInfinity;

        p = m;

        do
        {
            if (hx >= pool[p].X && pool[p].X >= mx && hx != pool[p].X &&
                PointInTriangle(
                    hy < my ? hx : qx, hy,
                    mx, my,
                    hy < my ? qx : hx, hy,
                    pool[p].X, pool[p].Y))
            {
                double tan = Math.Abs(hy - pool[p].Y) / (hx - pool[p].X);

                if (LocallyInside(p, hole, ref pool) &&
                    (tan < tanMin ||
                     (tan == tanMin &&
                      (pool[p].X > pool[m].X ||
                       (pool[p].X == pool[m].X && SectorContainsSector(m, p, ref pool))))))
                {
                    m = p;
                    tanMin = tan;
                }
            }

            p = pool[p].Next;
        }
        while (p != stop);

        return m;
    }

    // ───────────────────── z-order spatial indexing ─────────────────────────

    /// <summary>Interlinks polygon nodes in z-order.</summary>
    private static void IndexCurve(
        int start,
        double minX,
        double minY,
        double invSize,
        ref NodePool pool)
    {
        int p = start;

        do
        {
            if (pool[p].Z == 0)
            {
                pool[p].Z = ZOrder(pool[p].X, pool[p].Y, minX, minY, invSize);
            }

            pool[p].PrevZ = pool[p].Prev;
            pool[p].NextZ = pool[p].Next;
            p = pool[p].Next;
        }
        while (p != start);

        int startPrevZ = pool[start].PrevZ;
        pool[startPrevZ].NextZ = -1;
        pool[start].PrevZ = -1;

        SortLinked(start, ref pool);
    }

    /// <summary>
    /// Simon Tatham's linked-list merge sort.
    /// <see href="http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html"/>
    /// </summary>
    private static int SortLinked(int list, ref NodePool pool)
    {
        int numMerges;
        int inSize = 1;

        do
        {
            int p = list;
            list = -1;
            int tail = -1;
            numMerges = 0;

            while (p != -1)
            {
                numMerges++;
                int q = p;
                int pSize = 0;

                for (int i = 0; i < inSize; i++)
                {
                    pSize++;
                    q = pool[q].NextZ;
                    if (q == -1)
                    {
                        break;
                    }
                }

                int qSize = inSize;

                while (pSize > 0 || (qSize > 0 && q != -1))
                {
                    int e;

                    if (pSize != 0 &&
                        (qSize == 0 || q == -1 || pool[p].Z <= pool[q].Z))
                    {
                        e = p;
                        p = pool[p].NextZ;
                        pSize--;
                    }
                    else
                    {
                        e = q;
                        q = pool[q].NextZ;
                        qSize--;
                    }

                    if (tail != -1)
                    {
                        pool[tail].NextZ = e;
                    }
                    else
                    {
                        list = e;
                    }

                    pool[e].PrevZ = tail;
                    tail = e;
                }

                p = q;
            }

            pool[tail].NextZ = -1;
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
    private static double Area(int p, int q, int r, ref NodePool pool) =>
        (pool[q].Y - pool[p].Y) * (pool[r].X - pool[q].X) -
        (pool[q].X - pool[p].X) * (pool[r].Y - pool[q].Y);

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
    private static bool SectorContainsSector(int m, int p, ref NodePool pool) =>
        Area(pool[m].Prev, m, pool[p].Prev, ref pool) < 0.0 &&
        Area(pool[p].Next, m, pool[m].Next, ref pool) < 0.0;

    private static bool IsValidDiagonal(int a, int b, ref NodePool pool) =>
        pool[pool[a].Next].I != pool[b].I &&
        pool[pool[a].Prev].I != pool[b].I &&
        !IntersectsPolygon(a, b, ref pool) &&
        (LocallyInside(a, b, ref pool) && LocallyInside(b, a, ref pool) && MiddleInside(a, b, ref pool) &&
         (Area(pool[a].Prev, a, pool[b].Prev, ref pool) != 0.0 || Area(a, pool[b].Prev, b, ref pool) != 0.0) ||
         (pool[a].X == pool[b].X && pool[a].Y == pool[b].Y) &&
         Area(pool[a].Prev, a, pool[a].Next, ref pool) > 0.0 &&
         Area(pool[b].Prev, b, pool[b].Next, ref pool) > 0.0);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool Intersects(int p1, int q1, int p2, int q2, ref NodePool pool)
    {
        int o1 = Math.Sign(Area(p1, q1, p2, ref pool));
        int o2 = Math.Sign(Area(p1, q1, q2, ref pool));
        int o3 = Math.Sign(Area(p2, q2, p1, ref pool));
        int o4 = Math.Sign(Area(p2, q2, q1, ref pool));

        return (o1 != o2 && o3 != o4) ||
               (o1 == 0 && OnSegment(p1, p2, q1, ref pool)) ||
               (o2 == 0 && OnSegment(p1, q2, q1, ref pool)) ||
               (o3 == 0 && OnSegment(p2, p1, q2, ref pool)) ||
               (o4 == 0 && OnSegment(p2, q1, q2, ref pool));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool OnSegment(int p, int q, int r, ref NodePool pool) =>
        pool[q].X <= Math.Max(pool[p].X, pool[r].X) && pool[q].X >= Math.Min(pool[p].X, pool[r].X) &&
        pool[q].Y <= Math.Max(pool[p].Y, pool[r].Y) && pool[q].Y >= Math.Min(pool[p].Y, pool[r].Y);

    private static bool IntersectsPolygon(int a, int b, ref NodePool pool)
    {
        int p = a;

        do
        {
            int pNext = pool[p].Next;
            if (pool[p].I != pool[a].I && pool[pNext].I != pool[a].I &&
                pool[p].I != pool[b].I && pool[pNext].I != pool[b].I &&
                Intersects(p, pNext, a, b, ref pool))
            {
                return true;
            }

            p = pNext;
        }
        while (p != a);

        return false;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool LocallyInside(int a, int b, ref NodePool pool) =>
        Area(pool[a].Prev, a, pool[a].Next, ref pool) < 0.0
            ? Area(a, b, pool[a].Next, ref pool) >= 0.0 && Area(a, pool[a].Prev, b, ref pool) >= 0.0
            : Area(a, b, pool[a].Prev, ref pool) < 0.0 || Area(a, pool[a].Next, b, ref pool) < 0.0;

    private static bool MiddleInside(int a, int b, ref NodePool pool)
    {
        int p = a;
        bool inside = false;
        double px = (pool[a].X + pool[b].X) * 0.5;
        double py = (pool[a].Y + pool[b].Y) * 0.5;

        do
        {
            int pNext = pool[p].Next;
            if (((pool[p].Y > py) != (pool[pNext].Y > py)) &&
                pool[pNext].Y != pool[p].Y &&
                px < (pool[pNext].X - pool[p].X) * (py - pool[p].Y) / (pool[pNext].Y - pool[p].Y) + pool[p].X)
            {
                inside = !inside;
            }

            p = pNext;
        }
        while (p != a);

        return inside;
    }

    private static int GetLeftmost(int start, ref NodePool pool)
    {
        int p = start;
        int leftmost = start;

        do
        {
            if (pool[p].X < pool[leftmost].X ||
                (pool[p].X == pool[leftmost].X && pool[p].Y < pool[leftmost].Y))
            {
                leftmost = p;
            }

            p = pool[p].Next;
        }
        while (p != start);

        return leftmost;
    }

    // ──────────────────── linked-list manipulation ─────────────────────────

    private static int SplitPolygon(int a, int b, ref NodePool pool)
    {
        int a2 = pool.Alloc(pool[a].I, pool[a].X, pool[a].Y);
        int b2 = pool.Alloc(pool[b].I, pool[b].X, pool[b].Y);
        int an = pool[a].Next;
        int bp = pool[b].Prev;

        pool[a].Next = b;
        pool[b].Prev = a;

        pool[a2].Next = an;
        pool[an].Prev = a2;

        pool[b2].Next = a2;
        pool[a2].Prev = b2;

        pool[bp].Next = b2;
        pool[b2].Prev = bp;

        return b2;
    }

    private static int InsertNode(int i, double x, double y, int last, ref NodePool pool)
    {
        int p = pool.Alloc(i, x, y);

        if (last == -1)
        {
            pool[p].Prev = p;
            pool[p].Next = p;
        }
        else
        {
            int lastNext = pool[last].Next;
            pool[p].Next = lastNext;
            pool[p].Prev = last;
            pool[lastNext].Prev = p;
            pool[last].Next = p;
        }

        return p;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void RemoveNode(int p, ref NodePool pool)
    {
        pool[pool[p].Next].Prev = pool[p].Prev;
        pool[pool[p].Prev].Next = pool[p].Next;
        if (pool[p].PrevZ != -1) pool[pool[p].PrevZ].NextZ = pool[p].NextZ;
        if (pool[p].NextZ != -1) pool[pool[p].NextZ].PrevZ = pool[p].PrevZ;
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

    /// <summary>
    /// Maximum hole count for which the hole queue is stack-allocated.
    /// Larger counts fall back to <see cref="ArrayPool{T}"/>.
    /// </summary>
    private const int StackAllocHoleThreshold = 64;

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

    // ──────────────────────────── node types ───────────────────────────────

    /// <summary>
    /// Flat value-type representation of a vertex node.
    /// All link fields use -1 as the null sentinel.
    /// </summary>
    /// <remarks>
    /// Fields are ordered with the two <see langword="double"/> coordinates
    /// first to avoid 4 bytes of implicit padding before them.
    /// <see cref="LayoutKind.Sequential"/> with <c>Size = 64</c> pads the
    /// struct to exactly one CPU cache line (64 bytes), which lets the JIT
    /// replace the array-stride multiply <c>idx * 64</c> with a cheap bit
    /// shift <c>idx &lt;&lt; 6</c> in every hot-loop access.
    /// </remarks>
    [StructLayout(LayoutKind.Sequential, Size = 64)]
    private struct NodeData
    {
        public double X;        // X coordinate
        public double Y;        // Y coordinate
        public int    I;        // vertex index in the coordinate array
        public int    Prev;     // index of previous vertex in ring  (-1 = none)
        public int    Next;     // index of next vertex in ring      (-1 = none)
        public int    Z;        // z-order curve value
        public int    PrevZ;    // index of previous node in z-order (-1 = none)
        public int    NextZ;    // index of next node in z-order     (-1 = none)
        public bool   Steiner;  // is this a Steiner point?
        // 23 bytes implicit padding (Size = 64)
    }

    /// <summary>
    /// Pool of <see cref="NodeData"/> values backed by a rented array.
    /// All nodes for one <see cref="Triangulate"/> call live in a single
    /// contiguous allocation.
    /// </summary>
    private struct NodePool : IDisposable
    {
        private NodeData[] _nodes;
        private int        _count;

        public NodePool(int capacity)
        {
            _nodes = ArrayPool<NodeData>.Shared.Rent(capacity);
            _count = 0;
        }

        /// <summary>Allocates one node and returns its index.</summary>
        public int Alloc(int i, double x, double y)
        {
            int idx = _count++;
            ref NodeData n = ref _nodes[idx];
            n.I       = i;
            n.X       = x;
            n.Y       = y;
            n.Prev    = -1;
            n.Next    = -1;
            n.Z       = 0;
            n.PrevZ   = -1;
            n.NextZ   = -1;
            n.Steiner = false;
            return idx;
        }

        /// <summary>Returns a ref to the node at the given index for in-place mutation.</summary>
        public ref NodeData this[int idx] => ref _nodes[idx];

        public void Dispose() => ArrayPool<NodeData>.Shared.Return(_nodes);
    }
}
