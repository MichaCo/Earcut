// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Buffers;
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

        // Every hole-bridge split adds 2 cloned nodes (SplitPolygon); extra headroom.
        int nodeCapacity = vertexCount + holeIndices.Length * 2 + 2;
        var pool = new NodePool(nodeCapacity);
        try
        {
            int outerNode = BuildLinkedList(data, 0, outerLen, dim, clockwise: true, ref pool);

            if (outerNode == -1 || pool.Next[outerNode] == pool.Prev[outerNode])
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
            int lastNext = pool.Next[last];
            if (pool.X[last] == pool.X[lastNext] && pool.Y[last] == pool.Y[lastNext])
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

        double[] x = pool.X, y = pool.Y;
        int[]    prev = pool.Prev, next = pool.Next;
        bool[]   steiner = pool.Steiner;

        int p = start;
        bool again;

        do
        {
            again = false;
            int pNext = next[p];

            if (!steiner[p] &&
                (x[pNext] == x[p] && y[pNext] == y[p] ||
                 Area(x, y, prev[p], p, pNext) == 0.0))
            {
                int pPrev = prev[p];
                RemoveNode(p, ref pool);
                p = end = pPrev;
                if (p == next[p])
                {
                    break;
                }

                again = true;
            }
            else
            {
                p = pNext;
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

        int[]    prev = pool.Prev, next = pool.Next;
        int[]    nodeI = pool.I;

        int stop = ear;

        while (true)
        {
            int prev_ = prev[ear];
            int next_ = next[ear];

            if (prev_ == next_)
            {
                break;
            }

            if (invSize != 0.0
                ? IsEarHashed(ear, minX, minY, invSize, ref pool)
                : IsEar(ear, ref pool))
            {
                // Emit triangle.
                triangles.Add(nodeI[prev_], nodeI[ear], nodeI[next_]);

                RemoveNode(ear, ref pool);

                // Skip the next vertex — produces fewer sliver triangles.
                ear = next[next_];
                stop = next[next_];
                continue;
            }

            ear = next_;

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
        double[] x = pool.X, y = pool.Y;
        int[]    prev = pool.Prev, next = pool.Next;

        int a = prev[ear];
        int c = next[ear];

        if (Area(x, y, a, ear, c) >= 0.0)
        {
            return false;
        }

        double ax = x[a], ay = y[a];
        double bx = x[ear], by = y[ear];
        double cx = x[c], cy = y[c];

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int p = next[c];

        while (p != a)
        {
            double px = x[p], py = y[p];

            if (px >= x0 && px <= x1 &&
                py >= y0 && py <= y1 &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, px, py) &&
                Area(x, y, prev[p], p, next[p]) >= 0.0)
            {
                return false;
            }

            p = next[p];
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
        // Cache array refs as locals so the JIT can keep them in registers.
        double[] x     = pool.X,    y     = pool.Y;
        int[]    prev  = pool.Prev, next  = pool.Next;
        int[]    z     = pool.Z,    prevZ = pool.PrevZ, nextZ = pool.NextZ;

        int a = prev[ear];
        int c = next[ear];

        if (Area(x, y, a, ear, c) >= 0.0)
        {
            return false;
        }

        double ax = x[a], ay = y[a];
        double bx = x[ear], by = y[ear];
        double cx = x[c], cy = y[c];

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int minZ = ZOrder(x0, y0, minX, minY, invSize);
        int maxZ = ZOrder(x1, y1, minX, minY, invSize);

        int p = prevZ[ear];
        int n = nextZ[ear];

        while (p != -1 && z[p] >= minZ &&
               n != -1 && z[n] <= maxZ)
        {
            double px = x[p], py = y[p];
            if (px >= x0 && px <= x1 && py >= y0 && py <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, px, py) &&
                Area(x, y, prev[p], p, next[p]) >= 0.0)
            {
                return false;
            }

            p = prevZ[p];

            double nx = x[n], ny = y[n];
            if (nx >= x0 && nx <= x1 && ny >= y0 && ny <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, nx, ny) &&
                Area(x, y, prev[n], n, next[n]) >= 0.0)
            {
                return false;
            }

            n = nextZ[n];
        }

        while (p != -1 && z[p] >= minZ)
        {
            double px = x[p], py = y[p];
            if (px >= x0 && px <= x1 && py >= y0 && py <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, px, py) &&
                Area(x, y, prev[p], p, next[p]) >= 0.0)
            {
                return false;
            }

            p = prevZ[p];
        }

        while (n != -1 && z[n] <= maxZ)
        {
            double nx = x[n], ny = y[n];
            if (nx >= x0 && nx <= x1 && ny >= y0 && ny <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, nx, ny) &&
                Area(x, y, prev[n], n, next[n]) >= 0.0)
            {
                return false;
            }

            n = nextZ[n];
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
        double[] x = pool.X, y = pool.Y;
        int[]    nodeI = pool.I, prev = pool.Prev, next = pool.Next;

        int p = start;

        do
        {
            int a     = prev[p];
            int pNext = next[p];
            int b     = next[pNext];

            if (!(x[a] == x[b] && y[a] == y[b]) &&
                Intersects(a, p, pNext, b, ref pool) &&
                LocallyInside(a, b, ref pool) &&
                LocallyInside(b, a, ref pool))
            {
                triangles.Add(nodeI[a], nodeI[p], nodeI[b]);

                RemoveNode(p, ref pool);
                RemoveNode(pNext, ref pool);

                p = start = b;
            }

            p = next[p];
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
        int[] prev = pool.Prev, next = pool.Next;

        int a = start;

        do
        {
            int b = next[next[a]];

            while (b != prev[a])
            {
                if (pool.I[a] != pool.I[b] && IsValidDiagonal(a, b, ref pool))
                {
                    int c = SplitPolygon(a, b, ref pool);

                    a = FilterPoints(a, ref pool, next[a]);
                    c = FilterPoints(c, ref pool, next[c]);

                    EarcutLinked(a, triangles, dim, minX, minY, invSize, pass: 0, ref pool);
                    EarcutLinked(c, triangles, dim, minX, minY, invSize, pass: 0, ref pool);
                    return;
                }

                b = next[b];
            }

            a = next[a];
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
        int[]? rentedQueue = queueCount > StackAllocHoleThreshold
            ? ArrayPool<int>.Shared.Rent(queueCount)
            : null;
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
                if (list == pool.Next[list])
                {
                    pool.Steiner[list] = true;
                }

                queue[actualQueueCount++] = GetLeftmost(list, ref pool);
            }
        }

        // Sort holes by leftmost point X (then Y, then slope) so we process
        // left-to-right. Capture array refs (reference types) rather than
        // copying the pool struct.
        double[] snapX = pool.X, snapY = pool.Y;
        int[]    snapNext = pool.Next;
        queue[..actualQueueCount].Sort((a, b) =>
        {
            int cmp = snapX[a].CompareTo(snapX[b]);
            if (cmp != 0)
            {
                return cmp;
            }

            cmp = snapY[a].CompareTo(snapY[b]);
            if (cmp != 0)
            {
                return cmp;
            }

            double aSlope = (snapY[snapNext[a]] - snapY[a]) / (snapX[snapNext[a]] - snapX[a]);
            double bSlope = (snapY[snapNext[b]] - snapY[b]) / (snapX[snapNext[b]] - snapX[b]);
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
        FilterPoints(bridgeReverse, ref pool, pool.Next[bridgeReverse]);
        return FilterPoints(bridge, ref pool, pool.Next[bridge]);
    }

    /// <summary>
    /// David Eberly's algorithm for finding a bridge between a hole and the
    /// outer polygon.
    /// </summary>
    private static int FindHoleBridge(int hole, int outerNode, ref NodePool pool)
    {
        double[] x = pool.X, y = pool.Y;
        int[]    next = pool.Next;

        int    p  = outerNode;
        double hx = x[hole];
        double hy = y[hole];
        double qx = double.NegativeInfinity;
        int    m  = -1;

        if (x[hole] == x[p] && y[hole] == y[p])
        {
            return p;
        }

        do
        {
            int pNext = next[p];

            if (x[hole] == x[pNext] && y[hole] == y[pNext])
            {
                return pNext;
            }

            if (hy <= y[p] && hy >= y[pNext] && y[pNext] != y[p])
            {
                double xi = x[p] + (hy - y[p]) * (x[pNext] - x[p]) / (y[pNext] - y[p]);

                if (xi <= hx && xi > qx)
                {
                    qx = xi;
                    m  = x[p] < x[pNext] ? p : pNext;
                    if (xi == hx)
                    {
                        return m;
                    }
                }
            }

            p = pNext;
        }
        while (p != outerNode);

        if (m == -1)
        {
            return -1;
        }

        int    stop   = m;
        double mx     = x[m];
        double my     = y[m];
        double tanMin = double.PositiveInfinity;

        p = m;

        do
        {
            double px = x[p], py = y[p];

            if (hx >= px && px >= mx && hx != px &&
                PointInTriangle(
                    hy < my ? hx : qx, hy,
                    mx, my,
                    hy < my ? qx : hx, hy,
                    px, py))
            {
                double tan = Math.Abs(hy - py) / (hx - px);

                if (LocallyInside(p, hole, ref pool) &&
                    (tan < tanMin ||
                     (tan == tanMin &&
                      (px > x[m] || (px == x[m] && SectorContainsSector(m, p, ref pool))))))
                {
                    m      = p;
                    tanMin = tan;
                }
            }

            p = next[p];
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
        double[] x     = pool.X,    y     = pool.Y;
        int[]    prev  = pool.Prev, next  = pool.Next;
        int[]    z     = pool.Z,    prevZ = pool.PrevZ, nextZ = pool.NextZ;

        int p = start;

        do
        {
            if (z[p] == 0)
            {
                z[p] = ZOrder(x[p], y[p], minX, minY, invSize);
            }

            prevZ[p] = prev[p];
            nextZ[p] = next[p];
            p = next[p];
        }
        while (p != start);

        int startPrevZ = prevZ[start];
        nextZ[startPrevZ] = -1;
        prevZ[start]      = -1;

        SortLinked(start, ref pool);
    }

    /// <summary>
    /// Simon Tatham's linked-list merge sort.
    /// <see href="http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html"/>
    /// </summary>
    private static int SortLinked(int list, ref NodePool pool)
    {
        int[] z     = pool.Z;
        int[] prevZ = pool.PrevZ;
        int[] nextZ = pool.NextZ;

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
                int q     = p;
                int pSize = 0;

                for (int i = 0; i < inSize; i++)
                {
                    pSize++;
                    q = nextZ[q];
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
                        (qSize == 0 || q == -1 || z[p] <= z[q]))
                    {
                        e = p;
                        p = nextZ[p];
                        pSize--;
                    }
                    else
                    {
                        e = q;
                        q = nextZ[q];
                        qSize--;
                    }

                    if (tail != -1)
                    {
                        nextZ[tail] = e;
                    }
                    else
                    {
                        list = e;
                    }

                    prevZ[e] = tail;
                    tail     = e;
                }

                p = q;
            }

            nextZ[tail] = -1;
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
    private static double Area(double[] x, double[] y, int p, int q, int r) =>
        (y[q] - y[p]) * (x[r] - x[q]) - (x[q] - x[p]) * (y[r] - y[q]);

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
    private static bool SectorContainsSector(int m, int p, ref NodePool pool)
    {
        double[] x = pool.X, y = pool.Y;
        int[]    prev = pool.Prev, next = pool.Next;
        return Area(x, y, prev[m], m, prev[p]) < 0.0 &&
               Area(x, y, next[p], m, next[m]) < 0.0;
    }

    private static bool IsValidDiagonal(int a, int b, ref NodePool pool)
    {
        double[] x = pool.X, y = pool.Y;
        int[]    nodeI = pool.I, prev = pool.Prev, next = pool.Next;
        return nodeI[next[a]] != nodeI[b] &&
               nodeI[prev[a]] != nodeI[b] &&
               !IntersectsPolygon(a, b, ref pool) &&
               (LocallyInside(a, b, ref pool) && LocallyInside(b, a, ref pool) && MiddleInside(a, b, ref pool) &&
                (Area(x, y, prev[a], a, prev[b]) != 0.0 || Area(x, y, a, prev[b], b) != 0.0) ||
                (x[a] == x[b] && y[a] == y[b]) &&
                Area(x, y, prev[a], a, next[a]) > 0.0 &&
                Area(x, y, prev[b], b, next[b]) > 0.0);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool Intersects(int p1, int q1, int p2, int q2, ref NodePool pool)
    {
        double[] x = pool.X, y = pool.Y;
        int o1 = Math.Sign(Area(x, y, p1, q1, p2));
        int o2 = Math.Sign(Area(x, y, p1, q1, q2));
        int o3 = Math.Sign(Area(x, y, p2, q2, p1));
        int o4 = Math.Sign(Area(x, y, p2, q2, q1));

        return (o1 != o2 && o3 != o4) ||
               (o1 == 0 && OnSegment(x, y, p1, p2, q1)) ||
               (o2 == 0 && OnSegment(x, y, p1, q2, q1)) ||
               (o3 == 0 && OnSegment(x, y, p2, p1, q2)) ||
               (o4 == 0 && OnSegment(x, y, p2, q1, q2));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool OnSegment(double[] x, double[] y, int p, int q, int r) =>
        x[q] <= Math.Max(x[p], x[r]) && x[q] >= Math.Min(x[p], x[r]) &&
        y[q] <= Math.Max(y[p], y[r]) && y[q] >= Math.Min(y[p], y[r]);

    private static bool IntersectsPolygon(int a, int b, ref NodePool pool)
    {
        int[]  nodeI = pool.I, next = pool.Next;

        int p = a;

        do
        {
            int pNext = next[p];
            if (nodeI[p] != nodeI[a] && nodeI[pNext] != nodeI[a] &&
                nodeI[p] != nodeI[b] && nodeI[pNext] != nodeI[b] &&
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
    private static bool LocallyInside(int a, int b, ref NodePool pool)
    {
        double[] x = pool.X, y = pool.Y;
        int[]    prev = pool.Prev, next = pool.Next;
        return Area(x, y, prev[a], a, next[a]) < 0.0
            ? Area(x, y, a, b, next[a]) >= 0.0 && Area(x, y, a, prev[a], b) >= 0.0
            : Area(x, y, a, b, prev[a]) < 0.0   || Area(x, y, a, next[a], b) < 0.0;
    }

    private static bool MiddleInside(int a, int b, ref NodePool pool)
    {
        double[] x = pool.X, y = pool.Y;
        int[]    next = pool.Next;

        int    p      = a;
        bool   inside = false;
        double px     = (x[a] + x[b]) * 0.5;
        double py     = (y[a] + y[b]) * 0.5;

        do
        {
            int pNext = next[p];
            if (((y[p] > py) != (y[pNext] > py)) &&
                y[pNext] != y[p] &&
                px < (x[pNext] - x[p]) * (py - y[p]) / (y[pNext] - y[p]) + x[p])
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
        double[] x = pool.X, y = pool.Y;
        int[]    next = pool.Next;

        int p        = start;
        int leftmost = start;

        do
        {
            if (x[p] < x[leftmost] || (x[p] == x[leftmost] && y[p] < y[leftmost]))
            {
                leftmost = p;
            }

            p = next[p];
        }
        while (p != start);

        return leftmost;
    }

    // ──────────────────── linked-list manipulation ─────────────────────────

    private static int SplitPolygon(int a, int b, ref NodePool pool)
    {
        int[] prev = pool.Prev, next = pool.Next;

        int a2 = pool.Alloc(pool.I[a], pool.X[a], pool.Y[a]);
        int b2 = pool.Alloc(pool.I[b], pool.X[b], pool.Y[b]);
        int an = next[a];
        int bp = prev[b];

        next[a]  = b;
        prev[b]  = a;

        next[a2] = an;
        prev[an] = a2;

        next[b2] = a2;
        prev[a2] = b2;

        next[bp] = b2;
        prev[b2] = bp;

        return b2;
    }

    private static int InsertNode(int i, double x, double y, int last, ref NodePool pool)
    {
        int p = pool.Alloc(i, x, y);

        if (last == -1)
        {
            pool.Prev[p] = p;
            pool.Next[p] = p;
        }
        else
        {
            int lastNext = pool.Next[last];
            pool.Next[p]        = lastNext;
            pool.Prev[p]        = last;
            pool.Prev[lastNext] = p;
            pool.Next[last]     = p;
        }

        return p;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void RemoveNode(int p, ref NodePool pool)
    {
        int[]  prev = pool.Prev, next = pool.Next;
        int[]  prevZ = pool.PrevZ, nextZ = pool.NextZ;

        int pPrev = prev[p];
        int pNext = next[p];
        next[pPrev] = pNext;
        prev[pNext] = pPrev;

        int pPrevZ = prevZ[p];
        int pNextZ = nextZ[p];
        if (pPrevZ != -1) nextZ[pPrevZ] = pNextZ;
        if (pNextZ != -1) prevZ[pNextZ] = pPrevZ;
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

            _buffer[_count]     = a;
            _buffer[_count + 1] = b;
            _buffer[_count + 2] = c;
            _count += 3;
        }

        public int[] ToArray() => _count == _buffer.Length ? _buffer : _buffer[.._count];
    }

    // ──────────────────────────── node pool ────────────────────────────────

    /// <summary>
    /// Structure-of-Arrays node pool for one <see cref="Triangulate"/> call.
    /// </summary>
    /// <remarks>
    /// Separating the immutable coordinate arrays (<see cref="X"/>, <see cref="Y"/>,
    /// <see cref="I"/>) from the mutable link arrays (<see cref="Prev"/>,
    /// <see cref="Next"/>, …) means that writes to link data during
    /// <c>RemoveNode</c> / <c>InsertNode</c> / <c>SplitPolygon</c> can never
    /// dirty the cache lines holding coordinate data.  The coordinate arrays
    /// remain hot in L1/L2 cache throughout the inner loops of
    /// <c>IsEar</c> / <c>IsEarHashed</c>, eliminating the false-sharing
    /// cache thrash that occurred with the previous Array-of-Structures layout.
    /// </remarks>
    private struct NodePool : IDisposable
    {
        // Immutable coordinate data — set once in Alloc, never written again.
        internal double[] X;
        internal double[] Y;
        internal int[]    I;

        // Mutable link data — frequently updated by linked-list operations.
        internal int[]  Prev;
        internal int[]  Next;
        internal int[]  Z;
        internal int[]  PrevZ;
        internal int[]  NextZ;
        internal bool[] Steiner;

        private int _count;

        public NodePool(int capacity)
        {
            X       = ArrayPool<double>.Shared.Rent(capacity);
            Y       = ArrayPool<double>.Shared.Rent(capacity);
            I       = ArrayPool<int>.Shared.Rent(capacity);
            Prev    = ArrayPool<int>.Shared.Rent(capacity);
            Next    = ArrayPool<int>.Shared.Rent(capacity);
            Z       = ArrayPool<int>.Shared.Rent(capacity);
            PrevZ   = ArrayPool<int>.Shared.Rent(capacity);
            NextZ   = ArrayPool<int>.Shared.Rent(capacity);
            Steiner = ArrayPool<bool>.Shared.Rent(capacity);
            _count  = 0;
        }

        /// <summary>Allocates one node and returns its index.</summary>
        public int Alloc(int i, double x, double y)
        {
            int idx   = _count++;
            X[idx]       = x;
            Y[idx]       = y;
            I[idx]       = i;
            Prev[idx]    = -1;
            Next[idx]    = -1;
            Z[idx]       = 0;
            PrevZ[idx]   = -1;
            NextZ[idx]   = -1;
            Steiner[idx] = false;
            return idx;
        }

        public void Dispose()
        {
            ArrayPool<double>.Shared.Return(X);
            ArrayPool<double>.Shared.Return(Y);
            ArrayPool<int>.Shared.Return(I);
            ArrayPool<int>.Shared.Return(Prev);
            ArrayPool<int>.Shared.Return(Next);
            ArrayPool<int>.Shared.Return(Z);
            ArrayPool<int>.Shared.Return(PrevZ);
            ArrayPool<int>.Shared.Return(NextZ);
            ArrayPool<bool>.Shared.Return(Steiner);
        }
    }
}
