// This is an automated Csharp port of https://github.com/mapbox/earcut.
// Copyright 2026 Michael Conrad.
// Licensed under the MIT License.
// See LICENSE file for details.

using System.Runtime.CompilerServices;
using System.Numerics;
using System.Buffers;

[module: SkipLocalsInit]

namespace Earcut;

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
    public static IReadOnlyList<int> Triangulate(
        ReadOnlySpan<double> data,
        ReadOnlySpan<int> holeIndices = default,
        int dim = 2)
    {
        bool hasHoles = !holeIndices.IsEmpty;
        int outerLen = hasHoles ? holeIndices[0] * dim : data.Length;

        Node? outerNode = BuildLinkedList(data, 0, outerLen, dim, clockwise: true);

        if (outerNode is null || outerNode.Next == outerNode.Prev)
        {
            return [];
        }

        // Pre-allocate: a simple polygon with V vertices produces at most V-2 triangles.
        int vertexCount = data.Length / dim;
        var triangles = new List<int>(Math.Max(0, (vertexCount - 2) * 3));

        if (hasHoles)
        {
            outerNode = EliminateHoles(data, holeIndices, outerNode, dim);
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

        EarcutLinked(outerNode, triangles, dim, minX, minY, invSize, pass: 0);

        return triangles;
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

        var vertices = new List<double>(totalCoords);
        var holes = new List<int>(data.Length - 1);
        int holeIndex = 0;
        int prevLen = 0;

        foreach (double[][] ring in data)
        {
            foreach (double[] point in ring)
            {
                for (int d = 0; d < dimensions; d++)
                {
                    vertices.Add(point[d]);
                }
            }

            if (prevLen > 0)
            {
                holeIndex += prevLen;
                holes.Add(holeIndex);
            }

            prevLen = ring.Length;
        }

        return (vertices.ToArray(), holes.ToArray(), dimensions);
    }

    // ──────────────────────── linked-list helpers ───────────────────────────

    /// <summary>
    /// Creates a circular doubly-linked list from polygon points in the
    /// specified winding order.
    /// </summary>
    private static Node? BuildLinkedList(
        ReadOnlySpan<double> data,
        int start,
        int end,
        int dim,
        bool clockwise)
    {
        Node? last = null;

        if (clockwise == (SignedArea(data, start, end, dim) > 0))
        {
            for (int i = start; i < end; i += dim)
            {
                last = InsertNode(i / dim, data[i], data[i + 1], last);
            }
        }
        else
        {
            for (int i = end - dim; i >= start; i -= dim)
            {
                last = InsertNode(i / dim, data[i], data[i + 1], last);
            }
        }

        if (last is not null && CoordsEqual(last, last.Next!))
        {
            RemoveNode(last);
            last = last.Next;
        }

        return last;
    }

    /// <summary>Eliminates collinear or duplicate points.</summary>
    private static Node? FilterPoints(Node? start, Node? end = null)
    {
        if (start is null)
        {
            return null;
        }

        end ??= start;

        Node p = start;
        bool again;

        do
        {
            again = false;

            if (!p.Steiner &&
                (CoordsEqual(p, p.Next!) || Area(p.Prev!, p, p.Next!) == 0.0))
            {
                RemoveNode(p);
                p = end = p.Prev!;
                if (p == p.Next)
                {
                    break;
                }

                again = true;
            }
            else
            {
                p = p.Next!;
            }
        }
        while (again || p != end);

        return end;
    }

    // ───────────────────────── ear-clipping core ───────────────────────────

    /// <summary>Main ear-slicing loop.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void EarcutLinked(
        Node? ear,
        List<int> triangles,
        int dim,
        double minX,
        double minY,
        double invSize,
        int pass)
    {
        if (ear is null)
        {
            return;
        }

        if (pass == 0 && invSize != 0.0)
        {
            IndexCurve(ear, minX, minY, invSize);
        }

        Node stop = ear;

        while (ear!.Prev != ear.Next)
        {
            Node prev = ear.Prev!;
            Node next = ear.Next!;

            if (invSize != 0.0
                ? IsEarHashed(ear, minX, minY, invSize)
                : IsEar(ear))
            {
                // Emit triangle.
                triangles.Add(prev.I);
                triangles.Add(ear.I);
                triangles.Add(next.I);

                RemoveNode(ear);

                // Skip the next vertex — produces fewer sliver triangles.
                ear = next.Next!;
                stop = next.Next!;
                continue;
            }

            ear = next;

            if (ear == stop)
            {
                switch (pass)
                {
                    case 0:
                        EarcutLinked(
                            FilterPoints(ear), triangles, dim,
                            minX, minY, invSize, pass: 1);
                        break;

                    case 1:
                        ear = CureLocalIntersections(FilterPoints(ear)!, triangles);
                        EarcutLinked(
                            ear, triangles, dim,
                            minX, minY, invSize, pass: 2);
                        break;

                    case 2:
                        SplitEarcut(ear, triangles, dim, minX, minY, invSize);
                        break;
                }

                break;
            }
        }
    }

    /// <summary>Checks whether a polygon node forms a valid ear.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool IsEar(Node ear)
    {
        Node a = ear.Prev!;
        Node b = ear;
        Node c = ear.Next!;

        if (Area(a, b, c) >= 0.0)
        {
            return false;
        }

        double ax = a.X, ay = a.Y;
        double bx = b.X, by = b.Y;
        double cx = c.X, cy = c.Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        Node p = c.Next!;

        while (p != a)
        {
            if (p.X >= x0 && p.X <= x1 &&
                p.Y >= y0 && p.Y <= y1 &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, p.X, p.Y) &&
                Area(p.Prev!, p, p.Next!) >= 0.0)
            {
                return false;
            }

            p = p.Next!;
        }

        return true;
    }

    /// <summary>Ear check optimised with z-order spatial indexing.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEarHashed(
        Node ear,
        double minX,
        double minY,
        double invSize)
    {
        Node a = ear.Prev!;
        Node b = ear;
        Node c = ear.Next!;

        if (Area(a, b, c) >= 0.0)
        {
            return false;
        }

        double ax = a.X, ay = a.Y;
        double bx = b.X, by = b.Y;
        double cx = c.X, cy = c.Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int minZ = ZOrder(x0, y0, minX, minY, invSize);
        int maxZ = ZOrder(x1, y1, minX, minY, invSize);

        Node? p = ear.PrevZ;
        Node? n = ear.NextZ;

        while (p is not null && p.Z >= minZ &&
               n is not null && n.Z <= maxZ)
        {
            if (p.X >= x0 && p.X <= x1 && p.Y >= y0 && p.Y <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, p.X, p.Y) &&
                Area(p.Prev!, p, p.Next!) >= 0.0)
            {
                return false;
            }

            p = p.PrevZ;

            if (n.X >= x0 && n.X <= x1 && n.Y >= y0 && n.Y <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, n.X, n.Y) &&
                Area(n.Prev!, n, n.Next!) >= 0.0)
            {
                return false;
            }

            n = n.NextZ;
        }

        while (p is not null && p.Z >= minZ)
        {
            if (p.X >= x0 && p.X <= x1 && p.Y >= y0 && p.Y <= y1 &&
                p != a && p != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, p.X, p.Y) &&
                Area(p.Prev!, p, p.Next!) >= 0.0)
            {
                return false;
            }

            p = p.PrevZ;
        }

        while (n is not null && n.Z <= maxZ)
        {
            if (n.X >= x0 && n.X <= x1 && n.Y >= y0 && n.Y <= y1 &&
                n != a && n != c &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, n.X, n.Y) &&
                Area(n.Prev!, n, n.Next!) >= 0.0)
            {
                return false;
            }

            n = n.NextZ;
        }

        return true;
    }

    // ───────────────────────── intersection cures ──────────────────────────

    /// <summary>
    /// Walks the polygon and fixes small local self-intersections by
    /// emitting a triangle at each crossing.
    /// </summary>
    private static Node CureLocalIntersections(Node start, List<int> triangles)
    {
        Node p = start;

        do
        {
            Node a = p.Prev!;
            Node b = p.Next!.Next!;

            if (!CoordsEqual(a, b) &&
                Intersects(a, p, p.Next!, b) &&
                LocallyInside(a, b) &&
                LocallyInside(b, a))
            {
                triangles.Add(a.I);
                triangles.Add(p.I);
                triangles.Add(b.I);

                RemoveNode(p);
                RemoveNode(p.Next!);

                p = start = b;
            }

            p = p.Next!;
        }
        while (p != start);

        return FilterPoints(p)!;
    }

    /// <summary>
    /// Last-resort: find a valid diagonal, split the polygon in two, and
    /// triangulate each half independently.
    /// </summary>
    private static void SplitEarcut(
        Node start,
        List<int> triangles,
        int dim,
        double minX,
        double minY,
        double invSize)
    {
        Node a = start;

        do
        {
            Node b = a.Next!.Next!;

            while (b != a.Prev)
            {
                if (a.I != b.I && IsValidDiagonal(a, b))
                {
                    Node c = SplitPolygon(a, b);

                    a = FilterPoints(a, a.Next)!;
                    c = FilterPoints(c, c.Next)!;

                    EarcutLinked(a, triangles, dim, minX, minY, invSize, pass: 0);
                    EarcutLinked(c, triangles, dim, minX, minY, invSize, pass: 0);
                    return;
                }

                b = b.Next!;
            }

            a = a.Next!;
        }
        while (a != start);
    }

    // ────────────────────────── hole elimination ───────────────────────────

    /// <summary>
    /// Links every hole into the outer loop, producing a single-ring polygon
    /// without holes.
    /// </summary>
    private static Node EliminateHoles(
        ReadOnlySpan<double> data,
        ReadOnlySpan<int> holeIndices,
        Node outerNode,
        int dim)
    {
        int holeCount = holeIndices.Length;
        
        // Use ArrayPool for the queue
        Node[] queue = ArrayPool<Node>.Shared.Rent(holeCount);

        try
        {
            int queueSize = 0;
            for (int i = 0; i < holeIndices.Length; i++)
            {
                int start = holeIndices[i] * dim;
                int end = i < holeIndices.Length - 1
                    ? holeIndices[i + 1] * dim
                    : data.Length;

                Node? list = BuildLinkedList(data, start, end, dim, clockwise: false);

                if (list is not null)
                {
                    if (list == list.Next)
                    {
                        list.Steiner = true;
                    }

                    queue[queueSize++] = GetLeftmost(list);
                }
            }

            // Sort only the valid portion
            Array.Sort(queue, 0, queueSize, Comparer<Node>.Create(static (a, b) =>
            {
                int cmp = a.X.CompareTo(b.X);
                if (cmp != 0)
                {
                    return cmp;
                }

                cmp = a.Y.CompareTo(b.Y);
                if (cmp != 0)
                {
                    return cmp;
                }

                double aSlope = (a.Next!.Y - a.Y) / (a.Next.X - a.X);
                double bSlope = (b.Next!.Y - b.Y) / (b.Next.X - b.X);
                return aSlope.CompareTo(bSlope);
            }));

            for (int i = 0; i < queueSize; i++)
            {
                outerNode = EliminateHole(queue[i], outerNode);
            }

            return outerNode;
        }
        finally
        {
            ArrayPool<Node>.Shared.Return(queue, clearArray: true);
        }
    }

    private static Node EliminateHole(Node hole, Node outerNode)
    {
        Node? bridge = FindHoleBridge(hole, outerNode);

        if (bridge is null)
        {
            return outerNode;
        }

        Node bridgeReverse = SplitPolygon(bridge, hole);
        FilterPoints(bridgeReverse, bridgeReverse.Next);
        return FilterPoints(bridge, bridge.Next)!;
    }

    /// <summary>
    /// David Eberly's algorithm for finding a bridge between a hole and the
    /// outer polygon.
    /// </summary>
    private static Node? FindHoleBridge(Node hole, Node outerNode)
    {
        Node p = outerNode;
        double hx = hole.X;
        double hy = hole.Y;
        double qx = double.NegativeInfinity;
        Node? m = null;

        if (CoordsEqual(hole, p))
        {
            return p;
        }

        do
        {
            if (CoordsEqual(hole, p.Next!))
            {
                return p.Next;
            }

            if (hy <= p.Y && hy >= p.Next!.Y && p.Next.Y != p.Y)
            {
                double x = p.X + (hy - p.Y) * (p.Next.X - p.X) / (p.Next.Y - p.Y);

                if (x <= hx && x > qx)
                {
                    qx = x;
                    m = p.X < p.Next.X ? p : p.Next;
                    if (x == hx)
                    {
                        return m;
                    }
                }
            }

            p = p.Next!;
        }
        while (p != outerNode);

        if (m is null)
        {
            return null;
        }

        Node stop = m;
        double mx = m.X;
        double my = m.Y;
        double tanMin = double.PositiveInfinity;

        p = m;

        do
        {
            if (hx >= p.X && p.X >= mx && hx != p.X &&
                PointInTriangle(
                    hy < my ? hx : qx, hy,
                    mx, my,
                    hy < my ? qx : hx, hy,
                    p.X, p.Y))
            {
                double tan = Math.Abs(hy - p.Y) / (hx - p.X);

                if (LocallyInside(p, hole) &&
                    (tan < tanMin ||
                     (tan == tanMin &&
                      (p.X > m.X || (p.X == m.X && SectorContainsSector(m, p))))))
                {
                    m = p;
                    tanMin = tan;
                }
            }

            p = p.Next!;
        }
        while (p != stop);

        return m;
    }

    // ───────────────────── z-order spatial indexing ─────────────────────────

    /// <summary>Interlinks polygon nodes in z-order.</summary>
    private static void IndexCurve(
        Node start,
        double minX,
        double minY,
        double invSize)
    {
        Node p = start;

        do
        {
            if (p.Z == 0)
            {
                p.Z = ZOrder(p.X, p.Y, minX, minY, invSize);
            }

            p.PrevZ = p.Prev;
            p.NextZ = p.Next;
            p = p.Next!;
        }
        while (p != start);

        p.PrevZ!.NextZ = null;
        p.PrevZ = null;

        SortLinked(p);
    }

    /// <summary>
    /// Simon Tatham's linked-list merge sort.
    /// <see href="http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html"/>
    /// </summary>
    private static Node SortLinked(Node list)
    {
        int numMerges;
        int inSize = 1;

        do
        {
            Node? p = list;
            list = null!;
            Node? tail = null;
            numMerges = 0;

            while (p is not null)
            {
                numMerges++;
                Node? q = p;
                int pSize = 0;

                for (int i = 0; i < inSize; i++)
                {
                    pSize++;
                    q = q!.NextZ;
                    if (q is null)
                    {
                        break;
                    }
                }

                int qSize = inSize;

                while (pSize > 0 || (qSize > 0 && q is not null))
                {
                    Node e;

                    if (pSize != 0 &&
                        (qSize == 0 || q is null || p!.Z <= q.Z))
                    {
                        e = p!;
                        p = p!.NextZ;
                        pSize--;
                    }
                    else
                    {
                        e = q!;
                        q = q!.NextZ;
                        qSize--;
                    }

                    if (tail is not null)
                    {
                        tail.NextZ = e;
                    }
                    else
                    {
                        list = e;
                    }

                    e.PrevZ = tail;
                    tail = e;
                }

                p = q;
            }

            tail!.NextZ = null;
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
        int ix = (int)((x - minX) * invSize);
        int iy = (int)((y - minY) * invSize);

        ix = (ix | (ix << 8)) & 0x00FF00FF;
        ix = (ix | (ix << 4)) & 0x0F0F0F0F;
        ix = (ix | (ix << 2)) & 0x33333333;
        ix = (ix | (ix << 1)) & 0x55555555;

        iy = (iy | (iy << 8)) & 0x00FF00FF;
        iy = (iy | (iy << 4)) & 0x0F0F0F0F;
        iy = (iy | (iy << 2)) & 0x33333333;
        iy = (iy | (iy << 1)) & 0x55555555;

        return ix | (iy << 1);
    }

    // ──────────────────────── geometry predicates ──────────────────────────

    /// <summary>Signed area of triangle (p → q → r).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Area(Node p, Node q, Node r) =>
        (q.Y - p.Y) * (r.X - q.X) - (q.X - p.X) * (r.Y - q.Y);

    /// <summary>Coordinate equality (not reference equality).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool CoordsEqual(Node a, Node b) =>
        a.X == b.X && a.Y == b.Y;

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
    private static bool SectorContainsSector(Node m, Node p) =>
        Area(m.Prev!, m, p.Prev!) < 0.0 && Area(p.Next!, m, m.Next!) < 0.0;

    private static bool IsValidDiagonal(Node a, Node b) =>
        a.Next!.I != b.I &&
        a.Prev!.I != b.I &&
        !IntersectsPolygon(a, b) &&
        (LocallyInside(a, b) && LocallyInside(b, a) && MiddleInside(a, b) &&
         (Area(a.Prev, a, b.Prev!) != 0.0 || Area(a, b.Prev!, b) != 0.0) ||
         CoordsEqual(a, b) &&
         Area(a.Prev, a, a.Next) > 0.0 &&
         Area(b.Prev!, b, b.Next!) > 0.0);

    private static bool Intersects(Node p1, Node q1, Node p2, Node q2)
    {
        int o1 = Math.Sign(Area(p1, q1, p2));
        int o2 = Math.Sign(Area(p1, q1, q2));
        int o3 = Math.Sign(Area(p2, q2, p1));
        int o4 = Math.Sign(Area(p2, q2, q1));

        if (o1 != o2 && o3 != o4)
        {
            return true;
        }

        if (o1 == 0 && OnSegment(p1, p2, q1))
        {
            return true;
        }

        if (o2 == 0 && OnSegment(p1, q2, q1))
        {
            return true;
        }

        if (o3 == 0 && OnSegment(p2, p1, q2))
        {
            return true;
        }

        if (o4 == 0 && OnSegment(p2, q1, q2))
        {
            return true;
        }

        return false;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool OnSegment(Node p, Node q, Node r) =>
        q.X <= Math.Max(p.X, r.X) && q.X >= Math.Min(p.X, r.X) &&
        q.Y <= Math.Max(p.Y, r.Y) && q.Y >= Math.Min(p.Y, r.Y);

    private static bool IntersectsPolygon(Node a, Node b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        Node p = a;

        do
        {
            if (p.I != a.I && p.Next!.I != a.I &&
                p.I != b.I && p.Next.I != b.I &&
                Intersects(p, p.Next, a, b))
            {
                return true;
            }

            p = p.Next;
        }
        while (p != a);

        return false;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool LocallyInside(Node a, Node b) =>
        Area(a.Prev!, a, a.Next!) < 0.0
            ? Area(a, b, a.Next!) >= 0.0 && Area(a, a.Prev!, b) >= 0.0
            : Area(a, b, a.Prev!) < 0.0 || Area(a, a.Next!, b) < 0.0;

    private static bool MiddleInside(Node a, Node b)
    {
        Node p = a;
        bool inside = false;
        double px = (a.X + b.X) * 0.5;
        double py = (a.Y + b.Y) * 0.5;

        do
        {
            if (((p.Y > py) != (p.Next!.Y > py)) &&
                p.Next.Y != p.Y &&
                px < (p.Next.X - p.X) * (py - p.Y) / (p.Next.Y - p.Y) + p.X)
            {
                inside = !inside;
            }

            p = p.Next;
        }
        while (p != a);

        return inside;
    }

    private static Node GetLeftmost(Node start)
    {
        Node p = start;
        Node leftmost = start;

        do
        {
            if (p.X < leftmost.X ||
                (p.X == leftmost.X && p.Y < leftmost.Y))
            {
                leftmost = p;
            }

            p = p.Next!;
        }
        while (p != start);

        return leftmost;
    }

    // ──────────────────── linked-list manipulation ─────────────────────────

    private static Node SplitPolygon(Node a, Node b)
    {
        var a2 = new Node(a.I, a.X, a.Y);
        var b2 = new Node(b.I, b.X, b.Y);
        Node an = a.Next!;
        Node bp = b.Prev!;

        a.Next = b;
        b.Prev = a;

        a2.Next = an;
        an.Prev = a2;

        b2.Next = a2;
        a2.Prev = b2;

        bp.Next = b2;
        b2.Prev = bp;

        return b2;
    }

    private static Node InsertNode(int i, double x, double y, Node? last)
    {
        var p = new Node(i, x, y);

        if (last is null)
        {
            p.Prev = p;
            p.Next = p;
        }
        else
        {
            p.Next = last.Next;
            p.Prev = last;
            last.Next!.Prev = p;
            last.Next = p;
        }

        return p;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void RemoveNode(Node p)
    {
        p.Next!.Prev = p.Prev;
        p.Prev!.Next = p.Next;

        if (p.PrevZ is not null)
        {
            p.PrevZ.NextZ = p.NextZ;
        }

        if (p.NextZ is not null)
        {
            p.NextZ.PrevZ = p.PrevZ;
        }
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

    // ──────────────────────────── node type ────────────────────────────────

    /// <summary>
    /// Vertex node in the circular doubly-linked polygon ring.
    /// </summary>
    private sealed class Node(int i, double x, double y)
    {
        public readonly int I = i;       // vertex index in the coordinate array
        public readonly double X = x;    // X coordinate
        public readonly double Y = y;    // Y coordinate

        public Node? Prev;               // previous vertex in ring
        public Node? Next;               // next vertex in ring

        public int Z;                    // z-order curve value
        public Node? PrevZ;              // previous node in z-order
        public Node? NextZ;              // next node in z-order

        public bool Steiner;             // is this a Steiner point?
    }
}
