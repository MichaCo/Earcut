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

        var state = new EarcutState();

        Node? outerNode = BuildLinkedList(data, 0, outerLen, dim, clockwise: true);

        if (outerNode is null || ReferenceEquals(outerNode.Next, outerNode.Prev))
        {
            return [];
        }

        // Pre-allocate: a simple polygon with V vertices produces at most V-2 triangles.
        int vertexCount = data.Length / dim;
        var triangles = new TriangleList(Math.Max(0, (vertexCount - 2) * 3));

        if (hasHoles)
        {
            outerNode = EliminateHoles(data, holeIndices, outerNode, dim, state);
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

        EarcutLinked(outerNode, triangles, dim, minX, minY, invSize, state);

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

    /// <summary>
    /// Refine a triangulation toward the constrained Delaunay triangulation by legalizing
    /// every interior edge with Lawson flips — maximizes minimum angle and removes most
    /// slivers. An optional post-pass for <see cref="Triangulate"/> output, or any manifold
    /// triangle-index array indexing into <paramref name="coords"/>.
    /// Adapted from delaunator's edge legalization.
    /// </summary>
    /// <param name="triangles">Triangle indices as returned by <see cref="Triangulate"/>; mutated in place.</param>
    /// <param name="coords">Flat vertex coordinates passed to <see cref="Triangulate"/>.</param>
    /// <param name="dim">Number of coordinates per vertex in <paramref name="coords"/> (default 2).</param>
    public static void Refine(int[] triangles, ReadOnlySpan<double> coords, int dim = 2)
    {
        int n = triangles.Length;
        if (n == 0)
        {
            return;
        }

        int v = coords.Length / dim;
        EnsureRefineScratch(n);

        uint gen = ++s_refineGen;  // bump generation to logically empty the hash

        Span<int> he = s_he!;
        he[..n].Fill(-1);

        Span<int> hTable = s_hTable!;
        Span<uint> hStamp = s_hStamp!;
        int hMask = s_hMask;

        // Build half-edge twins with an undirected-edge hash; consumed slots mark linked pairs.
        for (int e = 0; e < n; e++)
        {
            int a = triangles[e], b = triangles[NextHE(e)];
            int lo = a < b ? a : b;
            int hi = a < b ? b : a;
            int h = (int)(((uint)(hi * v + lo) * 0x9e3779b1u) & (uint)hMask);
            while (hStamp[h] == gen)
            {
                int s = hTable[h];
                // s == -1 marks a consumed slot (a pair already linked) — skip past it
                if (s != -1)
                {
                    int sa = triangles[s], sb = triangles[NextHE(s)];
                    if ((sa == lo && sb == hi) || (sa == hi && sb == lo))
                    {
                        he[e] = s;
                        he[s] = e;
                        hTable[h] = -1; // link, then consume the slot
                        break;
                    }
                }

                h = (h + 1) & hMask;
            }

            if (hStamp[h] != gen)
            {
                hTable[h] = e;
                hStamp[h] = gen;
            }
        }

        Span<int> edgeStack = s_edgeStack!;
        int stackTop = 0;

        for (int e0 = 0; e0 < n; e0++)
        {
            if (he[e0] == -1)
            {
                continue;
            }

            int a = e0;
            while (true)
            {
                int b = he[a];
                int a0 = a - a % 3;
                int ar = a0 + (a + 2) % 3;

                if (b == -1)
                {
                    if (stackTop == 0)
                    {
                        break;
                    }

                    a = edgeStack[--stackTop];
                    continue;
                }

                int b0 = b - b % 3;
                int al = a0 + (a + 1) % 3;
                int bl = b0 + (b + 2) % 3;

                int p0 = triangles[ar], pr = triangles[a], pl = triangles[al], p1 = triangles[bl];

                double x0 = coords[p0 * dim], y0 = coords[p0 * dim + 1];
                double xr = coords[pr * dim], yr = coords[pr * dim + 1];
                double xl = coords[pl * dim], yl = coords[pl * dim + 1];
                double x1 = coords[p1 * dim], y1 = coords[p1 * dim + 1];

                // Both triangles of the flipped diagonal p0-p1 must be CCW (quad must be convex).
                // Flipping a reflex quad would push a triangle outside the polygon.
                bool convex = Orient(x0, y0, xr, yr, x1, y1) > 0.0 &&
                              Orient(x0, y0, x1, y1, xl, yl) > 0.0;

                if (convex && !InCircle(x0, y0, xr, yr, xl, yl, x1, y1))
                {
                    triangles[a] = p1;
                    triangles[b] = p0;

                    int hbl = he[bl], har = he[ar];
                    he[a] = hbl;
                    if (hbl != -1) he[hbl] = a;
                    he[b] = har;
                    if (har != -1) he[har] = b;
                    he[ar] = bl;
                    he[bl] = ar;

                    if (stackTop < edgeStack.Length)
                    {
                        edgeStack[stackTop++] = b0 + (b + 1) % 3;
                    }
                }
                else
                {
                    if (stackTop == 0)
                    {
                        break;
                    }

                    a = edgeStack[--stackTop];
                }
            }
        }
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

        if (last is not null && last == last.Next!)
        {
            RemoveNode(last);
            last = last.Next;
        }

        return last;
    }

    /// <summary>
    /// Eliminates collinear or duplicate points.
    /// When called with the same start and end (or end omitted), sweeps the full ring.
    /// When called with an explicit end different from start, heals only the window around a
    /// bridge/diagonal cut — O(window) instead of O(ring).
    /// </summary>
    private static Node? FilterPoints(Node? start, Node? end, EarcutState state)
    {
        if (start is null)
        {
            return null;
        }

        end ??= start;
        bool full = ReferenceEquals(end, start);

        Node p = start;
        bool again;

        do
        {
            again = false;

            if (!ReferenceEquals(p, p.Next!) &&
                !state.IsSteiner(p) &&
                (p == p.Next! || Area(p.Prev!, p, p.Next!) == 0.0))
            {
                // Pull the stop bound back past the removal
                if (full || ReferenceEquals(p, end))
                {
                    end = p.Prev!;
                }

                state.FilteredOut = true;
                RemoveNode(p, state.BlockIdx);
                p = p.Prev!;  // re-check the predecessor
                again = true;
            }
            else if (full || !ReferenceEquals(p, end))
            {
                p = p.Next!;
                again = !full;  // local heal: keep looping until the sweep reaches end
            }
        }
        while (again || !ReferenceEquals(p, end));

        return end;
    }

    // ───────────────────────── ear-clipping core ───────────────────────────

    /// <summary>Main ear-slicing loop.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void EarcutLinked(
        Node? ear,
        TriangleList triangles,
        int dim,
        double minX,
        double minY,
        double invSize,
        EarcutState state)
    {
        if (ear is null)
        {
            return;
        }

        if (invSize != 0.0)
        {
            IndexCurve(ear, minX, minY, invSize);
        }

        Node stop = ear;
        bool cured = false;

        while (!ReferenceEquals(ear!.Prev, ear.Next))
        {
            Node prev = ear.Prev!;
            Node next = ear.Next!;

            if (Area(prev, ear, next) < 0.0 &&
                (invSize != 0.0
                    ? IsEarHashed(ear, minX, minY, invSize)
                    : IsEar(ear)))
            {
                // Emit triangle.
                triangles.Add(prev.I, ear.I, next.I);

                RemoveNode(ear, state.BlockIdx);

                ear = next;
                stop = next;
                continue;
            }

            ear = next;

            // If we looped through the whole remaining polygon and can't find any more ears
            if (ReferenceEquals(ear, stop))
            {
                // Try filtering collinear/coincident points and slicing again — repeat as long as
                // filtering actually removes nodes, since each removal can expose new ears
                state.FilteredOut = false;
                ear = FilterPoints(ear, null, state)!;
                if (state.FilteredOut)
                {
                    stop = ear;
                    continue;
                }

                // Filtering is exhausted: cure small local self-intersections once, then retry
                if (!cured)
                {
                    ear = CureLocalIntersections(ear, triangles, state);
                    stop = ear;
                    cured = true;
                    continue;
                }

                // As a last resort, try splitting the remaining polygon into two
                SplitEarcut(ear, triangles, dim, minX, minY, invSize, state);
                break;
            }
        }
    }

    /// <summary>Checks whether a polygon node forms a valid ear.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEar(Node ear)
    {
        // Reflex check (area(a, b, c) >= 0) is hoisted into the EarcutLinked caller.
        Node a = ear.Prev!;
        Node b = ear;
        Node c = ear.Next!;

        double ax = a.X, ay = a.Y;
        double bx = b.X, by = b.Y;
        double cx = c.X, cy = c.Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        Node p = c.Next!;

        while (!ReferenceEquals(p, a))
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
        // Reflex check is hoisted into the EarcutLinked caller (see IsEar).
        Node a = ear.Prev!;
        Node b = ear;
        Node c = ear.Next!;

        double ax = a.X, ay = a.Y;
        double bx = b.X, by = b.Y;
        double cx = c.X, cy = c.Y;

        double x0 = Min3(ax, bx, cx);
        double y0 = Min3(ay, by, cy);
        double x1 = Max3(ax, bx, cx);
        double y1 = Max3(ay, by, cy);

        int minZ = ZOrder(x0, y0, minX, minY, invSize);
        int maxZ = ZOrder(x1, y1, minX, minY, invSize);

        // Look for points inside the triangle in decreasing z-order
        Node? p = ear.PrevZ;
        while (p is not null && p.Z >= minZ)
        {
            if (p.X >= x0 && p.X <= x1 && p.Y >= y0 && p.Y <= y1 &&
                !ReferenceEquals(p, c) &&
                PointInTriangleExceptFirst(ax, ay, bx, by, cx, cy, p.X, p.Y) &&
                Area(p.Prev!, p, p.Next!) >= 0.0)
            {
                return false;
            }

            p = p.PrevZ;
        }

        // Look for points in increasing z-order
        Node? n = ear.NextZ;
        while (n is not null && n.Z <= maxZ)
        {
            if (n.X >= x0 && n.X <= x1 && n.Y >= y0 && n.Y <= y1 &&
                !ReferenceEquals(n, c) &&
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
    private static Node CureLocalIntersections(Node start, TriangleList triangles, EarcutState state)
    {
        Node p = start;
        bool cured = false;

        do
        {
            Node a = p.Prev!;
            Node b = p.Next!.Next!;

            if (a != b &&
                Intersects(a, p, p.Next!, b, includeBoundary: false) &&
                LocallyInside(a, b) &&
                LocallyInside(b, a))
            {
                triangles.Add(a.I, p.I, b.I);

                RemoveNode(p, state.BlockIdx);
                RemoveNode(p.Next!, state.BlockIdx);

                p = start = b;
                cured = true;
            }

            p = p.Next!;
        }
        while (!ReferenceEquals(p, start));

        return cured ? FilterPoints(p, null, state)! : p;
    }

    /// <summary>
    /// Last-resort: find a valid diagonal, split the polygon in two, and
    /// triangulate each half independently.
    /// </summary>
    private static void SplitEarcut(
        Node start,
        TriangleList triangles,
        int dim,
        double minX,
        double minY,
        double invSize,
        EarcutState state)
    {
        Node a = start;

        do
        {
            Node b = a.Next!.Next!;

            while (!ReferenceEquals(b, a.Prev))
            {
                if (a.I != b.I && IsValidDiagonal(a, b))
                {
                    Node c = SplitPolygon(a, b);

                    a = FilterPoints(a, a.Next, state)!;
                    c = FilterPoints(c, c.Next, state)!;

                    EarcutLinked(a, triangles, dim, minX, minY, invSize, state);
                    EarcutLinked(c, triangles, dim, minX, minY, invSize, state);
                    return;
                }

                b = b.Next!;
            }

            a = a.Next!;
        }
        while (!ReferenceEquals(a, start));
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
        int dim,
        EarcutState state)
    {
        var queue = new Node[holeIndices.Length];
        int queueCount = 0;

        for (int i = 0; i < holeIndices.Length; i++)
        {
            int start = holeIndices[i] * dim;
            int end = i < holeIndices.Length - 1
                ? holeIndices[i + 1] * dim
                : data.Length;

            Node? list = BuildLinkedList(data, start, end, dim, clockwise: false);

            if (list is not null)
            {
                if (ReferenceEquals(list, list.Next))
                {
                    list.Steiner = true;
                    state.AddSteiner(list);
                }

                queue[queueCount++] = GetLeftmost(list);
            }
        }

        queue.AsSpan(0, queueCount).Sort(static (a, b) =>
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
        });

        // Block-bbox index for FindHoleBridge, grown append-only as holes merge.
        // Seed it with the outer ring, then append each merged hole.
        var blockIdx = new BlockIndex(data.Length / dim, holeIndices.Length);
        state.BlockIdx = blockIdx;
        blockIdx.IndexSegment(outerNode, outerNode);

        for (int i = 0; i < queueCount; i++)
        {
            outerNode = EliminateHole(queue[i], outerNode, state);
        }

        state.BlockIdx = null;

        // Collapse collinear/coincident points across the whole merged ring once before clipping
        return FilterPoints(outerNode, null, state)!;
    }

    private static Node EliminateHole(Node hole, Node outerNode, EarcutState state)
    {
        Node? bridge = FindHoleBridge(hole, outerNode, state.BlockIdx!);

        if (bridge is null)
        {
            return outerNode;
        }

        Node bridgeReverse = SplitPolygon(bridge, hole);

        // Index the merged-in segment before filtering: in ring order the splice runs
        // bridge -> hole -> bridgeReverse -> bridge2 -> (bridge's old next), covering the
        // hole's edges and both new slit edges. FilterPoints below only drops collinear /
        // coincident points, so these bboxes stay valid (conservative) supersets.
        Node bridge2 = bridgeReverse.Next!;
        state.BlockIdx!.IndexSegment(bridge, bridge2.Next!);

        // Heal collinear/coincident points around the two new slit edges
        FilterPoints(bridgeReverse, bridgeReverse.Next, state);
        return FilterPoints(bridge, bridge.Next, state)!;
    }

    /// <summary>
    /// David Eberly's algorithm for finding a bridge between a hole and the
    /// outer polygon, accelerated with a block-bbox index.
    /// </summary>
    private static Node? FindHoleBridge(Node hole, Node outerNode, BlockIndex blockIdx)
    {
        double hx = hole.X;
        double hy = hole.Y;
        double qx = double.NegativeInfinity;
        Node? m = null;

        if (hole == outerNode)
        {
            return outerNode;
        }

        // Scan blocks; skip any whose bbox can't hold a crossing that beats qx and lies left of hx
        for (int b = 0, g = 0; b < blockIdx.NumBlocks; b++, g += 4)
        {
            double[] bbox = blockIdx.BBox;
            if (hy < bbox[g + 1] || hy > bbox[g + 3] || bbox[g] > hx || bbox[g + 2] <= qx)
            {
                continue;
            }

            // Ensure the walk's exclusive bound is live so we don't overrun into other blocks
            Node stop = blockIdx.LiveBlockStop(b);
            Node p = blockIdx.LiveBlockHead(b);

            do
            {
                if (ReferenceEquals(p.Prev!.Next, p)) // skip nodes removed by FilterPoints (stale in the index)
                {
                    if (hole == p.Next!)
                    {
                        return p.Next;
                    }
                    else if (hy <= p.Y && hy >= p.Next!.Y && p.Next.Y != p.Y)
                    {
                        double x = p.X + (hy - p.Y) * (p.Next.X - p.X) / (p.Next.Y - p.Y);
                        if (x <= hx && x > qx)
                        {
                            qx = x;
                            m = p.X < p.Next.X ? p : p.Next;
                            if (x == hx)
                            {
                                return m; // hole touches outer segment; pick leftmost endpoint
                            }
                        }
                    }
                }

                p = p.Next!;
            }
            while (!ReferenceEquals(p, stop));
        }

        if (m is null)
        {
            return null;
        }

        double mx = m.X;
        double my = m.Y;
        double tminY = Math.Min(hy, my);
        double tmaxY = Math.Max(hy, my);
        double tanMin = double.PositiveInfinity;

        // Scan the same blocks; skip any whose bbox can't overlap the triangle's [mx,hx]×[tminY,tmaxY] box
        for (int b = 0, g = 0; b < blockIdx.NumBlocks; b++, g += 4)
        {
            double[] bbox = blockIdx.BBox;
            if (bbox[g + 2] < mx || bbox[g] > hx || bbox[g + 3] < tminY || bbox[g + 1] > tmaxY)
            {
                continue;
            }

            Node stop = blockIdx.LiveBlockStop(b);
            Node p = blockIdx.LiveBlockHead(b);

            do
            {
                if (ReferenceEquals(p.Prev!.Next, p) && hx >= p.X && p.X >= mx && hx != p.X &&
                    PointInTriangle(
                        hy < my ? hx : qx, hy,
                        mx, my,
                        hy < my ? qx : hx, hy,
                        p.X, p.Y))
                {
                    double tan = Math.Abs(hy - p.Y) / (hx - p.X);

                    // If hole point sits on p's horizontal edge (T-junction touch): the bridge runs
                    // along that edge — LocallyInside rejects it as collinear, but it's valid
                    if ((LocallyInside(p, hole) || (p.Y == hy && p.Next!.Y == hy && p.Next.X > hx)) &&
                        (tan < tanMin || (tan == tanMin &&
                            (p.X > m.X || (p.X == m.X && SectorContainsSector(m, p))))))
                    {
                        m = p;
                        tanMin = tan;
                    }
                }

                p = p.Next!;
            }
            while (!ReferenceEquals(p, stop));
        }

        return m;
    }

    // ───────────────────── z-order spatial indexing ─────────────────────────

    /// <summary>
    /// Interlinks polygon nodes in z-order: collect into an array, sort by z, relink.
    /// </summary>
    private static void IndexCurve(
        Node start,
        double minX,
        double minY,
        double invSize)
    {
        // Count nodes and collect them into a flat array for sorting
        int n = 0;
        Node p = start;
        do
        {
            // Always (re)compute z: may still hold a block index left over from EliminateHoles
            p.Z = ZOrder(p.X, p.Y, minX, minY, invSize);
            n++;
            p = p.Next!;
        }
        while (!ReferenceEquals(p, start));

        // Ensure thread-local sort scratch
        Node[] sortArr = EnsureSortScratch(n);

        p = start;
        for (int i = 0; i < n; i++)
        {
            sortArr[i] = p;
            p = p.Next!;
        }

        SortNodes(sortArr, n);

        // Relink in sorted order
        Node? prev = null;
        for (int i = 0; i < n; i++)
        {
            Node node = sortArr[i];
            node.PrevZ = prev;
            if (prev is not null) prev.NextZ = node;
            prev = node;
        }

        prev!.NextZ = null;
        sortArr[0].PrevZ = null;
    }

    /// <summary>
    /// Sorts the first n nodes of the sort scratch array by z-order value.
    /// Uses insertion sort for small arrays (≤32) and LSD radix sort otherwise.
    /// </summary>
    private static void SortNodes(Node[] sortArr, int n)
    {
        if (n <= 32)
        {
            // Insertion sort — cheaper than radix for short arrays
            for (int i = 1; i < n; i++)
            {
                Node node = sortArr[i];
                int z = node.Z;
                int j = i - 1;
                while (j >= 0 && sortArr[j].Z > z)
                {
                    sortArr[j + 1] = sortArr[j];
                    j--;
                }

                sortArr[j + 1] = node;
            }

            return;
        }

        // LSD radix sort in four 8-bit passes (covering z's 30 bits)
        Node[] sortBuf = EnsureSortBuf(n);
        uint[] zArr = EnsureSortZArr(n);
        uint[] zBuf = EnsureSortZBuf(n);
        uint[] counts = GetRadixCounts();

        for (int i = 0; i < n; i++) zArr[i] = (uint)sortArr[i].Z;

        RadixPass(n, sortArr, zArr, sortBuf, zBuf, 0, counts);
        RadixPass(n, sortBuf, zBuf, sortArr, zArr, 8, counts);
        RadixPass(n, sortArr, zArr, sortBuf, zBuf, 16, counts);
        RadixPass(n, sortBuf, zBuf, sortArr, zArr, 24, counts);
    }

    /// <summary>
    /// One LSD radix pass: stably scatter the first n nodes from src to dst,
    /// bucketed by the 8-bit digit of z at the given bit shift.
    /// </summary>
    private static void RadixPass(
        int n,
        Node[] src,
        uint[] srcZ,
        Node[] dst,
        uint[] dstZ,
        int shift,
        uint[] counts)
    {
        Array.Clear(counts, 0, 256);

        for (int i = 0; i < n; i++)
        {
            counts[(srcZ[i] >> shift) & 0xFF]++;
        }

        // Turn per-bucket counts into start offsets (prefix sum)
        uint sum = 0;
        for (int b = 0; b < 256; b++)
        {
            uint c = counts[b];
            counts[b] = sum;
            sum += c;
        }

        for (int i = 0; i < n; i++)
        {
            uint z = srcZ[i];
            int pos = (int)counts[(z >> shift) & 0xFF]++;
            dst[pos] = src[i];
            dstZ[pos] = z;
        }
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
    private static double Area(Node p, Node q, Node r) =>
        (q.Y - p.Y) * (r.X - q.X) - (q.X - p.X) * (r.Y - q.Y);

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

    private static bool IsValidDiagonal(Node a, Node b)
    {
        // Degenerate case: zero-length diagonal between coincident vertices
        bool zeroLength = a == b &&
            Area(a.Prev!, a, a.Next!) > 0.0 &&
            Area(b.Prev!, b, b.Next!) > 0.0;

        return a.Next!.I != b.I &&
            (zeroLength ||
             (LocallyInside(a, b) && LocallyInside(b, a) &&
              (Area(a.Prev!, a, b.Prev!) != 0.0 || Area(a, b.Prev!, b) != 0.0))) &&
            !IntersectsPolygon(a, b) &&
            (zeroLength || MiddleInside(a, b));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool Intersects(Node p1, Node q1, Node p2, Node q2, bool includeBoundary = true)
    {
        double o1 = Area(p1, q1, p2);
        double o2 = Area(p1, q1, q2);
        double o3 = Area(p2, q2, p1);
        double o4 = Area(p2, q2, q1);

        if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
            ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0)))
        {
            return true;
        }

        if (!includeBoundary)
        {
            return false;
        }

        if (o1 == 0.0 && OnSegment(p1, p2, q1)) return true;
        if (o2 == 0.0 && OnSegment(p1, q2, q1)) return true;
        if (o3 == 0.0 && OnSegment(p2, p1, q2)) return true;
        if (o4 == 0.0 && OnSegment(p2, q1, q2)) return true;

        return false;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool OnSegment(Node p, Node q, Node r) =>
        q.X <= Math.Max(p.X, r.X) && q.X >= Math.Min(p.X, r.X) &&
        q.Y <= Math.Max(p.Y, r.Y) && q.Y >= Math.Min(p.Y, r.Y);

    private static bool IntersectsPolygon(Node a, Node b)
    {
        // Diagonal bbox; an edge whose bbox can't overlap it can't intersect it, so
        // skip the orientation test for those (the common case — the diagonal is short)
        double minX = Math.Min(a.X, b.X);
        double maxX = Math.Max(a.X, b.X);
        double minY = Math.Min(a.Y, b.Y);
        double maxY = Math.Max(a.Y, b.Y);

        Node p = a;

        do
        {
            Node nn = p.Next!;
            if ((p.X > maxX && nn.X > maxX) || (p.X < minX && nn.X < minX) ||
                (p.Y > maxY && nn.Y > maxY) || (p.Y < minY && nn.Y < minY))
            {
                p = nn;
                continue;
            }

            if (p.I != a.I && nn.I != a.I &&
                p.I != b.I && nn.I != b.I &&
                Intersects(p, nn, a, b))
            {
                return true;
            }

            p = nn;
        }
        while (!ReferenceEquals(p, a));

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
        while (!ReferenceEquals(p, a));

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
        while (!ReferenceEquals(p, start));

        return leftmost;
    }

    // ──────────────────────── refine geometry helpers ──────────────────────

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Orient(double ax, double ay, double bx, double by, double cx, double cy) =>
        (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);

    /// <summary>
    /// Whether p is inside the circumcircle of triangle (a, b, c).
    /// Sign is negated vs the usual predicate to match earcut's CCW winding.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool InCircle(
        double ax, double ay,
        double bx, double by,
        double cx, double cy,
        double px, double py)
    {
        double dx = ax - px, dy = ay - py;
        double ex = bx - px, ey = by - py;
        double fx = cx - px, fy = cy - py;
        double ap = dx * dx + dy * dy;
        double bp = ex * ex + ey * ey;
        double cp = fx * fx + fy * fy;
        return dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0.0;
    }

    /// <summary>Next half-edge within the same triangle.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int NextHE(int e) => e - e % 3 + (e + 1) % 3;

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
    private static void RemoveNode(Node p, BlockIndex? blockIdx = null)
    {
        p.Next!.Prev = p.Prev;
        p.Prev!.Next = p.Next;
        p.PrevZ?.NextZ = p.NextZ;
        p.NextZ?.PrevZ = p.PrevZ;

        // Keep the hole-bridge index's block bboxes covering the healed prev->next edge
        blockIdx?.GrowBlock(p.Prev!, p.Next!);
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

    // ─────────────────── thread-local sort scratch ─────────────────────────

    [ThreadStatic] private static Node[]? s_sortArr;
    [ThreadStatic] private static Node[]? s_sortBuf;
    [ThreadStatic] private static uint[]? s_sortZArr;
    [ThreadStatic] private static uint[]? s_sortZBuf;
    [ThreadStatic] private static uint[]? s_radixCounts;

    private static Node[] EnsureSortScratch(int n)
    {
        if (s_sortArr is null || s_sortArr.Length < n)
        {
            s_sortArr = new Node[Math.Max(n, 64)];
        }

        return s_sortArr;
    }

    private static Node[] EnsureSortBuf(int n)
    {
        if (s_sortBuf is null || s_sortBuf.Length < n)
        {
            s_sortBuf = new Node[Math.Max(n, 64)];
        }

        return s_sortBuf;
    }

    private static uint[] EnsureSortZArr(int n)
    {
        if (s_sortZArr is null || s_sortZArr.Length < n)
        {
            s_sortZArr = new uint[Math.Max(n, 64)];
        }

        return s_sortZArr;
    }

    private static uint[] EnsureSortZBuf(int n)
    {
        if (s_sortZBuf is null || s_sortZBuf.Length < n)
        {
            s_sortZBuf = new uint[Math.Max(n, 64)];
        }

        return s_sortZBuf;
    }

    private static uint[] GetRadixCounts()
    {
        return s_radixCounts ??= new uint[256];
    }

    // ─────────────────── thread-local refine scratch ───────────────────────

    [ThreadStatic] private static int[]? s_edgeStack;
    [ThreadStatic] private static int[]? s_he;
    [ThreadStatic] private static int[]? s_hTable;
    [ThreadStatic] private static uint[]? s_hStamp;
    [ThreadStatic] private static int s_hMask;
    [ThreadStatic] private static uint s_refineGen;

    private static void EnsureRefineScratch(int n)
    {
        s_edgeStack ??= new int[512];

        if (s_he is null || s_he.Length < n)
        {
            s_he = new int[n];
        }

        // Power-of-two hash table, load factor ≤ 0.25
        int size = 1;
        while (size < n * 4)
        {
            size <<= 1;
        }

        if (s_hTable is null || s_hTable.Length < size)
        {
            s_hTable = new int[size];
            s_hStamp = new uint[size];
        }

        s_hMask = size - 1;
    }

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

    // ─────────────────── per-triangulation mutable state ───────────────────

    /// <summary>
    /// Mutable state for a single triangulation call.
    /// </summary>
    private sealed class EarcutState
    {
        /// <summary>Set by <see cref="FilterPoints"/> when it removes at least one node.</summary>
        public bool FilteredOut;

        /// <summary>
        /// Block-bbox index active during hole elimination; null otherwise.
        /// </summary>
        public BlockIndex? BlockIdx;

        /// <summary>
        /// Steiner nodes (single-vertex holes) preserved through <see cref="FilterPoints"/>.
        /// Kept null when empty so non-steiner inputs pay nothing.
        /// </summary>
        private HashSet<Node>? _steiners;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void AddSteiner(Node n)
        {
            _steiners ??= new HashSet<Node>(ReferenceEqualityComparer.Instance);
            _steiners.Add(n);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool IsSteiner(Node n) =>
            _steiners is not null && _steiners.Contains(n);
    }

    // ─────────────────────── block-bbox index ──────────────────────────────

    /// <summary>
    /// Block-bbox index for <see cref="FindHoleBridge"/> (upstream issue #183):
    /// one [minX,minY,maxX,maxY] bbox per K consecutive ring edges stored in a flat
    /// array, so the leftward-ray scan can skip whole blocks in O(1) instead of
    /// walking the entire merged ring. Grown append-only as holes merge.
    /// </summary>
    private sealed class BlockIndex
    {
        private const int K = 16; // edges per block

        public double[] BBox;
        public int NumBlocks;
        private Node[] _head;
        private Node[] _stop;

        public BlockIndex(int maxNodes, int numHoles)
        {
            // Upper bound: every input node indexed once, +2 bridge nodes per hole, plus a partial
            // trailing block per appended segment (outer ring + one per hole)
            int maxBlocks = (int)Math.Ceiling((double)(maxNodes + 2 * numHoles) / K) + numHoles + 2;
            BBox = new double[maxBlocks * 4];
            _head = new Node[maxBlocks];
            _stop = new Node[maxBlocks];
            NumBlocks = 0;
        }

        /// <summary>
        /// Index the ring run head..stop (exclusive) as ⌈len/K⌉ blocks;
        /// head === stop means the whole ring.
        /// Each block's bbox covers both endpoints of every edge it owns.
        /// </summary>
        public void IndexSegment(Node head, Node stop)
        {
            Node p = head;
            do
            {
                int b = NumBlocks++;
                if (b >= _head.Length)
                {
                    // Grow if needed (safety)
                    int newSize = _head.Length * 2;
                    Array.Resize(ref _head, newSize);
                    Array.Resize(ref _stop, newSize);
                    Array.Resize(ref BBox, newSize * 4);
                }

                _head[b] = p;
                double bMinX = double.PositiveInfinity, bMinY = double.PositiveInfinity;
                double bMaxX = double.NegativeInfinity, bMaxY = double.NegativeInfinity;
                int k = 0;

                do
                {
                    Node c = p.Next!; // edge p->c; bbox must bound both endpoints
                    p.Z = b; // reuse z as the owning block during EliminateHoles (see GrowBlock)

                    if (p.X < bMinX) bMinX = p.X;
                    if (p.X > bMaxX) bMaxX = p.X;
                    if (p.Y < bMinY) bMinY = p.Y;
                    if (p.Y > bMaxY) bMaxY = p.Y;
                    if (c.X < bMinX) bMinX = c.X;
                    if (c.X > bMaxX) bMaxX = c.X;
                    if (c.Y < bMinY) bMinY = c.Y;
                    if (c.Y > bMaxY) bMaxY = c.Y;

                    p = c;
                }
                while (++k < K && !ReferenceEquals(p, stop));

                _stop[b] = p;
                int g = b * 4;
                BBox[g] = bMinX;
                BBox[g + 1] = bMinY;
                BBox[g + 2] = bMaxX;
                BBox[g + 3] = bMaxY;
            }
            while (!ReferenceEquals(p, stop));
        }

        /// <summary>
        /// When <see cref="FilterPoints"/> heals an edge head→tail (removing the collinear node
        /// between them), grow head's block bbox to cover tail so the leftward-ray prune can't
        /// false-skip it.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void GrowBlock(Node head, Node tail)
        {
            int g = head.Z * 4;
            if (g < 0 || g + 3 >= BBox.Length) return;
            if (tail.X < BBox[g]) BBox[g] = tail.X;
            if (tail.Y < BBox[g + 1]) BBox[g + 1] = tail.Y;
            if (tail.X > BBox[g + 2]) BBox[g + 2] = tail.X;
            if (tail.Y > BBox[g + 3]) BBox[g + 3] = tail.Y;
        }

        /// <summary>
        /// The block's stop node can be stale (removed by <see cref="FilterPoints"/>);
        /// advance to the next live node.
        /// </summary>
        public Node LiveBlockStop(int b)
        {
            Node stop = _stop[b];
            while (!ReferenceEquals(stop.Prev!.Next, stop)) stop = stop.Next!;
            _stop[b] = stop;
            return stop;
        }

        /// <summary>
        /// The block's head node can be removed by <see cref="FilterPoints"/> during merges;
        /// advance to the next live node.
        /// </summary>
        public Node LiveBlockHead(int b)
        {
            Node head = _head[b];
            while (!ReferenceEquals(head.Prev!.Next, head)) head = head.Next!;
            _head[b] = head;
            return head;
        }
    }

    // ──────────────────────────── node type ────────────────────────────────

    /// <summary>
    /// Vertex node in the circular doubly-linked polygon ring.
    /// </summary>
    private sealed class Node(int i, double x, double y) : IEquatable<Node>
    {
        public int I { get; } = i;          // vertex index in the coordinate array

        public double X { get; } = x;       // X coordinate

        public double Y { get; } = y;       // Y coordinate

        public Node? Prev { get; set; }     // previous vertex in ring

        public Node? Next { get; set; }     // next vertex in ring

        public int Z { get; set; }          // z-order curve value; doubles as block index during EliminateHoles

        public Node? PrevZ { get; set; }    // previous node in z-order

        public Node? NextZ { get; set; }    // next node in z-order

        public bool Steiner { get; set; }   // is this a Steiner point?

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Equals(Node? other) => other is not null && X == other.X && Y == other.Y;

        public override bool Equals(object? obj) => obj is Node other && Equals(other);

        public override int GetHashCode() => HashCode.Combine(X, Y);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Node? left, Node? right) =>
            left is null ? right is null : left.Equals(right);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Node? left, Node? right) => !(left == right);
    }
}
