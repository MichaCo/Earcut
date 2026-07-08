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

        Node? outerNode = BuildLinkedList(data, 0, outerLen, dim, clockwise: true);

        if (outerNode is null || ReferenceEquals(outerNode.Next, outerNode.Prev))
        {
            return [];
        }

        // Pre-allocate: a simple polygon with V vertices produces at most V-2 triangles.
        int vertexCount = data.Length / dim;
        var triangles = new TriangleList(Math.Max(0, (vertexCount - 2) * 3));
        var workspace = new TriangulationWorkspace();

        if (hasHoles)
        {
            outerNode = EliminateHoles(data, holeIndices, outerNode, dim, workspace);
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

        EarcutLinked(outerNode, triangles, minX, minY, invSize, workspace);

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
    /// Refine a triangulation toward the constrained Delaunay triangulation by
    /// legalizing every interior edge in place with Lawson flips — maximising the
    /// minimum angle and removing most slivers. Mutates <paramref name="triangles"/>
    /// in place. An optional post-pass for <see cref="Triangulate"/> output, or any
    /// manifold triangle-index array. Adapted from delaunator's edge legalization.
    /// </summary>
    /// <param name="triangles">Triangle indices as returned by <see cref="Triangulate"/>; mutated in place.</param>
    /// <param name="coords">The flat vertex coordinates passed to <see cref="Triangulate"/>.</param>
    /// <param name="dim">Number of coordinates per vertex in <paramref name="coords"/> (default 2).</param>
    public static void Refine(int[] triangles, ReadOnlySpan<double> coords, int dim = 2)
    {
        int n = triangles.Length;
        if (n < 6)
        {
            return;
        }

        var workspace = new RefineWorkspace(n);
        const uint generation = 1;

        int[] halfEdges = workspace.HalfEdges;
        halfEdges.AsSpan(0, n).Fill(-1);
        int[] hashTable = workspace.HashTable;
        uint[] hashStamp = workspace.HashStamp;
        byte[] edgeStamp = workspace.EdgeStamp;
        int[] edgeStack = workspace.EdgeStack;
        int hashMask = workspace.HashMask;

        // Build half-edge twins with an undirected-edge hash; seed interior edges onto the stack.
        int stackTop = 0;
        for (int e = 0; e < n; e++)
        {
            int a = triangles[e];
            int b = triangles[NextHE(e)];
            int lo = a < b ? a : b;
            int hi = a < b ? b : a;
            int h = (int)((unchecked((uint)lo * 0x9e3779b1u) ^ unchecked((uint)hi * 0x85ebca6bu)) & (uint)hashMask);

            while (hashStamp[h] == generation)
            {
                int s = hashTable[h];
                if (s != -1)
                {
                    int sa = triangles[s];
                    int sb = triangles[NextHE(s)];
                    if ((sa == lo && sb == hi) || (sa == hi && sb == lo))
                    {
                        halfEdges[e] = s;
                        halfEdges[s] = e;
                        hashTable[h] = -1;
                        edgeStamp[s] = 1;
                        edgeStack[stackTop++] = s;
                        break;
                    }
                }

                h = (h + 1) & hashMask;
            }

            if (hashStamp[h] != generation)
            {
                hashTable[h] = e;
                hashStamp[h] = generation;
            }
        }

        while (stackTop > 0)
        {
            int a = edgeStack[--stackTop];
            edgeStamp[a] = 0;
            int b = halfEdges[a];
            if (b == -1)
            {
                continue;
            }

            int a0 = a - a % 3;
            int b0 = b - b % 3;
            int ar = a0 + (a + 2) % 3;
            int al = a0 + (a + 1) % 3;
            int bl = b0 + (b + 2) % 3;
            int br = b0 + (b + 1) % 3;
            int p0 = triangles[ar];
            int pr = triangles[a];
            int pl = triangles[al];
            int p1 = triangles[bl];

            double x0 = coords[p0 * dim];
            double y0 = coords[p0 * dim + 1];
            double xr = coords[pr * dim];
            double yr = coords[pr * dim + 1];
            double xl = coords[pl * dim];
            double yl = coords[pl * dim + 1];
            double x1 = coords[p1 * dim];
            double y1 = coords[p1 * dim + 1];

            // Flip the shared edge when it would improve the Delaunay condition and the quad
            // is convex (both resulting triangles must remain CCW after the flip).
            if (!InCircle(x0, y0, xr, yr, xl, yl, x1, y1) &&
                Orient(x0, y0, xr, yr, x1, y1) > 0 &&
                Orient(x0, y0, x1, y1, xl, yl) > 0)
            {
                triangles[a] = p1;
                triangles[b] = p0;

                int hbl = halfEdges[bl];
                int har = halfEdges[ar];
                halfEdges[a] = hbl;
                if (hbl != -1)
                {
                    halfEdges[hbl] = a;
                }

                halfEdges[b] = har;
                if (har != -1)
                {
                    halfEdges[har] = b;
                }

                halfEdges[ar] = bl;
                halfEdges[bl] = ar;

                if (hbl != -1 && edgeStamp[a] == 0)
                {
                    edgeStamp[a] = 1;
                    edgeStack[stackTop++] = a;
                }

                if (har != -1 && edgeStamp[b] == 0)
                {
                    edgeStamp[b] = 1;
                    edgeStack[stackTop++] = b;
                }

                if (halfEdges[al] != -1 && edgeStamp[al] == 0)
                {
                    edgeStamp[al] = 1;
                    edgeStack[stackTop++] = al;
                }

                if (halfEdges[br] != -1 && edgeStamp[br] == 0)
                {
                    edgeStamp[br] = 1;
                    edgeStack[stackTop++] = br;
                }
            }
        }
    }

    private const int BlockSize = 16;

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

    /// <summary>Eliminates collinear or duplicate points.</summary>
    private static Node? FilterPoints(
        Node? start,
        Node? end = null,
        TriangulationWorkspace? workspace = null)
    {
        if (start is null)
        {
            return null;
        }

        bool full = end is null;
        end ??= start;

        Node p = start;
        bool again;

        do
        {
            again = false;

            if (!ReferenceEquals(p, p.Next) && !p.Steiner &&
                (p == p.Next! || Area(p.Prev!, p, p.Next!) == 0.0))
            {
                if (full || ReferenceEquals(p, end))
                {
                    end = p.Prev!;
                }
                if (workspace is not null)
                {
                    workspace.FilteredOut = true;
                }

                RemoveNode(p, workspace);
                p = p.Prev!;
                again = true;
            }
            else if (full || !ReferenceEquals(p, end))
            {
                p = p.Next!;
                again = !full;
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
        double minX,
        double minY,
        double invSize,
        TriangulationWorkspace workspace)
    {
        if (ear is null)
        {
            return;
        }

        if (invSize != 0.0)
        {
            IndexCurve(ear, minX, minY, invSize, workspace);
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

                RemoveNode(ear, workspace);

                ear = next;
                stop = next;
                continue;
            }

            ear = next;

            if (ReferenceEquals(ear, stop))
            {
                // Try filtering collinear/coincident points and re-clipping.
                workspace.FilteredOut = false;
                ear = FilterPoints(ear, workspace: workspace)!;
                if (workspace.FilteredOut)
                {
                    stop = ear;
                    continue;
                }

                // Cure small local self-intersections once, then retry.
                if (!cured)
                {
                    ear = CureLocalIntersections(ear, triangles, workspace);
                    stop = ear;
                    cured = true;
                    continue;
                }

                // Last resort: split the polygon.
                SplitEarcut(ear, triangles, minX, minY, invSize, workspace);
                break;
            }
        }
    }

    /// <summary>Checks whether a polygon node forms a valid ear.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static bool IsEar(Node ear)
    {
        // Reflex (area >= 0) check is hoisted into EarcutLinked caller.
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
        // Reflex (area >= 0) check is hoisted into EarcutLinked caller.
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

        // Look for points inside the triangle in decreasing z-order.
        Node? p = ear.PrevZ;
        while (p is not null && p.Z >= minZ)
        {
            if (p.X >= x0 && p.X <= x1 && p.Y >= y0 && p.Y <= y1 &&
                !ReferenceEquals(p, c) && !(ax == p.X && ay == p.Y) &&
                PointInTriangle(ax, ay, bx, by, cx, cy, p.X, p.Y) &&
                Area(p.Prev!, p, p.Next!) >= 0.0)
            {
                return false;
            }

            p = p.PrevZ;
        }

        // Look for points in increasing z-order.
        Node? n = ear.NextZ;
        while (n is not null && n.Z <= maxZ)
        {
            if (n.X >= x0 && n.X <= x1 && n.Y >= y0 && n.Y <= y1 &&
                !ReferenceEquals(n, c) && !(ax == n.X && ay == n.Y) &&
                PointInTriangle(ax, ay, bx, by, cx, cy, n.X, n.Y) &&
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
    private static Node CureLocalIntersections(
        Node start,
        TriangleList triangles,
        TriangulationWorkspace workspace)
    {
        Node p = start;
        bool cured = false;

        do
        {
            Node a = p.Prev!;
            Node b = p.Next!.Next!;

            if (!ReferenceEquals(a, b) &&
                Intersects(a, p, p.Next!, b, includeBoundary: false) &&
                LocallyInside(a, b) &&
                LocallyInside(b, a))
            {
                triangles.Add(a.I, p.I, b.I);

                RemoveNode(p, workspace);
                RemoveNode(p.Next!, workspace);

                p = start = b;
                cured = true;
            }

            p = p.Next!;
        }
        while (!ReferenceEquals(p, start));

        return cured ? FilterPoints(p, workspace: workspace)! : p;
    }

    /// <summary>
    /// Last-resort: find a valid diagonal, split the polygon in two, and
    /// triangulate each half independently.
    /// </summary>
    private static void SplitEarcut(
        Node start,
        TriangleList triangles,
        double minX,
        double minY,
        double invSize,
        TriangulationWorkspace workspace)
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

                    a = FilterPoints(a, a.Next, workspace)!;
                    c = FilterPoints(c, c.Next, workspace)!;

                    EarcutLinked(a, triangles, minX, minY, invSize, workspace);
                    EarcutLinked(c, triangles, minX, minY, invSize, workspace);
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
        TriangulationWorkspace workspace)
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
                }

                queue[queueCount++] = GetLeftmost(list);
            }
        }

        queue.AsSpan(0, queueCount).Sort(static (a, b) =>
        {
            double dx = a.X - b.X;
            if (dx != 0.0)
            {
                return dx < 0.0 ? -1 : 1;
            }

            double dy = a.Y - b.Y;
            if (dy != 0.0)
            {
                return dy < 0.0 ? -1 : 1;
            }

            double aSlope = (a.Next!.Y - a.Y) / (a.Next.X - a.X);
            double bSlope = (b.Next!.Y - b.Y) / (b.Next.X - b.X);
            return aSlope.CompareTo(bSlope);
        });

        // Build block-bbox index seeded with the outer ring.
        BuildBlockIndex(workspace, data.Length / dim, holeIndices.Length);
        IndexSegment(outerNode, outerNode, workspace);   // outerNode === outerNode → full ring

        // Process holes from left to right; indexActive keeps GrowBlock live.
        workspace.IndexActive = true;
        try
        {
            for (int i = 0; i < queueCount; i++)
            {
                outerNode = EliminateHole(queue[i], outerNode, workspace);
            }
        }
        finally
        {
            workspace.IndexActive = false;
        }

        // Collapse collinear/coincident points across the whole merged ring before clipping.
        return FilterPoints(outerNode, workspace: workspace)!;
    }

    private static Node EliminateHole(
        Node hole,
        Node outerNode,
        TriangulationWorkspace workspace)
    {
        Node? bridge = FindHoleBridge(hole, outerNode, workspace);

        if (bridge is null)
        {
            return outerNode;
        }

        Node bridgeReverse = SplitPolygon(bridge, hole);

        // Index the merged-in segment before filtering.
        Node bridge2 = bridgeReverse.Next!;
        IndexSegment(bridge, bridge2.Next!, workspace);

        // Heal collinear/coincident points around the two new slit edges.
        FilterPoints(bridgeReverse, bridgeReverse.Next, workspace);
        return FilterPoints(bridge, bridge.Next, workspace)!;
    }

    // ──────────────── block-bbox index for FindHoleBridge ────────────────

    private static void BuildBlockIndex(
        TriangulationWorkspace workspace,
        int maxNodes,
        int numHoles) =>
        workspace.EnsureBlockIndexCapacity(maxNodes, numHoles, BlockSize);

    /// <summary>
    /// Index the ring segment head..stop (exclusive; head==stop means the whole ring)
    /// as ceil(len/BlockSize) blocks, each with a [minX,minY,maxX,maxY] bbox.
    /// Reuses <see cref="Node.Z"/> as the owning block index during this phase.
    /// </summary>
    private static void IndexSegment(
        Node head,
        Node stop,
        TriangulationWorkspace workspace)
    {
        Node p = head;
        do
        {
            int blockIndex = workspace.BlockCount++;
            workspace.BlockHeads![blockIndex] = p;
            double minX = double.PositiveInfinity;
            double minY = double.PositiveInfinity;
            double maxX = double.NegativeInfinity;
            double maxY = double.NegativeInfinity;
            int edgeCount = 0;
            do
            {
                Node current = p.Next!;
                p.Z = blockIndex;
                if (p.X < minX)
                {
                    minX = p.X;
                }

                if (p.X > maxX)
                {
                    maxX = p.X;
                }

                if (p.Y < minY)
                {
                    minY = p.Y;
                }

                if (p.Y > maxY)
                {
                    maxY = p.Y;
                }

                if (current.X < minX)
                {
                    minX = current.X;
                }

                if (current.X > maxX)
                {
                    maxX = current.X;
                }

                if (current.Y < minY)
                {
                    minY = current.Y;
                }

                if (current.Y > maxY)
                {
                    maxY = current.Y;
                }

                p = current;
            }
            while (++edgeCount < BlockSize && !ReferenceEquals(p, stop));

            workspace.BlockStops![blockIndex] = p;
            int bboxIndex = blockIndex * 4;
            double[] blockBoundingBoxes = workspace.BlockBoundingBoxes!;
            blockBoundingBoxes[bboxIndex] = minX;
            blockBoundingBoxes[bboxIndex + 1] = minY;
            blockBoundingBoxes[bboxIndex + 2] = maxX;
            blockBoundingBoxes[bboxIndex + 3] = maxY;
        }
        while (!ReferenceEquals(p, stop));
    }

    /// <summary>
    /// When <see cref="FilterPoints"/> heals an edge head→tail, grow head's block bbox
    /// to cover tail so the leftward-ray prune cannot false-skip the new edge.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void GrowBlock(
        Node head,
        Node tail,
        TriangulationWorkspace workspace)
    {
        int bboxIndex = head.Z * 4;
        double[] blockBoundingBoxes = workspace.BlockBoundingBoxes!;
        if (tail.X < blockBoundingBoxes[bboxIndex])
        {
            blockBoundingBoxes[bboxIndex] = tail.X;
        }

        if (tail.Y < blockBoundingBoxes[bboxIndex + 1])
        {
            blockBoundingBoxes[bboxIndex + 1] = tail.Y;
        }

        if (tail.X > blockBoundingBoxes[bboxIndex + 2])
        {
            blockBoundingBoxes[bboxIndex + 2] = tail.X;
        }

        if (tail.Y > blockBoundingBoxes[bboxIndex + 3])
        {
            blockBoundingBoxes[bboxIndex + 3] = tail.Y;
        }
    }

    private static Node LiveBlockStop(int blockIndex, TriangulationWorkspace workspace)
    {
        Node stop = workspace.BlockStops![blockIndex];
        while (!ReferenceEquals(stop.Prev!.Next, stop))
        {
            stop = stop.Next!;
        }

        workspace.BlockStops[blockIndex] = stop;
        return stop;
    }

    private static Node LiveBlockHead(int blockIndex, TriangulationWorkspace workspace)
    {
        Node head = workspace.BlockHeads![blockIndex];
        while (!ReferenceEquals(head.Prev!.Next, head))
        {
            head = head.Next!;
        }

        workspace.BlockHeads[blockIndex] = head;
        return head;
    }

    /// <summary>
    /// David Eberly's algorithm for finding a bridge between a hole and the
    /// outer polygon — accelerated with the block-bbox index.
    /// </summary>
    private static Node? FindHoleBridge(
        Node hole,
        Node outerNode,
        TriangulationWorkspace workspace)
    {
        double hx = hole.X;
        double hy = hole.Y;
        double qx = double.NegativeInfinity;
        Node? m = null;

        if (hole == outerNode)
        {
            return outerNode;
        }

        double[] blockBoundingBoxes = workspace.BlockBoundingBoxes!;

        // First scan: find the leftward ray crossing to locate the best bridge candidate.
        for (int blockIndex = 0, bboxIndex = 0; blockIndex < workspace.BlockCount; blockIndex++, bboxIndex += 4)
        {
            if (hy < blockBoundingBoxes[bboxIndex + 1] ||
                hy > blockBoundingBoxes[bboxIndex + 3] ||
                blockBoundingBoxes[bboxIndex] > hx ||
                blockBoundingBoxes[bboxIndex + 2] <= qx)
            {
                continue;
            }

            Node blockStop = LiveBlockStop(blockIndex, workspace);
            Node p = LiveBlockHead(blockIndex, workspace);
            do
            {
                if (ReferenceEquals(p.Prev!.Next, p))   // skip removed nodes
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
                                return m;
                            }
                        }
                    }
                }
                p = p.Next!;
            }
            while (!ReferenceEquals(p, blockStop));
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

        // Second scan: find the point of minimum angle within the candidate triangle.
        for (int blockIndex = 0, bboxIndex = 0; blockIndex < workspace.BlockCount; blockIndex++, bboxIndex += 4)
        {
            if (blockBoundingBoxes[bboxIndex + 2] < mx ||
                blockBoundingBoxes[bboxIndex] > hx ||
                blockBoundingBoxes[bboxIndex + 3] < tminY ||
                blockBoundingBoxes[bboxIndex + 1] > tmaxY)
            {
                continue;
            }

            Node blockStop = LiveBlockStop(blockIndex, workspace);
            Node p = LiveBlockHead(blockIndex, workspace);
            do
            {
                if (ReferenceEquals(p.Prev!.Next, p) &&
                    hx >= p.X && p.X >= mx && hx != p.X &&
                    PointInTriangle(
                        hy < my ? hx : qx, hy,
                        mx, my,
                        hy < my ? qx : hx, hy,
                        p.X, p.Y))
                {
                    double tan = Math.Abs(hy - p.Y) / (hx - p.X);
                    if ((LocallyInside(p, hole) ||
                         (p.Y == hy && p.Next!.Y == hy && p.Next.X > hx)) &&
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
            while (!ReferenceEquals(p, blockStop));
        }

        return m;
    }

    // ───────────────────── z-order spatial indexing ─────────────────────────

    /// <summary>
    /// Interlinks polygon nodes in z-order using radix sort for large rings,
    /// insertion sort for small ones.
    /// </summary>
    private static void IndexCurve(
        Node start,
        double minX,
        double minY,
        double invSize,
        TriangulationWorkspace workspace)
    {
        // Count nodes.
        int nodeCount = 0;
        Node p = start;
        do
        {
            nodeCount++;
            p = p.Next!;
        }
        while (!ReferenceEquals(p, start));

        workspace.EnsureSortCapacity(nodeCount);

        // Collect nodes and (re-)compute z-order values; Z may hold a stale block index.
        int i = 0;
        p = start;
        do
        {
            p.Z = ZOrder(p.X, p.Y, minX, minY, invSize);
            workspace.SortNodes![i++] = p;
            p = p.Next!;
        } while (!ReferenceEquals(p, start));

        SortNodes(nodeCount, workspace);

        // Relink the z-order list.
        Node? prev = null;
        for (int j = 0; j < nodeCount; j++)
        {
            Node node = workspace.SortNodes![j];
            node.PrevZ = prev;
            if (prev is not null)
            {
                prev.NextZ = node;
            }

            prev = node;
        }
        prev!.NextZ = null;

        // Clear node references so the GC can collect them after this call.
        Array.Clear(workspace.SortNodes!, 0, nodeCount);
        Array.Clear(workspace.SortBuffer!, 0, nodeCount);
    }

    /// <summary>
    /// Sort the first <paramref name="n"/> entries of <see cref="TriangulationWorkspace.SortNodes"/> by z-order:
    /// insertion sort for small n, four-pass LSD radix sort otherwise.
    /// The result is always in <see cref="TriangulationWorkspace.SortNodes"/> (even number of passes).
    /// </summary>
    private static void SortNodes(int n, TriangulationWorkspace workspace)
    {
        if (n <= 32)
        {
            for (int i = 1; i < n; i++)
            {
                Node node = workspace.SortNodes![i];
                uint z = (uint)node.Z;
                int j = i - 1;
                while (j >= 0 && (uint)workspace.SortNodes[j].Z > z)
                {
                    workspace.SortNodes[j + 1] = workspace.SortNodes[j];
                    j--;
                }

                workspace.SortNodes[j + 1] = node;
            }
            return;
        }

        // Copy z values for the radix passes.
        for (int i = 0; i < n; i++)
        {
            workspace.SortKeys![i] = (uint)workspace.SortNodes![i].Z;
        }

        // Four 8-bit LSD passes; even pass count returns result in SortNodes.
        RadixPass(n, workspace.SortNodes!, workspace.SortKeys!, workspace.SortBuffer!, workspace.SortKeyBuffer!, 0, workspace);
        RadixPass(n, workspace.SortBuffer!, workspace.SortKeyBuffer!, workspace.SortNodes!, workspace.SortKeys!, 8, workspace);
        RadixPass(n, workspace.SortNodes!, workspace.SortKeys!, workspace.SortBuffer!, workspace.SortKeyBuffer!, 16, workspace);
        RadixPass(n, workspace.SortBuffer!, workspace.SortKeyBuffer!, workspace.SortNodes!, workspace.SortKeys!, 24, workspace);
    }

    /// <summary>
    /// One LSD radix pass: stably scatter the first <paramref name="n"/> nodes
    /// from <paramref name="src"/> to <paramref name="dst"/>, bucketed by the
    /// 8-bit digit of z at the given bit <paramref name="shift"/>.
    /// </summary>
    private static void RadixPass(
        int n,
        Node[] src,
        uint[] srcZ,
        Node[] dst,
        uint[] dstZ,
        int shift,
        TriangulationWorkspace workspace)
    {
        uint[] counts = workspace.Counts;
        Array.Clear(counts, 0, counts.Length);

        for (int i = 0; i < n; i++)
        {
            counts[(srcZ[i] >> shift) & 0xff]++;
        }

        uint sum = 0;
        for (int bucket = 0; bucket < counts.Length; bucket++)
        {
            uint count = counts[bucket];
            counts[bucket] = sum;
            sum += count;
        }

        for (int i = 0; i < n; i++)
        {
            uint z = srcZ[i];
            uint position = counts[(z >> shift) & 0xff]++;
            dst[position] = src[i];
            dstZ[position] = z;
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
        bool zeroLength = a == b &&
            Area(a.Prev!, a, a.Next!) > 0.0 &&
            Area(b.Prev!, b, b.Next!) > 0.0;

        return a.Next!.I != b.I &&
            (zeroLength ||
             LocallyInside(a, b) && LocallyInside(b, a) &&
             (Area(a.Prev!, a, b.Prev!) != 0.0 || Area(a, b.Prev!, b) != 0.0)) &&
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

        if (o1 == 0.0 && OnSegment(p1, p2, q1))
        {
            return true;
        }

        if (o2 == 0.0 && OnSegment(p1, q2, q1))
        {
            return true;
        }

        if (o3 == 0.0 && OnSegment(p2, p1, q2))
        {
            return true;
        }

        if (o4 == 0.0 && OnSegment(p2, q1, q2))
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
        double minX = Math.Min(a.X, b.X);
        double maxX = Math.Max(a.X, b.X);
        double minY = Math.Min(a.Y, b.Y);
        double maxY = Math.Max(a.Y, b.Y);

        Node p = a;
        do
        {
            Node n = p.Next!;
            // Fast bbox skip — if both endpoints are on the same outside of the diagonal
            // bbox, the edge can't cross the diagonal.
            if ((p.X > maxX && n.X > maxX) || (p.X < minX && n.X < minX) ||
                (p.Y > maxY && n.Y > maxY) || (p.Y < minY && n.Y < minY))
            {
                p = n;
                continue;
            }
            if (p.I != a.I && n.I != a.I && p.I != b.I && n.I != b.I &&
                Intersects(p, n, a, b))
            {
                return true;
            }

            p = n;
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
            Node n = p.Next!;
            if (((p.Y > py) != (n.Y > py)) &&
                px < (n.X - p.X) * (py - p.Y) / (n.Y - p.Y) + p.X)
            {
                inside = !inside;
            }
            p = n;
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
    private static void RemoveNode(Node p, TriangulationWorkspace? workspace = null)
    {
        Node prev = p.Prev!;
        Node next = p.Next!;
        next.Prev = prev;
        prev.Next = next;
        p.PrevZ?.NextZ = p.NextZ;
        p.NextZ?.PrevZ = p.PrevZ;
        if (workspace?.IndexActive == true)
        {
            GrowBlock(prev, next, workspace);
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

    // ───────────────────── Refine geometry helpers ─────────────────────────

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Orient(double ax, double ay, double bx, double by, double cx, double cy) =>
        (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);

    /// <summary>
    /// Returns true when point p is inside or on the circumcircle of (a,b,c).
    /// Sign convention matches earcut's CCW winding — the tolerance margin
    /// prevents infinite flip loops on near-cocircular quads.
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
        double s = ap + bp + cp;
        return dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) <= 1e-13 * s * s;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int NextHE(int e) => e - e % 3 + (e + 1) % 3;

    private sealed class TriangulationWorkspace
    {
        public bool FilteredOut { get; set; }

        public bool IndexActive { get; set; }

        public double[]? BlockBoundingBoxes { get; private set; }

        public int BlockCount { get; set; }

        public Node[]? BlockHeads { get; private set; }

        public Node[]? BlockStops { get; private set; }

        public Node[]? SortNodes { get; private set; }

        public Node[]? SortBuffer { get; private set; }

        public uint[]? SortKeys { get; private set; }

        public uint[]? SortKeyBuffer { get; private set; }

        public uint[] Counts { get; } = new uint[256];

        public void EnsureBlockIndexCapacity(int maxNodes, int numHoles, int blockSize)
        {
            int maxBlocks = (maxNodes + 2 * numHoles + blockSize - 1) / blockSize + numHoles + 2;
            int bboxLength = maxBlocks * 4;

            if (BlockBoundingBoxes is null || BlockBoundingBoxes.Length < bboxLength)
            {
                BlockBoundingBoxes = new double[bboxLength];
            }

            if (BlockHeads is null || BlockHeads.Length < maxBlocks)
            {
                BlockHeads = new Node[maxBlocks];
            }

            if (BlockStops is null || BlockStops.Length < maxBlocks)
            {
                BlockStops = new Node[maxBlocks];
            }

            BlockCount = 0;
        }

        public void EnsureSortCapacity(int nodeCount)
        {
            if (SortNodes is null || SortNodes.Length < nodeCount)
            {
                SortNodes = new Node[nodeCount];
                SortBuffer = new Node[nodeCount];
            }

            if (SortKeys is null || SortKeys.Length < nodeCount)
            {
                SortKeys = new uint[nodeCount];
                SortKeyBuffer = new uint[nodeCount];
            }
        }
    }

    private sealed class RefineWorkspace
    {
        public RefineWorkspace(int triangleIndexCount)
        {
            EdgeStack = new int[triangleIndexCount];
            HalfEdges = new int[triangleIndexCount];
            EdgeStamp = new byte[triangleIndexCount];

            int hashSize = 1;
            while (hashSize < triangleIndexCount * 4)
            {
                hashSize <<= 1;
            }

            HashTable = new int[hashSize];
            HashStamp = new uint[hashSize];
            HashMask = hashSize - 1;
        }

        public int[] EdgeStack { get; }

        public int[] HalfEdges { get; }

        public int[] HashTable { get; }

        public uint[] HashStamp { get; }

        public byte[] EdgeStamp { get; }

        public int HashMask { get; }
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

        public int Z { get; set; }          // z-order curve value; doubles as block index during eliminateHoles

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
