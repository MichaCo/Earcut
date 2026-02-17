using System;
using System.Collections.Generic;

namespace Earcut;

/// <summary>
/// Fast polygon triangulation algorithm using ear clipping technique.
/// Direct port from the JavaScript earcut library.
/// </summary>
public static class Earcut
{
    /// <summary>
    /// Triangulates a polygon given as a flat array of vertex coordinates.
    /// </summary>
    /// <param name="data">Flat array of vertex coordinates like [x0,y0, x1,y1, x2,y2, ...]</param>
    /// <param name="holeIndices">Array of hole indices if any (e.g. [5, 8] for a 12-vertex input would mean one hole with vertices 5–7 and another with 8–11)</param>
    /// <param name="dim">Number of coordinates per vertex in the input array (2 by default)</param>
    /// <returns>Array of triangulated indices, where each group of three indices forms a triangle</returns>
    public static int[] Triangulate(ReadOnlySpan<double> data, ReadOnlySpan<int> holeIndices = default, int dim = 2)
    {
        bool hasHoles = holeIndices.Length > 0;
        int outerLen = hasHoles ? holeIndices[0] * dim : data.Length;
        Node? outerNode = LinkedList(data, 0, outerLen, dim, true);
        var triangles = new List<int>();

        if (outerNode == null || outerNode.next == outerNode.prev) return triangles.ToArray();

        double minX = 0, minY = 0, maxX = 0, maxY = 0, x, y, invSize = 0;

        if (data.Length > 80 * dim)
        {
            minX = maxX = data[0];
            minY = maxY = data[1];

            // calculate polygon bbox
            for (int i = dim; i < outerLen; i += dim)
            {
                x = data[i];
                y = data[i + 1];
                if (x < minX) minX = x;
                if (y < minY) minY = y;
                if (x > maxX) maxX = x;
                if (y > maxY) maxY = y;
            }

            // minX, minY and invSize are later used to transform coords into integers for z-order calculation
            invSize = Math.Max(maxX - minX, maxY - minY);
            invSize = invSize != 0 ? 32767 / invSize : 0;
        }

        if (outerNode != null)
        {
            if (hasHoles) outerNode = EliminateHoles(data, holeIndices, outerNode, dim);

            // run earcut algorithm
            if (data.Length > 80 * dim)
            {
                EarcutLinked(outerNode, triangles, dim, minX, minY, invSize, 0);
            }
            else
            {
                EarcutLinked(outerNode, triangles, dim, 0, 0, 0, 0);
            }
        }

        return triangles.ToArray();
    }

    // create a circular doubly linked list from polygon points in the specified winding order
    private static Node? LinkedList(ReadOnlySpan<double> data, int start, int end, int dim, bool clockwise)
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

        if (last != null && Equals(last, last.next))
        {
            RemoveNode(last);
            last = last.next;
        }

        if (last != null)
        {
            last.next!.prev = last;
            last.prev!.next = last;
        }

        return last;
    }

    // eliminate colinear or duplicate points
    private static Node? FilterPoints(Node? start, Node? end = null)
    {
        if (start == null) return start;
        if (end == null) end = start;

        Node p = start;
        bool again;

        do
        {
            again = false;

            if (!p.steiner && (Equals(p, p.next) || Area(p.prev, p, p.next) == 0))
            {
                RemoveNode(p);
                p = end = p.prev!;
                if (p == p.next) break;
                again = true;
            }
            else
            {
                p = p.next!;
            }
        } while (again || p != end);

        return end;
    }

    // main ear slicing loop which triangulates a polygon (given as a linked list)
    private static void EarcutLinked(Node? ear, List<int> triangles, int dim, double minX, double minY, double invSize, int pass)
    {
        if (ear == null) return;

        // iterate through ears, slicing them one by one
        if (pass == 0 && invSize != 0) IndexCurve(ear, minX, minY, invSize);

        Node stop = ear;
        Node? prev, next;

        while (ear!.prev != ear.next)
        {
            prev = ear.prev;
            next = ear.next;

            bool isEar = invSize != 0 ? IsEarHashed(ear, minX, minY, invSize) : IsEar(ear);

            if (isEar)
            {
                // output previous, current and next vertices to form a triangle
                triangles.Add(prev!.i);
                triangles.Add(ear.i);
                triangles.Add(next!.i);

                RemoveNode(ear);

                // skip the next vertex
                ear = next.next;
                stop = next.next;

                continue;
            }

            ear = next;

            // if we looped through all possible ears and can't find any more triangles to cut
            if (ear == stop)
            {
                // try filtering points and slicing again
                if (pass == 0)
                {
                    EarcutLinked(FilterPoints(ear), triangles, dim, minX, minY, invSize, 1);
                }
                // if this didn't work, try curing all small self-intersections locally
                else if (pass == 1)
                {
                    ear = CureLocalIntersections(FilterPoints(ear), triangles);
                    EarcutLinked(ear, triangles, dim, minX, minY, invSize, 2);
                }
                // as a last resort, try splitting remaining polygon into two
                else if (pass == 2)
                {
                    SplitEarcut(ear, triangles, dim, minX, minY, invSize);
                }

                break;
            }
        }
    }

    // check whether a polygon node forms a valid ear with adjacent nodes
    private static bool IsEar(Node ear)
    {
        Node a = ear.prev!;
        Node b = ear;
        Node c = ear.next!;

        if (Area(a, b, c) >= 0) return false; // reflex, can't be an ear

        // now make sure we don't have other points inside the potential ear
        double ax = a.x, bx = b.x, cx = c.x, ay = a.y, by = b.y, cy = c.y;

        // triangle bbox
        double x0 = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx);
        double y0 = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy);
        double x1 = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx);
        double y1 = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy);

        Node? p = c.next;
        while (p != a)
        {
            if (p!.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1 &&
                PointInTriangle(ax, ay, bx, by, cx, cy, p.x, p.y) &&
                Area(p.prev, p, p.next) >= 0) return false;
            p = p.next;
        }

        return true;
    }

    private static bool IsEarHashed(Node ear, double minX, double minY, double invSize)
    {
        Node a = ear.prev!;
        Node b = ear;
        Node c = ear.next!;

        if (Area(a, b, c) >= 0) return false;

        double ax = a.x, bx = b.x, cx = c.x, ay = a.y, by = b.y, cy = c.y;

        // triangle bbox
        double x0 = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx);
        double y0 = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy);
        double x1 = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx);
        double y1 = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy);

        // z-order range for the current triangle bbox
        int minZ = ZOrder((int)((x0 - minX) * invSize), (int)((y0 - minY) * invSize));
        int maxZ = ZOrder((int)((x1 - minX) * invSize), (int)((y1 - minY) * invSize));

        Node? p = ear.prevZ;
        Node? n = ear.nextZ;

        // look for points inside the triangle in both directions
        while (p != null && p.z >= minZ && n != null && n.z <= maxZ)
        {
            if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1 && p != a && p != c &&
                PointInTriangle(ax, ay, bx, by, cx, cy, p.x, p.y) && Area(p.prev, p, p.next) >= 0) return false;

            p = p.prevZ;

            if (n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1 && n != a && n != c &&
                PointInTriangle(ax, ay, bx, by, cx, cy, n.x, n.y) && Area(n.prev, n, n.next) >= 0) return false;

            n = n.nextZ;
        }

        // look for remaining points in decreasing z-order
        while (p != null && p.z >= minZ)
        {
            if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1 && p != a && p != c &&
                PointInTriangle(ax, ay, bx, by, cx, cy, p.x, p.y) && Area(p.prev, p, p.next) >= 0) return false;
            p = p.prevZ;
        }

        // look for remaining points in increasing z-order
        while (n != null && n.z <= maxZ)
        {
            if (n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1 && n != a && n != c &&
                PointInTriangle(ax, ay, bx, by, cx, cy, n.x, n.y) && Area(n.prev, n, n.next) >= 0) return false;
            n = n.nextZ;
        }

        return true;
    }

    // go through all polygon nodes and cure small local self-intersections
    private static Node? CureLocalIntersections(Node? start, List<int> triangles)
    {
        Node? p = start;
        do
        {
            Node a = p!.prev!;
            Node b = p.next!.next!;

            if (!Equals(a, b) && Intersects(a, p, p.next, b) && LocallyInside(a, b) && LocallyInside(b, a))
            {
                triangles.Add(a.i);
                triangles.Add(p.i);
                triangles.Add(b.i);

                // remove two nodes involved
                RemoveNode(p);
                RemoveNode(p.next);

                p = start = b;
            }
            p = p.next;
        } while (p != start);

        return FilterPoints(p);
    }

    // try splitting polygon and triangulate them independently
    private static void SplitEarcut(Node start, List<int> triangles, int dim, double minX, double minY, double invSize)
    {
        // look for a valid diagonal that divides the polygon into two
        Node? a = start;
        do
        {
            Node? b = a!.next!.next;
            while (b != a.prev)
            {
                if (a.i != b!.i && IsValidDiagonal(a, b))
                {
                    // split the polygon in two by the diagonal
                    Node? c = SplitPolygon(a, b);

                    // filter colinear points around the cuts
                    a = FilterPoints(a, a.next);
                    c = FilterPoints(c, c!.next);

                    // run earcut on each half
                    EarcutLinked(a, triangles, dim, minX, minY, invSize, 0);
                    EarcutLinked(c, triangles, dim, minX, minY, invSize, 0);
                    return;
                }
                b = b.next;
            }
            a = a.next;
        } while (a != start);
    }

    // link every hole into the outer loop, producing a single-ring polygon without holes
    private static Node EliminateHoles(ReadOnlySpan<double> data, ReadOnlySpan<int> holeIndices, Node outerNode, int dim)
    {
        var queue = new List<Node?>();

        int len = holeIndices.Length;
        for (int i = 0; i < len; i++)
        {
            int start = holeIndices[i] * dim;
            int end = i < len - 1 ? holeIndices[i + 1] * dim : data.Length;
            Node? list = LinkedList(data, start, end, dim, false);
            if (list == list!.next) list.steiner = true;
            queue.Add(GetLeftmost(list));
        }

        queue.Sort((a, b) =>
        {
            if (a == null || b == null) return 0;
            return a.x.CompareTo(b.x);
        });

        // process holes from left to right
        for (int i = 0; i < queue.Count; i++)
        {
            if (queue[i] != null)
            {
                EliminateHole(queue[i]!, outerNode);
                outerNode = FilterPoints(outerNode, outerNode.next)!;
            }
        }

        return outerNode;
    }

    // find a bridge connecting given point and the polygon
    private static void EliminateHole(Node hole, Node outerNode)
    {
        Node? bridge = FindHoleBridge(hole, outerNode);
        if (bridge == null) return;

        Node b = SplitPolygon(bridge, hole)!;

        // filter out colinear points around cuts
        FilterPoints(FilterPoints(b, b.next), b.next);
    }

    // David Eberly's algorithm for finding a bridge between hole and outer polygon
    private static Node? FindHoleBridge(Node hole, Node outerNode)
    {
        Node? p = outerNode;
        double hx = hole.x;
        double hy = hole.y;
        double qx = double.NegativeInfinity;
        Node? m = null;

        // find a segment intersected by a ray from the hole's leftmost point to the left
        do
        {
            if (hy <= p!.y && hy >= p.next!.y && p.next.y != p.y)
            {
                double x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
                if (x <= hx && x > qx)
                {
                    qx = x;
                    m = p.x < p.next.x ? p : p.next;
                    if (x == hx) return m; // hole touches outer segment; pick leftmost endpoint
                }
            }
            p = p.next;
        } while (p != outerNode);

        if (m == null) return null;

        // look for points inside the triangle of hole point, segment intersection and endpoint;
        // if there are no points found, we have a valid connection;
        // otherwise choose the point of the minimum angle with the ray as connection point

        Node stop = m;
        double mx = m.x;
        double my = m.y;
        double tanMin = double.PositiveInfinity;
        double tan;

        p = m;

        do
        {
            if (hx >= p!.x && p.x >= mx && hx != p.x &&
                PointInTriangle(hx < mx ? hx : mx, hy, mx, my, hx < mx ? mx : hx, hy, p.x, p.y))
            {
                tan = Math.Abs(hy - p.y) / (hx - p.x); // tangential

                if (LocallyInside(p, hole) &&
                    (tan < tanMin || (tan == tanMin && (p.x > m.x || (p.x == m.x && SectorContainsSector(m, p))))))
                {
                    m = p;
                    tanMin = tan;
                }
            }

            p = p.next;
        } while (p != stop);

        return m;
    }

    // whether sector in node m contains sector in node p in the same coordinates
    private static bool SectorContainsSector(Node m, Node p)
    {
        return Area(m.prev, m, p.prev) < 0 && Area(p.next, m, m.next) < 0;
    }

    // interlink polygon nodes in z-order
    private static void IndexCurve(Node start, double minX, double minY, double invSize)
    {
        Node? p = start;
        do
        {
            if (p!.z == 0) p.z = ZOrder((int)((p.x - minX) * invSize), (int)((p.y - minY) * invSize));
            p.prevZ = p.prev;
            p.nextZ = p.next;
            p = p.next;
        } while (p != start);

        p.prevZ!.nextZ = null;
        p.prevZ = null;

        SortLinked(p);
    }

    // Simon Tatham's linked list merge sort algorithm
    private static Node SortLinked(Node? list)
    {
        int inSize = 1;
        Node? p, q, e, tail;
        int numMerges, pSize, qSize, i;

        do
        {
            p = list;
            list = null;
            tail = null;
            numMerges = 0;

            while (p != null)
            {
                numMerges++;
                q = p;
                pSize = 0;
                for (i = 0; i < inSize && q != null; i++)
                {
                    pSize++;
                    q = q.nextZ;
                }
                qSize = inSize;

                while (pSize > 0 || (qSize > 0 && q != null))
                {
                    if (pSize != 0 && (qSize == 0 || q == null || p!.z <= q.z))
                    {
                        e = p;
                        p = p!.nextZ;
                        pSize--;
                    }
                    else
                    {
                        e = q;
                        q = q!.nextZ;
                        qSize--;
                    }

                    if (tail != null) tail.nextZ = e;
                    else list = e;

                    e!.prevZ = tail;
                    tail = e;
                }

                p = q;
            }

            tail!.nextZ = null;
            inSize *= 2;

        } while (numMerges > 1);

        return list!;
    }

    // z-order of a point given coords and inverse of the longer side of data bbox
    private static int ZOrder(int x, int y)
    {
        x = (x | (x << 8)) & 0x00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F;
        x = (x | (x << 2)) & 0x33333333;
        x = (x | (x << 1)) & 0x55555555;

        y = (y | (y << 8)) & 0x00FF00FF;
        y = (y | (y << 4)) & 0x0F0F0F0F;
        y = (y | (y << 2)) & 0x33333333;
        y = (y | (y << 1)) & 0x55555555;

        return x | (y << 1);
    }

    // find the leftmost node of a polygon ring
    private static Node GetLeftmost(Node start)
    {
        Node? p = start;
        Node leftmost = start;

        do
        {
            if (p!.x < leftmost.x || (p.x == leftmost.x && p.y < leftmost.y)) leftmost = p;
            p = p.next;
        } while (p != start);

        return leftmost;
    }

    // check if a point lies within a convex triangle
    private static bool PointInTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py)
    {
        return (cx - px) * (ay - py) >= (ax - px) * (cy - py) &&
               (ax - px) * (by - py) >= (bx - px) * (ay - py) &&
               (bx - px) * (cy - py) >= (cx - px) * (by - py);
    }

    // check if a diagonal between two polygon nodes is valid (lies in polygon interior)
    private static bool IsValidDiagonal(Node a, Node b)
    {
        return a.next!.i != b.i && a.prev!.i != b.i && !IntersectsPolygon(a, b) && // doesn't intersect other edges
               (LocallyInside(a, b) && LocallyInside(b, a) && MiddleInside(a, b) && // locally visible
                (Area(a.prev, a, b.prev) != 0 || Area(a, b.prev, b) != 0) || // does not create opposite-facing sectors
                Equals(a, b) && Area(a.prev, a, a.next) > 0 && Area(b.prev, b, b.next) > 0); // special zero-length case
    }

    // signed area of a triangle
    private static double Area(Node? p, Node? q, Node? r)
    {
        return (q!.y - p!.y) * (r!.x - q.x) - (q.x - p.x) * (r.y - q.y);
    }

    // check if two points are equal
    private static bool Equals(Node? p1, Node? p2)
    {
        return p1!.x == p2!.x && p1.y == p2.y;
    }

    // check if two segments intersect
    private static bool Intersects(Node p1, Node? q1, Node? p2, Node q2)
    {
        int o1 = Sign(Area(p1, q1, p2));
        int o2 = Sign(Area(p1, q1, q2));
        int o3 = Sign(Area(p2, q2, p1));
        int o4 = Sign(Area(p2, q2, q1));

        if (o1 != o2 && o3 != o4) return true; // general case

        if (o1 == 0 && OnSegment(p1, p2!, q1!)) return true; // p1, q1 and p2 are collinear and p2 lies on p1q1
        if (o2 == 0 && OnSegment(p1, q2, q1!)) return true; // p1, q1 and q2 are collinear and q2 lies on p1q1
        if (o3 == 0 && OnSegment(p2!, p1, q2)) return true; // p2, q2 and p1 are collinear and p1 lies on p2q2
        if (o4 == 0 && OnSegment(p2!, q1!, q2)) return true; // p2, q2 and q1 are collinear and q1 lies on p2q2

        return false;
    }

    // for collinear points p, q, r, check if point q lies on segment pr
    private static bool OnSegment(Node p, Node q, Node r)
    {
        return q.x <= Math.Max(p.x, r.x) && q.x >= Math.Min(p.x, r.x) &&
               q.y <= Math.Max(p.y, r.y) && q.y >= Math.Min(p.y, r.y);
    }

    private static int Sign(double val)
    {
        return val > 0 ? 1 : val < 0 ? -1 : 0;
    }

    // check if a polygon diagonal intersects with any existing edges
    private static bool IntersectsPolygon(Node a, Node b)
    {
        Node? p = a;
        do
        {
            if (p!.i != a.i && p.next!.i != a.i && p.i != b.i && p.next.i != b.i &&
                Intersects(p, p.next, a, b)) return true;
            p = p.next;
        } while (p != a);

        return false;
    }

    // check if a polygon diagonal is locally inside the polygon
    private static bool LocallyInside(Node a, Node b)
    {
        return Area(a.prev, a, a.next) < 0 ?
            Area(a, b, a.next) >= 0 && Area(a, a.prev, b) >= 0 :
            Area(a, b, a.prev) < 0 || Area(a, a.next, b) < 0;
    }

    // check if the middle point of a polygon diagonal is inside the polygon
    private static bool MiddleInside(Node a, Node b)
    {
        Node? p = a;
        bool inside = false;
        double px = (a.x + b.x) / 2;
        double py = (a.y + b.y) / 2;

        do
        {
            double sx = p!.x;
            double sy = p.y;
            double ex = p.next!.x;
            double ey = p.next.y;

            if (((sy > py) != (ey > py)) && (px < (ex - sx) * (py - sy) / (ey - sy) + sx)) inside = !inside;
            p = p.next;
        } while (p != a);

        return inside;
    }

    // link two polygon vertices with a bridge
    private static Node? SplitPolygon(Node a, Node b)
    {
        Node a2 = new Node(a.i, a.x, a.y);
        Node b2 = new Node(b.i, b.x, b.y);
        Node an = a.next!;
        Node bp = b.prev!;

        a.next = b;
        b.prev = a;

        a2.next = an;
        an.prev = a2;

        b2.next = a2;
        a2.prev = b2;

        bp.next = b2;
        b2.prev = bp;

        return b2;
    }

    // create a node and optionally link it with a previous one (in a circular doubly linked list)
    private static Node InsertNode(int i, double x, double y, Node? last)
    {
        Node p = new Node(i, x, y);

        if (last == null)
        {
            p.prev = p;
            p.next = p;
        }
        else
        {
            p.next = last.next;
            p.prev = last;
            last.next!.prev = p;
            last.next = p;
        }
        return p;
    }

    private static void RemoveNode(Node p)
    {
        p.next!.prev = p.prev;
        p.prev!.next = p.next;

        if (p.prevZ != null) p.prevZ.nextZ = p.nextZ;
        if (p.nextZ != null) p.nextZ.prevZ = p.prevZ;
    }

    private static double SignedArea(ReadOnlySpan<double> data, int start, int end, int dim)
    {
        double sum = 0;
        int j = end - dim;

        for (int i = start; i < end; i += dim)
        {
            sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
            j = i;
        }

        return sum;
    }

    /// <summary>
    /// Flattens a multi-dimensional array (e.g. GeoJSON Polygon) into a format expected by Earcut.
    /// </summary>
    /// <param name="data">Multi-dimensional array of polygons, where each polygon is an array of rings, and each ring is an array of [x, y] coordinates</param>
    /// <returns>Flattened data with vertices, holes, and dimensions</returns>
    public static (double[] vertices, int[] holes, int dimensions) Flatten(double[][][] data)
    {
        int dim = data[0][0].Length;
        var vertices = new List<double>();
        var holes = new List<int>();
        int holeIndex = 0;

        for (int i = 0; i < data.Length; i++)
        {
            for (int j = 0; j < data[i].Length; j++)
            {
                for (int d = 0; d < dim; d++)
                {
                    vertices.Add(data[i][j][d]);
                }
            }
            if (i > 0)
            {
                holeIndex += data[i - 1].Length;
                holes.Add(holeIndex);
            }
        }

        return (vertices.ToArray(), holes.ToArray(), dim);
    }

    /// <summary>
    /// Returns the relative difference between the total area of triangles and the area of the input polygon.
    /// </summary>
    /// <param name="data">Flat array of vertex coordinates</param>
    /// <param name="holeIndices">Array of hole indices</param>
    /// <param name="dim">Number of coordinates per vertex</param>
    /// <param name="triangles">Array of triangle indices</param>
    /// <returns>Deviation value where 0 means the triangulation is fully correct</returns>
    public static double Deviation(ReadOnlySpan<double> data, ReadOnlySpan<int> holeIndices, int dim, ReadOnlySpan<int> triangles)
    {
        bool hasHoles = holeIndices.Length > 0;
        int outerLen = hasHoles ? holeIndices[0] * dim : data.Length;

        double polygonArea = Math.Abs(SignedArea(data, 0, outerLen, dim));
        if (hasHoles)
        {
            for (int i = 0; i < holeIndices.Length; i++)
            {
                int start = holeIndices[i] * dim;
                int end = i < holeIndices.Length - 1 ? holeIndices[i + 1] * dim : data.Length;
                polygonArea -= Math.Abs(SignedArea(data, start, end, dim));
            }
        }

        double trianglesArea = 0;
        for (int i = 0; i < triangles.Length; i += 3)
        {
            int a = triangles[i] * dim;
            int b = triangles[i + 1] * dim;
            int c = triangles[i + 2] * dim;
            trianglesArea += Math.Abs(
                (data[a] - data[c]) * (data[b + 1] - data[a + 1]) -
                (data[a] - data[b]) * (data[c + 1] - data[a + 1]));
        }

        return polygonArea == 0 && trianglesArea == 0 ? 0 :
            Math.Abs((trianglesArea - polygonArea) / polygonArea);
    }

    private class Node
    {
        public int i; // vertex index in original input array
        public double x; // x coordinate
        public double y; // y coordinate
        public Node? prev; // previous vertex node in a polygon ring
        public Node? next; // next vertex node in a polygon ring
        public int z; // z-order curve value
        public Node? prevZ; // previous vertex node in z-order
        public Node? nextZ; // next vertex node in z-order
        public bool steiner; // indicates whether this is a steiner point

        public Node(int i, double x, double y)
        {
            this.i = i;
            this.x = x;
            this.y = y;
            this.prev = null;
            this.next = null;
            this.z = 0;
            this.prevZ = null;
            this.nextZ = null;
            this.steiner = false;
        }
    }
}
