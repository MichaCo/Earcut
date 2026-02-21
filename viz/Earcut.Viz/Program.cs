// Port of https://github.com/mapbox/earcut/blob/main/viz/viz.js
// Generates a PNG visualization of earcut triangulation.

using System.Diagnostics;
using System.Text.Json;
using EarcutDotNet;
using SixLabors.ImageSharp;
using SixLabors.ImageSharp.Drawing;
using SixLabors.ImageSharp.Drawing.Processing;
using SixLabors.ImageSharp.PixelFormats;
using SixLabors.ImageSharp.Processing;
using IOPath = System.IO.Path;

string fixtureName = args.Length > 0 ? args[0] : "water";
int rotation = args.Length > 1 ? int.Parse(args[1]) : 0;

// Load fixture
string fixturesDir = IOPath.Combine(AppContext.BaseDirectory, "fixtures");
string fixturePath = IOPath.Combine(fixturesDir, $"{fixtureName}.json");
if (!File.Exists(fixturePath))
{
    Console.Error.WriteLine($"Fixture not found: {fixturePath}");
    return 1;
}

double[][][] testPoints = JsonSerializer.Deserialize<double[][][]>(File.ReadAllText(fixturePath))!;

// Apply rotation
double theta = rotation * Math.PI / 180;
bool isAxis = rotation % 90 == 0;
double xx = isAxis ? Math.Round(Math.Cos(theta)) : Math.Cos(theta);
double xy = isAxis ? Math.Round(-Math.Sin(theta)) : -Math.Sin(theta);
double yx = isAxis ? Math.Round(Math.Sin(theta)) : Math.Sin(theta);
double yy = isAxis ? Math.Round(Math.Cos(theta)) : Math.Cos(theta);

foreach (double[][] ring in testPoints)
{
    foreach (double[] coord in ring)
    {
        double x = coord[0];
        double y = coord[1];
        coord[0] = xx * x + xy * y;
        coord[1] = yx * x + yy * y;
    }
}

// Compute bounds from the outer ring
double minX = double.PositiveInfinity, maxX = double.NegativeInfinity;
double minY = double.PositiveInfinity, maxY = double.NegativeInfinity;

for (int i = 0; i < testPoints[0].Length; i++)
{
    minX = Math.Min(minX, testPoints[0][i][0]);
    maxX = Math.Max(maxX, testPoints[0][i][0]);
    minY = Math.Min(minY, testPoints[0][i][1]);
    maxY = Math.Max(maxY, testPoints[0][i][1]);
}

double polyWidth = maxX - minX;
double polyHeight = maxY - minY;

int canvasWidth = 512;
int canvasHeight = (int)(canvasWidth * polyHeight / polyWidth + 10);
double ratio = (canvasWidth - 10.0) / polyWidth;

// Flatten and triangulate
var (vertices, holes, dimensions) = Earcut.Flatten(testPoints);

var sw = Stopwatch.StartNew();
int[] result = Earcut.Triangulate(vertices, holes, dimensions);
sw.Stop();

Console.WriteLine($"earcut: {sw.ElapsedMilliseconds}ms");
Console.WriteLine($"deviation: {Earcut.Deviation(vertices, holes, dimensions, result)}");

// Build triangle vertex list (each entry is a [x,y] point in polygon coordinates)
var triangleVerts = new double[result.Length][];
for (int i = 0; i < result.Length; i++)
{
    int idx = result[i] * dimensions;
    triangleVerts[i] = [vertices[idx], vertices[idx + 1]];
}

// Render
using var image = new Image<Rgba32>(canvasWidth, canvasHeight, Color.White);

image.Mutate(ctx =>
{
    var yellowFill = new Rgba32(255, 255, 0, 51);   // rgba(255,255,0,0.2)
    var redStroke = new Rgba32(255, 0, 0, 51);       // rgba(255,0,0,0.2)

    var fillOpts = new DrawingOptions
    {
        ShapeOptions = new ShapeOptions { IntersectionRule = IntersectionRule.EvenOdd },
        GraphicsOptions = new GraphicsOptions { Antialias = true }
    };

    // Draw each triangle: fill yellow, stroke red
    for (int i = 0; i < triangleVerts.Length; i += 3)
    {
        var pts = new PointF[]
        {
            ToCanvas(triangleVerts[i],     minX, minY, ratio),
            ToCanvas(triangleVerts[i + 1], minX, minY, ratio),
            ToCanvas(triangleVerts[i + 2], minX, minY, ratio),
        };
        IPath triPath = new Polygon(pts);
        ctx.Fill(fillOpts, yellowFill, triPath);
        ctx.Draw(new DrawingOptions(), redStroke, 1f, triPath);
    }

    // Draw polygon outline (all rings) in black
    foreach (double[][] ring in testPoints)
    {
        var pts = ring.Select(pt => ToCanvas(pt, minX, minY, ratio)).ToArray();
        ctx.Draw(new DrawingOptions(), Color.Black, 1f, new Polygon(pts));
    }
});

string outputPath = $"{fixtureName}.png";
await image.SaveAsPngAsync(outputPath).ConfigureAwait(false);

Console.WriteLine($"Saved: {IOPath.GetFullPath(outputPath)}");

// Attempt to open the image with the system viewer
try
{
    Process.Start(new ProcessStartInfo(IOPath.GetFullPath(outputPath)) { UseShellExecute = true });
}
catch
{
    // Not all environments support shell execute — that's fine.
}

return 0;

static PointF ToCanvas(double[] pt, double minX, double minY, double ratio)
    => new((float)((pt[0] - minX) * ratio + 5), (float)((pt[1] - minY) * ratio + 5));
