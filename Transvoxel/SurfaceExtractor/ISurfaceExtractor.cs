using Transvoxel.Geometry;
using Transvoxel.Math;

namespace Transvoxel.SurfaceExtractor
{
	public interface ISurfaceExtractor
	{
		Mesh ExtractMesh(Vector3i offsetPosition, ExtractionSettings settings);
	}
}