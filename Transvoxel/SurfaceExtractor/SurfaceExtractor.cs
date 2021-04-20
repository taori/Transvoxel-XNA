using Transvoxel.Geometry;
using Transvoxel.VolumeData;
using Transvoxel.Math;
using Transvoxel.Lengyel;
using System.Diagnostics;

namespace Transvoxel.SurfaceExtractor
{
	public class TransvoxelExtractor : ISurfaceExtractor
    {
        public bool UseCache { get; set; }
        private readonly IVolumeData<sbyte> _volumeData;
        private readonly RegularCellCache _cache;

        public TransvoxelExtractor(IVolumeData<sbyte> volumeDataData)
        {
            _volumeData = volumeDataData;
            _cache = new RegularCellCache(_volumeData.Size.SideLength * 10);
            UseCache = true;
        }

        public Mesh ExtractMesh(Vector3i offsetPosition, ExtractionSettings settings)
        {
            var mesh = new Mesh();
            for (var x = 0; x < settings.MeshLength ; x++)
            {
                for (var y = 0; y < settings.MeshLength; y ++)
                {
                    for (var z = 0; z < settings.MeshLength; z ++)
                    {
                        var position = new Vector3i(x, y, z);
                        PolygonizeCell(offsetPosition, position, ref mesh, settings.LevelOfDetail);
                    }
                }
            }
            return mesh;
        }


        private void PolygonizeCell(Vector3i offsetPos, Vector3i position, ref Mesh mesh, int levelOfDetail)
        {
            Debug.Assert(levelOfDetail >= 1, "Level of Detail must be greater than 1");
            offsetPos += position * levelOfDetail;

            var directionMask = (byte)((position.X > 0 ? 1 : 0) | ((position.Z > 0 ? 1 : 0) << 1) | ((position.Y > 0 ? 1 : 0) << 2));

            var density = new sbyte[8];

            for (var i = 0; i < density.Length; i++)
            {
                density[i] = _volumeData[offsetPos + Tables.CornerIndex[i] * levelOfDetail];
            }

            var caseCode = GetCaseCode(density);
            if ((caseCode ^ ((density[7] >> 7) & 0xFF)) == 0) //for this cases there is no triangulation
                return;

            var cornerNormals = new Vector3f[8];
            for (var i = 0; i < 8; i++)
            {
                var p = offsetPos + Tables.CornerIndex[i] * levelOfDetail;
                var nx = (_volumeData[p + Vector3i.UnitX] - _volumeData[p - Vector3i.UnitX]) * 0.5f;
                var ny = (_volumeData[p + Vector3i.UnitY] - _volumeData[p - Vector3i.UnitY]) * 0.5f;
                var nz = (_volumeData[p + Vector3i.UnitZ] - _volumeData[p - Vector3i.UnitZ]) * 0.5f;
                //cornerNormals[i] = new Vector3f(nx, ny, nz);

                cornerNormals[i].X = nx;
                cornerNormals[i].Y = ny;
                cornerNormals[i].Z = nz;
                cornerNormals[i].Normalize();
            }

            var regularCellClass = Tables.RegularCellClass[caseCode];
            var vertexLocations = Tables.RegularVertexData[caseCode];

            var c = Tables.RegularCellData[regularCellClass];
            var vertexCount = c.GetVertexCount();
            var triangleCount = c.GetTriangleCount();
            var indexOffset = c.Indizes(); //index offsets for current cell
            var mappedIndizes = new ushort[indexOffset.Length]; //array with real indizes for current cell

            for (var i = 0; i < vertexCount; i++)
            {
                var edge = (byte)(vertexLocations[i] >> 8);
                var reuseIndex = (byte)(edge & 0xF); //Vertex id which should be created or reused 1,2 or 3
                var rDir = (byte)(edge >> 4); //the direction to go to reach a previous cell for reusing 

                var v1 = (byte)((vertexLocations[i]) & 0x0F); //Second Corner Index
                var v0 = (byte)((vertexLocations[i] >> 4) & 0x0F); //First Corner Index

                var d0 = density[v0];
                var d1 = density[v1];

                //Vector3f n0 = cornerNormals[v0];
                //Vector3f n1 = cornerNormals[v1];

                Debug.Assert(v1 > v0);

                var t = (d1 << 8) / (d1 - d0);
                var u = 0x0100 - t;
                var t0 = t / 256f;
                var t1 = u / 256f;

                var index = -1;

                if (UseCache && v1 != 7 && (rDir & directionMask) == rDir)
                {
                    Debug.Assert(reuseIndex != 0);
                    var cell = _cache.GetReusedIndex(position, rDir);
                    index = cell.Verts[reuseIndex];
                }

                if (index == -1)
                {
                    var normal = cornerNormals[v0] * t0 + cornerNormals[v1] * t1;
                    GenerateVertex(ref offsetPos, ref position, mesh, levelOfDetail, t, ref v0, ref v1, ref d0, ref d1, normal);
                    index = mesh.LatestAddedVertIndex();
                }

                if ((rDir & 8) != 0)
                {
                    _cache.SetReusableIndex(position, reuseIndex, mesh.LatestAddedVertIndex());
                }

                mappedIndizes[i] = (ushort)index;
            }

            for (var t = 0; t < triangleCount; t++)
            {
                for (var i = 0; i < 3; i++)
                {
                    mesh.AddIndex(mappedIndizes[c.Indizes()[t * 3 + i]]);
                }
            }
        }

        private void GenerateVertex(ref Vector3i offsetPos, ref Vector3i pos, Mesh mesh, int lod, long t, ref byte v0, ref byte v1, ref sbyte d0, ref sbyte d1, Vector3f normal)
        {
            var iP0 = (offsetPos + Tables.CornerIndex[v0] * lod);
            Vector3f P0;// = new Vector3f(iP0.X, iP0.Y, iP0.Z);
            P0.X = iP0.X;
            P0.Y = iP0.Y;
            P0.Z = iP0.Z;

            var iP1 = (offsetPos + Tables.CornerIndex[v1] * lod);
            Vector3f P1;// = new Vector3f(iP1.X, iP1.Y, iP1.Z);
            P1.X = iP1.X;
            P1.Y = iP1.Y;
            P1.Z = iP1.Z;

            //EliminateLodPositionShift(lod, ref d0, ref d1, ref t, ref iP0, ref P0, ref iP1, ref P1);


            var Q = InterpolateVoxelVector(t, P0, P1);

            mesh.AddVertex(Q, normal);
        }

        private void EliminateLodPositionShift(int lod, ref sbyte d0, ref sbyte d1, ref long t, ref Vector3i iP0, ref Vector3f P0, ref Vector3i iP1, ref Vector3f P1)
        {
	        for (var k = 0; k < lod - 1; k++)
            {
                var vm = (P0 + P1) / 2.0f;
                var pm = (iP0 + iP1) / 2;
                var sm = _volumeData[pm];

                if ((d0 & 0x8F) != (d1 & 0x8F))
                {
                    P1 = vm;
                    iP1 = pm;
                    d1 = sm;
                }
                else
                {
                    P0 = vm;
                    iP0 = pm;
                    d0 = sm;
                }
            }

            if (d1 == d0) // ?????????????
                return;
            t = (d1 << 8) / (d1 - d0); // recalc
        }

        /*private static ReuseCell getReuseCell(ReuseCell[, ,] cells, int lod, byte rDir, Vector3i pos)
        {
            int rx = rDir & 0x01;
            int rz = (rDir >> 1) & 0x01;
            int ry = (rDir >> 2) & 0x01;

            int dx = pos.X / lod - rx;
            int dy = pos.Y / lod - ry;
            int dz = pos.Z / lod - rz;

            ReuseCell ccc = cells[dx, dy, dz];
            return ccc;
        }*/

        private static Vector3f InterpolateVoxelVector(long t, Vector3f P0, Vector3f P1)
        {
            var u = 0x0100 - t; //256 - t
            var s = 1.0f / 256.0f;
            var Q = P0 * t + P1 * u; //Density Interpolation
            Q *= s; // shift to shader ! 
            return Q;
        }

        /* internal void mapIndizes2Vertice(int vertexNr, ushort index, ushort[] mappedIndizes, byte[] indexOffset)
         {
             for (int j = 0; j < mappedIndizes.Length; j++)
             {
                 if (vertexNr == indexOffset[j])
                 {
                     mappedIndizes[j] = index;
                 }
             }
         }*/

        private static byte GetCaseCode(sbyte[] density)
        {
            byte code = 0;
            byte konj = 0x01;
            for (var i = 0; i < density.Length; i++)
            {
                code |= (byte)((density[i] >> (density.Length - 1 - i)) & konj);
                konj <<= 1;
            }

            return code;
        }
    }
}
