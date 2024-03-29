// DendroAPI.cpp : Defines the exported functions for the DLL application.
//
#include "stdafx.h"
#include "DendroAPI.h"

#include"DendroParticle.h"
#include"DendroMesh.h"
#include <openvdb/util/Util.h>
#include <vector>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/openvdb.h>

// grid class constructors
DENDRO_API DendroGrid* DendroCreate()
{
	DendroGrid *grid = new DendroGrid();
	return grid;
}

DENDRO_API void DendroDelete(DendroGrid * grid)
{
	if (grid != NULL) {
		delete grid;
		grid = NULL;
	}
}

DENDRO_API DendroGrid* DendroDuplicate(DendroGrid * grid)
{
	DendroGrid *dup = new DendroGrid(grid);
	return dup;
}

DENDRO_API bool DendroRead(DendroGrid * grid, const char * filename)
{
	return grid->Read(filename);
}

DENDRO_API bool DendroWrite(DendroGrid * grid, const char * filename)
{
	return grid->Write(filename);
}


// grid conversion methods
DENDRO_API bool DendroFromPoints(DendroGrid * grid, double *vPoints, int pCount, double *vRadius, int rCount, double voxelSize, double bandwidth)
{
	double inverseVoxelSize = 1.0 / voxelSize;

	std::vector<openvdb::Vec3R> particleList;

	int i = 0;
	while (i < pCount)
	{
		double x = vPoints[i]; // *inverseVoxelSize;
		double y = vPoints[i + 1]; // * inverseVoxelSize;
		double z = vPoints[i + 2]; // * inverseVoxelSize;

		particleList.push_back(openvdb::Vec3R(x, y, z));

		i += 3;
	}

	DendroParticle ps;
	ps.clear();

	if (particleList.size() == rCount)
	{

		int i = 0;
		for (auto it = particleList.begin(); it != particleList.end(); ++it)
		{
			ps.add((*it), openvdb::Real(vRadius[i]));
			i++;
		}
	}
	else
	{

		double average = 0.0;
		for (int i = 0; i < rCount; i++)
		{
			average += vRadius[i];
		}
		average /= rCount;
		openvdb::Real radius = openvdb::Real(average);

		for (auto it = particleList.begin(); it != particleList.end(); ++it)
		{
			ps.add((*it), radius);
		}
	}

	return grid->CreateFromPoints(ps, voxelSize, bandwidth);

}

DENDRO_API bool DendroFromMesh(DendroGrid * grid, float* vPoints, int vCount, int * vFaces, int fCount, double voxelSize, double bandwidth)
{
	double inverseVoxelSize = 1.0 / voxelSize;

	DendroMesh vMesh;
	vMesh.Clear();

	int i = 0;
	while (i < vCount) {

		openvdb::Vec3s vertex(vPoints[i], vPoints[i + 1], vPoints[i + 2]);

		vertex *= inverseVoxelSize;

		vMesh.AddVertice(vertex);

		i += 3;
	}

	i = 0;
	while (i<fCount) {
		openvdb::Vec4I face(vFaces[i], vFaces[i + 1], vFaces[i + 2], openvdb::util::INVALID_IDX);

		vMesh.AddFace(face);
		i += 3;
	}

	return grid->CreateFromMesh(vMesh, voxelSize, bandwidth);

}


// grid render methods
DENDRO_API void DendroToMesh(DendroGrid * grid)
{
	grid->UpdateDisplay();
}

DENDRO_API void DendroToMeshSettings(DendroGrid * grid, double isovalue, double adaptivity)
{
	grid->UpdateDisplay(isovalue, adaptivity);
}

DENDRO_API float* DendroVertexBuffer(DendroGrid * grid, int * size)
{
	float *verticeArray = grid->GetMeshVertices();

	*size = grid->GetVertexCount();

	return verticeArray;
}

DENDRO_API int * DendroFaceBuffer(DendroGrid * grid, int * size)
{
	int *faceArray = grid->GetMeshFaces();

	*size = grid->GetFaceCount();

	return faceArray;
}


// grid transformation methods
DENDRO_API bool DendroTransform(DendroGrid *grid, double *matrix, int mCount)
{
	if (mCount != 16) {
		return false;
	}

	openvdb::math::Mat4d xform = openvdb::math::Mat4d(matrix[0], matrix[1], matrix[2], matrix[3],
		matrix[4], matrix[5], matrix[6], matrix[7],
		matrix[8], matrix[9], matrix[10], matrix[11],
		matrix[12], matrix[13], matrix[14], matrix[15]);

	grid->Transform(xform);

	return true;
}


// grid csg methods
DENDRO_API void DendroUnion(DendroGrid * grid, DendroGrid * csgGrid)
{
	grid->BooleanUnion(*csgGrid);
}

DENDRO_API void DendroDifference(DendroGrid * grid, DendroGrid * csgGrid)
{
	grid->BooleanDifference(*csgGrid);
}

DENDRO_API void DendroIntersection(DendroGrid * grid, DendroGrid * csgGrid)
{
	grid->BooleanIntersection(*csgGrid);
}


// grid filter methods
DENDRO_API void DendroOffset(DendroGrid * grid, double amount)
{
	grid->Offset(amount);
}

DENDRO_API void DendroOffsetMask(DendroGrid * grid, double amount, DendroGrid * mask, double min, double max, bool invert)
{
	grid->Offset(amount, *mask, min, max, invert);
}

DENDRO_API void DendroSmooth(DendroGrid * grid, int type, int iterations, int width)
{
	grid->Smooth(type, iterations, width);
}

DENDRO_API void DendroSmoothMask(DendroGrid * grid, int type, int iterations, int width, DendroGrid * mask, double min, double max, bool invert)
{
	grid->Smooth(type, iterations, width, *mask, min, max, invert);
}

DENDRO_API void DendroBlend(DendroGrid * bGrid, DendroGrid * eGrid, double bPosition, double bEnd)
{
	bGrid->Blend(*eGrid, bPosition, bEnd);
}

DENDRO_API void DendroBlendMask(DendroGrid * bGrid, DendroGrid * eGrid, double bPosition, double bEnd, DendroGrid * mask, double min, double max, bool invert)
{
	bGrid->Blend(*eGrid, bPosition, bEnd, *mask, min, max, invert);
}

//DENDRO_API bool CreateSphere(DendroGrid* grid, double* vPoints, int pCount, double* vRadius, int rCount, double voxelSize, double bandwidth)
//{
//
//		double inverseVoxelSize = 1.0 / voxelSize;
//
//		std::vector<openvdb::Vec3R> particleList;
//
//		int i = 0;
//		while (i < pCount)
//		{
//			double x = vPoints[i]; // *inverseVoxelSize;
//			double y = vPoints[i + 1]; // * inverseVoxelSize;
//			double z = vPoints[i + 2]; // * inverseVoxelSize;
//
//			particleList.push_back(openvdb::Vec3R(x, y, z));
//
//			i += 3;
//		}
//
//		DendroParticle ps;
//		ps.clear();
//
//		if (particleList.size() == rCount)
//		{
//
//			int i = 0;
//			for (auto it = particleList.begin(); it != particleList.end(); ++it)
//			{
//				ps.add((*it), openvdb::Real(vRadius[i]));
//				i++;
//			}
//		}
//		else
//		{
//
//			double average = 0.0;
//			for (int i = 0; i < rCount; i++)
//			{
//				average += vRadius[i];
//			}
//			average /= rCount;
//			openvdb::Real radius = openvdb::Real(average);
//
//			for (auto it = particleList.begin(); it != particleList.end(); ++it)
//			{
//				ps.add((*it), radius);
//			}
//		}
//
//		return grid->CreateSphereFromPoints(ps, voxelSize, bandwidth);

//}

//Method 2 :
//DENDRO_API bool CreateSphereFromOnePoint(DendroGrid* grid, double center_x,double center_y,double center_z,double radius)
//{
//
//	return grid->CreateSphereFromPoint(center_x,center_y,center_z,radius);
//
//}


//Method 3 ：
DENDRO_API void CreateSphereFromOnePoint(DendroGrid* grid, double center_x, double center_y, double center_z, double radius)
{



		
	const openvdb::Vec3f &c=openvdb::Vec3f(center_x, center_y, center_z);
	
			
	

	return grid->makeSphere(radius,  c);


}

DENDRO_API void CreateSphereFromTwoPoint(DendroGrid* grid, double center_x1, double center_y1, double center_z1,double radius)
{
	const openvdb::Vec3f& c = openvdb::Vec3f(center_x1, center_y1, center_z1);


	return grid->addSphere(radius, c);
}


//CYLINDER
DENDRO_API void GenerateCylinder(DendroGrid* grid, double radius, double x1, double y1, double z1, double x2, double y2, double z2)
{
	


	return grid->cylinder(radius, x1,y1,z1,x2,y2,z2);
}

//A series of CYLINDER
DENDRO_API void GenerateCylinders(DendroGrid* grid, double length_mm, double width_mm, double height_mm,  double num)
{



	return grid->RandomCylinder(length_mm, width_mm, height_mm,  num);
}

//generate fibers 2021.8.4 after graduation
DENDRO_API void GenerateFibers (DendroGrid* grid, double aspectRatioMean, double aspectRatioDe, double diameterMean, double diameterDe, double orientationMean, double orientationDe, double volumeFraction)
{
	return grid->RandomFibers(aspectRatioMean, aspectRatioDe, diameterMean, diameterDe, orientationMean, orientationDe, volumeFraction);
}

// volume utilities
DENDRO_API float* DendroClosestPoint(DendroGrid* grid, float* vPoints, int vCount, int* rSize)
{
	std::vector<openvdb::Vec3R> points;
	std::vector<float> distances;

	int i = 0;
	while (i < vCount) {

		openvdb::Vec3R vertex(vPoints[i], vPoints[i + 1], vPoints[i + 2]);

		points.push_back(vertex);

		i += 3;
	}

	grid->ClosestPoint(points, distances);

	*rSize = points.size() * 3;

	float* pArray = reinterpret_cast<float*>(malloc(*rSize * sizeof(float)));

	i = 0;
	for (auto it = points.begin(); it != points.end(); ++it) {
		pArray[i] = it->x();
		pArray[i + 1] = it->y();
		pArray[i + 2] = it->z();
		i += 3;
	}

	return pArray;
}