#pragma once
#pragma once

#ifndef __DENDROAPI_H__
#define __DENDROAPI_H__

#include "DendroGrid.h"

#ifdef DENDROAPI_EXPORTS
#define DENDRO_API __declspec(dllexport)
#else
#define DENDRO_API __declspec(dllimport)
#endif

#ifdef __cplusplus
extern "C" {
#endif
	// ovdb volume class constructors
	extern DENDRO_API DendroGrid* DendroCreate();
	extern DENDRO_API void DendroDelete(DendroGrid* grid);
	extern DENDRO_API DendroGrid* DendroDuplicate(DendroGrid * grid);

	extern DENDRO_API bool DendroRead(DendroGrid * grid, const char * filename);
	extern DENDRO_API bool DendroWrite(DendroGrid * grid, const char * filename);

	// volume conversion methods
	extern DENDRO_API bool DendroFromPoints(DendroGrid * grid, double *vPoints, int pCount, double *vRadius, int rCount, double voxelSize, double bandwidth);
	extern DENDRO_API bool DendroFromMesh(DendroGrid * grid, float* vPoints, int vCount, int * vFaces, int fCount, double voxelSize, double bandwidth);

	// volume render methods
	extern DENDRO_API void DendroToMesh(DendroGrid * grid);
	extern DENDRO_API void DendroToMeshSettings(DendroGrid * grid, double isovalue, double adaptivity);

	extern DENDRO_API float* DendroVertexBuffer(DendroGrid * grid, int* size);
	extern DENDRO_API int* DendroFaceBuffer(DendroGrid * grid, int* size);

	// volume transformation methods
	extern DENDRO_API bool DendroTransform(DendroGrid * grid, double* matrix, int mCount);

	// volume csg methods
	extern DENDRO_API void DendroUnion(DendroGrid * grid, DendroGrid * csgGrid);
	extern DENDRO_API void DendroDifference(DendroGrid * grid, DendroGrid * csgGrid);
	extern DENDRO_API void DendroIntersection(DendroGrid * grid, DendroGrid * csgGrid);

	// volume filter methods
	extern DENDRO_API void DendroOffset(DendroGrid * grid, double amount);
	extern DENDRO_API void DendroOffsetMask(DendroGrid * grid, double amount, DendroGrid * mask, double min, double max, bool invert);
	extern DENDRO_API void DendroSmooth(DendroGrid * grid, int type, int iterations, int width);
	extern DENDRO_API void DendroSmoothMask(DendroGrid * grid, int type, int iterations, int width, DendroGrid * mask, double min, double max, bool invert);
	extern DENDRO_API void DendroBlend(DendroGrid * bGrid, DendroGrid * eGrid, double bPosition, double bEnd);
	extern DENDRO_API void DendroBlendMask(DendroGrid * bGrid, DendroGrid * eGrid, double bPosition, double bEnd, DendroGrid * mask, double min, double max, bool invert);

	/*extern DENDRO_API bool CreateSphere(DendroGrid* grid, double* vPoints, int pCount, double* vRadius, int rCount, double voxelSize, double bandwidth);*/

	//extern DENDRO_API bool CreateSphereFromOnePoint(DendroGrid* grid, double center_x,double center_y,double center_z,double radius);

	//Method 3:
	extern DENDRO_API void CreateSphereFromOnePoint(DendroGrid * grid, double center_x, double center_y, double center_z, double radius);

	extern DENDRO_API void CreateSphereFromTwoPoint(DendroGrid* grid, double center_x1, double center_y1, double center_z1, double radius);

	//cylinder
	extern DENDRO_API void GenerateCylinder(DendroGrid* grid, double radius, double x1, double y1, double z1, double x2, double y2, double z2);

	// a series of cylinders
	extern DENDRO_API void GenerateCylinders(DendroGrid* grid, double length_mm, double width_mm, double height_mm, double num);

	//fibers 2021.8.4 after graduation
	extern DENDRO_API void GenerateFibers(DendroGrid* grid, double aspectRatioMean, double aspectRatioDe, double diameterMean, double diameterDe, double orientationMean, double orientationDe, double volumeFraction);

	// utilities and analysis
	extern DENDRO_API float* DendroClosestPoint(DendroGrid* grid, float* vPoints, int vCount, int* rSize);

#ifdef __cplusplus
}
#endif

#endif // __DENDROAPI_H__