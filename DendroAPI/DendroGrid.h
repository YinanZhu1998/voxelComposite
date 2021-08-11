#pragma once

#ifndef __DENDROGRID_H__
#define __DENDROGRID_H__

#include "DendroParticle.h"
#include "DendroMesh.h"
#include <openvdb/openvdb.h>
#include <vector>
#include <string>





class DendroGrid
{
public:
	DendroGrid();
	DendroGrid(DendroGrid * grid);
	~DendroGrid();

	openvdb::FloatGrid::Ptr Grid();

	bool Read(const char *vFile);
	bool Write(const char *vFile);

	bool CreateFromMesh(DendroMesh vMesh, double voxelSize, double bandwidth);
	bool CreateFromPoints(DendroParticle vPoints, double voxelSize, double bandwidth);

	//bool CreateSphereFromPoints(DendroParticle vPoints, double voxelSize, double bandwidth);method-1

	//bool CreateSphereFromPoint(double center_x,double center_y, double center_z, double radius);method-2

	
		// Method 3:
		void makeSphere(double radius, const openvdb::Vec3f& c);

		// ADD a sphere
		void addSphere(double radius, const openvdb::Vec3f& c);

    //create a cylinder
		void cylinder(double radius, double x1, double y1, double z1, double x2, double y2, double z2);

    //create a series of random cylinders
		void RandomCylinder(double length_mm, double width_mm, double height_mm,  double num);

	void Transform(openvdb::math::Mat4d xform);

	void BooleanUnion(DendroGrid vAdd);
	void BooleanIntersection(DendroGrid vIntersect);
	void BooleanDifference(DendroGrid vSubtract);

	void Offset(double amount);
	void Offset(double amount, DendroGrid vMask, double min, double max, bool invert);

	void Smooth(int type, int iterations, int width);
	void Smooth(int type, int iterations, int width, DendroGrid vMask, double min, double max, bool invert);

	void Blend(DendroGrid bGrid, double bPosition, double bEnd);
	void Blend(DendroGrid bGrid, double bPosition, double bEnd, DendroGrid vMask, double min, double max, bool invert);

	void ClosestPoint(std::vector<openvdb::Vec3R>& points, std::vector<float>& distances);

	DendroMesh Display();

	void UpdateDisplay();
	void UpdateDisplay(double isovalue, double adaptivity);

	float * GetMeshVertices();
	int * GetMeshFaces();
	int GetVertexCount();
	int GetFaceCount();

	void WriteGeometricModelToAbaqus(const std::string file_name, std::vector<double> isovalue, std::vector<double> elasticProperties, std::vector<double> poissionratio, std::vector<double> density, openvdb::CoordBBox box);

	void WriteConstraint(int node1, int node2, int dof, double c1, double c2, double c3, std::string RP, std::ofstream& file, std::string nodeset1, std::string nodeset2);

	bool checkIntersection(double length, double width, double height, double radius, double x1, double y1, double z1, double x2, double y2, double z2, int dim, double fiberLength, openvdb::FloatGrid::Ptr grid);

	void RandomFibers(double aspectRatioMean, double aspectRatioDe, double diameterMean, double diameterDe, double orientationMean, double orientationDe, double volumeFraction);

private:
	openvdb::FloatGrid::Ptr mGrid;
	DendroMesh mDisplay;
	int mFaceCount;
	int mVertexCount;
};

#endif // __DENDROGRID_H__