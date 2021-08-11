#include "stdafx.h"
#include "DendroGrid.h"

#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/LevelSetMorph.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/Types.h>
#include <openvdb/tools/VolumeToSpheres.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/ChangeBackground.h>



#include <cmath>
#include<random>
#include <algorithm>


#include <iostream>  
#include<fstream>

using namespace std;



DendroGrid::DendroGrid()
{
	openvdb::initialize();
}

DendroGrid::DendroGrid(DendroGrid * grid)
{
	openvdb::initialize();
	mGrid = grid->Grid()->deepCopy();
	mDisplay = grid->Display().Duplicate();
}

DendroGrid::~DendroGrid()
{
}

openvdb::FloatGrid::Ptr DendroGrid::Grid()
{
	return mGrid;
}

bool DendroGrid::Read(const char * vFile)
{
	openvdb::io::File file(vFile);

	file.open();

	openvdb::io::File::NameIterator nameIter = file.beginName();
	if (nameIter == file.endName()) {
		return false;
	}

	mGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(file.readGrid(nameIter.gridName()));

	return true;
}

bool DendroGrid::Write(const char * vFile)
{
	openvdb::GridPtrVec grids;
	grids.push_back(mGrid);

	openvdb::io::File file(vFile);
	file.write(grids);
	file.close();

	return true;
}

bool DendroGrid::CreateFromMesh(DendroMesh vMesh, double voxelSize, double bandwidth)
{
	if (!vMesh.IsValid()) {
		return false;
	}

	openvdb::math::Transform xform;
	xform.preScale(voxelSize);

	auto vertices = vMesh.Vertices();
	auto faces = vMesh.Faces();

	openvdb::tools::QuadAndTriangleDataAdapter<openvdb::Vec3s, openvdb::Vec4I> mesh(vertices, faces);
	mGrid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh, xform, static_cast<float>(bandwidth), static_cast<float>(bandwidth), 0, NULL);

	mDisplay = vMesh;

	return true;
}

bool DendroGrid::CreateFromPoints(DendroParticle vPoints, double voxelSize, double bandwidth)
{
	if (!vPoints.IsValid()) {
		return false;
	}

	mGrid = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, bandwidth);
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*mGrid);

	openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(voxelSize);
	mGrid->setTransform(xform);

	raster.setGrainSize(1);
	raster.rasterizeSpheres(vPoints);
	raster.finalize();

	return true;
}

//bool DendroGrid::CreateSphereFromPoint(double center_x, double center_y, double center_z, double radius)
//{
//	if (radius > 0) {
//
//
//
//		mGrid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
//			/*radius=*/radius, /*center=*/openvdb::Vec3f(center_x,center_y,center_z) ,
//			/*voxel size=*/0.5, /*width=*/4.0);
//
//
//		return true;
//	}
//	else
//		return false;
//


void DendroGrid::makeSphere(double radius, const openvdb::Vec3f& c)//define function "makeShphere"
{
	

		openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/2.0);
		openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.5);
		grid->setTransform(xform);


		// Distance value for the constant region exterior to the narrow band  
		const float outside = 2;
		// Distance value for the constant region interior to the narrow band
		// (by convention, the signed distance is negative in the interior of
		// a level set)
		const float inside = -outside;
		// Use the background value as the width in voxels of the narrow band.
		// (The narrow band is centered on the surface of the sphere, which
		// has distance 0.)
		int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
		// The bounding box of the narrow band is 2*dim voxels on a side.
		int dim = int(radius + padding);
		// Get a voxel accessor.
		openvdb::FloatGrid::Accessor  accessor = grid->getAccessor();
		// Compute the signed distance from the surface of the sphere of each
		// voxel within the bounding box and insert the value into the grid
		// if it is smaller in magnitude than the background value.
		openvdb::Coord ijk;
		int& i = ijk[0], & j = ijk[1], & k = ijk[2];
		for (i = c[0] - dim; i < c[0] + dim; ++i) {
			const float x2 = openvdb::math::Pow2(i - c[0]);
			for (j = c[1] - dim; j < c[1] + dim; ++j) {
				const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
				for (k = c[2] - dim; k < c[2] + dim; ++k) {
					// The distance from the sphere surface in voxels
					const float dist = openvdb::math::Sqrt(x2y2
						+ openvdb::math::Pow2(k - c[2])) - radius;

					// Only insert distances that are smaller in magnitude than
					// the background value.
					if (dist < inside || outside < dist) continue;
					// Set the distance for voxel (i,j,k).
					accessor.setValue(ijk, dist);
				}
			}
		}


		mGrid = grid;

	
	
	
}

void DendroGrid::addSphere(double radius, const openvdb::Vec3f& c)
{
	 if (mGrid == NULL)
	 {
		 //create a NEW grid 
		 mGrid = openvdb::FloatGrid::create(/*background value=*/2.0);
		 openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.5);
		 mGrid->setTransform(xform);
	 }
	 //settings
	 const float outside = 2;
	 const float inside = -outside;
	 int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
	
	 //generate the first sphere
	 int dim = int(radius + padding);
	 openvdb::FloatGrid::Accessor  accessor = mGrid->getAccessor();

	 openvdb::Coord ijk;
	 int& i = ijk[0], & j = ijk[1], & k = ijk[2];
	 for (i = c[0] - dim; i < c[0] + dim; ++i) {
		 const float x2 = openvdb::math::Pow2(i - c[0]);
		 for (j = c[1] - dim; j < c[1] + dim; ++j) {
			 const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
			 for (k = c[2] - dim; k < c[2] + dim; ++k) {
				 // The distance from the sphere surface in voxels
				 const float dist = openvdb::math::Sqrt(x2y2
					 + openvdb::math::Pow2(k - c[2])) - radius;

				 // Only insert distances that are smaller in magnitude than
				 // the background value.
				 if (dist < 0 || 0 < dist) continue;
				 // Set the distance for voxel (i,j,k).
				 accessor.setValue(ijk, dist);
			 }
		 }
	 }
	
}

void DendroGrid::cylinder(double radius, double x1, double y1, double z1, double x2, double y2, double z2)
{
	//create a NEW grid 
	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/2.0);
	openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.2);
	grid->setTransform(xform);

	//settings
	const float outside = 2;
	const float inside = -outside;
	int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));

	//create a cylinder
	openvdb::FloatGrid::Accessor  accessor = grid->getAccessor();

	//distance from Pa to Pb
	double distance_line = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));

	//mid point from Pa to Pb
	double center_x = 0.5 * (x1 + x2);
	double center_y = 0.5 * (y1 + y2);
	double center_z = 0.5 * (z1 + z2);

	//two types of shape of cylinder
	int dim;
	if (radius > 0.5 * distance_line)
		dim = (radius + padding);

	else
		dim = (0.5 * distance_line + padding);

	//create coord ijk
	openvdb::Coord ijk;
	int& i = ijk[0], & j = ijk[1], & k = ijk[2];

	//create cylinder
	for (i = center_x - dim; i < center_x + dim; ++i) {
		for (j = center_y - dim; j < center_y + dim; ++j) {
			for (k = center_z - dim; k < center_z + dim; ++k) {
				//coefficient U
				double U = (((i - x1) * (x2 - x1)) + ((j - y1) * (y2 - y1)) + ((k - z1) * (z2 - z1))) / (distance_line * distance_line);

				//coord of intersection point on the line Pa--Pb
				double intersection_x = x1 + U * (x2 - x1);
				double intersection_y = y1 + U * (y2 - y1);
				double intersection_z = z1 + U * (z2 - z1);

				//distance from intersection point to point(i,j,k)
				const float distance = openvdb::math::Sqrt((i - intersection_x) * (i - intersection_x) + (j - intersection_y) * (j - intersection_y) + (k - intersection_z) * (k - intersection_z))-radius;

				//check if point(i,j,k) is within cylinder region
				double condition1 = ((i - x1) * (x2 - x1) + (j - y1) * (y2 - y1) + (k - z1) * (z2 - z1)) / distance_line;
				double condition2= ((i - x2) * (x2 - x1) + (j - y2) * (y2 - y1) + (k - z2) * (z2 - z1)) / distance_line;
				if (condition1 < 0 || condition2>0 || distance > 0) continue;
				accessor.setValue(ijk, distance);
			}
		}
	}
	mGrid = grid;
}

bool DendroGrid::checkIntersection(double length, double width, double height, double radius, double x1, double y1, double z1, double x2, double y2, double z2, int dim, double fiberLength, openvdb::FloatGrid::Ptr grid)
{
	auto accessor = grid->getAccessor();

	double center_x = 0.5 * (x1 + x2);
	double center_y = 0.5 * (y1 + y2);
	double center_z = 0.5 * (z1 + z2);

	openvdb::Coord ijk;
	int& i = ijk[0], & j = ijk[1], & k = ijk[2];

	openvdb::Coord abc;
	int& a = abc[0], & b = abc[1], & c = abc[2];

	for (i = center_x - dim; i < center_x + dim; ++i) {
		for (j = center_y - dim; j < center_y + dim; ++j) {
			for (k = center_z - dim; k < center_z + dim; ++k) {

				//coefficient U
				double U = (((i - x1) * (x2 - x1)) + ((j - y1) * (y2 - y1)) + ((k - z1) * (z2 - z1))) / (fiberLength * fiberLength);

				//coord of intersection point on the line Pa--Pb
				double intersection_x = x1 + U * (x2 - x1);
				double intersection_y = y1 + U * (y2 - y1);
				double intersection_z = z1 + U * (z2 - z1);

				//distance from intersection point to point(i,j,k)
				const float distance = openvdb::math::Sqrt((i - intersection_x) * (i - intersection_x) + (j - intersection_y) * (j - intersection_y) + (k - intersection_z) * (k - intersection_z)) - radius;

				//check if point(i,j,k) is within cylinder region
				double condition1 = ((i - x1) * (x2 - x1) + (j - y1) * (y2 - y1) + (k - z1) * (z2 - z1)) / fiberLength;
				double condition2 = ((i - x2) * (x2 - x1) + (j - y2) * (y2 - y1) + (k - z2) * (z2 - z1)) / fiberLength;

				//identify valid points in its own cylinder region
				if (condition1 >= 0 && condition2 <= 0 && distance <= 0)
				{
					//check if intersect
					if (accessor.getValue(ijk) < 0)
					{
						return false;
					}
					else
					{
						if (i > length)
						{
							a = i - length;
							b = j;
							c = k;
							if (j > width)
							{
								b = j - width;
							}
							else if (j < 0)
							{
								b = j + width;
							}
							else if (k > height)
							{
								c = k - height;
							}
							else if (k < 0)
							{
								c = k + height;
							}
							accessor.setValue(abc, distance);
						}
						else if (i < 0)
						{
							a = i + length;
							b = j;
							c = k;
							if (j > width)
							{
								b = j - width;
							}
							else if (j < 0)
							{
								b = j + width;
							}
							else if (k > height)
							{
								c = k - height;
							}
							else if (k < 0)
							{
								c = k + height;
							}
							accessor.setValue(abc, distance);
						}
						else if (j >width)
						{
							a = i;
							b = j - width;
							c = k;
							if (k > height)
							{
								c = k - height;
							}
						    else if (k < 0)
						    {
							    c = k + height;
						    }
							accessor.setValue(abc, distance);
						}
						else if (j < 0)
						{
							a = i;
							b = j + width;
							c = k;
							if (k > height)
							{
								c = k - height;
							}
							else if (k < 0)
							{
								c = k + height;
							}
							accessor.setValue(abc, distance);
						}
						else if (k > height)
						{
							a = i;
							b = j;
							c = k - height;
							accessor.setValue(abc, distance);
						}
						else if (k < 0)
						{
							a = i;
							b = k;
							c = k + height;
							accessor.setValue(abc, distance);
						}
						else
					        accessor.setValue(ijk, distance);
					}
				}
			}
		}
	}
	//set value 0 to matrix region
	/*for (i = 0; i < length; ++i) {
		for (j = 0; j < width; ++j) {
			for (k = 0; k < height; ++k) {
				if (accessor.getValue(ijk) == NULL)
				{
					accessor.setValue(ijk, 0);
				}
			}
		}
	}
	*/
	return true;
}
	
void DendroGrid::RandomCylinder(double length_mm, double width_mm, double height_mm, double num)
{
	//transform global coord to index coord
	double length = length_mm / 0.2;
	double width = width_mm / 0.2;
	double height = height_mm / 0.2;

	//settings
	const float outside = 2;
	const float inside = -outside;
	int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));

	//MAIN "for loop"
	//srand((unsigned)time(NULL));
	std::random_device rd{};
	std::mt19937 gen{ rd() };

	//random length
	std::normal_distribution<double> FIBERLENGTH{ 1 / 0.2, 0.01 / 0.2 };

	// random diameter 
	std::normal_distribution<double> diameter{ 0.36*2 / 0.2, 0.005 / 0.2 };

	// generate the random number
	std::normal_distribution<double> THETA{ 3.14 / 2, 3.14 / 20 };
	std::normal_distribution<double> PHI{ 3.14 / 2, 3.14 / 20 };

	//generate random x,y,z coord of these random points
	std::uniform_real_distribution<double> startX{ 0,length };
	std::uniform_real_distribution<double> startY{ 0,width };
	std::uniform_real_distribution<double> startZ{ 0,height };

	//radnom direction between (0,1)
	std::uniform_real_distribution<double> inverse{ 0,1 };

	int fiberNum;
	fiberNum = 0;
	while( fiberNum<num)
	{

		if (mGrid == NULL)
		{
			mGrid = openvdb::FloatGrid::create(/*background value=*/2.0);
			openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.2);
			mGrid->setTransform(xform);
		}
	

		//random fiber length
		double fiberLength = FIBERLENGTH(gen);

	    //random radius
		double radius = (diameter(gen) / 2);

		//random start point
		double x1 = startX(gen);
		double y1 = startY(gen);
		double z1 = startZ(gen);

		//generate two angles phi and theta within range of (0,2*pi)
		
		double theta = THETA(gen);
		double phi = PHI(gen);
		if (inverse(gen) < 0.5)
		{
			theta = theta + 3.14;
			phi = phi + 3.14;
		}

		//calculate the x,y,z coord of the end point
		double x2 = x1 + fiberLength * sin(theta) * cos(phi);
		double y2 = y1 + fiberLength * sin(theta) * sin(phi);
		double z2 = z1 + fiberLength * cos(theta);

		//mid point from Pa to Pb
		double center_x = 0.5 * (x1 + x2);
		double center_y = 0.5 * (y1 + y2);
		double center_z = 0.5 * (z1 + z2);
		//two types of shape of cylinder
		int dim, coreRegion;
		if (radius > 0.5 * fiberLength)
		{
			dim = (radius + padding);
			coreRegion = 0.5 * fiberLength;

		}
		else
		{
			dim = (0.5 * fiberLength + padding);
			coreRegion = radius;
		}

		openvdb::FloatGrid::Ptr duplicateGrid = mGrid->deepCopy();

		if (checkIntersection( length, width, height, radius, x1, y1, z1, x2, y2, z2, dim, fiberLength, duplicateGrid))
		{
			//no intersection
			mGrid->clear();
			mGrid = duplicateGrid;
			fiberNum = fiberNum + 1;
		}

		else
		{
			//intersection detected
			duplicateGrid->clear();
		}
		
	}	
}

//generate fibers in 20*20*10 box 2021.8.4 after graduation
void DendroGrid::RandomFibers(double aspectRatioMean, double aspectRatioDe, double diameterMean, double diameterDe, double orientationMean, double orientationDe, double volumeFraction)
{
	//transform global coord to index coord
	double length = 20 / 0.2;
	double width = 20 / 0.2;
	double height = 10 / 0.2;

	//settings
	const float outside = 2;
	const float inside = -outside;
	int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));

	//srand((unsigned)time(NULL));
	std::random_device rd{};
	std::mt19937 gen{ rd() };

	//random aspect ratio
	std::normal_distribution<double> aspectRatio{ aspectRatioMean , aspectRatioDe };

	// random diameter 
	std::normal_distribution<double> diameter{ diameterMean / 0.2, diameterDe / 0.2 };

	// generate the random orientation
	std::normal_distribution<double> orientation{ orientationMean,orientationDe };

	//generate random x,y,z coord of these random points
	std::uniform_real_distribution<double> startX{ 0,length };
	std::uniform_real_distribution<double> startY{ 0,width };
	std::uniform_real_distribution<double> startZ{ 0,height };

	//radnom direction between (0,1)
	std::uniform_real_distribution<double> inverse{ 0,1 };

	int fiberNum = 0;
	double vf = 0;

	while (vf < volumeFraction)
	{
		if (mGrid == NULL)
		{
			mGrid = openvdb::FloatGrid::create(/*background value=*/2.0);
			openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/0.2);
			mGrid->setTransform(xform);
		}

		//random fiber length
		double fiberLength = aspectRatio(gen) * diameter(gen);

		//random radius
		double radius = (diameter(gen) / 2);

		//random start point
		double x1 = startX(gen);
		double y1 = startY(gen);
		double z1 = startZ(gen);

		//generate angle theta in degree
		double theta = orientation(gen);

		//calculate the x,y,z coord of the end point
		double x2 = x1 + fiberLength * cos(theta * 3.14 / 180);
		double y2 = y1 + fiberLength * sin(theta * 3.14 / 180);
		double z2 = z1;

		//mid point from Pa to Pb
		double center_x = 0.5 * (x1 + x2);
		double center_y = 0.5 * (y1 + y2);
		double center_z = 0.5 * (z1 + z2);
		//two types of shape of cylinder
		int dim, coreRegion;
		if (radius > 0.5 * fiberLength)
		{
			dim = (radius + padding);
			coreRegion = 0.5 * fiberLength;

		}
		else
		{
			dim = (0.5 * fiberLength + padding);
			coreRegion = radius;
		}

		openvdb::FloatGrid::Ptr duplicateGrid = mGrid->deepCopy();

		if (checkIntersection(length, width, height, radius, x1, y1, z1, x2, y2, z2, dim, fiberLength, duplicateGrid))
		{
			//no intersection
			mGrid->clear();
			mGrid = duplicateGrid;
			fiberNum = fiberNum + 1;
			double fiberVolume = 3.14 * radius * radius * fiberLength;
			double currentvf = fiberVolume / (length * width * height);
			vf = vf + currentvf;
		}
		else
		{
			duplicateGrid->clear();
		}

	}
	
	std::ofstream myfile;
	myfile.open("G:\\voxelFibers.txt");
	std::string heading = "*Heading\n**Job name : Job - 1 Model name : Model - 1\n**Generated by : CADForAM - 1\n*Preprint, echo = YES, model = YES, history = YES, contact = YES\n";
	myfile << heading;

	//Construct node and elements
	openvdb::CoordBBox boundbox(1, 1,1,length, width,height);
	int dimension_x = boundbox.dim().x();
	int dimension_y = boundbox.dim().y();
	int dimension_z = boundbox.dim().z();
	int size_of_node = (dimension_x + 1) * (dimension_y + 1) * (dimension_z + 1);
	openvdb::Vec3d voxelsize_mm = mGrid->transform().voxelSize();

	// Construct element sets
	std::vector<std::vector<int>> listOfElementsets;
	for (int i = 0; i < 2; i++)
	{
		std::vector<int> elementset = std::vector<int>();
		listOfElementsets.push_back(elementset);
	}
	
	// write node
	std::string nodetitle = "**\n** PARTS\n**\n*Part, name = beam1\n*Node";
	myfile << nodetitle << std::endl;
	int node_id = 1;
	openvdb::Vec3d basenode = mGrid->indexToWorld(boundbox.min()) - 0.5 * voxelsize_mm;
	openvdb::Vec3d cnode;
	myfile << basenode.x() << "," << basenode.y() << "," << basenode.z() << std::endl;
	for (int i = 0; i < 100 + 1; i++)
	{
		for (int j = 0; j < 100 + 1; j++)
		{
			for (int k = 0; k < 50 + 1; k++)
			{
				cnode = basenode + openvdb::Vec3d((i * voxelsize_mm.x()), (j * voxelsize_mm.y()), (k * voxelsize_mm.z()));
				myfile << node_id << "," << cnode.x() << "," << cnode.y() << "," << cnode.z() << std::endl;
				node_id++;
			}
		}
	}

	// write element
	std::string element = "*Element, type=C3D8R\n";
	myfile << element;
	openvdb::Coord ijk;
	int elementid = 1;
	auto accessor = mGrid->getAccessor();
	for (int i = 0; i < dimension_x; i++)
	{
		for (int j = 0; j < dimension_y; j++)
		{
			for (int k = 0; k<dimension_z; k++)
			{
				int node0 = 1 + i * ((dimension_y + 1) * (dimension_z + 1)) + j * (dimension_z + 1) + k;
				int node1 = node0 + 1;
				int node2 = node0 + (dimension_z + 1);
				int node3 = node2 + 1;
				int node4 = node0 + ((dimension_z + 1) * (dimension_y + 1));
				int node5 = node4 + 1;
				int node6 = node4 + (dimension_z + 1);
				int node7 = node6 + 1;
				ijk[0] = i + boundbox.min().x();
				ijk[1] = j + boundbox.min().y();
				ijk[2] = k + boundbox.min().z();
				double value = accessor.getValue(ijk);
				if (value < 0)
				{
					listOfElementsets[0].push_back(elementid);
				}

				if (value > 0)
				{
					listOfElementsets[1].push_back(elementid);
				}

				myfile << elementid << "," << node0 << "," << node4 << "," << node6 << "," << node2 << "," << node1 << "," << node5 << "," << node7 << "," << node3 << std::endl;
				elementid++;
			}
		}
	}

	// write set
	for (int indexofesets = 0; indexofesets < 2; indexofesets++)
	{
		std::string elementsetTitle = "*Elset, elset = Set-" + std::to_string(indexofesets + 1);
		myfile << elementsetTitle << std::endl;
		int changeline = 0;
		for (int elementIndex = 0; elementIndex < listOfElementsets[indexofesets].size(); elementIndex++)
		{
			myfile << listOfElementsets[indexofesets][elementIndex];
			changeline++;
			if (changeline == 15)
			{
				myfile << "\n";
				changeline = 0;
			}
			else
			{
				myfile << ",";
			}

		}
		myfile << std::endl;
	}

	// write section
	for (int i = 0; i < 2; i++)
	{
		std::string sectiontitle = "**Section: Section-" + std::to_string(i + 1);
		myfile << sectiontitle << std::endl;
		myfile << "*Solid Section, elset=Set-" << i + 1 << "," << "material=material-" << i + 1 << std::endl;
	}
	myfile << "*End Part" << std::endl;
	myfile << "**" << "\n" << "**" << "\n" << "*Assembly, name = Assembly" << "\n" << "* *\n" << "*Instance, name = BEAM1 - 1, part = beam1\n" << "*End Instance" << "\n" << "**" << std::endl;
	
	//assembly codes start
	// Define boundary condition
	// define reference point
	std::vector<int> referenceNode;
	for (int i = 1; i <= 6; i++)
	{
		myfile << "*Node" << endl;
		double x = 201.003;
		myfile << node_id << " ," << x + i << "," << "  2.5,  " << "2.5  " << endl;
		referenceNode.push_back(node_id);
		node_id++;
	}
	// Reference NSET
	for (int i = 1; i <= 6; i++)
	{
		myfile << "*Nset, nset=RP" << i << endl;
		myfile << referenceNode[i - 1] << endl;
	}

	// Generate Nset for BACK BC
	myfile << "*Nset, nset=backbc, instance=BEAM1-1, generate" << endl;
	myfile << 1 << "," << (dimension_y + 1) * (dimension_z + 1) << "," << endl;
	for (int i = 1; i <= (dimension_y + 1) * (dimension_z + 1); i++)
	{
		myfile << "*Nset, nset=backbc" << i << ", instance=BEAM1-1\n" << i << ", " << endl;
	}

	std::vector<int> backs;
	int linenumber = 0;
	myfile << "*Nset, nset=backs, instance=BEAM1-1" << endl;
	for (int i = 2; i <= dimension_y; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = j + (i - 1) * (dimension_z + 1);
			backs.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < backs.size(); i++)
	{
		myfile << "*Nset, nset = backs" << backs[i] << ", instance = BEAM1 - 1" << endl;
		myfile << backs[i] << endl;
	}

	// Back left edge
	myfile << "\n";
	myfile << "*Nset, nset=bledge, instance=BEAM1-1," << "generate" << endl;
	myfile << (dimension_z + 1) * 2 << "," << (dimension_z + 1) * dimension_y << "," << (dimension_z + 1) << endl;
	std::vector<int> bledge;
	for (int i = 2; i <= dimension_y; i++)
	{
		myfile << "*Nset, nset = bledge" << (dimension_z + 1) * i << ", instance = BEAM1 - 1" << endl;
		myfile << (dimension_z + 1) * i << endl;
		bledge.push_back((dimension_z + 1) * i);
	}
	// Back right edge
	myfile << "\n";
	myfile << "*Nset, nset=bredge, instance=BEAM1-1," << "generate" << endl;
	myfile << (dimension_z + 1) + 1 << "," << (dimension_z + 1) * (dimension_y - 1) + 1 << "," << (dimension_z + 1) << endl;
	std::vector<int> bredge;
	for (int i = 1; i <= dimension_y - 1; i++)
	{
		myfile << "*Nset, nset = bredge" << (dimension_z + 1) * i + 1 << ", instance = BEAM1 - 1" << endl;
		myfile << (dimension_z + 1) * i + 1 << endl;
		bredge.push_back((dimension_z + 1) * i + 1);
	}
	// Back Bottom edge
	myfile << "\n";
	myfile << "*Nset, nset=bbedge, instance=BEAM1-1" << ", generate" << endl;
	myfile << "2," << dimension_z << "," << 1 << endl;
	std::vector<int> bbedge;
	for (int i = 2; i <= dimension_z; i++)
	{
		myfile << "*Nset, nset = bbedge" << i << ", instance = BEAM1 - 1" << endl;
		myfile << i << endl;
		bbedge.push_back(i);
	}
	// Back Top edge
	myfile << "\n";
	myfile << "*Nset, nset=btedge, instance=BEAM1-1" << ", generate" << endl;
	myfile << (dimension_z + 1) * dimension_y + 2 << "," << (dimension_z + 1) * (dimension_y + 1) - 1 << "," << 1 << endl;
	std::vector<int> btedge;
	for (int i = 2; i <= dimension_z; i++)
	{
		myfile << "*Nset, nset = btedge" << i + (dimension_z + 1) * dimension_y << ", instance = BEAM1 - 1" << endl;
		myfile << i + (dimension_z + 1) * dimension_y << endl;
		btedge.push_back(i + (dimension_z + 1) * dimension_y);
	}

	// rightbc
	myfile << "*Nset, nset=rightbc, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> rightbc;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_y + 1; j++)
		{
			int number = 1 + (dimension_z + 1) * (j - 1) + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			rightbc.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < rightbc.size(); i++)
	{
		myfile << "*Nset, nset = rightbc" << rightbc[i] << ", instance = BEAM1 - 1" << endl;
		myfile << rightbc[i] << endl;
	}

	// leftbc
	myfile << "*Nset, nset=leftbc, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> leftbc;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_y + 1; j++)
		{
			int number = (dimension_z + 1) * j + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			leftbc.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < leftbc.size(); i++)
	{
		myfile << "*Nset, nset = leftbc" << leftbc[i] << ", instance = BEAM1 - 1" << endl;
		myfile << leftbc[i] << endl;
	}
	// lefts
	myfile << "*Nset, nset=lefts, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> lefts;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_y; j++)
		{
			int number = (dimension_z + 1) * j + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			lefts.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < lefts.size(); i++)
	{
		myfile << "*Nset, nset = lefts" << lefts[i] << ", instance = BEAM1 - 1" << endl;
		myfile << lefts[i] << endl;
	}
	// rights
	myfile << "*Nset, nset=rights, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> rights;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_y; j++)
		{
			int number = 1 + (dimension_z + 1) * (j - 1) + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			rights.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < rights.size(); i++)
	{
		myfile << "*Nset, nset = rights" << rights[i] << ", instance = BEAM1 - 1" << endl;
		myfile << rights[i] << endl;
	}
	// leftopedge
	myfile << "*Nset, nset=ltedge, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> ltedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * i;
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
		ltedge.push_back(number);
	}
	myfile << "\n";
	for (int i = 0; i < ltedge.size(); i++)
	{
		myfile << "*Nset, nset = ltedge" << ltedge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << ltedge[i] << endl;
	}
	// Leftbottomedge
	myfile << "*Nset, nset=lbedge, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> lbedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * (i - 1) + dimension_z + 1;
		lbedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < lbedge.size(); i++)
	{
		myfile << "*Nset, nset = lbedge" << lbedge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << lbedge[i] << endl;
	}
	// righttop
	myfile << "*Nset, nset=rtedge, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> rtedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * (i - 1) + (dimension_z + 1) * dimension_y + 1;
		rtedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < rtedge.size(); i++)
	{
		myfile << "*Nset, nset = rtedge" << rtedge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << rtedge[i] << endl;
	}
	// rightbottom
	myfile << "*Nset, nset=rtedge, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> rbedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * (i - 1) + 1;
		rbedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < rbedge.size(); i++)
	{
		myfile << "*Nset, nset = rbedge" << rbedge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << rbedge[i] << endl;
	}
	// frontbc
	myfile << "*Nset, nset=frontbc, instance=BEAM1-1,generate" << endl;
	myfile << (dimension_y + 1) * (dimension_z + 1) * dimension_x + 1 << "," << (dimension_y + 1) * (dimension_z + 1) * (dimension_x + 1) << "," << 1 << endl;

	for (int i = 0; i < (dimension_y + 1) * (dimension_z + 1); i++)
	{
		myfile << "*Nset, nset = frontbc" << (dimension_y + 1) * (dimension_z + 1) * dimension_x + 1 + i << ", instance = BEAM1 - 1" << endl;
		myfile << (dimension_y + 1) * (dimension_z + 1) * dimension_x + 1 + i << endl;
	}
	// fronts
	myfile << "*Nset, nset=fronts, instance=BEAM1-1" << endl;
	linenumber = 0;
	std::vector<int> fronts;
	for (int i = 2; i <= dimension_y; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * dimension_x + j + (i - 1) * (dimension_z + 1);
			myfile << number;
			fronts.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < fronts.size(); i++)
	{
		myfile << "*Nset, nset = fronts" << fronts[i] << ", instance = BEAM1 - 1" << endl;
		myfile << fronts[i] << endl;
	}
	//fronttop
	std::vector<int> ftedge;
	myfile << "*Nset, nset=ftedge, instance=BEAM1-1" << endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_z; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + (dimension_z + 1) * dimension_y + i;
		ftedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < ftedge.size(); i++)
	{
		myfile << "*Nset, nset = ftedge" << ftedge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << ftedge[i] << endl;
	}
	//frontbottom
	std::vector<int> fbedge;
	myfile << "*Nset, nset=fbedge, instance=BEAM1-1" << endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_z; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + i;
		fbedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < fbedge.size(); i++)
	{
		myfile << "*Nset, nset = fbedge" << fbedge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << fbedge[i] << endl;
	}
	//frontleft
	std::vector<int> fledge;
	myfile << "*Nset, nset=fledge, instance=BEAM1-1" << endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_y; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + i * (dimension_z + 1);
		fledge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < fledge.size(); i++)
	{
		myfile << "*Nset, nset = fledge" << fledge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << fledge[i] << endl;
	}
	//frontright
	std::vector<int> fredge;
	myfile << "*Nset, nset=fredge, instance=BEAM1-1" << endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_y; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + (i - 1) * (dimension_z + 1) + 1;
		fredge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < fredge.size(); i++)
	{
		myfile << "*Nset, nset = fredge" << fredge[i] << ", instance = BEAM1 - 1" << endl;
		myfile << fredge[i] << endl;
	}
	// tops
	std::vector<int> tops;
	linenumber = 0;
	myfile << "*Nset, nset=tops, instance=BEAM1-1" << endl;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + (dimension_z + 1) * dimension_y + j;
			tops.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < tops.size(); i++)
	{
		myfile << "*Nset, nset = tops" << tops[i] << ", instance = BEAM1 - 1" << endl;
		myfile << tops[i] << endl;
	}
	// Bottom
	std::vector<int> bots;
	linenumber = 0;
	myfile << "*Nset, nset=bots, instance=BEAM1-1" << endl;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + j;
			bots.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < bots.size(); i++)
	{
		myfile << "*Nset, nset = bots" << bots[i] << ", instance = BEAM1 - 1" << endl;
		myfile << bots[i] << endl;
	}
	// topbc
	std::vector<int> topbc;
	linenumber = 0;
	myfile << "*Nset, nset=topbc, instance=BEAM1-1" << endl;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_z + 1; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + (dimension_z + 1) * dimension_y + j;
			topbc.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < topbc.size(); i++)
	{
		myfile << "*Nset, nset = topbc" << topbc[i] << ", instance = BEAM1 - 1" << endl;
		myfile << topbc[i] << endl;
	}
	// Bottom
	std::vector<int> botbc;
	linenumber = 0;
	myfile << "*Nset, nset=botbc, instance=BEAM1-1" << endl;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_z + 1; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + j;
			botbc.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < botbc.size(); i++)
	{
		myfile << "*Nset, nset = botbc" << botbc[i] << ", instance = BEAM1 - 1" << endl;
		myfile << botbc[i] << endl;
	}


	myfile << "*Nset, nset = c1" << ", instance = BEAM1 - 1" << endl;
	myfile << (dimension_x + 1) * (dimension_z + 1) * (dimension_y + 1) << endl;
	myfile << "*Nset, nset = c2" << ", instance = BEAM1 - 1" << endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) << endl;
	myfile << "*Nset, nset = c3" << ", instance = BEAM1 - 1" << endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) - dimension_z << endl;
	myfile << "*Nset, nset = c4" << ", instance = BEAM1 - 1" << endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) * (dimension_x + 1) - dimension_z << endl;
	myfile << "*Nset, nset = c5" << ", instance = BEAM1 - 1" << endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) * (dimension_x)+dimension_z + 1 << endl;
	myfile << "*Nset, nset = c6" << ", instance = BEAM1 - 1" << endl;
	myfile << dimension_z + 1 << endl;
	myfile << "*Nset, nset = c7" << ", instance = BEAM1 - 1" << endl;
	myfile << 1 << endl;
	myfile << "*Nset, nset = c8" << ", instance = BEAM1 - 1" << endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) * (dimension_x)+1 << endl;


	// Define Equation
	// Tops Bots
	for (int i = 0; i < tops.size(); i++)
	{
		WriteConstraint(tops[i], bots[i], 1, 1.0, -1.0, 0.0, "na", myfile, "tops", "bots");
		WriteConstraint(tops[i], bots[i], 2, 1.0, -1.0, -1.0, "RP5", myfile, "tops", "bots");
		WriteConstraint(tops[i], bots[i], 3, 1.0, -1.0, 0.0, "na", myfile, "tops", "bots");
	}

	// Fronts Backs
	for (int i = 0; i < fronts.size(); i++)
	{
		WriteConstraint(fronts[i], backs[i], 1, 1.0, -1.0, -1.0, "RP4", myfile, "fronts", "backs");
		WriteConstraint(fronts[i], backs[i], 2, 1.0, -1.0, 0.0, "na", myfile, "fronts", "backs");
		WriteConstraint(fronts[i], backs[i], 3, 1.0, -1.0, 0.0, "na", myfile, "fronts", "backs");
	}
	// Left Right
	for (int i = 0; i < lefts.size(); i++)
	{
		WriteConstraint(lefts[i], rights[i], 1, 1.0, -1.0, 0.0, "na", myfile, "lefts", "rights");
		WriteConstraint(lefts[i], rights[i], 2, 1.0, -1.0, 0.0, "na", myfile, "lefts", "rights");
		WriteConstraint(lefts[i], rights[i], 3, 1.0, -1.0, -1.0, "RP6", myfile, "lefts", "rights");

	}
	// ft bt
	for (int i = 0; i < ftedge.size(); i++)
	{
		WriteConstraint(ftedge[i], btedge[i], 1, 1.0, -1.0, -1.0, "RP4", myfile, "ftedge", "btedge");
		WriteConstraint(ftedge[i], btedge[i], 2, 1.0, -1.0, 0.0, "na", myfile, "ftedge", "btedge");
		WriteConstraint(ftedge[i], btedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "ftedge", "btedge");
	}
	// bt bb
	for (int i = 0; i < btedge.size(); i++)
	{
		WriteConstraint(btedge[i], bbedge[i], 1, 1.0, -1.0, 0.0, "na", myfile, "btedge", "bbedge");
		WriteConstraint(btedge[i], bbedge[i], 2, 1.0, -1.0, -1.0, "RP5", myfile, "btedge", "bbedge");
		WriteConstraint(btedge[i], bbedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "btedge", "bbedge");
	}
	// BB fb
	for (int i = 0; i < bbedge.size(); i++)
	{
		WriteConstraint(bbedge[i], fbedge[i], 1, 1.0, -1.0, 1.0, "RP4", myfile, "bbedge", "fbedge");
		WriteConstraint(bbedge[i], fbedge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "bbedge", "fbedge");
		WriteConstraint(bbedge[i], fbedge[i], 3, 1.0, -1.0, 0.0, "NA", myfile, "bbedge", "fbedge");

	}
	// FL BL
	for (int i = 0; i < fledge.size(); i++)
	{
		WriteConstraint(fledge[i], bledge[i], 1, 1.0, -1.0, -1.0, "RP4", myfile, "fledge", "bledge");
		WriteConstraint(fledge[i], bledge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "fledge", "bledge");
		WriteConstraint(fledge[i], bledge[i], 3, 1.0, -1.0, 0.0, "NA", myfile, "fledge", "bledge");

	}
	//BL BR
	for (int i = 0; i < bledge.size(); i++)
	{
		WriteConstraint(bledge[i], bredge[i], 1, 1.0, -1.0, 0.0, "na", myfile, "bledge", "bredge");
		WriteConstraint(bledge[i], bredge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "bledge", "bredge");
		WriteConstraint(bledge[i], bredge[i], 3, 1.0, -1.0, -1.0, "RP6", myfile, "bledge", "bredge");

	}
	//br fr
	for (int i = 0; i < bredge.size(); i++)
	{
		WriteConstraint(bredge[i], fredge[i], 1, 1.0, -1.0, 1.0, "RP4", myfile, "bredge", "fredge");
		WriteConstraint(bredge[i], fredge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "bredge", "fredge");
		WriteConstraint(bredge[i], fredge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "bredge", "fredge");
	}
	//LT LB
	for (int i = 0; i < ltedge.size(); i++)
	{
		WriteConstraint(ltedge[i], lbedge[i], 1, 1.0, -1.0, 0.0, "na", myfile, "ltedge", "lbedge");
		WriteConstraint(ltedge[i], lbedge[i], 2, 1.0, -1.0, -1.0, "RP5", myfile, "ltedge", "lbedge");
		WriteConstraint(ltedge[i], lbedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "ltedge", "lbedge");
	}
	//LB RB
	for (int i = 0; i < lbedge.size(); i++)
	{
		WriteConstraint(lbedge[i], rbedge[i], 1, 1.0, -1.0, 0.0, "0", myfile, "lbedge", "rbedge");
		WriteConstraint(lbedge[i], rbedge[i], 2, 1.0, -1.0, 0.0, "0", myfile, "lbedge", "rbedge");
		WriteConstraint(lbedge[i], rbedge[i], 3, 1.0, -1.0, -1.0, "RP6", myfile, "lbedge", "rbedge");
	}
	//RB RT
	for (int i = 0; i < rbedge.size(); i++)
	{
		WriteConstraint(rbedge[i], rtedge[i], 1, 1.0, -1.0, 0.0, "0", myfile, "rbedge", "rtedge");
		WriteConstraint(rbedge[i], rtedge[i], 2, 1.0, -1.0, 1.0, "RP5", myfile, "rbedge", "rtedge");
		WriteConstraint(rbedge[i], rtedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "rbedge", "rtedge");
	}
	//62
	WriteConstraint(6, 2, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(6, 2, 2, 1.0, -1.0, 1.0, "RP5", myfile, "c", "c");
	WriteConstraint(6, 2, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//23
	WriteConstraint(2, 3, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(2, 3, 2, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(2, 3, 3, 1.0, -1.0, -1.0, "RP6", myfile, "c", "c");

	//34
	WriteConstraint(3, 4, 1, 1.0, -1.0, 1.0, "RP4", myfile, "c", "c");
	WriteConstraint(3, 4, 2, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(3, 4, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//48
	WriteConstraint(4, 8, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(4, 8, 2, 1.0, -1.0, -1.0, "RP5", myfile, "c", "c");
	WriteConstraint(4, 8, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//85
	WriteConstraint(8, 5, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(8, 5, 2, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(8, 5, 3, 1.0, -1.0, 1.0, "RP6", myfile, "c", "c");
	//51
	WriteConstraint(5, 1, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(5, 1, 2, 1.0, -1.0, 1.0, "RP5", myfile, "c", "c");
	WriteConstraint(5, 1, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//17
	WriteConstraint(1, 7, 1, 1.0, -1.0, -1.0, "RP4", myfile, "c", "c");
	WriteConstraint(1, 7, 2, 1.0, -1.0, -1.0, "RP5", myfile, "c", "c");
	WriteConstraint(1, 7, 3, 1.0, -1.0, -1.0, "RP6", myfile, "c", "c");

	//end assembly
	myfile << "*End Assembly" << std::endl;

	// write materials
	std::vector<double> elasticProperties = { 10,20 };
	std::vector<double> poissionratio = { 0.33,0.33 };
	std::string materialtitle = "**\n**MATERIALS\n**";
	myfile << materialtitle << std::endl;
	for (int indexofm = 0; indexofm < 2; indexofm++)
	{
		std::string material = "*Material, name = material-" + std::to_string(indexofm + 1) + "\n*Elastic\n" + std::to_string(elasticProperties[indexofm]) + "," + std::to_string(poissionratio[indexofm]);
		myfile << material << std::endl;
	}
	
	myfile.close();
}
				

//this is a basic step of how to generate random numbers,
//formula1 is r = rand()%(n - m + 1) + m;
//formula2 is value = min + (double)rand() * (max - min) / (double)RAND_MAX;
//notice: srand should be out of for loop

/*int n, m, x1
srand((unsigned)time(NULL));
for (int i = 0; i < 5; i++) {

		n = 100;
		m = 0;
		x1 = rand() % (n - m + 1) + m;
		cout << x1 << endl;

}*/  //formula1

/*int main()
{
	srand((unsigned)time(NULL));
	int min = 0;
	int max = 100;
	double value;
	for (int i = 0; i < 10; ++i)
	{
		value = min + (double)rand() * (max - min) / (double)RAND_MAX;
		cout << value << endl;
	}
}*/  //formula2


//this is a new C++ random engine
//#include<random>
//#include<iostream>
//#include<ctime>
//	using namespace std;
//
//	int main()
//	{
//		default_random_engine e(time(0));
//		uniform_real_distribution<double> u(-1.2, 3.5);//change ¡°real¡± to int£¬  <> filled "unsigned"£¬and it can generate integers
//		for (int i = 0; i < 10; ++i)
//			cout << u(e) << endl;
//		return 0;
//	}

void DendroGrid::Transform(openvdb::math::Mat4d xform)
{
	mGrid->transform().postMult(xform);
}

void DendroGrid::BooleanUnion(DendroGrid vAdd)
{
	auto csgGrid = vAdd.Grid();

	// store current tranforms of both csg volumes
	const openvdb::math::Transform
		&sourceXform = csgGrid->transform(),
		&targetXform = mGrid->transform();

	// create a copy of the source grid for resampling
	openvdb::FloatGrid::Ptr cGrid = openvdb::createLevelSet<openvdb::FloatGrid>(mGrid->voxelSize()[0]);
	cGrid->transform() = mGrid->transform();

	// compute a source grid to target grid transform
	openvdb::Mat4R xform =
		sourceXform.baseMap()->getAffineMap()->getMat4() *
		targetXform.baseMap()->getAffineMap()->getMat4().inverse();

	// create the transformer
	openvdb::tools::GridTransformer transformer(xform);

	// resample using trilinear interpolation 
	transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*csgGrid, *cGrid);

	// solve for the csg operation with result being stored in mGrid
	openvdb::tools::csgUnion(*mGrid, *cGrid, true);
}

void DendroGrid::BooleanIntersection(DendroGrid vIntersect)
{
	auto csgGrid = vIntersect.Grid();

	// store current tranforms of both csg volumes
	const openvdb::math::Transform
		&sourceXform = csgGrid->transform(),
		&targetXform = mGrid->transform();

	// create a copy of the source grid for resampling
	openvdb::FloatGrid::Ptr cGrid = openvdb::createLevelSet<openvdb::FloatGrid>(mGrid->voxelSize()[0]);
	cGrid->transform() = mGrid->transform();

	// compute a source grid to target grid transform
	openvdb::Mat4R xform =
		sourceXform.baseMap()->getAffineMap()->getMat4() *
		targetXform.baseMap()->getAffineMap()->getMat4().inverse();

	// create the transformer
	openvdb::tools::GridTransformer transformer(xform);

	// resample using trilinear interpolation 
	transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*csgGrid, *cGrid);

	// solve for the csg operation with result being stored in mGrid
	openvdb::tools::csgIntersection(*mGrid, *cGrid, true);
}

void DendroGrid::BooleanDifference(DendroGrid vSubtract)
{
	auto csgGrid = vSubtract.Grid();

	// store current tranforms of both csg volumes
	const openvdb::math::Transform
		&sourceXform = csgGrid->transform(),
		&targetXform = mGrid->transform();

	// create a copy of the source grid for resampling
	openvdb::FloatGrid::Ptr cGrid = openvdb::createLevelSet<openvdb::FloatGrid>(mGrid->voxelSize()[0]);
	cGrid->transform() = mGrid->transform();

	// compute a source grid to target grid transform
	openvdb::Mat4R xform =
		sourceXform.baseMap()->getAffineMap()->getMat4() *
		targetXform.baseMap()->getAffineMap()->getMat4().inverse();

	// create the transformer
	openvdb::tools::GridTransformer transformer(xform);

	// resample using trilinear interpolation 
	transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*csgGrid, *cGrid);

	// solve for the csg operation with result being stored in mGrid
	openvdb::tools::csgDifference(*mGrid, *cGrid, true);
}

void DendroGrid::Offset(double amount)
{
	// create a new filter to operate on grid with
	openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*mGrid);

	filter.setGrainSize(1);

	amount = amount * -1;

	// apply offset to grid of supplied amount
	filter.offset((float)amount);
}

void DendroGrid::Offset(double amount, DendroGrid vMask, double min, double max, bool invert)
{
	// create a new filter to operate on grid with
	openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*mGrid);

	filter.invertMask(invert);
	filter.setMaskRange((float)min, (float)max);
	filter.setGrainSize(1);

	// create filter mask
	openvdb::Grid<openvdb::FloatTree> mMask(*vMask.Grid());

	amount = amount * -1;

	// apply offset to grid of supplied amount
	filter.offset((float)amount, &mMask);
}

void DendroGrid::Smooth(int type, int iterations, int width)
{
	// create a new filter to operate on grid with
	openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*mGrid);
	filter.setGrainSize(1);

	// apply filter for the number iterations supplied
	for (int i = 0; i < iterations; i++) {

		// filter by desired type supplied
		switch (type) {
		case 0:
			filter.gaussian(width);
			break;
		case 1:
			filter.laplacian();
			break;
		case 2:
			filter.mean(width);
			break;
		case 3:
			filter.median(width);
			break;
		default:
			filter.laplacian();
			break;
		}
	}
}

void DendroGrid::Smooth(int type, int iterations, int width, DendroGrid vMask, double min, double max, bool invert)
{
	// create a new filter to operate on grid with
	openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*mGrid);

	filter.invertMask(invert);
	filter.setMaskRange((float)min, (float)max);
	filter.setGrainSize(1);

	// create filter mask
	openvdb::Grid<openvdb::FloatTree> mMask(*vMask.Grid());

	// apply filter for the number iterations supplied
	for (int i = 0; i < iterations; i++) {

		// filter by desired type supplied
		switch (type) {
		case 0:
			filter.gaussian(width, &mMask);
			break;
		case 1:
			filter.laplacian(&mMask);
			break;
		case 2:
			filter.mean(width, &mMask);
			break;
		case 3:
			filter.median(width, &mMask);
			break;
		default:
			filter.laplacian(&mMask);
			break;
		}
	}
}

void DendroGrid::Blend(DendroGrid bGrid, double bPosition, double bEnd)
{
	openvdb::tools::LevelSetMorphing<openvdb::FloatGrid> morph(*mGrid, *bGrid.Grid());
	morph.setSpatialScheme(openvdb::math::HJWENO5_BIAS);
	morph.setTemporalScheme(openvdb::math::TVD_RK3);
	morph.setTrackerSpatialScheme(openvdb::math::HJWENO5_BIAS);
	morph.setTrackerTemporalScheme(openvdb::math::TVD_RK2);
	morph.setGrainSize(1);

	double bStart = bPosition * bEnd;
	morph.advect(bStart, bEnd);
}

void DendroGrid::Blend(DendroGrid bGrid, double bPosition, double bEnd, DendroGrid vMask, double mMin, double mMax, bool invert)
{
	openvdb::tools::LevelSetMorphing<openvdb::FloatGrid> morph(*mGrid, *bGrid.Grid());
	morph.setSpatialScheme(openvdb::math::HJWENO5_BIAS);
	morph.setTemporalScheme(openvdb::math::TVD_RK3);
	morph.setTrackerSpatialScheme(openvdb::math::HJWENO5_BIAS);
	morph.setTrackerTemporalScheme(openvdb::math::TVD_RK2);

	morph.setAlphaMask(*vMask.Grid());
	morph.invertMask(invert);
	morph.setMaskRange((float)mMin, (float)mMax);
	morph.setGrainSize(1);

	double bStart = bPosition * bEnd;
	morph.advect(bStart, bEnd);
}

void DendroGrid::ClosestPoint(std::vector<openvdb::Vec3R>& points, std::vector<float>& distances)
{
	auto csp = openvdb::tools::ClosestSurfacePoint<openvdb::FloatGrid>::create(*mGrid);
	csp->searchAndReplace(points, distances);
}

DendroMesh DendroGrid::Display()
{
	return mDisplay;
}

void DendroGrid::UpdateDisplay()
{
	using openvdb::Index64;

	openvdb::tools::VolumeToMesh mesher(mGrid->getGridClass() == openvdb::GRID_LEVEL_SET ? 0.0 : 0.01);
	mesher(*mGrid);

	mDisplay.Clear();

	for (Index64 n = 0, i = 0, N = mesher.pointListSize(); n < N; ++n)
	{
		auto v = mesher.pointList()[n];
		mDisplay.AddVertice(v);
	}

	openvdb::tools::PolygonPoolList &polygonPoolList = mesher.polygonPoolList();

	for (Index64 n = 0, N = mesher.polygonPoolListSize(); n < N; ++n)
	{
		const openvdb::tools::PolygonPool &polygons = polygonPoolList[n];
		for (Index64 i = 0, I = polygons.numQuads(); i < I; ++i)
		{
			auto face = polygons.quad(i);
			mDisplay.AddFace(face);
		}
	}
}

void DendroGrid::UpdateDisplay(double isovalue, double adaptivity)
{
	isovalue /= mGrid->voxelSize().x();

	std::vector<openvdb::Vec3s> points;
	std::vector<openvdb::Vec4I> quads;
	std::vector<openvdb::Vec3I> triangles;

	openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*mGrid, points, triangles, quads, isovalue, adaptivity);

	mDisplay.Clear();

	mDisplay.AddVertice(points);

	auto begin = triangles.begin();
	auto end = triangles.end();

	for (auto it = begin; it != end; ++it) {
		int w = -1;
		int x = it->x();
		int y = it->y();
		int z = it->z();

		openvdb::Vec4I face(x,y,z,w);

		mDisplay.AddFace(face);
	}

	mDisplay.AddFace(quads);
}

float * DendroGrid::GetMeshVertices()
{
	auto vertices = mDisplay.Vertices();
	
	mVertexCount = vertices.size() * 3;

	float *verticeArray = reinterpret_cast<float*>(malloc(mVertexCount * sizeof(float)));

	int i = 0;
	for (auto it = vertices.begin(); it != vertices.end(); ++it) {
		verticeArray[i] = it->x();
		verticeArray[i + 1] = it->y();
		verticeArray[i + 2] = it->z();
		i += 3;
	}

	return verticeArray;
}

int * DendroGrid::GetMeshFaces()
{
	auto faces = mDisplay.Faces();

	mFaceCount = faces.size() * 4;

	int *faceArray = reinterpret_cast<int*>(malloc(mFaceCount * sizeof(int)));

	int i = 0;
	for (auto it = faces.begin(); it != faces.end(); ++it) {
		faceArray[i] = it->w();
		faceArray[i + 1] = it->x();
		faceArray[i + 2] = it->y();
		faceArray[i + 3] = it->z();
		i += 4;
	}

	return faceArray;
}

int DendroGrid::GetVertexCount()
{
	return mVertexCount;
}

int DendroGrid::GetFaceCount()
{
	return mFaceCount;
}



void DendroGrid::WriteConstraint(int node1, int node2, int dof, double c1, double c2, double c3, std::string RP, std::ofstream& file, std::string nodeset1, std::string nodeset2)
{
	
		file << "**Constraint: E - " << nodeset1 << " - " << nodeset2 << " " << node1 << std::endl;
		file << "*Equation" << std::endl;
		if (abs(c3) < 1e-5)
		{
			file << 2 << std::endl;
			file << nodeset1 << node1 << ", " << dof << ", " << c1 << std::endl;
			file << nodeset2 << node2 << ", " << dof << ", " << c2 << std::endl;
		}
		else
		{
			file << 3 << std::endl;
			file << nodeset1 << node1 << ", " << dof << ", " << c1 << std::endl;
			file << nodeset2 << node2 << ", " << dof << ", " << c2 << std::endl;
			file << RP << ", " << dof << ", " << c3 << std::endl;
		}
		
}




void DendroGrid::WriteGeometricModelToAbaqus(const std::string file_name, std::vector<double> isovalue, std::vector<double> elasticProperties, std::vector<double> poissionratio, std::vector<double> density, openvdb::CoordBBox box)
{
	auto geometry_grid = this->mGrid;
	std::ofstream myfile;
	myfile.open(file_name);
	std::string heading = "*Heading\n**Job name : Job - 1 Model name : Model - 1\n**Generated by : CADForAM - 1\n*Preprint, echo = YES, model = YES, history = YES, contact = YES\n";
	myfile << heading;
	// Construct node and element
	openvdb::CoordBBox boundbox = box;
	int dimension_x = boundbox.dim().x();
	int dimension_y = boundbox.dim().y();
	int dimension_z = boundbox.dim().z();
	int size_of_node = (dimension_x + 1) * (dimension_y + 1) * (dimension_z + 1);
	openvdb::Vec3d voxelsize_mm = geometry_grid->transform().voxelSize();

	// Construct element sets
	std::vector<std::vector<int>> listOfElementsets;
	for (int i = 0; i < density.size(); i++)
	{
		std::vector<int> elementset = std::vector<int>();
		listOfElementsets.push_back(elementset);
	}
	// write node
	std::string nodetitle = "**\n** PARTS\n**\n*Part, name = beam1\n*Node";
	myfile << nodetitle << std::endl;
	int node_id = 1;
	openvdb::Vec3d basenode = geometry_grid->indexToWorld(boundbox.min()) - 0.5 * voxelsize_mm;
	openvdb::Vec3d cnode;
	for (int i = 0; i < dimension_x + 1; i++)
	{
		for (int j = 0; j < dimension_y + 1; j++)
		{
			for (int k = 0; k < dimension_z + 1; k++)
			{
				cnode = basenode + openvdb::Vec3d((i * voxelsize_mm.x()), (j * voxelsize_mm.y()), (k * voxelsize_mm.z()));
				myfile << node_id << "," << cnode.x() << "," << cnode.y() << "," << cnode.z() << std::endl;
				node_id++;
			}
		}
	}
	// write element
	std::string element = "*Element, type=C3D8R\n";
	myfile << element;
	openvdb::Coord ijk;
	int elementid = 1;
	auto accessor = geometry_grid->getAccessor();
	for (int i = 0; i < dimension_x; i++)
	{
		for (int j = 0; j < dimension_y; j++)
		{
			for (int k = 0; k < dimension_z; k++)
			{
				int node0 = 1 + i * ((dimension_y + 1) * (dimension_z + 1)) + j * (dimension_z + 1) + k;
				int node1 = node0 + 1;
				int node2 = node0 + (dimension_z + 1);
				int node3 = node2 + 1;
				int node4 = node0 + ((dimension_z + 1) * (dimension_y + 1));
				int node5 = node4 + 1;
				int node6 = node4 + (dimension_z + 1);
				int node7 = node6 + 1;
				ijk[0] = i + boundbox.min().x();
				ijk[1] = j + boundbox.min().y();
				ijk[2] = k + boundbox.min().z();
				double value = accessor.getValue(ijk);
				for (int indexofIso = 0; indexofIso < isovalue.size(); indexofIso++)
				{
					if (indexofIso == 0)
					{
						if (value < isovalue[indexofIso])
						{
							listOfElementsets[indexofIso].push_back(elementid);
							break;
						}
					}
					if (indexofIso == isovalue.size() - 1)
					{
						if (value > isovalue[indexofIso])
						{
							listOfElementsets[indexofIso + 1].push_back(elementid);
							break;
						}
					}
					if (value < isovalue[indexofIso] && value >= isovalue[indexofIso - 1])
					{
						listOfElementsets[indexofIso].push_back(elementid);
						break;
					}




				}
				myfile << elementid << "," << node0 << "," << node4 << "," << node6 << "," << node2 << "," << node1 << "," << node5 << "," << node7 << "," << node3 << std::endl;
				elementid++;
			}
		}
	}

	// write set

	for (int indexofesets = 0; indexofesets < isovalue.size() + 1; indexofesets++)
	{
		std::string elementsetTitle = "*Elset, elset = Set-" + std::to_string(indexofesets + 1);
		myfile << elementsetTitle << std::endl;
		int changeline = 0;
		for (int elementIndex = 0; elementIndex < listOfElementsets[indexofesets].size(); elementIndex++)
		{
			myfile << listOfElementsets[indexofesets][elementIndex];
			changeline++;
			if (changeline == 15)
			{
				myfile << "\n";
				changeline = 0;
			}
			else
			{
				myfile << ",";
			}

		}
		myfile << std::endl;
	}
	// write section
	for (int i = 0; i < isovalue.size() + 1; i++)
	{
		std::string sectiontitle = "**Section: Section-" + std::to_string(i + 1);
		myfile << sectiontitle << std::endl;
		myfile << "*Solid Section, elset=Set-" << i + 1 << "," << "material=material-" << i + 1 << std::endl;
	}
	myfile << "*End Part" << std::endl;
	myfile << "**" << "\n" << "**" << "\n" << "*Assembly, name = Assembly" << "\n" << "* *\n" << "*Instance, name = BEAM1 - 1, part = beam1\n" << "*End Instance" << "\n" << "**" << std::endl;

	// Define boundary condition
	// define reference point
	std::vector<int> referenceNode;
	for (int i = 1; i <= 6; i++)
	{
		myfile << "*Node" << std::endl;
		double x = 201.003;
		myfile << node_id << " ," << x + i << "," << "  2.5,  " << "2.5  " << std::endl;
		referenceNode.push_back(node_id);
		node_id++;
	}
	// Reference NSET
	for (int i = 1; i <= 6; i++)
	{
		myfile << "*Nset, nset=RP" << i << std::endl;
		myfile << referenceNode[i - 1] << std::endl;
	}

	// Generate Nset for BACK BC
	myfile << "*Nset, nset=backbc, instance=BEAM1-1, generate" <<std:: endl;
	myfile << 1 << "," << (dimension_y + 1) * (dimension_z + 1) << "," << std::endl;
	for (int i = 1; i <= (dimension_y + 1) * (dimension_z + 1); i++)
	{
		myfile << "*Nset, nset=backbc" << i << ", instance=BEAM1-1\n" << i << ", " <<std:: endl;
	}

	std::vector<int> backs;
	int linenumber = 0;
	myfile << "*Nset, nset=backs, instance=BEAM1-1" <<std:: endl;
	for (int i = 2; i <= dimension_y; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = j + (i - 1) * (dimension_z + 1);
			backs.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < backs.size(); i++)
	{
		myfile << "*Nset, nset = backs" << backs[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << backs[i] << std::endl;
	}

	// Back left edge
	myfile << "\n";
	myfile << "*Nset, nset=bledge, instance=BEAM1-1," << "generate" << std::endl;
	myfile << (dimension_z + 1) * 2 << "," << (dimension_z + 1) * dimension_y << "," << (dimension_z + 1) << std::endl;
	std::vector<int> bledge;
	for (int i = 2; i <= dimension_y; i++)
	{
		myfile << "*Nset, nset = bledge" << (dimension_z + 1) * i << ", instance = BEAM1 - 1" << std::endl;
		myfile << (dimension_z + 1) * i << std::endl;
		bledge.push_back((dimension_z + 1) * i);
	}
	// Back right edge
	myfile << "\n";
	myfile << "*Nset, nset=bredge, instance=BEAM1-1," << "generate" << std::endl;
	myfile << (dimension_z + 1) + 1 << "," << (dimension_z + 1) * (dimension_y - 1) + 1 << "," << (dimension_z + 1) << std::endl;
	std::vector<int> bredge;
	for (int i = 1; i <= dimension_y - 1; i++)
	{
		myfile << "*Nset, nset = bredge" << (dimension_z + 1) * i + 1 << ", instance = BEAM1 - 1" <<std:: endl;
		myfile << (dimension_z + 1) * i + 1 << std::endl;
		bredge.push_back((dimension_z + 1) * i + 1);
	}
	// Back Bottom edge
	myfile << "\n";
	myfile << "*Nset, nset=bbedge, instance=BEAM1-1" << ", generate" <<std:: endl;
	myfile << "2," << dimension_z << "," << 1 <<std:: endl;
	std::vector<int> bbedge;
	for (int i = 2; i <= dimension_z; i++)
	{
		myfile << "*Nset, nset = bbedge" << i << ", instance = BEAM1 - 1" << std::endl;
		myfile << i << std::endl;
		bbedge.push_back(i);
	}
	// Back Top edge
	myfile << "\n";
	myfile << "*Nset, nset=btedge, instance=BEAM1-1" << ", generate" <<std:: endl;
	myfile << (dimension_z + 1) * dimension_y + 2 << "," << (dimension_z + 1) * (dimension_y + 1) - 1 << "," << 1 << std::endl;
	std::vector<int> btedge;
	for (int i = 2; i <= dimension_z; i++)
	{
		myfile << "*Nset, nset = btedge" << i + (dimension_z + 1) * dimension_y << ", instance = BEAM1 - 1" << std::endl;
		myfile << i + (dimension_z + 1) * dimension_y << std::endl;
		btedge.push_back(i + (dimension_z + 1) * dimension_y);
	}

	// rightbc
	myfile << "*Nset, nset=rightbc, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> rightbc;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_y + 1; j++)
		{
			int number = 1 + (dimension_z + 1) * (j - 1) + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			rightbc.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < rightbc.size(); i++)
	{
		myfile << "*Nset, nset = rightbc" << rightbc[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << rightbc[i] << std::endl;
	}

	// leftbc
	myfile << "*Nset, nset=leftbc, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> leftbc;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_y + 1; j++)
		{
			int number = (dimension_z + 1) * j + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			leftbc.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < leftbc.size(); i++)
	{
		myfile << "*Nset, nset = leftbc" << leftbc[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << leftbc[i] << std::endl;
	}
	// lefts
	myfile << "*Nset, nset=lefts, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> lefts;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_y; j++)
		{
			int number = (dimension_z + 1) * j + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			lefts.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < lefts.size(); i++)
	{
		myfile << "*Nset, nset = lefts" << lefts[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << lefts[i] << std::endl;
	}
	// rights
	myfile << "*Nset, nset=rights, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> rights;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_y; j++)
		{
			int number = 1 + (dimension_z + 1) * (j - 1) + (i - 1) * (dimension_z + 1) * (dimension_y + 1);
			myfile << number;
			rights.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < rights.size(); i++)
	{
		myfile << "*Nset, nset = rights" << rights[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << rights[i] << std::endl;
	}
	// leftopedge
	myfile << "*Nset, nset=ltedge, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> ltedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * i;
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
		ltedge.push_back(number);
	}
	myfile << "\n";
	for (int i = 0; i < ltedge.size(); i++)
	{
		myfile << "*Nset, nset = ltedge" << ltedge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << ltedge[i] << std::endl;
	}
	// Leftbottomedge
	myfile << "*Nset, nset=lbedge, instance=BEAM1-1" <<std:: endl;
	linenumber = 0;
	std::vector<int> lbedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * (i - 1) + dimension_z + 1;
		lbedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < lbedge.size(); i++)
	{
		myfile << "*Nset, nset = lbedge" << lbedge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << lbedge[i] << std::endl;
	}
	// righttop
	myfile << "*Nset, nset=rtedge, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> rtedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * (i - 1) + (dimension_z + 1) * dimension_y + 1;
		rtedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < rtedge.size(); i++)
	{
		myfile << "*Nset, nset = rtedge" << rtedge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << rtedge[i] << std::endl;
	}
	// rightbottom
	myfile << "*Nset, nset=rtedge, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> rbedge;
	for (int i = 2; i <= dimension_x; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * (i - 1) + 1;
		rbedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < rbedge.size(); i++)
	{
		myfile << "*Nset, nset = rbedge" << rbedge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << rbedge[i] << std::endl;
	}
	// frontbc
	myfile << "*Nset, nset=frontbc, instance=BEAM1-1,generate" << std::endl;
	myfile << (dimension_y + 1) * (dimension_z + 1) * dimension_x + 1 << "," << (dimension_y + 1) * (dimension_z + 1) * (dimension_x + 1) << "," << 1 << std::endl;

	for (int i = 0; i < (dimension_y + 1) * (dimension_z + 1); i++)
	{
		myfile << "*Nset, nset = frontbc" << (dimension_y + 1) * (dimension_z + 1) * dimension_x + 1 + i << ", instance = BEAM1 - 1" << std::endl;
		myfile << (dimension_y + 1) * (dimension_z + 1) * dimension_x + 1 + i << std::endl;
	}
	// fronts
	myfile << "*Nset, nset=fronts, instance=BEAM1-1" << std::endl;
	linenumber = 0;
	std::vector<int> fronts;
	for (int i = 2; i <= dimension_y; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * dimension_x + j + (i - 1) * (dimension_z + 1);
			myfile << number;
			fronts.push_back(number);
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < fronts.size(); i++)
	{
		myfile << "*Nset, nset = fronts" << fronts[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << fronts[i] << std::endl;
	}
	//fronttop
	std::vector<int> ftedge;
	myfile << "*Nset, nset=ftedge, instance=BEAM1-1" << std::endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_z; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + (dimension_z + 1) * dimension_y + i;
		ftedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < ftedge.size(); i++)
	{
		myfile << "*Nset, nset = ftedge" << ftedge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << ftedge[i] << std::endl;
	}
	//frontbottom
	std::vector<int> fbedge;
	myfile << "*Nset, nset=fbedge, instance=BEAM1-1" <<std:: endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_z; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + i;
		fbedge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < fbedge.size(); i++)
	{
		myfile << "*Nset, nset = fbedge" << fbedge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << fbedge[i] << std::endl;
	}
	//frontleft
	std::vector<int> fledge;
	myfile << "*Nset, nset=fledge, instance=BEAM1-1" <<std:: endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_y; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + i * (dimension_z + 1);
		fledge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < fledge.size(); i++)
	{
		myfile << "*Nset, nset = fledge" << fledge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << fledge[i] << std::endl;
	}
	//frontright
	std::vector<int> fredge;
	myfile << "*Nset, nset=fredge, instance=BEAM1-1" << std::endl;
	linenumber == 0;
	for (int i = 2; i <= dimension_y; i++)
	{
		int number = (dimension_z + 1) * (dimension_y + 1) * dimension_x + (i - 1) * (dimension_z + 1) + 1;
		fredge.push_back(number);
		myfile << number;
		linenumber++;
		if (linenumber != 15)
		{
			myfile << ",";
		}
		else
		{
			myfile << "\n";
			linenumber = 0;
		}
	}
	myfile << "\n";
	for (int i = 0; i < fredge.size(); i++)
	{
		myfile << "*Nset, nset = fredge" << fredge[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << fredge[i] << std::endl;
	}
	// tops
	std::vector<int> tops;
	linenumber = 0;
	myfile << "*Nset, nset=tops, instance=BEAM1-1" << std::endl;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + (dimension_z + 1) * dimension_y + j;
			tops.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < tops.size(); i++)
	{
		myfile << "*Nset, nset = tops" << tops[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << tops[i] << std::endl;
	}
	// Bottomf
	std::vector<int> bots;
	linenumber = 0;
	myfile << "*Nset, nset=bots, instance=BEAM1-1" << std::endl;
	for (int i = 2; i <= dimension_x; i++)
	{
		for (int j = 2; j <= dimension_z; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + j;
			bots.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < bots.size(); i++)
	{
		myfile << "*Nset, nset = bots" << bots[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << bots[i] << std::endl;
	}
	// topbc
	std::vector<int> topbc;
	linenumber = 0;
	myfile << "*Nset, nset=topbc, instance=BEAM1-1" << std::endl;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_z + 1; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + (dimension_z + 1) * dimension_y + j;
			topbc.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < topbc.size(); i++)
	{
		myfile << "*Nset, nset = topbc" << topbc[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << topbc[i] << std::endl;
	}
	// Bottom
	std::vector<int> botbc;
	linenumber = 0;
	myfile << "*Nset, nset=botbc, instance=BEAM1-1" <<std:: endl;
	for (int i = 1; i <= dimension_x + 1; i++)
	{
		for (int j = 1; j <= dimension_z + 1; j++)
		{
			int number = (dimension_y + 1) * (dimension_z + 1) * (i - 1) + j;
			botbc.push_back(number);
			myfile << number;
			linenumber++;
			if (linenumber != 15)
			{
				myfile << ",";
			}
			else
			{
				myfile << "\n";
				linenumber = 0;
			}
		}
	}
	myfile << "\n";
	for (int i = 0; i < botbc.size(); i++)
	{
		myfile << "*Nset, nset = botbc" << botbc[i] << ", instance = BEAM1 - 1" << std::endl;
		myfile << botbc[i] << std::endl;
	}


	myfile << "*Nset, nset = c1" << ", instance = BEAM1 - 1" << std::endl;
	myfile << (dimension_x + 1) * (dimension_z + 1) * (dimension_y + 1) << std::endl;
	myfile << "*Nset, nset = c2" << ", instance = BEAM1 - 1" << std::endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) << std::endl;
	myfile << "*Nset, nset = c3" << ", instance = BEAM1 - 1" << std::endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) - dimension_z << std::endl;
	myfile << "*Nset, nset = c4" << ", instance = BEAM1 - 1" << std::endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) * (dimension_x + 1) - dimension_z << std::endl;
	myfile << "*Nset, nset = c5" << ", instance = BEAM1 - 1" << std::endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) * (dimension_x)+dimension_z + 1 << std::endl;
	myfile << "*Nset, nset = c6" << ", instance = BEAM1 - 1" << std::endl;
	myfile << dimension_z + 1 << std::endl;
	myfile << "*Nset, nset = c7" << ", instance = BEAM1 - 1" << std::endl;
	myfile << 1 << std::endl;
	myfile << "*Nset, nset = c8" << ", instance = BEAM1 - 1" << std::endl;
	myfile << (dimension_z + 1) * (dimension_y + 1) * (dimension_x)+1 << std::endl;


	// Define Equation
	// Tops Bots
	for (int i = 0; i < tops.size(); i++)
	{
		WriteConstraint(tops[i], bots[i], 1, 1.0, -1.0, 0.0, "na", myfile, "tops", "bots");
		WriteConstraint(tops[i], bots[i], 2, 1.0, -1.0, -1.0, "RP5", myfile, "tops", "bots");
		WriteConstraint(tops[i], bots[i], 3, 1.0, -1.0, 0.0, "na", myfile, "tops", "bots");
	}

	// Fronts Backs
	for (int i = 0; i < fronts.size(); i++)
	{
		WriteConstraint(fronts[i], backs[i], 1, 1.0, -1.0, -1.0, "RP4", myfile, "fronts", "backs");
		WriteConstraint(fronts[i], backs[i], 2, 1.0, -1.0, 0.0, "na", myfile, "fronts", "backs");
		WriteConstraint(fronts[i], backs[i], 3, 1.0, -1.0, 0.0, "na", myfile, "fronts", "backs");
	}
	// Left Right
	for (int i = 0; i < lefts.size(); i++)
	{
		WriteConstraint(lefts[i], rights[i], 1, 1.0, -1.0, 0.0, "na", myfile, "lefts", "rights");
		WriteConstraint(lefts[i], rights[i], 2, 1.0, -1.0, 0.0, "na", myfile, "lefts", "rights");
		WriteConstraint(lefts[i], rights[i], 3, 1.0, -1.0, -1.0, "RP6", myfile, "lefts", "rights");

	}
	// ft bt
	for (int i = 0; i < ftedge.size(); i++)
	{
		WriteConstraint(ftedge[i], btedge[i], 1, 1.0, -1.0, -1.0, "RP4", myfile, "ftedge", "btedge");
		WriteConstraint(ftedge[i], btedge[i], 2, 1.0, -1.0, 0.0, "na", myfile, "ftedge", "btedge");
		WriteConstraint(ftedge[i], btedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "ftedge", "btedge");
	}
	// bt bb
	for (int i = 0; i < btedge.size(); i++)
	{
		WriteConstraint(btedge[i], bbedge[i], 1, 1.0, -1.0, 0.0, "na", myfile, "btedge", "bbedge");
		WriteConstraint(btedge[i], bbedge[i], 2, 1.0, -1.0, -1.0, "RP5", myfile, "btedge", "bbedge");
		WriteConstraint(btedge[i], bbedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "btedge", "bbedge");
	}
	// BB fb
	for (int i = 0; i < bbedge.size(); i++)
	{
		WriteConstraint(bbedge[i], fbedge[i], 1, 1.0, -1.0, 1.0, "RP4", myfile, "bbedge", "fbedge");
		WriteConstraint(bbedge[i], fbedge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "bbedge", "fbedge");
		WriteConstraint(bbedge[i], fbedge[i], 3, 1.0, -1.0, 0.0, "NA", myfile, "bbedge", "fbedge");

	}
	// FL BL
	for (int i = 0; i < fledge.size(); i++)
	{
		WriteConstraint(fledge[i], bledge[i], 1, 1.0, -1.0, -1.0, "RP4", myfile, "fledge", "bledge");
		WriteConstraint(fledge[i], bledge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "fledge", "bledge");
		WriteConstraint(fledge[i], bledge[i], 3, 1.0, -1.0, 0.0, "NA", myfile, "fledge", "bledge");

	}
	//BL BR
	for (int i = 0; i < bledge.size(); i++)
	{
		WriteConstraint(bledge[i], bredge[i], 1, 1.0, -1.0, 0.0, "na", myfile, "bledge", "bredge");
		WriteConstraint(bledge[i], bredge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "bledge", "bredge");
		WriteConstraint(bledge[i], bredge[i], 3, 1.0, -1.0, -1.0, "RP6", myfile, "bledge", "bredge");

	}
	//br fr
	for (int i = 0; i < bredge.size(); i++)
	{
		WriteConstraint(bredge[i], fredge[i], 1, 1.0, -1.0, 1.0, "RP4", myfile, "bredge", "fredge");
		WriteConstraint(bredge[i], fredge[i], 2, 1.0, -1.0, 0.0, "NA", myfile, "bredge", "fredge");
		WriteConstraint(bredge[i], fredge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "bredge", "fredge");
	}
	//LT LB
	for (int i = 0; i < ltedge.size(); i++)
	{
		WriteConstraint(ltedge[i], lbedge[i], 1, 1.0, -1.0, 0.0, "na", myfile, "ltedge", "lbedge");
		WriteConstraint(ltedge[i], lbedge[i], 2, 1.0, -1.0, -1.0, "RP5", myfile, "ltedge", "lbedge");
		WriteConstraint(ltedge[i], lbedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "ltedge", "lbedge");
	}
	//LB RB
	for (int i = 0; i < lbedge.size(); i++)
	{
		WriteConstraint(lbedge[i], rbedge[i], 1, 1.0, -1.0, 0.0, "0", myfile, "lbedge", "rbedge");
		WriteConstraint(lbedge[i], rbedge[i], 2, 1.0, -1.0, 0.0, "0", myfile, "lbedge", "rbedge");
		WriteConstraint(lbedge[i], rbedge[i], 3, 1.0, -1.0, -1.0, "RP6", myfile, "lbedge", "rbedge");
	}
	//RB RT
	for (int i = 0; i < rbedge.size(); i++)
	{
		WriteConstraint(rbedge[i], rtedge[i], 1, 1.0, -1.0, 0.0, "0", myfile, "rbedge", "rtedge");
		WriteConstraint(rbedge[i], rtedge[i], 2, 1.0, -1.0, 1.0, "RP5", myfile, "rbedge", "rtedge");
		WriteConstraint(rbedge[i], rtedge[i], 3, 1.0, -1.0, 0.0, "na", myfile, "rbedge", "rtedge");
	}
	//62
	WriteConstraint(6, 2, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(6, 2, 2, 1.0, -1.0, 1.0, "RP5", myfile, "c", "c");
	WriteConstraint(6, 2, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//23
	WriteConstraint(2, 3, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(2, 3, 2, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(2, 3, 3, 1.0, -1.0, -1.0, "RP6", myfile, "c", "c");

	//34
	WriteConstraint(3, 4, 1, 1.0, -1.0, 1.0, "RP4", myfile, "c", "c");
	WriteConstraint(3, 4, 2, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(3, 4, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//48
	WriteConstraint(4, 8, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(4, 8, 2, 1.0, -1.0, -1.0, "RP5", myfile, "c", "c");
	WriteConstraint(4, 8, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//85
	WriteConstraint(8, 5, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(8, 5, 2, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(8, 5, 3, 1.0, -1.0, 1.0, "RP6", myfile, "c", "c");
	//51
	WriteConstraint(5, 1, 1, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	WriteConstraint(5, 1, 2, 1.0, -1.0, 1.0, "RP5", myfile, "c", "c");
	WriteConstraint(5, 1, 3, 1.0, -1.0, 0.0, "0", myfile, "c", "c");
	//17
	WriteConstraint(1, 7, 1, 1.0, -1.0, -1.0, "RP4", myfile, "c", "c");
	WriteConstraint(1, 7, 2, 1.0, -1.0, -1.0, "RP5", myfile, "c", "c");
	WriteConstraint(1, 7, 3, 1.0, -1.0, -1.0, "RP6", myfile, "c", "c");


	//
	myfile << "*End Assembly" << std::endl;

	// write materials
	std::string materialtitle = "**\n**MATERIALS\n**";
	myfile << materialtitle << std::endl;
	for (int indexofm = 0; indexofm < elasticProperties.size(); indexofm++)
	{
		std::string material = "*Material, name = material-" + std::to_string(indexofm + 1) + "\n*Elastic\n" + std::to_string(elasticProperties[indexofm]) + "," + std::to_string(poissionratio[indexofm]);
		myfile << material << std::endl;
	}
	// write step
	myfile << "** ----------------------------------------------------------------" << std::endl;
	myfile << "**\n" << "** STEP:Step-1\n" << "** " << std::endl;
	myfile << "*Step, name = Step - 1, nlgeom=NO" << std::endl;
	myfile << "*Static\n" << "1., 1. , 1e-05,1.0" << std::endl;
	myfile << "**\n" << "** BOUNDARY CONDITION" << "**\n" << "** NAME:E11-1 Type:Displacement/Rotation" << std::endl;
	myfile << "*Boundary\n" << "RP6,3,3,0.02" << std::endl;
	myfile << "**\n" << "** OUTPUT REQUESTS\n" << "**\n" << "*Restart, write, frequency = 0\n" << "**\n" << "** FIELD OUTPUT : F - Output - 1\n" << "**" << std::endl;
	myfile << "*Output, field, variable = PRESELECT" << std::endl;
	myfile << "**\n" << "**HISTORY OUTPUT : H - Output - 2\n" << "**" << std::endl;
	myfile << "*Output, history\n" <<
		"*Node Output, nset = c1\n" << "RT" << std::endl;
	myfile << "**\n" << "**HISTORY OUTPUT : H - Output - 1\n" << "**" << std::endl;
	myfile << "*Output, history, variable = PRESELECT" << std::endl;
	myfile << "*End Step" << std::endl;

	myfile.close();

	return;
}