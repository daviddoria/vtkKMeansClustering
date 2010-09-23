#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cstdlib>

//#include "Cluster.h"
#include "vtkKMeansClustering.h"

vtkStandardNewMacro(vtkKMeansClustering);

std::vector<unsigned int> KMeans(const std::vector<OrientedPoint> &Points, const double k)
{

	std::vector<vgl_point_3d<double> > Centroids(k);

	std::vector<unsigned int> RandomIndices = Tools::UniqueRandomIndices(Points.size(), k);
	for(unsigned int i = 0; i < k; i++)
	{
		Centroids[i] = Points[RandomIndices[i]].getCoord();
	}
	/*
	//pick k random cluster centers between min and max
	double MIN = -5.0;
	double MAX = 5.0;
	for(unsigned int i = 0; i < k; i++)
	{
		Centroids[i] = vgl_point_3d<double> (Tools::RandomDouble(MIN, MAX),
									  Tools::RandomDouble(MIN, MAX),
									Tools::RandomDouble(MIN, MAX));
	}
	*/

	std::cout << "Initial cluster centers:" << std::endl;
	Tools::OutputVector(Centroids);

	std::vector<unsigned int> Labels(Points.size());
	std::vector<unsigned int> OldLabels;

	bool done = false;
	while(!done)
	{
		//save the old labels
		OldLabels = Labels;

		//assign each point to the closest cluster
		for(unsigned int point = 0; point < Points.size(); point++)
		{
			unsigned int ClosestCentroid = geom::ClosestPointIndex(Points[point].getCoord(), Centroids);
			Labels[point] = ClosestCentroid;
		}

		//count the number of points which did not change clusters
		unsigned int unchanged = 0;
		for(unsigned int i = 0; i < Points.size(); i++)
		{
			if(Labels[i] == OldLabels[i]) //if something changed
				unchanged++;
		}

		//quit if nothing changed
		if(unchanged == Labels.size())
			done = true;

		//find the centers of the new clusters
		if(!done)
			Centroids = EstimateClusterCenters(Points, Labels);
	}

	return Labels;
}

std::vector<vgl_point_3d<double> > EstimateClusterCenters(const std::vector<OrientedPoint> &Points, const std::vector<unsigned int> &Labels)
{
	unsigned int NumLabels = Tools::VectorMax(Labels) + 1;

	std::vector<vgl_point_3d<double> > Centers(NumLabels);
	for(unsigned int cluster = 0; cluster < NumLabels; cluster++)
	{
		std::vector<vgl_point_3d<double> > ClassPoints;
		for(unsigned int point = 0; point < Points.size(); point++)
		{
			if(Labels[point] == cluster)
				ClassPoints.push_back(Points[point].getCoord());
		}
		Centers[cluster] = geom::CenterOfMass(ClassPoints);
	}

	return Centers;
}
