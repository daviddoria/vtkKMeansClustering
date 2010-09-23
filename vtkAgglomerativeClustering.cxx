#include <iostream>
#include <vector>
//#include <list>
#include <set>
#include <map>
#include <cstdlib>

#include <assert.h>

#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_distance.h>

//#include <vnl/vnl_matrix.h>
#include <vbl/vbl_array_2d.h>

#include <Geometry/Geometry.h>

#include <Tools/Tools.h>
#include <ValidType/ValidType.h>

#include <VXLHelpers/VXLHelpers.h>

//#include "Cluster.h"
#include "vtkAgglomerativeClustering.h"

vtkStandardNewMacro(vtkAgglomerativeClustering);

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


///////////////////////////////////////////////////////////////////

/////////// Constructor ////////////
Agglomerative::Agglomerative(const std::vector<OrientedPoint> &points, const double maxdist) : Points(points), MaxDist(maxdist)
{
		//initialize the labels to their index
	PointLabels.resize(Points.size());
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		PointLabels[i] = i;
	}

		//initialize the cluster centers to the input points
	ClusterCenters.resize(Points.size());
	for(unsigned int i = 0; i < Points.size(); i++)
	{
		ClusterCenters[i] = ValidType<vgl_point_3d<double> >(Points[i].getCoord());
	}

		//create the distance matrix
	CreateDistanceMatrix();

}

void Agglomerative::PerformClustering()
{

	//while(CombineClosestClusters())
	std::cout << std::endl << "Clustering..." << std::endl;
	boost::progress_display show_progress(ClusterDistanceMatrix.rows());
	while(CombineSingleLinkage())
	{
		std::vector<unsigned int> Unique = Tools::UniqueElements(PointLabels);
		//cout << "There are currently " << Unique.size() << " clusters." << endl;
		++show_progress;
	}

	std::vector<unsigned int> Unique = Tools::UniqueElements(PointLabels);

	//throw away clusters that are too small (group them all together)
	unsigned int UnusedIndex = Tools::FirstUnusedNumber(Unique);
	unsigned int MinPoints = 100;

	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		std::vector<unsigned int> ClusterPoints = GetClusterPointIndices(i);
		if(ClusterPoints.size() < MinPoints)
		{
			for(unsigned int u = 0; u < ClusterPoints.size(); u++)
			{
				PointLabels[ClusterPoints[u]] = UnusedIndex;
			}
		}
	}

}

bool Agglomerative::CombineClosestClusters()
{
	//cluster distance function - distance between cluster centers

	//return true if two clusters were merged
	//return false if there was an error or the distance threshold was reached

	//cout << "There are " << Clusters.size() << " clusters." << endl;

	//cout << endl << "CombineClosestClusters" << endl << "--------------" << endl;

	//OutputDistanceMatrix();

	ValidType<double> MinVal;
	std::pair<unsigned int, unsigned int> MinLoc;
	MinLoc = MinValueLocation(ClusterDistanceMatrix, MinVal);

	if(!MinVal.Valid)
		return false;

	if(MinVal.Value > MaxDist)
		return false;

	//update the point labels
	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		if(PointLabels[i] == MinLoc.second)
			PointLabels[i] = MinLoc.first;
	}

	//update the cluster center of the new cluster
	std::vector<vgl_point_3d<double> > ClusterPoints;
	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		if(PointLabels[i] == MinLoc.first)
		{
			ClusterPoints.push_back(Points[i].getCoord());
		}
	}
	ClusterCenters[MinLoc.first] = geom::CenterOfMass(ClusterPoints);

	//invalidate the cluster center of the cluster that is being removed
	ClusterCenters[MinLoc.second].Valid = false;

	//cout << "Closest Distance: " << MinVal.Value << endl;

	//invalidate the column of the distance matrix corresponding to the cluster that is being removed
	for(unsigned int r = 0; r < Points.size(); r++)
	{
		ClusterDistanceMatrix(r, MinLoc.second).Value = 0.0;
		ClusterDistanceMatrix(r, MinLoc.second).Valid = false;
	}

	//invalidate the row of the distance matrix corresponding to the cluster that is being removed
	for(unsigned int c = 0; c < Points.size(); c++)
	{
		ClusterDistanceMatrix(MinLoc.second, c).Value = 0.0;
		ClusterDistanceMatrix(MinLoc.second, c).Valid = false;
	}

	//update the distances for the new cluster column
	for(unsigned int r = 0; r < Points.size(); r++)
	{
		if(ClusterDistanceMatrix(r,MinLoc.first).Valid == false)
			continue;

		if(r == MinLoc.first)
		{
			ClusterDistanceMatrix(r,MinLoc.first).Valid = false;
			continue;
		}

		ValidType<double> dist = vgl_distance(ClusterCenters[r].Value, ClusterCenters[MinLoc.first].Value);
		ClusterDistanceMatrix(r,MinLoc.first) = dist;
	}

	//update the distances for the new cluster row
	for(unsigned int c = 0; c < Points.size(); c++)
	{
		if(ClusterDistanceMatrix(MinLoc.first, c).Valid == false)
			continue;

		if(c == MinLoc.first)
		{
			ClusterDistanceMatrix(MinLoc.first, c).Valid = false;
			continue;
		}

		ValidType<double> dist = vgl_distance(ClusterCenters[MinLoc.first].Value, ClusterCenters[c].Value);
		ClusterDistanceMatrix(MinLoc.first, c) = dist;
	}
	return true;

}

std::vector<vgl_point_3d<double> > Agglomerative::GetClusterPoints(const unsigned int ClusterLabel) const
{
	vector<vgl_point_3d<double> > ClusterPoints;
	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		if(PointLabels[i] == ClusterLabel)
			ClusterPoints.push_back(Points[i].getCoord());
	}
	return ClusterPoints;
}

std::vector<unsigned int> Agglomerative::GetClusterPointIndices(const unsigned int ClusterLabel) const
{
	vector<unsigned int> Indices;
	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		if(PointLabels[i] == ClusterLabel)
			Indices.push_back(i);
	}
	return Indices;
}

bool Agglomerative::CombineSingleLinkage()
{
	//cluster distance function - "Single linkage clustering" - distance between clusters is the distance between the two closest elements

	//return true if two clusters were merged
	//return false if there was an error or the distance threshold was reached

	ValidType<double> MinVal;
	std::pair<unsigned int, unsigned int> MinLoc;
	MinLoc = MinValueLocation(ClusterDistanceMatrix, MinVal);

	if(!MinVal.Valid)
		return false;

	if(MinVal.Value > MaxDist)
		return false;

	//update the point labels
	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		if(PointLabels[i] == MinLoc.second)
			PointLabels[i] = MinLoc.first;
	}

	//invalidate the column of the distance matrix corresponding to the cluster that is being removed
	for(unsigned int r = 0; r < Points.size(); r++)
	{
		ClusterDistanceMatrix(r, MinLoc.second).Value = 0.0;
		ClusterDistanceMatrix(r, MinLoc.second).Valid = false;
	}

	//invalidate the row of the distance matrix corresponding to the cluster that is being removed
	for(unsigned int c = 0; c < Points.size(); c++)
	{
		ClusterDistanceMatrix(MinLoc.second, c).Value = 0.0;
		ClusterDistanceMatrix(MinLoc.second, c).Valid = false;
	}

	vector<unsigned int> ThisClusterPointIndices = GetClusterPointIndices(MinLoc.first);

	//update the distances for the new cluster column
//#pragma omp parallel for //may have to change something
	for(unsigned int r = 0; r < Points.size(); r++)
	{
		//if the distance is already invalid, skip it
		if(ClusterDistanceMatrix(r,MinLoc.first).Valid == false)
			continue;

		if(r == MinLoc.first) //don't compare a cluster with itself
			continue;

		std::vector<unsigned int> CurrentClusterPointIndices = GetClusterPointIndices(r);
		std::vector<std::pair<unsigned int, unsigned int> > Pairs = CreatePairs(ThisClusterPointIndices, CurrentClusterPointIndices);
		std::pair<unsigned int, unsigned int> MinPair = GetMinPair(Pairs);

		double dist = PointDistanceMatrix(MinPair.first, MinPair.second);

		ClusterDistanceMatrix(r,MinLoc.first) = ValidType<double>(dist);
	}

	//update the distances for the new cluster row
	for(unsigned int c = 0; c < Points.size(); c++)
	{
		//if the distance is already invalid, skip it
		if(ClusterDistanceMatrix(MinLoc.first, c).Valid == false)
			continue;

		if(c == MinLoc.first) //don't compare a cluster with itself
			continue;

		std::vector<unsigned int> CurrentClusterPointIndices = GetClusterPointIndices(c);
		std::vector<std::pair<unsigned int, unsigned int> > Pairs = CreatePairs(ThisClusterPointIndices, CurrentClusterPointIndices);
		std::pair<unsigned int, unsigned int> MinPair = GetMinPair(Pairs);

		double dist = PointDistanceMatrix(MinPair.first, MinPair.second);

		ClusterDistanceMatrix(MinLoc.first, c) = ValidType<double>(dist);
	}
	return true;

}

void Agglomerative::CreateDistanceMatrix()
{
	ValidType<double> InvalidValue(0.0);
	InvalidValue.Valid = false;
	ClusterDistanceMatrix = vbl_array_2d<ValidType<double> > (ClusterCenters.size(), ClusterCenters.size(), InvalidValue);
	PointDistanceMatrix = vbl_array_2d<double> (ClusterCenters.size(), ClusterCenters.size());

	//ClusterDistanceMatrix = vbl_array_2d<ValidType<float> > (ClusterCenters.size(), ClusterCenters.size(), InvalidValue);
	//PointDistanceMatrix = vbl_array_2d<float> (ClusterCenters.size(), ClusterCenters.size());

	//ClusterDistanceMatrix = vbl_array_2d<ValidType<double> > (1000,1000, InvalidValue); //to test if not enough memory

	std::cout << std::endl << "Calculating Distance Matrix..." << std::endl;
	boost::progress_display show_progress(ClusterDistanceMatrix.rows() * ClusterDistanceMatrix.cols() / 2);

	//create the upper triangular distance matrix (it is symmetric)
	for(unsigned int r = 0; r < ClusterDistanceMatrix.rows(); r++)
	{
		for(unsigned int c = r; c < ClusterDistanceMatrix.cols(); c++)
		{
			if(r == c)
			{
				PointDistanceMatrix(r,c) = 0.0;
				continue;
			}

			double dist = vgl_distance(ClusterCenters[r].Value, ClusterCenters[c].Value);
			PointDistanceMatrix(r,c) = dist;
			PointDistanceMatrix(c,r) = dist;

			ValidType<double> ValidDist(dist);
			ClusterDistanceMatrix(r,c) = ValidDist;

			++show_progress;
		}
	}

}

void Agglomerative::OutputDistanceMatrix() const
{
	for(unsigned int r = 0; r < ClusterDistanceMatrix.rows(); r++)
	{
		for(unsigned int c = 0; c < ClusterDistanceMatrix.cols(); c++)
		{
			std::cout << ClusterDistanceMatrix(r,c) << " ";
		}
		std::cout << std::endl;
	}

}

pair<unsigned int, unsigned int> Agglomerative::GetMinPair(const vector<pair<unsigned int, unsigned int> > &Pairs) const
{
	unsigned int LowestPair = 0;
	double LowestValue = PointDistanceMatrix(Pairs[0].first, Pairs[0].second);
	for(unsigned int i = 0; i < Pairs.size(); i++)
	{
		if(PointDistanceMatrix(Pairs[i].first, Pairs[i].second) < LowestValue)
		{
			LowestPair = i;
			LowestValue = PointDistanceMatrix(Pairs[i].first, Pairs[i].second);
		}
	}
	return Pairs[LowestPair];
}

vector<pair<unsigned int, unsigned int> > Agglomerative::CreatePairs(const vector<unsigned int> &V1, const vector<unsigned int> &V2) const
{
	vector<pair<unsigned int, unsigned int> > Pairs;

	for(unsigned int i = 0; i < V1.size(); i++)
	{
		for(unsigned int j = 0; j < V2.size(); j++)
		{
			Pairs.push_back(pair<unsigned int, unsigned int>(V1[i], V2[j]));
		}
	}

	return Pairs;
}