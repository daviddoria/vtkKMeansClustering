#include "Cluster.h"

#include <Geometry/Helpers.h>
#include <Geometry/Geometry.h>

std::ostream& operator<<(std::ostream& output, const Cluster &C)
{
	output << "Cluster" << std::endl << "-----------" << std::endl
			<< C.NumPoints() << std::endl
			<< C.getCenter() << std::endl << std::endl;

	return output;
}

void Cluster::setPoints(const std::vector<OrientedPoint> &Points)
{
	Points_ = Points;
	Center_ = geom::CenterOfMass(GetOPCoords(Points_));
}


void Cluster::AppendCluster(const Cluster &C)
{
	for(unsigned int i = 0; i < C.NumPoints(); i++)
	{
		this->Points_.push_back(C.getPoint(i));
	}

	Center_ = geom::CenterOfMass(GetOPCoords(Points_));
}