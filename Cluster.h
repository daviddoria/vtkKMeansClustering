#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <Geometry/OrientedPoint.h>

#include <vgl/vgl_point_3d.h>

class Cluster
{
	std::vector<OrientedPoint> Points_;
	vgl_point_3d<double> Center_;

	public:

	Cluster(){}
	Cluster(const std::vector<OrientedPoint> &Points)
	{
		setPoints(Points);
	}
	
	vgl_point_3d<double> getCenter() const {return Center_;}

	OrientedPoint getPoint(unsigned int i) const {return Points_[i];}

	unsigned int NumPoints() const {return Points_.size();}

	void setPoints(const std::vector<OrientedPoint> &Points);

	void AppendCluster(const Cluster &C);
};

std::ostream& operator<<(std::ostream& output, const Cluster &C);

#endif