/*
KMeans clustering
*/

#ifndef __vtkKMeansClustering_h
#define __vtkKMeansClustering_h

#include "Cluster.h"

#include "vtkPolyDataAlgorithm.h"

class vtkKMeansClustering : public vtkPolyDataAlgorithm
{
  public:
  vtkTypeMacro(vtkKMeansClustering,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkKMeansClustering *New();

  std::vector<unsigned int> KMeans(const std::vector<OrientedPoint> &Points, const double k);
  std::vector<vgl_point_3d<double> > EstimateClusterCenters(const std::vector<OrientedPoint> &Points, const std::vector<unsigned int> &Labels);

private:

};

#endif