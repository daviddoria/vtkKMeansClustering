/*
KMeans clustering
*/

#ifndef __vtkKMeansClustering_h
#define __vtkKMeansClustering_h

#include "vtkPolyDataAlgorithm.h"

class vtkKMeansClustering : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkKMeansClustering,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkKMeansClustering *New();

  vtkKMeansClustering();

  vtkSetMacro(K, int);
  vtkGetMacro(K, int);

  void SetInitTypeToRandom(){this->InitType = RANDOM;}
  void SetInitTypeToRandomIndex(){this->InitType = RANDOMINDEX;}

protected:

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  std::vector<unsigned int> UniqueRandomIndices(const unsigned int MAX, const unsigned int Number);
  unsigned int ClosestPointIndex(vtkPoints* points, double queryPoint[3]);
  void EstimateClusterCenters(vtkPoints* data, vtkPoints* points, const std::vector<unsigned int> labels);
  void CenterOfMass(vtkPoints* points, double center[3]);
  void AssignLabels(vtkPoints* points, vtkPoints* clusterCenters, std::vector<unsigned int> &labels);
  bool CheckChanged(std::vector<unsigned int> labels, std::vector<unsigned int> oldLabels);

  // Member data
  int K;

  enum InitTypeEnum {RANDOM, RANDOMINDEX};
  int InitType;

private:
  vtkKMeansClustering(const vtkKMeansClustering&);  // Not implemented.
  void operator =(const vtkKMeansClustering&);  // Not implemented.
};

#endif