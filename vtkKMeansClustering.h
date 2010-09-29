/*
KMeans clustering
*/

#ifndef __vtkKMeansClustering_h
#define __vtkKMeansClustering_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkKMeansClustering : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkKMeansClustering,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkKMeansClustering *New();

  vtkKMeansClustering();

  vtkSetMacro(K, unsigned int);
  vtkGetMacro(K, unsigned int);

  std::vector<unsigned int> GetIndicesWithLabel(unsigned int label);
  void GetPointsWithLabel(unsigned int label, vtkPoints* points);
  
protected:

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  unsigned int ClosestPointIndex(vtkPoints* points, double queryPoint[3]);
  double ClosestPointDistance(vtkPoints* points, double queryPoint[3]);
  void EstimateClusterCenters(vtkPoints* data, vtkPoints* points);
  void CenterOfMass(vtkPoints* points, double center[3]);
  void AssignLabels(vtkPoints* points, vtkPoints* clusterCenters);
  bool CheckChanged(std::vector<unsigned int> labels, std::vector<unsigned int> oldLabels);
  void GetRandomPointInBounds(vtkPoints* data, double p[3]);
  unsigned int SelectWeightedIndex(std::vector<double> &weights);
  
  enum InitMethodEnum{RANDOM, KMEANSPP};
  
  // Member data
  int InitMethod;
  unsigned int K;
  
  std::vector<unsigned int> Labels;

  // These are stored in the class so that the GetPointsWithLabel function will work
  vtkSmartPointer<vtkPoints> Points;
  
private:
  vtkKMeansClustering(const vtkKMeansClustering&);  // Not implemented.
  void operator =(const vtkKMeansClustering&);  // Not implemented.
};

#endif