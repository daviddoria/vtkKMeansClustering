/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/*
KMeans clustering is a method in which to form K (known) clusters of points from
an unorganized set of input points.
*/

#ifndef __vtkKMeansClustering_h
#define __vtkKMeansClustering_h

#include <vtkKdTreePointLocator.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkKMeansClustering : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkKMeansClustering,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkKMeansClustering *New();

  vtkKMeansClustering();

  // The number of clusters to find
  vtkSetMacro(K, unsigned int);
  vtkGetMacro(K, unsigned int);

  std::vector<unsigned int> GetIndicesWithLabel(unsigned int label);
  void GetPointsWithLabel(unsigned int label, vtkPoints* points);

  /*
   * If this function is called, the randomness
   * is removed for repeatability for testing
   */
  vtkSetMacro(Random,bool);

  // Set which initialization method to use.
  vtkSetMacro(InitMethod,int);

  // Choices of initialization methods
  enum InitMethodEnum{RANDOM, KMEANSPP};

protected:

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // Randomly initialize cluster centers
  void RandomInit(vtkPoints* clusterCenters);

  // Initialize cluster centers using the KMeans++ algorithm
  void MeansPPInit(vtkPoints* clusterCenters);
  
  /* Create an array of colors, one for each cluster */
  void CreateLabeledPointsColorArray(vtkUnsignedCharArray* colors);

  void CreateCentersPolyData(vtkPoints* centers, vtkPolyData* polyData);

  void CreateOutputPointsPolyData(vtkPolyData* polyData);
  
  unsigned int ClosestPointIndex(vtkPoints* points, double queryPoint[3]);
  double ClosestPointDistance(vtkPoints* points, double queryPoint[3]);
  void EstimateClusterCenters(vtkPoints* data, vtkPoints* points);
  void CenterOfMass(vtkPoints* points, double center[3]);
  void AssignLabels(vtkPoints* points, vtkPoints* clusterCenters);
  bool CheckChanged(std::vector<unsigned int> labels, std::vector<unsigned int> oldLabels);
  void GetRandomPointInBounds(vtkPoints* data, double p[3]);
  unsigned int SelectWeightedIndex(std::vector<double> &weights);
  
  // The initialization method to use
  int InitMethod;

  // The number of clusters to find
  unsigned int K;

  // The label of each point
  std::vector<unsigned int> Labels;

  // When RequestData is called, store the input points so that they don't have to be passed to each member function that needs it.
  vtkSmartPointer<vtkPoints> Points;

  // Should use a KDTree to speed up
  vtkSmartPointer<vtkKdTreePointLocator> KDTree;
  
  // When RequestData is called, create a color lookup table so it doesn't have to be passed to each member function that needs it.
  vtkSmartPointer<vtkLookupTable> ColorLookupTable;

  /* Should the computation be random? If false, then it is repeatable (for testing)*/
  bool Random;

private:
  vtkKMeansClustering(const vtkKMeansClustering&);  // Not implemented.
  void operator =(const vtkKMeansClustering&);  // Not implemented.
};

#endif
