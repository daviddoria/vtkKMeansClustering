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

#include "vtkKMeansClustering.h"

#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkVertexGlyphFilter.h>

#include <set>
#include <numeric>
#include <limits>

vtkStandardNewMacro(vtkKMeansClustering);

vtkKMeansClustering::vtkKMeansClustering()
{
  this->SetNumberOfOutputPorts(2);
  this->K = 3;
  this->InitMethod = KMEANSPP;

  this->Points = vtkSmartPointer<vtkPoints>::New();
  this->ColorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    
  this->Random = true;

  this->KDTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
}

int vtkKMeansClustering::RequestData(vtkInformation *vtkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfoColoredPoints = outputVector->GetInformationObject(0);
  vtkInformation *outInfoClusterCenters = outputVector->GetInformationObject(1);

  // Get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));

  this->Points->ShallowCopy(input->GetPoints());

  vtkPolyData *outputColoredPoints = vtkPolyData::SafeDownCast(
      outInfoColoredPoints->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *outputClusterCenters = vtkPolyData::SafeDownCast(
    outInfoClusterCenters->Get(vtkDataObject::DATA_OBJECT()));

  // Create a KDTree of the points
  this->KDTree->SetDataSet(input);
  this->KDTree->BuildLocator();

  // Seed a random number generator
  if(this->Random)
    {
    vtkMath::RandomSeed(time(NULL));
    }
  else
    {
    vtkMath::RandomSeed(0);
    }

  // Initialize the structure in which to store the cluster centers
  vtkSmartPointer<vtkPoints> clusterCenters =
    vtkSmartPointer<vtkPoints>::New();
  clusterCenters->SetNumberOfPoints(this->K);

  if(this->InitMethod == RANDOM)
    {
    RandomInit(clusterCenters);
    }
  else if(this->InitMethod == KMEANSPP) // http://en.wikipedia.org/wiki/K-means%2B%2B
    {
    MeansPPInit(clusterCenters);
    }
  else
    {
    std::cerr << "An invalid initialization method has been specified!" << std::endl;
    exit(-1);
    }

  /*
  // Output cluster centers
  std::cout << "Initial cluster centers: " << std::endl;
  for(unsigned int i = 0; i < clusterCenters->GetNumberOfPoints(); i++)
    {
    double p[3];
    clusterCenters->GetPoint(i, p);
    std::cout << "Cluster center " << i << " : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  */

  // We must store the labels at the previous iteration to determine whether any labels changed at each iteration.
  std::vector<unsigned int> oldLabels(input->GetNumberOfPoints(), 0); // initialize to all zeros

  // Initialize the labels array
  this->Labels.resize(input->GetNumberOfPoints());

  // The current iteration number
  int iter = 0;

  // Track whether any labels changed in the last iteration
  bool changed = true;
  do
    {
    AssignLabels(input->GetPoints(), clusterCenters);

    EstimateClusterCenters(input->GetPoints(), clusterCenters);

    changed = CheckChanged(this->Labels, oldLabels);

    // Save the old labels
    oldLabels = this->Labels;
    iter++;
    }while(changed);
    //}while(iter < 100); // You could use this stopping criteria to make kmeans run for a specified number of iterations

  std::cout << "KMeans finished in " << iter << " iterations." << std::endl;

  // Create the color map
  this->ColorLookupTable->SetTableRange(0, this->K);
  this->ColorLookupTable->Build();

  CreateCentersPolyData(clusterCenters, outputClusterCenters);
  CreateOutputPointsPolyData(outputColoredPoints);

  return 1;
}

std::vector<unsigned int> vtkKMeansClustering::GetIndicesWithLabel(unsigned int label)
{
  std::vector<unsigned int> pointsWithLabel;
  for(unsigned int i = 0; i < this->Labels.size(); i++)
    {
    if(this->Labels[i] == label)
      {
      pointsWithLabel.push_back(i);
      }
    }

  return pointsWithLabel;
}

void vtkKMeansClustering::GetPointsWithLabel(unsigned int label, vtkPoints* points)
{
  points->Reset();
  points->Squeeze();

  std::vector<unsigned int> indicesWithLabel = GetIndicesWithLabel(label);

  for(unsigned int i = 0; i < indicesWithLabel.size(); i++)
    {
    double p[3];
    this->Points->GetPoint(indicesWithLabel[i], p);
    points->InsertNextPoint(p);
    }
}

unsigned int vtkKMeansClustering::SelectWeightedIndex(std::vector<double> &inputWeights)
{
  // Ensure all weights are positive
  for(unsigned int i = 0; i < inputWeights.size(); i++)
    {
    if(inputWeights[i] < 0)
      {
      std::cout << "inputWeights[" << i << "] is " << inputWeights[i] << " (must be positive!)" << std::endl;
      exit(-1);
      }
    }

  std::vector<double> weights = inputWeights;

  // Normalize
  double sum = std::accumulate(weights.begin(), weights.end(), 0.0f);
  //std::cout << "sum: " << sum << std::endl;
  if(sum <= 0)
    {
    std::cerr << "Sum must be positive!" << std::endl;
    exit(-1);
    }

  for(unsigned int i = 0; i < weights.size(); i++)
    {
    weights[i] /= sum;
    }

  double randomValue = vtkMath::Random(0, 1);

  double runningTotal = 0.0;
  for(unsigned int i = 0; i < weights.size(); i++)
    {
    runningTotal += weights[i];
    if(randomValue < runningTotal)
      {
      return i;
      }
    }

  std::cerr << "vtkKMeansClustering::SelectWeightedIndex reached end, we should never get here" << std::endl;
  std::cerr << "runningTotal: " << runningTotal << std::endl;
  std::cerr << "randomValue: " << randomValue << std::endl;
  exit(-1);
  
  return 0;
}

void vtkKMeansClustering::GetRandomPointInBounds(vtkPoints* data, double p[3])
{
  double bounds[6];
  data->GetBounds(bounds);
  //std::cout << "Bounds: " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << std::endl;

  p[0] = vtkMath::Random(bounds[0], bounds[1]);
  p[1] = vtkMath::Random(bounds[2], bounds[3]);
  p[2] = vtkMath::Random(bounds[4], bounds[5]);
}

bool vtkKMeansClustering::CheckChanged(std::vector<unsigned int> labels, std::vector<unsigned int> oldLabels)
{
  bool changed = false;
  for(unsigned int i = 0; i < labels.size(); i++)
    {
    if(labels[i] != oldLabels[i]) //if something changed
      {
      changed = true;
      break;
      }
    }
  return changed;
}

void vtkKMeansClustering::AssignLabels(vtkPoints* points, vtkPoints* clusterCenters)
{
  // Assign each point to the closest cluster
  for(unsigned int point = 0; point < points->GetNumberOfPoints(); point++)
    {
    unsigned int closestCluster = ClosestPointIndex(clusterCenters, points->GetPoint(point));
    this->Labels[point] = closestCluster;
    }
}

void vtkKMeansClustering::EstimateClusterCenters(vtkPoints* data, vtkPoints* clusterCenters)
{
  vtkSmartPointer<vtkPoints> oldCenters =
    vtkSmartPointer<vtkPoints>::New();
  oldCenters->ShallowCopy(clusterCenters);

  clusterCenters->Reset();
  clusterCenters->Squeeze();

  for(unsigned int cluster = 0; cluster < this->K; cluster++)
    {
    vtkSmartPointer<vtkPoints> classPoints =
      vtkSmartPointer<vtkPoints>::New();
    for(unsigned int point = 0; point < data->GetNumberOfPoints(); point++)
      {
      if(this->Labels[point] == cluster)
        {
        classPoints->InsertNextPoint(data->GetPoint(point));
        }
      }
    double center[3];
    if(classPoints->GetNumberOfPoints() == 0)
      {
      //GetRandomPoint(data, center);
      oldCenters->GetPoint(cluster, center);
      }
    else
      {
      CenterOfMass(classPoints, center);
      }
    clusterCenters->InsertNextPoint(center);
    }

}

unsigned int vtkKMeansClustering::ClosestPointIndex(vtkPoints* points, double queryPoint[3])
{
  // Should use the KDTree here!
  unsigned int closestPoint = 0;
  double minDist = std::numeric_limits<double>::max();
  for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
    //double dist = sqrt(vtkMath::Distance2BetweenPoints(points->GetPoint(i), queryPoint));
    double dist = vtkMath::Distance2BetweenPoints(points->GetPoint(i), queryPoint);
    if(dist < minDist)
      {
      minDist = dist;
      closestPoint = i;
      }
    }

  return closestPoint;
}

double vtkKMeansClustering::ClosestPointDistance(vtkPoints* points, double queryPoint[3])
{
  unsigned int closestPoint = 0;
  double minDist = std::numeric_limits<double>::infinity();
  for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
    double dist = sqrt(vtkMath::Distance2BetweenPoints(points->GetPoint(i), queryPoint));
    //double dist = vtkMath::Distance2BetweenPoints(points->GetPoint(i), queryPoint);
    if(dist < minDist)
      {
      minDist = dist;
      closestPoint = i;
      }
    }

  if(minDist == 0)
    {
    //std::cout << "minDist is 0 for " << queryPoint[0] << " "  <<queryPoint[1] << " " << queryPoint[2] << std::endl;
    }

  return minDist;
}

void vtkKMeansClustering::CenterOfMass(vtkPoints* points, double* center)
{
  if(points->GetNumberOfPoints() == 0)
    {
    std::cout << "Warning: tried to compute center of 0 points!" << std::endl;
    return;
    }
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;

  for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
    double point[3];
    points->GetPoint(i, point);
    vtkMath::Add(center, point, center);
    }

  double numberOfPoints = static_cast<double>(points->GetNumberOfPoints());
  vtkMath::MultiplyScalar(center, 1./static_cast<double>(numberOfPoints));

}

void vtkKMeansClustering::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "K: " << this->K << std::endl;
}

void vtkKMeansClustering::CreateLabeledPointsColorArray(vtkUnsignedCharArray* colors)
{
  // Generate the colors for each point based on the color map
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  for(int i = 0; i < this->Labels.size(); i++)
    {
    double dcolor[3]; // A double array of a color
    this->ColorLookupTable->GetColor(this->Labels[i], dcolor);

    unsigned char color[3];
    for(unsigned int j = 0; j < 3; j++)
      {
      color[j] = 255 * dcolor[j]/1.0;
      }

    colors->InsertNextTupleValue(color);
    }
}

void vtkKMeansClustering::CreateOutputPointsPolyData(vtkPolyData* polyData)
{
  // Create the colors array
  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  CreateLabeledPointsColorArray(colors);

  // We have to add the points to a polydata before we can glyph them
  vtkSmartPointer<vtkPolyData> tempPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  tempPolyData->SetPoints(this->Points);
  
  // Add a vertex to every point so the output is easily visualized
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexGlyphFilter->SetInputData(tempPolyData);
  vertexGlyphFilter->Update();
  polyData->ShallowCopy(vertexGlyphFilter->GetOutput());

  // Create the clusterId array
  vtkSmartPointer<vtkUnsignedIntArray> pointClusterId =
    vtkSmartPointer<vtkUnsignedIntArray>::New();
  pointClusterId->SetName("ClusterId");
  pointClusterId->SetNumberOfComponents(1);
  for(unsigned int i = 0; i < this->Labels.size(); ++i)
    {
    pointClusterId->InsertNextValue(this->Labels[i]);
    }

  polyData->GetPointData()->SetScalars(colors);
  polyData->GetPointData()->AddArray(pointClusterId);
}

void vtkKMeansClustering::CreateCentersPolyData(vtkPoints* clusterCenters, vtkPolyData* polyData)
{
  polyData->SetPoints(clusterCenters);

  vtkSmartPointer<vtkUnsignedCharArray> centerColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  centerColors->SetName("Colors");
  centerColors->SetNumberOfComponents(3);

  vtkSmartPointer<vtkUnsignedIntArray> clusterId =
    vtkSmartPointer<vtkUnsignedIntArray>::New();
  clusterId->SetName("ClusterId");
  clusterId->SetNumberOfComponents(1);
  
  for(int i = 0; i < polyData->GetNumberOfPoints(); i++)
    {
    double dcolor[3];
    this->ColorLookupTable->GetColor(i, dcolor);

    unsigned char color[3];
    for(unsigned int j = 0; j < 3; j++)
      {
      color[j] = 255 * dcolor[j]/1.0;
      }

    centerColors->InsertNextTupleValue(color);
    clusterId->InsertNextValue(i);
    }

  polyData->GetPointData()->AddArray(clusterId);
  polyData->GetPointData()->SetScalars(centerColors);
  
}

void vtkKMeansClustering::RandomInit(vtkPoints* clusterCenters)
{
  // Completely randomly choose initial cluster centers
  for(unsigned int i = 0; i < this->K; i++)
    {
    double p[3];
    GetRandomPointInBounds(this->Points, p);

    //std::cout << "p(" << i << ") = " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    clusterCenters->SetPoint(i, p);
    }
}

void vtkKMeansClustering::MeansPPInit(vtkPoints* clusterCenters)
{
  // Assign one center at random
  double p[3];
  this->Points->GetPoint(rand() % this->Points->GetNumberOfPoints(), p);
  clusterCenters->SetPoint(0, p);
  //std::cout << "Initial center set to " << p[0] << " " << p[1] << " " << p[2] << std::endl;

  // Assign the rest of the initial centers using a weighted probability of the distance to the nearest center
  std::vector<double> weights(this->Points->GetNumberOfPoints());
  for(unsigned int cluster = 1; cluster < this->K; cluster++)
    {
    // Create weight vector
    for(vtkIdType i = 0; i < this->Points->GetNumberOfPoints(); i++)
      {
      double currentPoint[3];
      this->Points->GetPoint(i,currentPoint);
      weights[i] = ClosestPointDistance(clusterCenters, currentPoint);
      }

    unsigned int selectedPoint = SelectWeightedIndex(weights);
    this->Points->GetPoint(selectedPoint, p);
    clusterCenters->SetPoint(cluster, p);
    //std::cout << "Center " << cluster << "set to " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
}
