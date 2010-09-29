#include "vtkKMeansClustering.h"

#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>

#include <vtkLookupTable.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
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

  vtkSmartPointer<vtkPoints> clusterCenters =
    vtkSmartPointer<vtkPoints>::New();
  clusterCenters->SetNumberOfPoints(this->K);

  if(this->InitMethod == RANDOM)
    {
    // Completely randomly choose initial cluster centers
    vtkMath::RandomSeed(time(NULL));
    for(unsigned int i = 0; i < this->K; i++)
      {
      double p[3];
      GetRandomPointInBounds(input->GetPoints(), p);

      std::cout << "p(" << i << ") = " << p[0] << " " << p[1] << " " << p[2] << std::endl;
      clusterCenters->SetPoint(i, p);
      }
    }
  else if(this->InitMethod == KMEANSPP) // http://en.wikipedia.org/wiki/K-means%2B%2B
    {
    // Assign one center at random
    double p[3];
    input->GetPoint(rand() % input->GetNumberOfPoints(), p);
    clusterCenters->SetPoint(0, p);
    // Assign the rest of the initial centers using a weighted probability of the distance to the nearest center
    std::vector<double> weights(input->GetNumberOfPoints());
    for(unsigned int cluster = 1; cluster < this->K; cluster++)
      {
      // Create weight vector
      for(vtkIdType i = 0; i < input->GetNumberOfPoints(); i++)
        {
        double currentPoint[3];
        input->GetPoint(i,currentPoint);
        weights[i] = ClosestPointDistance(clusterCenters, currentPoint);
        }

      unsigned int selectedPoint = SelectWeightedIndex(weights);
      input->GetPoint(selectedPoint, p);
      clusterCenters->SetPoint(cluster, p);
      }
    }
  else
    {
    std::cerr << "An invalid initialization method has been specified!" << std::endl;
    exit(-1);
    }
  /*
  std::cout << "Initial cluster centers: " << std::endl;
  for(unsigned int i = 0; i < clusterCenters->GetNumberOfPoints(); i++)
    {
    double p[3];
    clusterCenters->GetPoint(i, p);
    std::cout << "Cluster center " << i << " : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  */
  std::vector<unsigned int> oldLabels(input->GetNumberOfPoints(), 0); // initialize to all zeros
  this->Labels.resize(input->GetNumberOfPoints());

  int iter = 0;
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
    //}while(iter < 100);

  std::cout << "KMeans finished in " << iter << " iterations." << std::endl;

  // Create the color map
  vtkSmartPointer<vtkLookupTable> colorLookupTable =
    vtkSmartPointer<vtkLookupTable>::New();
  colorLookupTable->SetTableRange(0, this->K);
  colorLookupTable->Build();

  // Generate the colors for each point based on the color map
  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  for(int i = 0; i < input->GetNumberOfPoints(); i++)
    {
    double dcolor[3];
    colorLookupTable->GetColor(this->Labels[i], dcolor);

    unsigned char color[3];
    for(unsigned int j = 0; j < 3; j++)
      {
      color[j] = 255 * dcolor[j]/1.0;
      }

    colors->InsertNextTupleValue(color);
    }

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(input->GetProducerPort());
  glyphFilter->Update();

  outputColoredPoints->ShallowCopy(glyphFilter->GetOutput());
  outputColoredPoints->GetPointData()->SetScalars(colors);

  // Setup second output (cluster centers)
  vtkSmartPointer<vtkPolyData> clusterCentersPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  clusterCentersPolyData->SetPoints(clusterCenters);

  vtkSmartPointer<vtkUnsignedCharArray> centerColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  centerColors->SetName("Colors");
  centerColors->SetNumberOfComponents(3);

  for(int i = 0; i < clusterCentersPolyData->GetNumberOfPoints(); i++)
    {
    double dcolor[3];
    colorLookupTable->GetColor(i, dcolor);

    unsigned char color[3];
    for(unsigned int j = 0; j < 3; j++)
      {
      color[j] = 255 * dcolor[j]/1.0;
      }

    centerColors->InsertNextTupleValue(color);
    }

  clusterCentersPolyData->GetPointData()->SetScalars(centerColors);

  outputClusterCenters->ShallowCopy(clusterCentersPolyData);

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
  double sum = std::accumulate(weights.begin(), weights.end(), 0);
  std::cout << "sum: " << sum << std::endl;

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
  unsigned int closestPoint = 0;
  double minDist = 100000;
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