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
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>

#include <set>

vtkStandardNewMacro(vtkKMeansClustering);

vtkKMeansClustering::vtkKMeansClustering()
{
  this->InitType = RANDOMINDEX;
  this->SetNumberOfOutputPorts(2);
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

  vtkPolyData *outputColoredPoints = vtkPolyData::SafeDownCast(
      outInfoColoredPoints->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *outputClusterCenters = vtkPolyData::SafeDownCast(
    outInfoClusterCenters->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkPoints> clusterCenters =
    vtkSmartPointer<vtkPoints>::New();
  clusterCenters->SetNumberOfPoints(this->K);

  if(this->InitType == RANDOMINDEX)
    {
    // Initialize the cluster centers to be random existing points
    std::vector<unsigned int> randomIndices = UniqueRandomIndices(input->GetNumberOfPoints(), this->K);
    for(unsigned int i = 0; i < this->K; i++)
      {
      clusterCenters->SetPoint(i, input->GetPoint(randomIndices[i]));
      }
    }
  else if(this->InitType == RANDOM)
    {
    // Completely randomly choose initial cluster centers
    double bounds[6];
    input->GetBounds(bounds);
    for(unsigned int i = 0; i < this->K; i++)
      {
      double p[3];
      p[0] = vtkMath::Random(bounds[0], bounds[1]);
      p[1] = vtkMath::Random(bounds[2], bounds[3]);
      p[2] = vtkMath::Random(bounds[4], bounds[5]);
      clusterCenters->SetPoint(i, p);
      }
    }


  std::vector<unsigned int> oldLabels(input->GetNumberOfPoints(), 0); // initialize to all zeros
  std::vector<unsigned int> labels(input->GetNumberOfPoints());

  bool changed = true;
  do
    {
    AssignLabels(input->GetPoints(), clusterCenters, labels);

    EstimateClusterCenters(input->GetPoints(), clusterCenters, labels);

    changed = CheckChanged(labels, oldLabels);

    // Save the old labels
    oldLabels = labels;

    }while(changed);

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
    colorLookupTable->GetColor(labels[i], dcolor);

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

  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetRadius(0.1);
  sphereSource->Update();

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

  vtkSmartPointer<vtkGlyph3D> glyph3D =
    vtkSmartPointer<vtkGlyph3D>::New();
  glyph3D->SetSource(sphereSource->GetOutput());
  glyph3D->SetInput(clusterCentersPolyData);
  glyph3D->SetColorModeToColorByScalar();
  glyph3D->ScalingOff();
  glyph3D->Update();

  outputClusterCenters->ShallowCopy(glyph3D->GetOutput());
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

void vtkKMeansClustering::AssignLabels(vtkPoints* points, vtkPoints* clusterCenters, std::vector<unsigned int> &labels)
{
  // Assign each point to the closest cluster

  for(unsigned int point = 0; point < points->GetNumberOfPoints(); point++)
    {
    unsigned int closestCluster = ClosestPointIndex(clusterCenters, points->GetPoint(point));
    labels[point] = closestCluster;
    }

}

void vtkKMeansClustering::EstimateClusterCenters(vtkPoints* data, vtkPoints* clusterCenters, const std::vector<unsigned int> labels)
{
  clusterCenters->Reset();
  clusterCenters->Squeeze();

  for(unsigned int cluster = 0; cluster < this->K; cluster++)
    {
    vtkSmartPointer<vtkPoints> classPoints =
      vtkSmartPointer<vtkPoints>::New();
    for(unsigned int point = 0; point < data->GetNumberOfPoints(); point++)
      {
      if(labels[point] == cluster)
        {
        classPoints->InsertNextPoint(data->GetPoint(point));
        }
      }
    double center[3];
    CenterOfMass(classPoints, center);
    clusterCenters->InsertNextPoint(center);
    }

}

std::vector<unsigned int> vtkKMeansClustering::UniqueRandomIndices(const unsigned int MAX, const unsigned int Number)
{
  // Generate Number unique random indices from 0 to MAX
  vtkMath::RandomSeed(time(NULL));

  std::set<unsigned int> S;
  while(S.size() < Number)
    {
    S.insert(vtkMath::Random(0, MAX));
    }

  std::vector<unsigned int> indices;
  for(std::set<unsigned int>::iterator iter = S.begin(); iter != S.end(); iter++)
    {
    indices.push_back(*iter);
    }

  return indices;
}

unsigned int vtkKMeansClustering::ClosestPointIndex(vtkPoints* points, double queryPoint[3])
{
  unsigned int closestPoint = 0;
  double minDist = 100000;
  for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
    double dist = vtkMath::Distance2BetweenPoints(points->GetPoint(i), queryPoint);
    if(dist < minDist)
      {
      minDist = dist;
      closestPoint = i;
      }
    }

  return closestPoint;
}

void vtkKMeansClustering::CenterOfMass(vtkPoints* points, double* center)
{
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