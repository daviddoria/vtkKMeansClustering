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

#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>

// This function generates the very distinct clusters with 'numberOfPointsPerCluster' in each cluster.
// The information about which cluster each point belongs to is not stored, so that the kmeans algorithm
// can not cheat!
void GenerateData(vtkPolyData* polydata, int numberOfPointsPerCluster);

int main(int, char *[])
{
  vtkSmartPointer<vtkPolyData> input =
    vtkSmartPointer<vtkPolyData>::New();
  GenerateData(input, 5);

  vtkSmartPointer<vtkKMeansClustering> kmeans =
    vtkSmartPointer<vtkKMeansClustering>::New();
  kmeans->SetK(3);
  kmeans->SetInputConnection(input->GetProducerPort());
  kmeans->SetRandom(false);
  kmeans->Update();

  if(kmeans->GetOutput(0)->GetNumberOfPoints() != 15) // There should be 15 labeled points
    {
    std::cout << "The number of labeled points is: " << kmeans->GetOutput(0)->GetNumberOfPoints()
              << "but it should be " << 15 << "!" << std::endl;
    return EXIT_FAILURE;
    }
  
  if(kmeans->GetOutput(1)->GetNumberOfPoints() != 3) // There should be 3 clusters
    {
    std::cout << "The number of clusters is: " << kmeans->GetOutput(1)->GetNumberOfPoints()
              << "but it should be " << 3 << "!" << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}

void GenerateData(vtkPolyData* polydata, int numberOfPointsPerCluster)
{
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  vtkMath::RandomSeed(0);
  for(unsigned int i = 0; i < numberOfPointsPerCluster; i++)
    {
    points->InsertNextPoint(1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1));
    points->InsertNextPoint(-1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1));
    points->InsertNextPoint(1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1));
    }

  polydata->SetPoints(points);
}
