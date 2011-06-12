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

#include <vtkActor.h>
#include <vtkGlyph3D.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

void GenerateData(vtkPolyData* polydata, int n);

int main(int, char *[])
{
  vtkSmartPointer<vtkPolyData> input =
    vtkSmartPointer<vtkPolyData>::New();
  GenerateData(input, 5);

  std::cout << "There are " << input->GetNumberOfPoints() << " input points." << std::endl;

  vtkSmartPointer<vtkKMeansClustering> kmeans =
    vtkSmartPointer<vtkKMeansClustering>::New();
  kmeans->SetK(3);
  kmeans->SetInputConnection(input->GetProducerPort());
  kmeans->SetInitMethod(vtkKMeansClustering::KMEANSPP);
  kmeans->SetRandom(false); // for repeatable results
  //kmeans->SetRandom(true); // for real, random results
  kmeans->Update();

  std::cout << "There are " << kmeans->GetOutput(0)->GetNumberOfPoints() << " output points." << std::endl;

  // Write output files
  {
  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(kmeans->GetOutputPort(1));
  glyphFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> centersWriter =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  centersWriter->SetInputConnection(glyphFilter->GetOutputPort());
  centersWriter->SetFileName("Centers.vtp");
  centersWriter->Write();
  }

  {
  vtkSmartPointer<vtkXMLPolyDataWriter> labeledPointsWriter =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  labeledPointsWriter->SetInputConnection(kmeans->GetOutputPort(0));
  labeledPointsWriter->SetFileName("LabeledPoints.vtp");
  labeledPointsWriter->Write();
  }
  
  // Visualize result
  vtkSmartPointer<vtkPolyDataMapper> pointsMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  pointsMapper->SetInputConnection(kmeans->GetOutputPort(0));

  vtkSmartPointer<vtkActor> pointsActor =
    vtkSmartPointer<vtkActor>::New();
  pointsActor->SetMapper(pointsMapper);
  pointsActor->GetProperty()->SetPointSize(5);
  
  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetRadius(0.05);
  sphereSource->Update();
  
  vtkSmartPointer<vtkGlyph3D> glyph3D =
    vtkSmartPointer<vtkGlyph3D>::New();
  glyph3D->SetSource(sphereSource->GetOutput());
  glyph3D->SetInput(kmeans->GetOutput(1));
  glyph3D->SetColorModeToColorByScalar();
  glyph3D->ScalingOff();
  glyph3D->Update();

  vtkSmartPointer<vtkPolyDataMapper> centersMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  centersMapper->SetInputConnection(glyph3D->GetOutputPort());

  vtkSmartPointer<vtkActor> centersActor =
    vtkSmartPointer<vtkActor>::New();
  centersActor->SetMapper(centersMapper);

  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(pointsActor);
  renderer->AddActor(centersActor);
  renderer->SetBackground(1,1,1); // Background color white

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

void GenerateData(vtkPolyData* polydata, int n)
{
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  //vtkMath::RandomSeed(time(NULL));
  vtkMath::RandomSeed(0); // for repeatable results
  for(unsigned int i = 0; i < n; i++)
    {
    points->InsertNextPoint(1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1));
    points->InsertNextPoint(-1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1));
    points->InsertNextPoint(1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1));
    }

  polydata->SetPoints(points);
}
