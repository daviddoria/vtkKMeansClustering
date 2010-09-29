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
  kmeans->Update();

  std::cout << "There are " << kmeans->GetOutput(0)->GetNumberOfPoints() << " output points." << std::endl;

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

  vtkMath::RandomSeed(time(NULL));
  for(unsigned int i = 0; i < n; i++)
    {
    points->InsertNextPoint(1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1));
    points->InsertNextPoint(-1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1));
    points->InsertNextPoint(1.0 + vtkMath::Random(-.1, .1), -1.0 + vtkMath::Random(-.1, .1), 1.0 + vtkMath::Random(-.1, .1));
    }

  polydata->SetPoints(points);
}
