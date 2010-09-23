#include <iostream>

#include <Tools/Tools.h>

#include <ModelFile/ModelFile.h>

#include <Geometry/Color.h>
#include <Geometry/OrientedPoint.h>

#include <vgl/vgl_point_3d.h>

#include <vul/vul_timer.h>

#include "Clustering.h"

#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_vector_io.h>
#include <vnl/io/vnl_io_matrix.h>
#include <vnl/io/vnl_io_vector.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

void TestKMeans();
void TestAgglomerativeClusterFakeData();

void TestAgglomerativeClusterData(ModelFile &Model);

vector<OrientedPoint> GenerateData();

void WriteData(Agglomerative &A, ModelFile &Model);

int main(int argc, char *argv[])
{

	//Tools::AssertNumArgs(argc, 2);
	std::string InputFile = argv[1];
	//string OutputFile = argv[2];

	//create the model
	//ModelFile Model(InputFile);

	TestKMeans();
	
	//TestAgglomerativeClusterFakeData();

	/*
	vul_timer timer;
	TestAgglomerativeClusterData(Model);
	cout << "Time: " << timer.real() << " ms" << endl;
	cout << "Time: " << (timer.real()/1000.)/60./60. << " hours" << endl;
	*/
	return 0;
}

void TestKMeans()
{
	std::vector<OrientedPoint> Points = GenerateData();

	std::vector<unsigned int> Clusters = KMeans(Points, 3);

	Tools::OutputVector(Clusters);
	
}

#include <vsl/vsl_vector_io.txx>
#include <vbl/io/vbl_io_array_2d.txx>
VSL_VECTOR_IO_INSTANTIATE(vbl_array_2d<double>);
VSL_VECTOR_IO_INSTANTIATE(vbl_array_2d<bool>);

//VSL_VECTOR_IO_INSTANTIATE(vbl_array_2d<ValidType<double> >);
	    
void TestAgglomerativeClusterFakeData()
{
	std::vector<OrientedPoint> Points = GenerateData();
	//Tools::OutputVector(Points);
	
	Agglomerative A(Points, 0.2);
	//A.OutputDistanceMatrix();
	
	A.PerformClustering();
	
	std::vector<unsigned int> PointLabels = A.getPointLabels();
	
	//WriteData(A);
	
	
	//cout << "Final Clusters" << endl << "------------" << endl;
	//cout << "There are " << Clusters.size() << " clusters." << endl;
	
	Tools::OutputVector(PointLabels);
	
	
}

void TestAgglomerativeClusterData(ModelFile &Model)
{
	std::vector<OrientedPoint> ModelPoints = Model.getPoints();
	Agglomerative A(ModelPoints, 0.7);
	A.PerformClustering();
	
	WriteData(A, Model);
	
	cout << "There are : " << Tools::UniqueElements(A.getPointLabels()).size() << " clusters." << endl;

}

vector<OrientedPoint> GenerateData()
{

	OrientedPoint P0(vgl_point_3d<double> (1.0, 1.0, 1.0));
	OrientedPoint P1(vgl_point_3d<double> (1.01, 1.01, 1.01));
	OrientedPoint P2(vgl_point_3d<double> (1.02, 1.02, 1.02));

	OrientedPoint P3(vgl_point_3d<double> (-1.0, -1.0, -1.0));
	OrientedPoint P4(vgl_point_3d<double> (-1.01, -1.01, -1.01));
	OrientedPoint P5(vgl_point_3d<double> (-1.02, -1.02, -1.02));

	OrientedPoint P6(vgl_point_3d<double> (1.0, -1.0, 1.0));
	OrientedPoint P7(vgl_point_3d<double> (1.01, -1.01, 1.02));
	OrientedPoint P8(vgl_point_3d<double> (1.02, -1.02, 1.02));

	std::vector<OrientedPoint> Points;
	Points.push_back(P0);
	Points.push_back(P1);
	Points.push_back(P2);
	Points.push_back(P3);
	Points.push_back(P4);
	Points.push_back(P5);
	Points.push_back(P6);
	Points.push_back(P7);
	Points.push_back(P8);

	return Points;
}

void WriteData(Agglomerative &A, ModelFile &Model)
{

	std::vector<unsigned int> PointLabels = A.getPointLabels();
	
	std::vector<unsigned int> UniqueLabels = Tools::UniqueElements(PointLabels);
	std::vector<Color<unsigned char> > Colors = Spectrum(UniqueLabels.size());

	for(unsigned int label = 0; label < UniqueLabels.size(); label++)
	{
		for(unsigned int p = 0; p < Model.NumPoints(); p++)
		{
			if(PointLabels[p] == UniqueLabels[label])
				Model.Points_[p].setColor(Colors[label]);
		}
	}
	
	Model.Write("Clusters.vtp");
	
	vsl_b_ofstream outputPDM("PointDistanceMatrix.bin");
	vsl_b_write(outputPDM, A.PointDistanceMatrix);
	outputPDM.close();
	
	vbl_array_2d<double> CDMValue(PointLabels.size(), PointLabels.size());
	vbl_array_2d<bool> CDMValid(PointLabels.size(), PointLabels.size());
	for(unsigned int i = 0; i < PointLabels.size(); i++)
	{
		for(unsigned int j = 0; j < PointLabels.size(); j++)
		{
			CDMValue(i,j) = A.ClusterDistanceMatrix(i,j).Value;
			CDMValid(i,j) = A.ClusterDistanceMatrix(i,j).Valid;
		}
	}
	vsl_b_ofstream outputCDMValue("ClusterDistanceMatrixValue.bin");
	vsl_b_write(outputCDMValue, CDMValue);
	outputCDMValue.close();
	
	vsl_b_ofstream outputCDMValid("ClusterDistanceMatrixValid.bin");
	vsl_b_write(outputCDMValid, CDMValid);
	outputCDMValid.close();
	/*
	vsl_b_ofstream outputCDM("ClusterDistanceMatrix.bin");
	vsl_b_write(outputCDM, A.ClusterDistanceMatrix);
	outputCDM.close();
	*/
	
	vsl_b_ofstream outputPL("PointLabels.bin");
	vsl_b_write(outputPL, PointLabels);
	outputPL.close();
}