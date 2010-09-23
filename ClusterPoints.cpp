#include <iostream>
#include <string>
#include <vector>

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

std::vector<OrientedPoint> GenerateData();

void WriteData(Agglomerative &A, ModelFile &Model);

void CheckRequiredArgs(const po::variables_map vm);

int main(int argc, char *argv[])
{
	std::string InputFile;
	double maxsize;
	
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "Help message.")
			("input", po::value<std::string>(&InputFile), "Set input file")
			//("output", po::value<std::string>(&OutputFile), "Set output file")
			//("x", po::value<double>(&X)->default_value(0.0), "Center X")
			//("y", po::value<double>(&Y)->default_value(0.0), "Center Y")
			//("z", po::value<double>(&Z)->default_value(0.0), "Center Z")
			("maxsize", po::value<double>(&maxsize), "Max size between clusters.")
			;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	CheckRequiredArgs(vm);

	//create the model
	ModelFile Model(InputFile);

	TestAgglomerativeClusterData(Model);

	return 0;
}


void CheckRequiredArgs(const po::variables_map vm)
{
	if(!vm.count("input"))
	{
		std::cout << "input is required!" << std::endl;
		exit(-1);
	}
/*	
	if(!vm.count("output"))
	{
		std::cout << "output is required!" << std::endl;
		exit(-1);
	}
*/
}

#include <vsl/vsl_vector_io.txx>
#include <vbl/io/vbl_io_array_2d.txx>
				 VSL_VECTOR_IO_INSTANTIATE(vbl_array_2d<double>);
		 VSL_VECTOR_IO_INSTANTIATE(vbl_array_2d<bool>);

//VSL_VECTOR_IO_INSTANTIATE(vbl_array_2d<ValidType<double> >);
	    

void TestAgglomerativeClusterData(ModelFile &Model)
{
	std::vector<OrientedPoint> ModelPoints = Model.getPoints();
	Agglomerative A(ModelPoints, 0.7);
	A.PerformClustering();
	
	WriteData(A, Model);
	
	std::cout << "There are : " << Tools::UniqueElements(A.getPointLabels()).size() << " clusters." << std::endl;

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