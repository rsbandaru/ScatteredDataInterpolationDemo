#pragma once



//Testing create random image
#include <itkRandomImageSource.h>

//Testing bspline
#include <cstdlib> 
#include <ctime> 

//Testing Create image from buffer
#include "itkImage.h"
#include "itkImportImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"

#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkOrientImageFilter.h"

#include "itkVectorMagnitudeImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <itkImageLinearConstIteratorWithIndex.h>

#include <vector>
#include "sort.h"

using namespace std;

/************************************************************************/
/*                        function declarations                         */
/************************************************************************/
double ModuloDouble(double a, double b);
template<typename Ttype> int ReadTXTToVector(std::vector<Ttype> &vectorData, std::string filename);
int WriteVectorToTXT(std::vector<double> vectorData, std::string filename);

/************************************************************************/
/* Setup constants and variables that will be defined by header file.   */
/************************************************************************/
const double pi = 3.14159265;
const unsigned int inputDimension = 3; // input image dimension
/************************************************************************/
/* Start type definitions.
/************************************************************************/
/** Define image type	*/
typedef unsigned char UCharType;
typedef short ShortType;
typedef itk::Image< UCharType, 1 >	UCharImageType1D;
typedef itk::Image< UCharType, 2 >	UCharImageType2D;
typedef itk::Image< UCharType, 3 >	UCharImageType3D;
typedef itk::Image< UCharType, 4 >	UCharImageType4D;
typedef itk::Image< UCharType, 6 >	UCharImageType6D;
typedef itk::Image< ShortType, 1 >	ShortImageType1D;
typedef itk::Image< ShortType, 2 >	ShortImageType2D;
typedef itk::Image< ShortType, 3 >	ShortImageType3D;
typedef itk::Image< ShortType, 4 >	ShortImageType4D;

typedef float RealType;
typedef itk::Vector<RealType, 1> VectorType1D;
typedef itk::Image< VectorType1D, 1 > VectorImageType1D;
typedef itk::Image< VectorType1D, 2 > VectorImageType2D;
typedef itk::Image< VectorType1D, 3 > VectorImageType3D;
typedef itk::Image< VectorType1D, 4 > VectorImageType4D;
typedef itk::PointSet<VectorImageType2D::PixelType, 1> PointSetType1D;
typedef itk::PointSet<VectorImageType2D::PixelType, 2> PointSetType2D;
typedef itk::PointSet<VectorImageType3D::PixelType, 3> PointSetType3D;
typedef itk::PointSet<VectorImageType4D::PixelType, 4> PointSetType4D;
PointSetType2D::Pointer pointSet2D = PointSetType2D::New();
PointSetType3D::Pointer pointSet3D = PointSetType3D::New();
PointSetType4D::Pointer pointSet4D = PointSetType4D::New();


/** Define image reader, name generator, writer and image filter */
typedef itk::ImageSeriesReader< ShortImageType2D >	ShortSeriesReaderType2D;
typedef itk::ImageFileReader< ShortImageType2D > ShortReaderType2D;
typedef itk::ImageSeriesReader< ShortImageType3D >	ShortSeriesReaderType3D;
typedef itk::ImageFileReader< ShortImageType3D > ShortReaderType3D;
typedef itk::ImageSeriesReader< ShortImageType4D >	ShortSeriesReaderType4D;
typedef itk::ImageFileReader< ShortImageType4D > ShortReaderType4D;
typedef itk::ImageSeriesReader< UCharImageType2D >	UCharSeriesReaderType2D;
typedef itk::ImageFileReader< UCharImageType2D > UCharReaderType2D;
typedef itk::ImageSeriesReader< UCharImageType3D >	UCharSeriesReaderType3D;
typedef itk::ImageFileReader< UCharImageType3D > UCharReaderType3D;
typedef itk::ImageSeriesReader< UCharImageType4D >	UCharSeriesReaderType4D;
typedef itk::ImageFileReader< UCharImageType4D > UCharReaderType4D;
typedef itk::NumericSeriesFileNames	NameGeneratorType;
typedef itk::ImageSeriesWriter< UCharImageType3D, UCharImageType2D >	UCharSeriesWriterType3D;
typedef itk::ImageSeriesWriter< UCharImageType4D, UCharImageType3D >	UCharSeriesWriterType4D;
typedef itk::ImageSeriesWriter< ShortImageType4D, ShortImageType3D >	ShortSeriesWriterType4D;
typedef itk::ImageSeriesWriter< ShortImageType3D, ShortImageType2D >	ShortSeriesWriterType3D;

typedef itk::ImageFileWriter< UCharImageType2D > UCharWriterType2D;
typedef itk::ImageFileWriter< UCharImageType3D > UCharWriterType3D;
typedef itk::ImageFileWriter< ShortImageType2D > ShortWriterType2D;
typedef itk::ImageFileWriter< ShortImageType3D > ShortWriterType3D;

typedef itk::ImageSliceConstIteratorWithIndex< ShortImageType3D> ShortSliceIteratorType;

typedef itk::ExtractImageFilter<ShortImageType4D, ShortImageType3D> ExtractFilterType3D;
typedef itk::ExtractImageFilter<ShortImageType3D, ShortImageType2D> ExtractFilterType2D;
typedef itk::ExtractImageFilter<ShortImageType2D, ShortImageType1D> ExtractFilterType1D;
typedef itk::ExtractImageFilter<UCharImageType4D, UCharImageType3D> UCharExtractFilterType3D;
typedef itk::ExtractImageFilter<UCharImageType3D, UCharImageType2D> UCharExtractFilterType2D;
typedef itk::ExtractImageFilter<UCharImageType2D, UCharImageType1D> UCharExtractFilterType1D;
typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType3D, VectorImageType3D> BSplineFilterType3D;
typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType2D, VectorImageType2D> BSplineFilterType2D;

typedef itk::ImportImageFilter< UCharType, 2 >   ImportFilterTypeUChar2D;
typedef itk::ImportImageFilter< UCharType, 3 >   ImportFilterTypeUChar3D;
typedef itk::ImportImageFilter< UCharType, 4 >   ImportFilterTypeUChar4D;
typedef itk::ImportImageFilter< ShortType, 2 >   ImportFilterTypeShort2D;
typedef itk::ImportImageFilter< ShortType, 3 >   ImportFilterTypeShort3D;
typedef itk::ImportImageFilter< ShortType, 4 >   ImportFilterTypeShort4D;

typedef itk::ImageLinearConstIteratorWithIndex<ShortImageType1D> ShortImageType1DLinearIterator;
typedef itk::ImageLinearConstIteratorWithIndex<ShortImageType2D> ShortImageType2DLinearIterator;
typedef itk::ImageLinearConstIteratorWithIndex<ShortImageType3D> ShortImageType3DLinearIterator;
typedef itk::ImageLinearConstIteratorWithIndex<ShortImageType4D> ShortImageType4DLinearIterator;

typedef itk::ImageRegionIterator<ShortImageType1D> ShortImageType1DRegionIterator;
typedef itk::ImageRegionIterator<ShortImageType2D> ShortImageType2DRegionIterator;
typedef itk::ImageRegionIterator<ShortImageType3D> ShortImageType3DRegionIterator;
typedef itk::ImageRegionIterator<ShortImageType4D> ShortImageType4DRegionIterator;

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}


/** convert x, y from 3D Cartesian to cylindrical coordinates */
void Cartesian2Cylindrical(double x, double y, double &theta, double &radius)
{
	double pi = 3.14159265;
	double precision = 10000.0;
	double tempTheta = floor((std::atan2(y, x)*180.0 / pi) * precision + 0.5); //Raja -180 to 180
	tempTheta = tempTheta / precision;
	double tempRadius = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
	//Raja shift to 0 to 360
	if (tempTheta < 0.0) {
		theta = 360.0 + tempTheta;
		radius = tempRadius;
	}
	else {
		theta = tempTheta;
		radius = tempRadius;
	}

}

/** convert cylindrical coordinates to x, y from 3D Cartesian*/
void Cylindrical2Cartesian(double theta, double radius, double &x, double &y)
{
	double pi = 3.14159265;
	double precision = 10000.0;
	if (theta > 180)
		theta = theta - 360;//shift from (0 to 360) to (-180 to 180)
	x = radius*std::cos(theta*pi / 180.0);
	y = radius*std::sin(theta*pi / 180.0);
}

//Compare 2d array by first column
struct FirstColumnOnlyCmp
{
	bool operator()(const std::vector<double>& lhs,
		const std::vector<double>& rhs) const
	{
		return lhs[0] < rhs[0];
	}
};

ShortImageType3D::Pointer extract3DImageSlice(ShortImageType4D::Pointer reader, int plane, int slice) {
	ExtractFilterType3D::Pointer filter = ExtractFilterType3D::New();

	ShortImageType4D::RegionType inputRegion = reader->GetLargestPossibleRegion();

	ShortImageType4D::SizeType size = inputRegion.GetSize();
	size[plane] = 0;

	ShortImageType4D::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = slice;
	start[plane] = sliceNumber;

	ShortImageType4D::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);

	filter->SetInput(reader);
	filter->SetDirectionCollapseToIdentity();

	ShortImageType3D::Pointer img = filter->GetOutput();
	img->Update();

	return img;
}

ShortImageType2D::Pointer extract2DImageSlice(ShortImageType3D::Pointer reader, int plane, int slice) {
	ExtractFilterType2D::Pointer filter = ExtractFilterType2D::New();

	ShortImageType3D::RegionType inputRegion = reader->GetLargestPossibleRegion();

	ShortImageType3D::SizeType size = inputRegion.GetSize();
	size[plane] = 0;

	ShortImageType3D::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = slice;
	start[plane] = sliceNumber;

	ShortImageType3D::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);

	filter->SetInput(reader);
	filter->SetDirectionCollapseToIdentity();

	ShortImageType2D::Pointer img = filter->GetOutput();
	img->Update();

	return img;
}

UCharImageType2D::Pointer extract2DImageSlice(UCharImageType3D::Pointer reader, int plane, int slice) {
	UCharExtractFilterType2D::Pointer filter = UCharExtractFilterType2D::New();

	UCharImageType3D::RegionType inputRegion = reader->GetLargestPossibleRegion();

	UCharImageType3D::SizeType size = inputRegion.GetSize();
	size[plane] = 0;

	UCharImageType3D::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = slice;
	start[plane] = sliceNumber;

	UCharImageType3D::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);

	filter->SetInput(reader);
	filter->SetDirectionCollapseToIdentity();

	UCharImageType2D::Pointer img = filter->GetOutput();
	img->Update();

	return img;
}


ShortImageType1D::Pointer extract1DImageSlice(ShortImageType2D::Pointer reader, int plane, int slice) {
	ExtractFilterType1D::Pointer filter = ExtractFilterType1D::New();

	ShortImageType2D::RegionType inputRegion = reader->GetLargestPossibleRegion();

	ShortImageType2D::SizeType size = inputRegion.GetSize();
	size[plane] = 0;

	ShortImageType2D::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = slice;
	start[plane] = sliceNumber;

	ShortImageType2D::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);

	filter->SetInput(reader);
	filter->SetDirectionCollapseToIdentity();

	ShortImageType1D::Pointer img = filter->GetOutput();
	img->Update();

	return img;
}

UCharImageType1D::Pointer extract1DImageSlice(UCharImageType2D::Pointer reader, int plane, int slice) {
	UCharExtractFilterType1D::Pointer filter = UCharExtractFilterType1D::New();

	UCharImageType2D::RegionType inputRegion = reader->GetLargestPossibleRegion();

	UCharImageType2D::SizeType size = inputRegion.GetSize();
	size[plane] = 0;

	UCharImageType2D::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = slice;
	start[plane] = sliceNumber;

	UCharImageType2D::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	filter->SetExtractionRegion(desiredRegion);

	filter->SetInput(reader);
	filter->SetDirectionCollapseToIdentity();

	UCharImageType1D::Pointer img = filter->GetOutput();
	img->Update();

	return img;
}

/************************************************************************/
/* Function Definitions.                                                */
/************************************************************************/
double ModuloDouble(double a, double b)
{
	int result = static_cast<int>(a / b);
	double mod = a - static_cast<double>(result) * b;
	return mod;
}

template<typename Ttype> int ReadTXTToVector(std::vector<Ttype> &vectorData, std::string filename)
{
	std::ifstream textfile(filename);
	std::string tempString;
	int errorCode;
	double temp;
	if (textfile.is_open())
	{
		//std::cout << "Header file opened. Reading data." << std::endl;
		while (!textfile.eof())
		{
			textfile >> temp;
			vectorData.push_back(temp);
		}

		errorCode = 1;
	}
	else
	{
		std::cout << "Header file could not be created. Exit with error code -1." << std::endl;
		errorCode = -1;
	}

	textfile.close();
	return errorCode;
}
int WriteVectorToTXT(std::vector<double> vectorData, std::string filename)
{
	std::ofstream textfile(filename);
	int errorCode;

	if (textfile.is_open())
	{
		//std::cout << "Header file created. Writing data." << std::endl;
		std::vector<double>::iterator vectorIt = vectorData.begin();
		int counter = 0;

		for (vectorIt; vectorIt < vectorData.end(); ++vectorIt)
		{
			//std::cout << counter << std::endl;
			if (vectorIt == vectorData.end() - 1)
			{
				textfile << *vectorIt;
			}
			else
			{
				textfile << *vectorIt << std::endl;
			}
			counter++;
		}
		errorCode = 1;
	}
	else
	{
		//std::cout << "Header file could not be created. Exit with error code -1." << std::endl;
		errorCode = -1;
	}

	textfile.close();
	return errorCode;
}