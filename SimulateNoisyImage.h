#pragma once

#include "common.h"


int SampleImagePolar2D(const char *inputImageFile, const char *outputImageFile, unsigned int numberOfAngles, int percentNoise)
{

	// We select the data type used to represent the image pixels. We
	// assume that the external block of memory uses the same data type to
	// represent the pixels.
	const unsigned int ParametricDimension = 2;
	const unsigned int DataDimension = 1;

	//Read image
	UCharReaderType2D::Pointer reader = UCharReaderType2D::New();
	reader->SetFileName(inputImageFile);

	try
	{
		reader->Update();
	}
	catch (...)
	{
		std::cerr << "TestBSplineScatteredInterpolation: reader exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}


	UCharImageType2D::Pointer inputImage = reader->GetOutput();
	UCharImageType2D::RegionType intputRegion = inputImage->GetLargestPossibleRegion();
	UCharImageType2D::SizeType inputSize = intputRegion.GetSize();
	UCharImageType2D::SpacingType  inputSpacing ;
	UCharImageType2D::IndexType inputIndex = inputImage->GetLargestPossibleRegion().GetIndex();
	UCharImageType2D::PointType inputOrigin;
	inputSpacing[0] = 1.0;
	inputSpacing[1] = 1.0;
	inputImage->SetSpacing(inputSpacing);
	inputOrigin[0] = 0; // inputIndex[0] - inputSpacing[0] * (inputSize[0] - 1) / 2;//x  																			  
	inputOrigin[1] = 0;//y
	inputImage->SetOrigin(inputOrigin);

	//Create Output image
	ImportFilterTypeUChar2D::Pointer importFilter = ImportFilterTypeUChar2D::New();
	ImportFilterTypeUChar2D::SizeType  outputSize = inputSize;
	ImportFilterTypeUChar2D::IndexType start;
	start.Fill(0);
	ImportFilterTypeUChar2D::RegionType outputRegion;
	outputRegion.SetIndex(start);
	outputRegion.SetSize(outputSize);
	importFilter->SetRegion(intputRegion);

	const itk::SpacePrecisionType outputOrigin[ParametricDimension] = { 0.0, 0.0 };
	importFilter->SetOrigin(outputOrigin);
	const itk::SpacePrecisionType  outputSpacing[ParametricDimension] = { 1.0, 1.0 };
	importFilter->SetSpacing(outputSpacing);

	const unsigned int numberOfPixels = outputSize[0] * outputSize[1];
	UCharType * localBuffer = new UCharType[numberOfPixels];
	memset(localBuffer, 0, sizeof(UCharType)*numberOfPixels);
	
	// Here we fill up the buffer with a binary sphere.
	UCharType * it = localBuffer;
	

	//Reset random pixels locations to zeros
	//Extract and random angles
	srand((unsigned)time(0));
	unsigned int lowest_point_a = 0, highest_point_a = 180;
	unsigned int range_point_a = (highest_point_a - lowest_point_a) + 1;
	const double radius = min(inputSize[0], inputSize[1])/2;
	UCharImageType2D::IndexType pixelIndex;
	for (unsigned int a = 0; a < numberOfAngles; a++)
	{
		unsigned int angle = lowest_point_a + int(range_point_a*rand() / (RAND_MAX + 1.0));

		for (int r = -radius; r < radius; r++)
		{
			unsigned int angle360;
			if (r > 0) {
				angle360 = fmod(angle + 180.0, 360.0);
			}
			else {
				angle360 = angle;
			}

			double x, y;
			Cylindrical2Cartesian(angle360, abs(r), x, y);

			//Regularize voxel location to
			x = x + inputSpacing[0] * (inputSize[0] ) / 2;//shift origin
			y = y + inputSpacing[0] * (inputSize[0] ) / 2;
			unsigned int X, Y;
			X = (unsigned int)x / outputSpacing[0];
			Y = (unsigned int)y / outputSpacing[1];			
			X = min(X, (unsigned int)inputSize[0]);
			Y = min(Y, (unsigned int)inputSize[1]);
			unsigned int pixelIndexOutput = Y * outputSize[0] + X;
			pixelIndex[0] = X;
			pixelIndex[1] = Y;
			

			it[pixelIndexOutput] = inputImage->GetPixel(pixelIndex);
		}
	}

	const bool importImageFilterWillOwnTheBuffer = true;
	importFilter->SetImportPointer(localBuffer, numberOfPixels, importImageFilterWillOwnTheBuffer);

	try
	{
		importFilter->Update();
	}
	catch (itk::ExceptionObject & exp)
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << exp << std::endl;
		return EXIT_FAILURE;
	}



	// Finally, we can connect the output of this filter to a pipeline.
	// For simplicity we just use a writer here, but it could be any other filter.
	typedef itk::ImageFileWriter< UCharImageType2D > WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName(outputImageFile);
	writer->SetInput(importFilter->GetOutput());
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & exp)
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << exp << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


