#pragma once
#include "common.h"

int TestBSplineScatteredInterpolation2D(const char *inputImageFile, const char *outputImageFile, unsigned int splineOrder = 3, unsigned int numberOfControlPoints = 4) {

	typedef unsigned char   PixelType;
	const unsigned int ParametricDimension = 2;
	const unsigned int DataDimension = 1;

	typedef itk::Image< PixelType, ParametricDimension > ImageType;

	typedef float RealType;
	typedef itk::Vector<RealType, DataDimension> VectorType;
	typedef itk::Image<VectorType, ParametricDimension> VectorImageType;
	typedef itk::PointSet
		<VectorImageType::PixelType, ParametricDimension> PointSetType;
	PointSetType::Pointer pointSet = PointSetType::New();

	//Read image
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
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

	//QuickView viewer;
	//viewer.AddImage<ImageType>(reader->GetOutput());
	//viewer.Visualize();

	//Sample image
	itk::ImageRegionIteratorWithIndex<ImageType> It(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

	// Iterate through the input image which consists.

	for (It.GoToBegin(); !It.IsAtEnd(); ++It)
	{
		ImageType::IndexType idx = It.GetIndex();
		PointSetType::PixelType V(DataDimension);
		V[0] = static_cast<RealType>(It.Get());

		if (It.Get() != itk::NumericTraits<PixelType>::Zero)
		{
			// We extract both the 2-D location of the point 
			// and the pixel value of that point.  

			ImageType::IndexType idx = It.GetIndex();
			PointSetType::PointType point;
			reader->GetOutput()->TransformIndexToPhysicalPoint(It.GetIndex(), point);

			unsigned long i = pointSet->GetNumberOfPoints();
			pointSet->SetPoint(i, point);

			PointSetType::PixelType V(DataDimension);
			V[0] = static_cast<RealType>(It.Get());
			pointSet->SetPointData(i, V);
		}
		++It;
	}


	// Instantiate the B-spline filter and set the desired parameters.
	typedef itk::BSplineScatteredDataPointSetToImageFilter	<PointSetType, VectorImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetSplineOrder(splineOrder);
	
	FilterType::ArrayType ncps;
	ncps.Fill(numberOfControlPoints);
	filter->SetNumberOfControlPoints(ncps);
	filter->SetNumberOfLevels(3);

	// Define the parametric domain.
	filter->SetOrigin(reader->GetOutput()->GetOrigin());
	filter->SetSpacing(reader->GetOutput()->GetSpacing());
	filter->SetSize(reader->GetOutput()->GetLargestPossibleRegion().GetSize());

	filter->SetInput(pointSet);

	try
	{
		filter->Update();
	}
	catch (...)
	{
		std::cerr << "TestBSplineScatteredInterpolation: itkBSplineScatteredDataImageFilter exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}

	// Write the output to an image.	
	typedef itk::Image< unsigned char, 2 >   UnsignedCharImage2D;

	typedef itk::VectorMagnitudeImageFilter<VectorImageType2D, UnsignedCharImage2D>  VectorMagnitudeFilterType;
	VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
	magnitudeFilter->SetInput(filter->GetOutput());
	// To write the magnitude image file, we should rescale the gradient values
	// to a reasonable range
	typedef itk::RescaleIntensityImageFilter<UnsignedCharImage2D, UnsignedCharImage2D> rescaleFilterType;
	rescaleFilterType::Pointer rescaler = rescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(magnitudeFilter->GetOutput());

	typedef itk::ImageFileWriter<UnsignedCharImage2D> WriterType;
	WriterType::Pointer writeroutput = WriterType::New();
	writeroutput->SetInput(rescaler->GetOutput());
	writeroutput->SetFileName(outputImageFile);
	try
	{
		writeroutput->Update();
	}
	catch (...)
	{
		std::cerr << "Writer exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int TestBSplineScatteredInterpolation3D(const char *inputImageFile, const char *outputImageFile, unsigned int splineOrder = 3, unsigned int numberOfControlPoints = 4) {


	PointSetType3D::Pointer pointSet = PointSetType3D::New();

	//Read image
	typedef itk::ImageFileReader<UCharImageType3D> ReaderType;
	UCharReaderType3D::Pointer reader = ReaderType::New();
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

	//QuickView viewer;
	//viewer.AddImage<ImageType>(reader->GetOutput());
	//viewer.Visualize();

	//Sample image
	itk::ImageRegionIteratorWithIndex<UCharImageType3D> It(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

	// Iterate through the input image which consists.
	for (It.GoToBegin(); !It.IsAtEnd(); ++It)
	{
		UCharImageType3D::IndexType idx = It.GetIndex();
		PointSetType3D::PixelType V(1);
		V[0] = static_cast<RealType>(It.Get());

		if (It.Get() != itk::NumericTraits<UCharType>::Zero)
		{
			// We extract both the 3-D location of the point 
			// and the pixel value of that point.  

			UCharImageType3D::IndexType idx = It.GetIndex();
			PointSetType3D::PointType point;
			reader->GetOutput()->TransformIndexToPhysicalPoint(It.GetIndex(), point);

			unsigned long i = pointSet->GetNumberOfPoints();
			pointSet->SetPoint(i, point);

			PointSetType3D::PixelType V(1);
			V[0] = static_cast<RealType>(It.Get());
			pointSet->SetPointData(i, V);
		}
		++It;
	}


	// Instantiate the B-spline filter and set the desired parameters.
	typedef itk::BSplineScatteredDataPointSetToImageFilter	<PointSetType3D, VectorImageType3D> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetSplineOrder(splineOrder);

	FilterType::ArrayType ncps;
	ncps.Fill(numberOfControlPoints);
	filter->SetNumberOfControlPoints(ncps);
	filter->SetNumberOfLevels(3);

	// Define the parametric domain.
	filter->SetOrigin(reader->GetOutput()->GetOrigin());
	filter->SetSpacing(reader->GetOutput()->GetSpacing());
	filter->SetSize(reader->GetOutput()->GetLargestPossibleRegion().GetSize());

	filter->SetInput(pointSet);

	try
	{
		filter->Update();
	}
	catch (...)
	{
		std::cerr << "TestBSplineScatteredInterpolation: itkBSplineScatteredDataImageFilter exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}

	// Write the output to an image.	
	typedef itk::VectorMagnitudeImageFilter<VectorImageType3D, UCharImageType3D>  VectorMagnitudeFilterType;
	VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
	magnitudeFilter->SetInput(filter->GetOutput());
	// To write the magnitude image file, we should rescale the gradient values
	// to a reasonable range
	typedef itk::RescaleIntensityImageFilter<UCharImageType3D, UCharImageType3D> rescaleFilterType;
	rescaleFilterType::Pointer rescaler = rescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(magnitudeFilter->GetOutput());

	typedef itk::ImageFileWriter<UCharImageType3D> WriterType;
	WriterType::Pointer writeroutput = WriterType::New();
	writeroutput->SetInput(rescaler->GetOutput());
	writeroutput->SetFileName(outputImageFile);
	try
	{
		writeroutput->Update();
	}
	catch (...)
	{
		std::cerr << "Writer exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int Test2DplusTime()
{
	const unsigned int ParametricDimension = 3;
	const unsigned int DataDimension = 1;

	typedef unsigned char PixelType;
	typedef float RealType;
	typedef itk::Vector<RealType, DataDimension> VectorType;
	typedef itk::Image<PixelType, ParametricDimension> ImageType;

	typedef itk::PointSet<VectorType, ParametricDimension> PointSetType;
	PointSetType::Pointer pointSet = PointSetType::New();
	unsigned int r_min = 70;
	unsigned int r_max = 100;
	double PI = 3.14159;
	ImageType::PointType inputOrigin;
	inputOrigin[0] = 0 - (200 - 1) / 2;
	inputOrigin[1] = 0 - (200 - 1) / 2;
	inputOrigin[2] = 0;
	
	//filterd output
	ImageType::SizeType outputSize;
	outputSize[0] = 200;
	outputSize[1] = 200;
	outputSize[2] = 15;
	ImageType::SpacingType outputSpacing;
	outputSpacing[0] = 1.0;
	outputSpacing[1] = 1.0;
	outputSpacing[2] = 1.0/ outputSize[2];
	ImageType::PointType outputOrigin;
	outputOrigin[0] = 0;
	outputOrigin[1] = 0;
	outputOrigin[2] = 0;

	
	int i = 0;
	for (float t = 0; t < 1.0; t += 0.2) {
		//temp save input points to image
		ImportFilterTypeUChar2D::Pointer importFilter = ImportFilterTypeUChar2D::New();
		ImportFilterTypeUChar2D::SizeType  importSize;
		importSize[0] = 200;  // size along X
		importSize[1] = 200;  // size along Y
		ImportFilterTypeUChar2D::IndexType importStart;
		importStart.Fill(0);
		ImportFilterTypeUChar2D::RegionType importRegion;
		importRegion.SetIndex(importStart);
		importRegion.SetSize(importSize);
		importFilter->SetRegion(importRegion);
		const itk::SpacePrecisionType importOrigin[2] = { 0, 0 };
		importFilter->SetOrigin(importOrigin);
		const itk::SpacePrecisionType  importSpacing[2] = { 1.0, 1.0 };
		importFilter->SetSpacing(importSpacing);
		const unsigned int numberOfPixels = importSize[0] * importSize[1];
		UCharType * localBuffer = new UCharType[numberOfPixels];
		memset(localBuffer, 0, sizeof(UCharType)*numberOfPixels);
		UCharType * it = localBuffer;
		float diagonal = sqrt(100*100 + 100*100);
		//circle filled with points
		unsigned int r = r_min + 30 * t;
		for (unsigned int y = 0; y < outputSize[1]; y += 2)
		{
			const double dy = static_cast<double>(y) - static_cast<double>(outputSize[1]) / 2.0;
			for (unsigned int x = 0; x < outputSize[0]; x += 2)
			{
				const double dx = static_cast<double>(x) - static_cast<double>(outputSize[0]) / 2.0;
				const double d = sqrt(dx*dx + dy*dy);
				
				if (d < r) 
				{
					double value = 255 * cos(pi*d / (2 * r));
					it[y*outputSize[0] + x] = value;

					PointSetType::PointType point;
					point[0] = x * importSpacing[0];
					point[1] = y * importSpacing[1];
					point[2] = t;
					unsigned long idx = pointSet->GetNumberOfPoints();
					pointSet->SetPoint(idx, point);
					VectorType V;
					V[0] = value;
					pointSet->SetPointData(idx, V);
				}
			}
		}

		//circle with points on circumference
		//for (unsigned int angle = 0; angle < 360; angle += 10){
		//	unsigned int x1 = (outputSize[0] - 1)/2 + r * cos(angle * PI / 180);
		//	unsigned int y1 = (outputSize[1] - 1)/2 + r * sin(angle * PI / 180);
		//	unsigned long idx = pointSet->GetNumberOfPoints();

		//	PointSetType::PointType point;
		//	point[0] = x1;
		//	point[1] = y1;
		//	point[2] = t;
		//	pointSet->SetPoint(idx, point);
		//	VectorType V;
		//	V[0] = 255;
		//	pointSet->SetPointData(idx, V);

		//	//temp
		//	unsigned int pixelIndex = y1 * importSize[0] + x1;
		//	it[pixelIndex] = 255;
		//}
		//temp
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
		char buffer[256]; sprintf(buffer, "%03d", i);
		i++;
		std::string path = std::string("test/test2D_") + std::string(buffer) + std::string(".png");
		writer->SetFileName(path);
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
	}

	
	// Instantiate the filter and set the parameters
	typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, VectorImageType3D>  FilterType;
	FilterType::Pointer filter = FilterType::New();

	// Define the parametric domain	
	filter->SetSize(outputSize);	
	filter->SetSpacing(outputSpacing);	
	filter->SetOrigin(outputOrigin);	
	filter->SetInput(pointSet);

	filter->SetSplineOrder(3);
	FilterType::ArrayType ncps;
	ncps.Fill(4);
	filter->SetNumberOfControlPoints(ncps);
	filter->SetNumberOfLevels(3);
	//filter->SetNumberOfThreads(1);
	//filter->SetGenerateOutputImage(false);

	try
	{
		filter->Update();
	}
	catch (...)
	{
		std::cerr << "TestBSplineScatteredInterpolation: itkBSplineScatteredDataImageFilter exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}
	

	// Write the output to an image.	
	typedef itk::VectorMagnitudeImageFilter<VectorImageType3D, UCharImageType3D>  VectorMagnitudeFilterType;
	VectorMagnitudeFilterType::Pointer magnitudeFilter = VectorMagnitudeFilterType::New();
	magnitudeFilter->SetInput(filter->GetOutput());
	// To write the magnitude image file, we should rescale the gradient values
	// to a reasonable range
	typedef itk::RescaleIntensityImageFilter<UCharImageType3D, UCharImageType3D> rescaleFilterType;
	rescaleFilterType::Pointer rescaler = rescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(magnitudeFilter->GetOutput());

	typedef itk::ImageFileWriter<UCharImageType3D> WriterType;
	WriterType::Pointer writeroutput = WriterType::New();
	writeroutput->SetInput(rescaler->GetOutput());
	writeroutput->SetFileName("test/test.mhd");
	try
	{
		writeroutput->Update();
	}
	catch (...)
	{
		std::cerr << "Writer exception thrown"
			<< std::endl;
		return EXIT_FAILURE;
	}


	//Series 3D writer
	UCharSeriesWriterType3D::Pointer writer = UCharSeriesWriterType3D::New();
	std::cout << "Writing 2D images..." << endl;
	NameGeneratorType::Pointer OutputNameGenerator = NameGeneratorType::New();
	OutputNameGenerator->SetSeriesFormat(std::string("test/Reconstructed2D_") + "%03d.png");
	OutputNameGenerator->SetStartIndex(0);
	OutputNameGenerator->SetEndIndex(outputSize[2] - 1);
	OutputNameGenerator->SetIncrementIndex(1);

	writer->SetFileNames(OutputNameGenerator->GetFileNames());
	writer->SetInput(rescaler->GetOutput());
	//writer->SetUseCompression(true);
	//writer->UseCompressionOn();
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
};















