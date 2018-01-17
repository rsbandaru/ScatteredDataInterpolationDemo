
#pragma once
//Generate input image
#include "SimulateNoisyImage.h"
#include "TestBSplineScatteredInterpolation.h"


#include <ctime>

#include <iostream>

int main() {
	clock_t begin = std::clock();

	int ret;

	//Test 2d lena simulation
	ret = SampleImagePolar2D("D:/OneDrive/Development/c++/ScatteredInterpolationDemo/lenaTest3.png", "D:/OneDrive/Development/c++/ScatteredInterpolationDemo/LinaSampled.png",60, 200);
	ret = TestBSplineScatteredInterpolation2D("D:/OneDrive/Development/c++/ScatteredInterpolationDemo/LinaSampled.png", "D:/OneDrive/Development/c++/ScatteredInterpolationDemo/LinaSampledOutput.png", 3, 64);

	
	clock_t end = std::clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Elapsed Seconds = " << elapsed_secs;

	system("pause");


	return 0;
}

