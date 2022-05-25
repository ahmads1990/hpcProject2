#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once
#include<mpi.h>

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{

	int* input;
	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}

void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}
static string ext = ".jpg";

//this function returns M*n^2 array m number of images
int** getImages(int n, int* w, int* h) //put the size of image in w & h
{
	int** image_array = new int* [n];
	System::String^ imagePath;
	string img, image_num, image;
	int width = 0, height = 0;
	for (int i = 1; i <= n; i++)
	{
		if (i <= 9) { image_num = "00" + to_string(i); }
		else if (i <= 99) { image_num = "0" + to_string(i); }
		else { image_num = to_string(i); }
		image = image_num + ext;
		img = "..//Data//Input//in000" + image;

		imagePath = marshal_as<System::String^>(img);
		int* imageData = inputImage(&*w, &*h, imagePath);

		image_array[i - 1] = new int[*w * *h];
		image_array[i - 1] = imageData;
	}
	return image_array;
}

int main()
{

	int ImageWidth = 0, ImageHeight = 0;
	int world_size;
	int world_rank;
	int ImageSize = 0;
	int start_s, stop_s, TotalTime = 0;
	// Initialize the MPI env
	MPI_Init(NULL, NULL);
	// Get processes number
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	// Get processor rank
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	int** TotalImageArray = 0;
	int* thresholdImage = 0;
	int IMAGE_NUMBER;
	int threshold;

	if (world_rank == 0)
	{
		//assign data
		IMAGE_NUMBER = 495;
		//selected threshbhold
		threshold = 50;
		// Sequential part
		TotalImageArray = getImages(IMAGE_NUMBER, &ImageWidth, &ImageHeight);

		//read threshold image
		string img = "..//Data//Input//in000158" + ext;
		System::String^ imagePath = marshal_as<System::String^>(img);
		thresholdImage = inputImage(&ImageWidth, &ImageHeight, imagePath);

		ImageSize = ImageWidth * ImageHeight;
	}
	//send number of images , full size of image, threshold number to all proccesors
	MPI_Bcast(&IMAGE_NUMBER, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ImageSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&threshold, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Buffer to store image data  inside each proccessor
	int** localArray = new int* [IMAGE_NUMBER];
	//store final result of the program
	int* backgroundAverage = new int[ImageSize];
	int* foregroundSubtraction = new int[ImageSize];
	//split the imagesize( columns) into parts for each proccessor
	int localArraySize = (ImageSize) / world_size;
	start_s = clock();
	//for each image
	for (int i = 0; i < IMAGE_NUMBER; i++)
	{
		//master proccessor
		if (world_rank == 0)
		{
			//for each image send (first n*size) columns for each n proccessor 
			for (int j = 1; j < world_size; j++)
			{
				MPI_Send(&TotalImageArray[i][j * localArraySize], localArraySize, MPI_INT, j, 0, MPI_COMM_WORLD);
			}

			//if processor 0 took all the first parts from all the images start summing them and divide it by the number of images
			int sum = 0;
			for (int x = 0; x < localArraySize; x++) {
				for (int k = 0; k < IMAGE_NUMBER; k++) {

					sum += TotalImageArray[k][x];

				}
				sum /= IMAGE_NUMBER;
				backgroundAverage[x] = sum;
				int diff = abs(backgroundAverage[x] - thresholdImage[x]);
				if (diff > threshold)
				{
					foregroundSubtraction[x] = diff;
				}

				else
				{
					foregroundSubtraction[x] = 0;
				}


			}
			//final step where master collects image data from other proccesors
			if (i == IMAGE_NUMBER - 1) {
				//for each proccessor send threshold image and receive its part of back,forground
				for (int k = 1; k < world_size; k++)
				{
					MPI_Status status;
					MPI_Send(&thresholdImage[k * localArraySize], localArraySize, MPI_INT, k, 0, MPI_COMM_WORLD);
					MPI_Recv(&backgroundAverage[k * localArraySize], localArraySize, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&foregroundSubtraction[k * localArraySize], localArraySize, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
				}
				//when done create those images
				createImage(backgroundAverage, ImageWidth, ImageHeight, 0);
				createImage(foregroundSubtraction, ImageWidth, ImageHeight, 1);
			}
		}
		else if (world_rank <= world_size)
		{
			//each processor recieve part of the image and save it in local image array
			localArray[i] = new int[localArraySize];

			MPI_Status status;
			MPI_Recv(localArray[i], localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			//last iteration has all data start proc
			if (i == IMAGE_NUMBER - 1)
			{
				//get threshold
				int* local_threshold = new int[localArraySize];
				MPI_Recv(local_threshold, localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

				int sum = 0;
				int* localbackgroundimage = new int[localArraySize];
				int* local_forground = new int[localArraySize];

				for (int column = 0; column < localArraySize; column++) {
					for (int k = 0; k < IMAGE_NUMBER; k++) {

						sum += localArray[k][column];

					}
					sum /= IMAGE_NUMBER;
					localbackgroundimage[column] = sum;
					int diff = abs(localbackgroundimage[column] - local_threshold[column]);
					if (diff > threshold) { local_forground[column] = diff; }
					else { local_forground[column] = 0; }
				}
				MPI_Send(localbackgroundimage, localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(local_forground, localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}

	}

	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	cout << "time: " << TotalTime << endl;
	MPI_Finalize();
	system("pause");
	return 0;
}