#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once
#include <mpi.h>
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
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

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
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}


int main()
{
	MPI_Init(NULL, NULL);
	//data
	int ImageWidth = 4, ImageHeight = 4;
	int ImageSize=4;
	//selecting data
	int threshold;
	int chosenThresholdImageNum = 158;
	//image containters
	int** totalImageArray=NULL;//M*n^2 holds all the data
	int* chosenThresholdImage=NULL;//hold the chosen image
	//495 pictures
	const int M = 495;
	//image path
	System::String^ imagePath;
	std::string imagePathFirstPart;
	imagePathFirstPart = "..//Data//Input//in000";
	string imageExtenstion = ".jpg";
	//mpi data
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//for the master proccessor
	if (rank == 0)
	{
		threshold = 20;
		//read selected thresh hold image
		//create path
		string NewPath = imagePathFirstPart + to_string(chosenThresholdImageNum) + imageExtenstion;
		imagePath = marshal_as<System::String^>(NewPath);
		//read it
		chosenThresholdImage = inputImage(&ImageWidth, &ImageHeight, imagePath);
		//assign image size
		ImageSize = ImageWidth * ImageHeight;
		//
		totalImageArray = new int* [M];
		//read all data
		for (int i = 1; i <= M; i++)
		{
			// fill index to be 3 digit (1 => 001 , 99 => 099)
			// get image path = orignal+index+extention
			string s;
			if (i <= 9) { s = "00" + to_string(i); }
			else if (i <= 99) { s = "0" + to_string(i); }
			else { s = to_string(i); }
			string NewPath = imagePathFirstPart + s + imageExtenstion;
			cout << NewPath << endl;
			// load it
			imagePath = marshal_as<System::String^>(NewPath);
			totalImageArray[i - 1] = inputImage(&ImageWidth, &ImageHeight, imagePath);
		}
	}
	cout << "thresh old " << threshold;
	int* imageArrayAverage = new int[ImageSize];
	int* imageArraySubtraction = new int[ImageSize];
	MPI_Bcast(&threshold, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&chosenThresholdImage, ImageSize, MPI_INT, 0, MPI_COMM_WORLD);
	//here start distrubuting the data
	//each proccessor will take portions(columns)
	int localArraySize = ImageSize / size;
	int** localArray = new int*[M];
	//for each image
	for (int i = 0; i < M; ++i)
	{    
		//master proccessor
		if (rank == 0)
		{
			//send first n columns in image i
			for (int procRank = 1; procRank < size; procRank++)
				MPI_Send(&totalImageArray[i][localArraySize * procRank], localArraySize, MPI_INT, procRank, 0, MPI_COMM_WORLD);
		
			if (i==M-1)
			{
				for (int k = 1; k < size; k++)
				{
					MPI_Status status;
					MPI_Send(&chosenThresholdImage[localArraySize * k], localArraySize, MPI_INT, k, 0, MPI_COMM_WORLD);
					MPI_Recv(&imageArrayAverage[k * localArraySize], localArraySize, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&imageArraySubtraction[k * localArraySize], localArraySize, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
				}
				createImage(imageArrayAverage, ImageWidth, ImageHeight, 0);
				createImage(imageArraySubtraction, ImageWidth, ImageHeight, 1);
			}
		}
		if (rank < size)
		{
			MPI_Status status;
			localArray[i] = new int[localArraySize];
			//receive n columnns in image i
			MPI_Recv(localArray[i], localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			//when i receive all the data in M images
			//i = 495  m=496-1=>495 final image
			if (i == M - 1)
			{
				int localSum = 0;
				//create 2 arrays 1-backgroundMean B 2-resultOfSubtraction
				int* localBackGroundMean = new int[M];
				int* localImageSubtraction = new int[M];
				//get the value of threshold
				int localThreshold=threshold;
				//MPI_Recv(&localThreshold, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				//get the chosen image
				int* localChosenImage = new int[ImageSize];
				MPI_Recv(&localChosenImage, ImageSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				//here sum all of them
				//for each pixel (column)
				for (int localColumn = 0; localColumn < localArraySize; localColumn++)
				{
					//for each image (row)
					for (int localRow = 0; localRow < M; localRow++)
					{
						localSum += localArray[localRow][localColumn];
					}
					//got the sum of pixel (in column) then get mean
					localSum /= M;
					localBackGroundMean[localColumn] = localSum;
					localSum = 0;
					//get the subtraction
					int subResult = abs(localBackGroundMean[localColumn] - localChosenImage[localColumn]);
					if (subResult > localThreshold)
					{
						localImageSubtraction[localColumn] = 0;
					}
					else
					{
						localImageSubtraction[localColumn] = subResult;
					}
				}
				MPI_Send(localBackGroundMean, localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(localChosenImage, localArraySize, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}
	}
	cout << "finished";
	MPI_Finalize();
	system("pause");
	return 0;
}



