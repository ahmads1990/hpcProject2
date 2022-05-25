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
static string ext = ".png";
//this function should take n images from the folder and then put them in a dynamic array
//the big array  is the number of the images and the small arrays is the size of each image 
int** inputVideoframe(int n, int* w, int* h) //put the size of image in w & h
{
	int** image_array = new int* [n];
	System::String^ imagePath;
	string img, image_num, image;
	int width = 0, hight = 0;
	for (int i = 1; i <= n; i++)
	{
		string s;
		if (i <= 9)
		{
			s = "00" + to_string(i);
		}
		else if (i <= 99)
		{
			s = "0" + to_string(i);
		}
		else
		{
			s = to_string(i);
		}
		image_num = s;
		image = image_num + ext;
		img = "..//Data//Input//in000" + s;
		string NewPath = img  + ".jpg";
		imagePath = marshal_as<System::String^>(NewPath);
		int* imageData = inputImage(&*w, &*h, imagePath);

		image_array[i - 1] = new int[*w * *h];
		image_array[i - 1] = imageData;


	}

	return image_array;
}
int main()
{

	MPI_Init(NULL, NULL);
	//getchar();
	int w = 0, h = 0;

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int full_size_of_img = 0;


	int** image_array = 0;
	int* threshold_img = 0;
	int imagenumber;
	int threshold;

	if (world_rank == 0)
	{
		//to input the image 
		imagenumber = 495;
		threshold = 20;
		image_array = inputVideoframe(imagenumber, &w, &h);

		//threshold pic
		string img = "..//Data//BackGround//in000160" + ext;
		System::String^ imagePath = marshal_as<System::String^>(img);
		threshold_img = inputImage(&w, &h, imagePath);

		full_size_of_img = w * h;

	}
	MPI_Bcast(&imagenumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&full_size_of_img, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&threshold, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int** local_image_array = new int* [imagenumber];
	int* Resultbackgroundimage = new int[full_size_of_img];
	int* Resultforeground = new int[full_size_of_img];
	int local_arr_size = (full_size_of_img) / world_size;

	for (int i = 0; i < imagenumber; i++)
	{
		if (world_rank == 0)
		{

			for (int j = 1; j < world_size; j++)
			{
				MPI_Send(&image_array[i][j * local_arr_size], local_arr_size, MPI_INT, j, 0, MPI_COMM_WORLD);
			}

			//if processor 0 took all the first parts from all the images start summing them and divide it by the number of images
			int sum = 0;
			for (int x = 0; x < local_arr_size; x++) {
				for (int k = 0; k < imagenumber; k++) {

					sum += image_array[k][x];

				}
				sum /= imagenumber;
				Resultbackgroundimage[x] = sum;
				int diff = abs(Resultbackgroundimage[x] - threshold_img[x]);
				if (diff > threshold)
				{
					Resultforeground[x] = diff;
				}

				else
				{
					Resultforeground[x] = 0;
				}


			}
			//if its the last iteration processor 0 will start collecting all the parts of the image from other processors and create it
			if (i == imagenumber - 1) {
				/*
				for (int j = 1; j < world_size; j++)
				{
					MPI_Send(&threshold_img[j * local_arr_size], local_arr_size, MPI_INT, j, 0, MPI_COMM_WORLD);

				}
				*/
				for (int k = 1; k < world_size; k++)
				{
					MPI_Status status;
					MPI_Send(&threshold_img[k * local_arr_size], local_arr_size, MPI_INT, k, 0, MPI_COMM_WORLD);
					MPI_Recv(&Resultbackgroundimage[k * local_arr_size], local_arr_size, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&Resultforeground[k * local_arr_size], local_arr_size, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
				}



				createImage(Resultbackgroundimage, w, h, 0);
				createImage(Resultforeground, w, h, 1);

			}

		}
		else if (world_rank <= world_size)
		{

			//each processor recieve part of the image and save it in local image array
			local_image_array[i] = new int[local_arr_size];

			MPI_Status status;
			MPI_Recv(local_image_array[i], local_arr_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);


			//if it recived all the parts of the image it will start working on them
			if (i == imagenumber - 1)
			{
				int* local_threshold = new int[local_arr_size];
				MPI_Recv(local_threshold, local_arr_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

				int sum = 0;
				int* localbackgroundimage = new int[local_arr_size];
				int* local_forground = new int[local_arr_size];

				for (int x = 0; x < local_arr_size; x++) {
					for (int k = 0; k < imagenumber; k++) {

						sum += local_image_array[k][x];

					}
					sum /= imagenumber;
					localbackgroundimage[x] = sum;
					int diff = abs(localbackgroundimage[x] - local_threshold[x]);
					if (diff > threshold)
					{
						local_forground[x] = diff;

					}
					else
					{
						local_forground[x] = 0;
					}


				}
				MPI_Send(localbackgroundimage, local_arr_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(local_forground, local_arr_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}

	}


	MPI_Finalize();
	system("pause");
	return 0;

}