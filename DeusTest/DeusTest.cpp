#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include "AP3P.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>


using namespace std;

void print_mat(cv::Mat A, int n, int m)
{
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			cout << A.at<double>(i, j) << " ";
		}
		cout << endl;
	}
}

void print_vec(vector<double> vec)
{
	cout << endl << "vec = ";
	for (size_t j = 0; j < vec.size(); j++)
	{
		cout << vec[j] << " ";
	}
	cout << endl;
}

vector<vector<double>> ToVec(vector<cv::Point3d> point)
{
	vector<vector<double>> vec;

	

	for (size_t i = 0; i < point.size(); i++)
	{
		vec.push_back({ point[i].x,point[i].y,point[i].z });
	}

	return vec;
}

vector<vector<double>> ToVec(vector<cv::Point2d> point)
{
	vector<vector<double>> vec;



	for (size_t i = 0; i < point.size(); i++)
	{
		vec.push_back({ point[i].x,point[i].y });
	}

	return vec;
}

int main()
{
	//подготовим входные данные.
	 std::vector<cv::Point3d> list_points3d;
	 std::vector<cv::Point2d> list_points2d;
	 int flags = cv::SOLVEPNP_AP3P;
	 
	 //координаты наших 3 точек в локальной системе координат О
	 list_points3d.push_back(cv::Point3d(0,100,0));
	 list_points3d.push_back(cv::Point3d(0, 100, 100));
	 list_points3d.push_back(cv::Point3d(0, 0,  0));

	 //координаты трех точек на "фотографии"
	 list_points2d.push_back(cv::Point2d(10, 10));
	 list_points2d.push_back(cv::Point2d(10, 150));
	 list_points2d.push_back(cv::Point2d(100, 10));

	 //преобразуем входные данные для использования в методе CalculateAP3P
	 std::vector<vector<double>> vec_points3d = ToVec(list_points3d);
	 std::vector<vector<double>> vec_points2d = ToVec(list_points2d);

	 //здесь будет хранится результат CalculateAP3P. То есть матрица поворота R и вектор трансляции(сдвига) Т
	 vector<vector<vector<double>>> MY_ANSER_R;
	 vector<vector<double>>		    MY_ANSER_T;
	 //Запуск метода AP3P
	 CalculateAP3P(vec_points3d, vec_points2d,	&MY_ANSER_R ,&MY_ANSER_T);


//настроим окончательно opencv параметры
cv::Mat distCoeffs = cv::Mat::zeros(4, 1, CV_64FC1);
std::vector<cv::Mat> rvec;
std::vector<cv::Mat> tvec;

// настройка параметров камеры. моя реализация не учитывает ее. 
// поэтому оставим данную матрицу единичной 
cv::Mat _A_matrix = cv::Mat::zeros(3, 3, CV_64FC1);   // intrinsic camera parameters
_A_matrix.at<double>(0, 0) = 1;       //      [ fx   0  cx ]
_A_matrix.at<double>(1, 1) = 1;       //      [  0  fy  cy ]
_A_matrix.at<double>(0, 2) = 0;       //      [  0   0   1 ]
_A_matrix.at<double>(1, 2) = 0;
_A_matrix.at<double>(2, 2) = 1;
/** The computed rotation matrix */
cv::Mat _R_matrix = cv::Mat::zeros(3, 3, CV_64FC1);
/** The computed translation matrix */
cv::Mat _t_matrix = cv::Mat::zeros(3, 1, CV_64FC1);

// Запуск метода opencv
bool correspondence = cv::solveP3P(list_points3d, list_points2d, _A_matrix, distCoeffs, rvec, tvec, flags);
	
// Transforms Rotation Vector to Matrix


std::cout << "number of solutions = " << tvec.size() <<std::endl;

cout << "4 Translation vectors:" << endl;
for (size_t j = 0; j < tvec.size(); j++)
{
	cout << endl << j << ") " << endl;
	cout << "my solution-> ";
	for (size_t i = 0; i < 3; i++)
	{
		std::cout << MY_ANSER_T[j][i] << "  ";
	}
	cout << "OpenCV solution-> ";
	for (size_t i = 0; i < 3; i++)
	{
		std::cout << tvec[j].at<double>(i, 0) << "  ";
	}
	
}

cout << endl  << endl << "4 Rotation matrixes:" << endl;
for (size_t j = 0; j < tvec.size(); j++)
{
	cout << endl << j << ") " << endl;	
		
	cout << "my solution-> " << endl;
	print_mat("", MY_ANSER_R[j]);
	cout << "OpenCV solution-> " << endl;
	Rodrigues(rvec[j], _R_matrix);
	print_mat(_R_matrix, 3, 3);	

}

char a;
std::cin >> a;
    return 0;
}

