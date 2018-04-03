// DeusTest.cpp: определяет точку входа для консольного приложения.
//
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
#include "Test.h";

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



	 std::vector<cv::Point3d> list_points3d;
	 std::vector<cv::Point2d> list_points2d;
	 int flags = cv::SOLVEPNP_AP3P;
	 

	 list_points3d.push_back(cv::Point3d(0,100,0));
	 list_points3d.push_back(cv::Point3d(0, 100, 100));
	 list_points3d.push_back(cv::Point3d(0, 0,  0));

	 list_points2d.push_back(cv::Point2d(10, 10));
	 list_points2d.push_back(cv::Point2d(10, 100));
	 list_points2d.push_back(cv::Point2d(100, 10));


	 std::vector<vector<double>> vec_points3d = ToVec(list_points3d);
	 std::vector<vector<double>> vec_points2d = ToVec(list_points2d);


	 ap3p method(1,1,0,0);
	 double RR[4][3][3];
	 double TT[4][3];
	 method.solve(RR, TT, list_points2d[0].x, list_points2d[0].y,
		 list_points3d[0].x, list_points3d[0].y, list_points3d[0].z,
		 list_points2d[1].x, list_points2d[1].y,
		 list_points3d[1].x, list_points3d[1].y, list_points3d[1].z,
		 list_points2d[2].x, list_points2d[2].y,
		 list_points3d[2].x, list_points3d[2].y, list_points3d[2].z);

	 vector<vector<vector<double>>> R_my;
	 vector<vector<double>> T_my;

	 CalculateAP3P(vec_points3d, vec_points2d,&R_my,&T_my);
	 print_vec(T_my[0]);




cv::Mat distCoeffs = cv::Mat::zeros(4, 1, CV_64FC1);
std::vector<cv::Mat> rvec;
std::vector<cv::Mat> tvec;// = cv::Mat::zeros(3, 1, CV_64FC1);

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
bool useExtrinsicGuess = false;

// Pose estimation
bool correspondence = cv::solveP3P(list_points3d, list_points2d, _A_matrix, distCoeffs, rvec, tvec, flags);
	
// Transforms Rotation Vector to Matrix

_t_matrix = tvec[0];

std::cout << "sovled = " << correspondence << std::endl;
std::cout << "number of solutions = " << tvec.size() <<std::endl;

for (size_t j = 0; j < tvec.size(); j++)
{
	cout << endl << j << endl;

	for (size_t i = 0; i < 3; i++)
	{
		std::cout << tvec[j].at<double>(i, 0) << "  ";
	}
	
}

cout << "rotation matrix" << endl;
for (size_t j = 0; j < tvec.size(); j++)
{
	cout << endl << j << endl;	
		Rodrigues(rvec[j], _R_matrix);
		print_mat(_R_matrix, 3, 3);
	

}

char a;
std::cin >> a;
    return 0;
}

