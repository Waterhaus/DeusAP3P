#pragma once
#include "stdafx.h"
#include <vector>
#include <complex>

using namespace std;

void CalculateAP3P(vector<vector<double>> points3d,
	vector<vector<double>> points2d,
	vector<vector<vector<double>>>  *ANSER_R,
	vector<vector<double>>		    *ANSER_T);