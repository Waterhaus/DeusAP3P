#include "stdafx.h"
#include "AP3P.h"
#include <iostream>



 vector<double> vec3_product(vector<double> a,vector<double> b) 
{
	 vector<double> c(3);
	 c[0] = a[1] * b[2] - a[2] * b[1];
	 c[1] = -(a[0] * b[2] - a[2] * b[0]);
	 c[2] = a[0] * b[1] - a[1] * b[0];

	 return c;
}

 double vec_scal(vector<double> a, vector<double> b)
 {
	 double S = 0.0;
	 for (size_t i = 0; i < a.size(); i++)
	 {
		 S = S + a[i] * b[i];
	 }

	 return S;
 }

	 

 double vec3_norm(vector<double> a) 
 {
	 return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
 }

 vector<double> vec3_mul( double s, vector<double> a) {
	 
	 vector<double> b(3);
	 
	 for (size_t i = 0; i < b.size(); i++)
	 {
		 b[i] = a[i] * s;
	 }

	 return b;
 }


 vector<double> vec3_sub(vector<double> a, vector<double> b) {

	 vector<double> c(3);

	 for (size_t i = 0; i < c.size(); i++)
	 {
		 c[i] = a[i] - b[i];
	 }

	 return c;
 }


 vector<double> vec3_div(double s, vector<double> a) {

	 return vec3_mul(1.0/s,a);
 }

 vector<double> mat3_mul_vec3(vector<vector<double>> A, vector<double> x)
 {
	 int n = 3;
	 vector<double> b(n);
	 for (size_t i = 0; i < n; i++)
	 {
		 for (size_t j = 0; j < n; j++)
		 {
			 b[i] += A[i][j] * x[j];
		 }
	 }
	 return b;
 }

 void solveQuartic(const double *factors, vector<double> *realRoots) {
	 const double &a4 = factors[0];
	 const double &a3 = factors[1];
	 const double &a2 = factors[2];
	 const double &a1 = factors[3];
	 const double &a0 = factors[4];

	 double a4_2 = a4 * a4;
	 double a3_2 = a3 * a3;
	 double a4_3 = a4_2 * a4;
	 double a2a4 = a2 * a4;

	 double p4 = (8 * a2a4 - 3 * a3_2) / (8 * a4_2);
	 double q4 = (a3_2 * a3 - 4 * a2a4 * a3 + 8 * a1 * a4_2) / (8 * a4_3);
	 double r4 = (256 * a0 * a4_3 - 3 * (a3_2 * a3_2) - 64 * a1 * a3 * a4_2 + 16 * a2a4 * a3_2) / (256 * (a4_3 * a4));

	 double p3 = ((p4 * p4) / 12 + r4) / 3; // /=-3
	 double q3 = (72 * r4 * p4 - 2 * p4 * p4 * p4 - 27 * q4 * q4) / 432; // /=2

	 double t; // *=2
	 complex<double> w;
	 if (q3 >= 0)
		 w = -sqrt(static_cast<complex<double> >(q3 * q3 - p3 * p3 * p3)) - q3;
	 else
		 w = sqrt(static_cast<complex<double> >(q3 * q3 - p3 * p3 * p3)) - q3;
	 if (w.imag() == 0.0) {
		 w.real(cbrt(w.real()));
		 t = 2.0 * (w.real() + p3 / w.real());
	 }
	 else {
		 w = pow(w, 1.0 / 3);
		 t = 4.0 * w.real();
	 }

	 complex<double> sqrt_2m = sqrt(static_cast<complex<double> >(-2 * p4 / 3 + t));
	 double B_4A = -a3 / (4 * a4);
	 double complex1 = 4 * p4 / 3 + t;

	 complex<double> complex2 = 2 * q4 / sqrt_2m;

	 double sqrt_2m_rh = sqrt_2m.real() / 2;
	 double sqrt1 = sqrt(-(complex1 + complex2)).real() / 2;
	 (*realRoots)[0] = B_4A + sqrt_2m_rh + sqrt1;
	 (*realRoots)[1] = B_4A + sqrt_2m_rh - sqrt1;
	 double sqrt2 = sqrt(-(complex1 - complex2)).real() / 2;
	 (*realRoots)[2] = B_4A - sqrt_2m_rh + sqrt2;
	 (*realRoots)[3] = B_4A - sqrt_2m_rh - sqrt2;
 }

 void polishQuarticRoots(double  *coeffs, vector < double > *roots) {
	 const int iterations = 2;
	 for (int i = 0; i < iterations; ++i) {
		 for (int j = 0; j < 4; ++j) {
			 double error =
				 (((coeffs[0] * (*roots)[j] + coeffs[1]) * (*roots)[j] + coeffs[2]) * (*roots)[j] + coeffs[3]) * (*roots)[j] +
				 coeffs[4];
			 double
				 derivative =
				 ((4 * coeffs[0] * (*roots)[j] + 3 * coeffs[1]) * (*roots)[j] + 2 * coeffs[2]) * (*roots)[j] + coeffs[3];
			 (*roots)[j] -= error / derivative;
		 }
	 }
 }

 double sign(double x) 
 {
	 if (x >= 0) return 1.0; else return -1.0;
 }

 void print_vec(string str,vector<double> vec)
 {
	 cout << endl << str << " ";
	 for (size_t i = 0; i < vec.size(); i++)
	 {
		 cout << vec[i] << " ";
	 }
	 cout << endl;
 }

 void print_mat(string str, vector<vector<double>> mat)
 {
	 cout << endl << str << " ";
	 for (size_t i = 0; i < mat.size(); i++)
	 {
		 for (size_t j = 0; j < mat.size(); j++)
		 {
			 cout << mat[i][j] << " ";
		 }
		 cout << endl;
	 }
	 cout << endl;
 }

 vector<vector<double>> mat3_mul(vector<vector<double>> A, vector<vector<double>> B)
 {
	 vector<vector<double>> C(3);
	 int n = 3;
	 for (size_t i = 0; i < 3; i++)
	 {
		 C[i].resize(3);
	 }

	 for (int i = 0; i < n; i++)
		 for (int j = 0; j < n; j++)
			 for (int r = 0; r < n; r++)
				 C[i][j] += A[i][ r] * B[r][j];

	 return C;
 }

 vector<vector<double>> transpose_mat3(vector<vector<double>> A) 
 {
	 vector<vector<double>> B(A.size());

	 for (size_t i = 0; i < A.size(); i++)
	 {
		 B[i].resize(A.size());
	 }

	 for (size_t i = 0; i < A.size(); i++)
	 {
		 for (size_t j = 0; j < A.size(); j++)
		 {
			 B[i][j] = A[j][i];
		 }
	 }

	 return B;
 }

 vector<vector<double>> GetBFrom2d( vector<double> X0, vector<double> u0,
		    vector<double> X1, vector<double> u1,
			vector<double> X2, vector<double> u2) 
 {
	 double mk0, mk1, mk2;
	 double norm;
	 double inv_fx = 1;
	 double cx_fx = 0;
	 double inv_fy = 1;
	 double cy_fy = 0;

	 u0[0] = inv_fx * u0[0] - cx_fx;
	 u0[1] = inv_fy * u0[1] - cy_fy;
	 norm = sqrt(u0[0] * u0[0] + u0[1] * u0[1] + 1);
	 mk0 = 1. / norm;
	 u0[0] *= mk0;
	 u0[1] *= mk0;

	 u1[0] = inv_fx * u1[0] - cx_fx;
	 u1[1] = inv_fy * u1[1] - cy_fy;
	 norm = sqrt(u1[0] * u1[0] + u1[1] * u1[1] + 1);
	 mk1 = 1. / norm;
	 u1[0] *= mk1;
	 u1[1] *= mk1;

	 u2[0] = inv_fx * u2[0] - cx_fx;
	 u2[1] = inv_fy * u2[1] - cy_fy;
	 norm = sqrt(u2[0] * u2[0] + u2[1] * u2[1] + 1);
	 mk2 = 1. / norm;
	 u2[0] *= mk2;
	 u2[1] *= mk2;

	

	 vector<vector<double>> featureVec =	{    { u0[0],u0[1],mk0 },
													 { u1[0],u1[1],mk1 },
													 { u2[0],u2[1],mk2 } };
	 
	 

	 return featureVec;
 }


 int AP3P_method( vector<vector<double>> B, 
					vector<vector<double>> P,
					vector<vector<vector<double>>>  *ANSER_R,
					vector<vector<double>>		    *ANSER_T)
 {
	 //1) k1,k2,k3,b
	 //k1
	 vector<double> k1;
	 vector<double> temp = vec3_sub(P[0], P[1]);
	 double norm = vec3_norm(temp);

	 k1 = vec3_div(norm, temp);
	 
	 //k3
	 vector<double> k3;
	 temp = vec3_product(B[0], B[1]);
	 norm = vec3_norm(temp);

	 k3 = vec3_div(norm, temp);
	 
	 //k2
	 vector<double> k2;
	 temp = vec3_product(k1,k3);
	 norm = vec3_norm(temp);

	 k2 = vec3_div(norm, temp);
	 
	 //2) считаем ui, vi; i=1,2
	 vector<double> u1 = vec3_sub(P[0],P[2]);
	 vector<double> u2 = vec3_sub(P[1], P[2]);
	 

	 vector<double> v1 = vec3_product(B[0], B[2]);
	 vector<double> v2 = vec3_product(B[1], B[2]);
	 
	 
	 //3) вычисл€ем дельта и K3

	 temp = vec3_product(u1, k1);
	 double delta = vec3_norm( temp );
	 vector<double> K3 = vec3_div(delta, temp);

	 //4) вычисл€ем fij; i=1,2; j = 0..4;
	 vector<double> f1(5);
	 vector<double> f2(5);

	 double k3b3 = vec_scal(k3, B[2]);
	 double norm_b1Mb2 = vec3_norm( vec3_product(B[0],B[1]) );
	 double u2k1 = vec_scal(u2, k1);
	 double b1b2 = vec_scal(B[0], B[1]);

	 f1[0] = delta*k3b3;
	 f2[0] = delta*b1b2*k3b3;

	 f2[1] = delta*k3b3*norm_b1Mb2;
	 f1[2] = delta*vec_scal(v1, k3);

	 f2[2] = delta*vec_scal(v2, k3);
	 f2[3] = u2k1*k3b3*norm_b1Mb2;

	 f1[4] = -vec_scal(u1, k1)*k3b3;
	 f2[4] = -u2k1*b1b2*k3b3;
	 
	 //5) alphai; i=0..4
	 double g1 = f1[2] * f2[1];
	 double g2 = f1[2] * f2[4] - f1[4] * f2[2];
	 double g3 = f1[0] * f2[2] - f1[2] * f2[0];
	 double g4 = -f1[2] * f2[3];
	 double g5 = f1[0] * f2[1];
	 double g6 = f1[0] * f2[4] - f1[4] * f2[0];
	 double g7 = -f1[4] * f2[3];

	 double alpha[5] = { g5 * g5 + g1 * g1 + g3 * g3,
		 2 * (g5 * g6 + g1 * g2 + g3 * g4),
		 g6 * g6 + 2 * g5 * g7 + g2 * g2 + g4 * g4 - g1 * g1 - g3 * g3,
		 2 * (g6 * g7 - g1 * g2 - g3 * g4),
		 g7 * g7 - g2 * g2 - g4 * g4 };

	 //6) найдем корни многочлена
	 vector<double> roots(4);
	 solveQuartic(alpha, &roots);
	 polishQuarticRoots(alpha, &roots);

	 double cos_th = 0;
	 double sin_th = 0;
	 //7) for
	 //посчитаем некоторые матрицы заранее 

	 temp = vec3_product(k1,K3);

	 vector<vector<double>> Ck1K3 =
	 { { k1[0], K3[0], temp[0] },
	 { k1[1], K3[1], temp[1] },
	 { k1[2], K3[2], temp[2] } };


	 temp = vec3_product(B[0], k3);

	 vector<vector<double>> Cb1k3 =
	 { { B[0][0], B[0][1], B[0][2] },
	 { k3[0], k3[1], k3[2] },
	 { temp[0], temp[1], temp[2] } };

	 
	 

	 for (size_t i = 0; i < 4; i++)
	 {
		 //8) считаем синус
		 cos_th = roots[i];
		 if (abs(cos_th) > 1.0)
			 continue;
		 sin_th = sign(k3b3)*sqrt(1.0 - cos_th*cos_th);
		 //9) найдем cos_th3 и sin_th3

		 double cos_th3 = g1 * cos_th + g2;
		 double sin_th3 = g3 * cos_th + g4;
		 double coef = sin_th / ((g5 * cos_th + g6) * cos_th + g7);

		 cos_th3 *= coef;
		 sin_th3 *= coef;

		 //10) —читаем матрицу поворота 

		 

		 vector<vector<double>>  C_cs =
		 { { cos_th3,            0,         -sin_th3 },
		 { sin_th * sin_th3, cos_th,  sin_th * cos_th3 },
		 { cos_th * sin_th3, -sin_th, cos_th * cos_th3 } };

		 

		 vector<vector<double>> temp_matrix;
		 vector<vector<double>> R;

		 temp_matrix = mat3_mul(Ck1K3, C_cs);
		 R =  mat3_mul(temp_matrix, Cb1k3);

		
		 R = transpose_mat3(R);
		 
		 (*ANSER_R).push_back(R);
		 //11) Ќайдем вектор перехода T

		 //vector<double> rb3 = mat3_mul_vec3(R, B[2]);

		 vector<double> b3p = vec3_mul((delta / k3b3), B[2]);

		 temp = vec3_mul(sin_th, b3p);
		 

		 vector<double> T = vec3_sub(temp, mat3_mul_vec3(R,P[2]));
		
		 (*ANSER_T).push_back(T);
	 }

	 return (*ANSER_T).size();
 }



 void CalculateAP3P(vector<vector<double>> points3d,
					vector<vector<double>> points2d,
					vector<vector<vector<double>>>  *ANSER_R,
					vector<vector<double>>		    *ANSER_T)
	 
 {
	 vector<vector<double>> B = GetBFrom2d(points3d[0],points2d[0], points3d[1], points2d[1],points3d[2], points2d[2]);
	
	 AP3P_method(B, points3d, ANSER_R, ANSER_T);
 }
