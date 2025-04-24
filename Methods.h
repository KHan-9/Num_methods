#pragma once
#include<iostream>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<random>
#define fi 0.61803398
#define E 2.71828182845904523536;

struct point {
	long double x;
	long double y;
	point() : x(0), y(0) {};
	point(long double x_v, long double y_v) : x(x_v), y(y_v) {}
};

void read(point*& arr, int& n, long double& x0, std::string file_name);	//Function for reading point coordinates from txt file

double Lagrange(point arr[], int n, long double x0); //Lagrange interpolation

double newton(point* arr, int n, double x0); //Newton interpolation

//---Numerical integration methods

double rectangle(double a, double b, int n, double(*f)(double));

double trapezoid(double a, double b, int n, double(*f)(double));

double Simpson(double a, double b, int n, double(*f)(double));

double Monte_Carlo(double a, double b, int n, double(*f)(double));

//-----------


//---Non-linear equation solving

double bisec(double(*f)(double), double a, double b, double prec);

double NR(double(*f)(double), double(*deriv_f)(double), double a, double b, double prec);	//Newtona-Raphson method

//------------

//---Differencial equation solving

int euler(double*& res, double x0, double y0, double h, double stop, double(*ptr)(double x, double y));

int RK4(double*& res, double x0, double y0, double h, double stop, double(*ptr)(double x, double y));

int Huen(double*& res, double x0, double y0, double h, double stop, double(*ptr)(double x, double y));

//---------


double gold_prop(double(*f)(double), double a, double b, double prec); //golden ratio optimization method