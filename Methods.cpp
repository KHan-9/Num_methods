#include"Methods.h"


void read(point*& arr, int& n, long double& x0, std::string file_name) {
	std::ifstream input(file_name);
	input >> n;
	input >> x0;
	arr = new point[n];
	for (point* i = arr; i < arr + n; i++) {
		input >> i->x;
		input >> i->y;
	}
	input.close();
}

double Lagrange(point arr[], int n, long double x0) {
	long double tmp, sum = 0;
	for (int j = 0; j < n; j++) {
		tmp = 1;
		for (int i = 0; i < n; i++) {
			if (i == j)continue;
			tmp *= (x0 - arr[i].x) / (arr[j].x - arr[i].x);
		}
		sum += tmp * arr[j].y;
	}
	return sum;
}

double newton(point* arr, int n, double x0) {
	double** tri_matrix = nullptr;
	double res = arr[0].y;
	tri_matrix = new double* [n - 1];
	for (int i = 0; i < n - 1; i++) {
		tri_matrix[i] = new double[n - 1 - i];
	}

	for (int i = 0; i < n - 1; i++) {
		tri_matrix[0][i] = (arr[i + 1].y - arr[i].y) / (arr[i + 1].x - arr[i].x);
	}

	for (int i = 1; i < n - 1; ++i) {
		for (int j = 0; j < n - 1 - i; j++) {
			tri_matrix[i][j] = (tri_matrix[i - 1][j + 1] - tri_matrix[i - 1][j]) / (arr[j + i + 1].x - arr[j].x);
		}
	}

	double multip = 1;
	for (int i = 0; i < n - 1; i++) {

		multip *= (x0 - arr[i].x);
		res += (tri_matrix[i][0] * multip);

	}

	return res;
}

double rectangle(double a, double b, int n, double(*f)(double)) {
	double dx = (b - a) / n;
	double integral = 0;
	for (double i = a; a < b; a += dx) {
		integral += f(a + dx) * dx;
	}

	return integral;
}

double trapezoid(double a, double b, int n, double(*f)(double)) {
	double dx = (b - a) / n;
	double integral = 0;
	for (double i = a; a < b; a += dx) {
		integral += (f(a) + f(a + dx));
	}
	return integral * 0.5 * dx;
}

double Simpson(double a, double b, int n, double(*f)(double)) {
	n += (n % 2);
	double h = (b - a) / n;
	double integral = 0;
	for (int i = 0; i <= n - 2; i += 2) {
		integral += (f(a + (i * h)) + 4 * f(a + ((i + 1) * h)) + f(a + ((i + 2) * h))) * h / 3;
	}
	return integral;
}

double Monte_Carlo(double a, double b, int n, double(*f)(double)) {

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> distribution(a, b);
	double* arr = new double[n];
	for (int i = 0; i < n; i++) {
		arr[i] = distribution(gen);
	}
	double integral = 0;

	for (int i = 0; i < n; i++) {
		integral += f(arr[i]);
	}
	integral /= n;
	integral *= fabs(b - a);

	return integral;
}

double bisec(double(*f)(double), double a, double b, double prec) {

	double lim = (a + b) / 2;
	if (fabs(f(lim)) < prec)return lim;
	else {
		if (f(lim) * f(a) < 0)b = lim;
		if (f(lim) * f(b) < 0)a = lim;

		return bisec(f, a, b, prec);
	}

}

double NR(double(*f)(double), double(*deriv_f)(double), double a, double b, double prec) {
	if (fabs(f(b)) < prec)return b;
	else {
		b = b - (f(b) / deriv_f(b));
		return NR(f, deriv_f, a, b, prec);
	}

}

int euler(double*& res, double x0, double y0, double h, double stop, double(*ptr)(double x, double y)) {
	int n = (stop - x0) / h;
	res = new double[n + 1];
	res[0] = y0;
	for (int i = 1; i <= n; i++) {
		res[i] = res[i - 1] + h * (*ptr)(x0 + ((i - 1) * h), res[i - 1]);
	}
	return n + 1;
}

int RK4(double*& res, double x0, double y0, double h, double stop, double(*ptr)(double x, double y)) {
	int n = (stop - x0) / h;
	double k[4];
	res = new double[n + 1];
	res[0] = y0;
	for (int i = 0; i <= n;) {
		k[0] = h * (*ptr)(x0 + (i * h), res[i]);
		k[1] = h * (*ptr)(x0 + (i * h) + 0.5 * h, res[i] + 0.5 * k[0]);
		k[2] = h * (*ptr)(0 + (i * h) + 0.5 * h, res[i] + 0.5 * k[1]);
		k[3] = h * (*ptr)(x0 + (i * h) + h, res[i] + k[2]);
		res[++i] = res[i - 1] + (k[0] + 2 * k[1] + 2 * k[2] + k[3]) / 6;
	}


	return n + 1;
}

int Huen(double*& res, double x0, double y0, double h, double stop, double(*ptr)(double x, double y)) {
	int n = (stop - x0) / h;
	res = new double[n + 1];
	res[0] = y0;
	for (int i = 1; i <= n; i++) {
		res[i] = res[i - 1] + (h / 2) * ((*ptr)(x0 + ((i - 1) * h), res[i - 1]) + (*ptr)(x0 + ((i - 1) * h) + h, res[i - 1] + h * (*ptr)(x0 + ((i - 1)) * h, res[i - 1])));
	}


	return n + 1;
}

double gold_prop(double(*f)(double), double a, double b, double prec) {

	if (fabs(b - a) < prec)return (a + b) / 2;
	double l, p;
	l = b - fi * (b - a);
	p = a + fi * (b - a);

	if (f(l) < f(p))b = p;
	if (f(l) > f(p))a = l;

	return gold_prop(f, a, b, prec);

}