double min(double a, double b){
	if (a < b) return a;
	else return b;
}

double Sphere(double y[], int task_id) //Sphere 
{
	double z[50], tmp[50];
	double LB = -100, UB =100;
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i];

	double s = 0;
	for (int i = 0; i < 50; i++)
		s += z[i] * z[i];
	return s;
}

double Sphere2(double y[], int task_id) //Sphere 
{
	double z[50], tmp[50];
	double LB = -100, UB = 100;
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i] - 80;

	double s = 0;
	for (int i = 0; i < 50; i++)
		s += z[i] * z[i];

	return s;
}

double Sphere3(double y[], int task_id) //Sphere 
{
	double z[50], tmp[50];
	double LB = -100, UB = 100;
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i] + 80;

	double s = 0;
	for (int i = 0; i < 50; i++)
		s += z[i] * z[i];

	return s;
}

double Rosenbrock(double y[], int task_id) //Rosenbrock
{
	double s = 0;
	double LB = -50, UB = 50;
	double z[50];
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i] + 1;
	for (int i = 0; i < 50 - 1; i++)
		s += 100 * (z[i] * z[i] - z[i + 1]) * (z[i] * z[i] - z[i + 1]) + (z[i] - 1) * (z[i] - 1);

	return s;
}


double Ackley(double y[], int task_id) //Ackley
{
	double LB = -50, UB = 50;
	double z[50], tmp[50];
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i] - 40;

	double s_a = 0, s_b = 0;
	for (int i = 0; i < 50; i++) {
		s_a += z[i] * z[i];
		s_b += cos(2.0 * PI * z[i]);
	}

	double s = 0.0;
	s = -20.0 * exp(-0.2 * sqrt(s_a / double(50))) - exp(s_b / double(50)) + 20.0 + exp(1.);

	return s;
}

double Rastrgin(double y[], int task_id) //Rastrgin
{
	double z[50], tmp[50];
	double LB = -50, UB = 50;
	for (int i = 0; i < 50; i++){
		if (i < 25) z[i] = LB + (UB - LB) * y[i] - 40;
		else z[i] = LB + (UB - LB) * y[i] + 40;
	}

	double s = 0.0;
	for (int i = 0; i < 50; i++)
		s += z[i] * z[i] - 10 * cos(2.0 * PI * z[i]) + 10;

	return s;
}

double Griewank(double y[], int task_id) //Griewank
{
	double z[50], tmp[50];
	double LB = -100, UB = 100;
	for (int i = 0; i < 50; i++){
		if (i < 25) z[i] = LB + (UB - LB) * y[i] + 80;
		else z[i] = LB + (UB - LB) * y[i] - 80;
	}

	double s_a = 0, s_b = 1.0;
	for (int i = 0; i < 50; i++) {
		s_a += z[i] * z[i];
		s_b *= cos(z[i] / sqrt(i + 1.0));
	}

	double s = 0;
	s = 1 + s_a / 4000.0 - s_b;

	return s;
}

double power(double y, int k)
{
	double s = 1;
	for (int i = 1; i <= k; i++)
		s *= y;
	return s;
}

double Weierstrass25D(double y[], int task_id) //Weierstrass with 25 dimensions
{
	
	double z[25], tmp[25];
	double LB = -0.5, UB = 0.5;
	for (int i = 0; i < 25; i++)
		z[i] = LB + (UB - LB) * y[i] + 0.4;

	double s = 0;
	int k_max = 20;
	double a = 0.5, b = 3;
	for (int i = 0; i < 25; i++)
		for (int j = 0; j <= k_max; j++)
			s += power(a, j) * cos(2 * PI * power(b, j) * (z[i] + 0.5));
	for (int i = 0; i <= k_max; i++)
		s -= double(25) * power(a, i) * cos(2 * PI * power(b, i) * 0.5);
	
	return s;
}

double Weierstrass50D(double y[], int task_id) //Weierstrass with 50 dimensions
{
	double z[50], tmp[50];
	double LB = -0.5, UB = 0.5;
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i] + 0.4;

	double s = 0;
	int k_max = 20;
	double a = 0.5, b = 3;
	for (int i = 0; i < 50; i++)
		for (int j = 0; j <= k_max; j++)
			s += power(a, j) * cos(2 * PI * power(b, j) * (z[i] + 0.5));
	for (int i = 0; i <= k_max; i++)
		s -= double(50) * power(a, i) * cos(2 * PI * power(b, i) * 0.5);

	return s;
}

double Schwefel(double y[], int task_id) //Schwefel 
{
	double s = 418.9829 * double(50);
	double LB = -500, UB = 500;
	double z[50];
	for (int i = 0; i < 50; i++)
		z[i] = LB + (UB - LB) * y[i];
	for (int i = 0; i < 50; i++)
		s -= z[i] * sin(sqrt(fabs(z[i])));

	return s;
}
