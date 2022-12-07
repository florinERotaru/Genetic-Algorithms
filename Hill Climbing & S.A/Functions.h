#pragma once

#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

float Rastrigin(const std::vector<float> &args)
{
	int dimension = args.size();
	float sum = 0;
	for (auto x : args)
	{
		sum += ( x*x - 10 * cos(2 * M_PI * x));
	}
	return (10 * dimension + sum);
}
float DeJong(const std::vector<float>& args)
{
	float sum = 0;
	for (auto x : args)
	{
		sum += x * x;
	}
	return sum;
}

float Schwefel(const std::vector<float>& args)
{
	float sum = 0;
	for (auto x : args)
	{
		sum += -1 * x * sin(sqrt(abs(x)));
	}
	return sum;
}

float  Michalewicz(const std::vector<float>& args)
{
	float sum = 0;
	int i = 0;
	for (auto x : args)
	{
		sum += sin(x) * pow((sin((i * x * x) / M_PI)), 20);
		i++;
	}
	return -1 * sum;
}