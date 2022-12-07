#pragma once
#include <math.h>
#include <random>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "Functions.h"
#include "Solution.h"
#include <time.h>
#include <fstream>
typedef std::vector< boost::dynamic_bitset<> > bitset_vector;
//MIND THE DIFFERENCE BETWEEN BIT VECTOR AND BITSET VECTOR !!!!, (QUADRATIC)

class Problem
{
	double a;
	double b;
	int dimension;
	int log10eps;
	int nr_of_intervals;
	int one_argument_length;
	int full_length;
public:
	Solution current_solution;
	Solution final_solution;
	float (*function_used)(const std::vector<float>& args);

	//ctrs
	Problem(double a, double b, int dimension, int log10eps)
		: a(a), b(b),
		dimension(dimension),
		log10eps(log10eps)
	{
		ComputeSpecs();
		current_solution.bits_representation.resize(full_length);
	}
	boost::dynamic_bitset<> getRandomSequence()
	{
		std::random_device device;
		std::uniform_int_distribution<int> distribution(0, 1);
		boost::dynamic_bitset<> res(full_length);
		for (int i = 0; i < full_length; i++)
		{
			res[i] = distribution(device);
		}
		return res;
	}
	void GenerateSolution()
	{
		current_solution.bits_representation = getRandomSequence();
		current_solution.solution_args = ConvertBitsetToVector(current_solution.bits_representation);
		current_solution.solution_value = function_used(current_solution.solution_args);
	}
	//setters
	void SetFunction(float (*arg_function)(const std::vector<float>& args))
	{
		function_used = arg_function;
	}
	void SetBoundaries(int a, int b)
	{
		this->a = a;
		this->b = b;
	}
	void SetDimension(int dim)
	{
		dimension = dim;
	}
	void SetEpsExp(int epsexp)
	{
		log10eps = epsexp;
	}
	void ComputeSpecs()
	{
		nr_of_intervals = (b - a) * pow(10, log10eps);
		one_argument_length = abs(ceil(log2(nr_of_intervals)));
		full_length = one_argument_length * dimension;
		std::cout << "Specificatiile situatiei:\n";
		std::cout << "numar intervale: " << nr_of_intervals << "\n";
		std::cout << "lungimea in biti a unei valori: " << one_argument_length <<" "<<pow(2, one_argument_length)<<"\n";
		std::cout << "lungimea in biti a vectorului solutii: " << full_length << "\n";
	}
	unsigned long NumberByIndex(int index, const boost::dynamic_bitset<> &bit_vector)
	{
		if (index >= dimension || index < 0)
		{
			std::cout << "Invalid Index";
			exit(2);
		}
		boost::dynamic_bitset<> atomic_number(one_argument_length);
		for (int i = 0; i<one_argument_length; i++)
		{
			atomic_number[one_argument_length-1-i] = bit_vector[full_length - 1 - (index * one_argument_length)-i];
		}
		
		return atomic_number.to_ulong();
	}
	std::vector<float> ConvertBitsetToVector(const boost::dynamic_bitset<> &bit_vector)//from binary to ready arguments
	{
		std::vector<float> float_vector;
		boost::dynamic_bitset<> atomic_number(one_argument_length);
		for (int i = 0; i < dimension; i++)
			float_vector.push_back(NumberByIndex(i, bit_vector));
		NormalizeSolutions(float_vector);
		return float_vector;

	}
	void NormalizeSolutions(std::vector<float>& ulong_vector) //gets them to interval
	{
		for (auto& x : ulong_vector)
		{
			x = (x / (pow(2, (one_argument_length)) -1 ));
			x *= (b - a);
			x += a;
			//std::cout << x << "    ";// PRINT TO DEBUG
		}
	}
	
	void GenerateAllNeighbors()
	{
		bitset_vector neighbors;
		std::vector<Solution> neighbors_vector;
		for (int i = 0; i < full_length; i++)
		{
			Solution new_neighbor;
			//negating a bit
			neighbors.push_back(current_solution.bits_representation);
			neighbors[i][i] = ~neighbors[i][i];

			new_neighbor.bits_representation = neighbors[i];
			
			auto bitset_converted_to_normal_args = ConvertBitsetToVector(neighbors[i]);
			new_neighbor.solution_args = bitset_converted_to_normal_args;
			
			new_neighbor.solution_value = function_used(new_neighbor.solution_args);

			neighbors_vector.push_back(new_neighbor);
			//std::cout << " function value:" << Rastrigin(bitset_converted_to_normal_args) << std::endl; //print to debug
		}
		//std::cout << "IN TOTAL AVEM " << float_neighbors.size() << std::endl;
		current_solution.all_neighbors = neighbors_vector;

	}
	Solution BestImprovement()
	{

		if (function_used == nullptr)
		{
			std::cout << "function not set!";
			exit(1);
		}
		Solution min = current_solution;
		for (auto sol : current_solution.all_neighbors)
		{
			if (sol.solution_value < min.solution_value)
				min = sol;
		}
		return min;

	}
	Solution FirstImprovement()
	{
		if (function_used == nullptr)
		{
			std::cout << "function not set!";
			exit(1);
		}

		std::random_device rd;
		std::uniform_int_distribution<int> dist(0, full_length-1);
		int index;
		for (int i=0; i<full_length; i++)
		{
			index = dist(rd);
			if (current_solution.all_neighbors[index].solution_value < current_solution.solution_value)
				return  current_solution.all_neighbors[index];
		}
		return current_solution;
	}
	Solution WorstImprovement()
	{
		if (function_used == nullptr)
		{
			std::cout << "function not set!";
			exit(1);
		}

		Solution res=current_solution;
		int j = -1;
		for (int i = 0; i < current_solution.all_neighbors.size(); i++)
		{
			if (current_solution.all_neighbors.at(i).solution_value < current_solution.solution_value)
			{
				res = current_solution.all_neighbors.at(i);
				j = i;
				break;
			}
		}
		if (j == -1)
			return res;
		for (int i = j ; i < current_solution.all_neighbors.size(); i++)
		{
			if (res.solution_value < current_solution.all_neighbors.at(i).solution_value && current_solution.all_neighbors.at(i).solution_value < current_solution.solution_value)
				res = current_solution.all_neighbors.at(i);
		}
		return res;
		
	}
	

	void SingleHillClimber(short improv)
	{
		while (1)
		{
			Solution sol;
			switch (improv)
			{
				case 0:  sol = FirstImprovement();
						break;
				case 1:  sol = BestImprovement();
						break;
				case 2:  sol = WorstImprovement();
						break;
				default: return;
			}
			if (sol.solution_value < current_solution.solution_value)
			{
				current_solution = sol;
				GenerateAllNeighbors(); //UPDATES THE NEIGHBORS
			}
			else
				return;
		}
	}
	

	void IterateHillClimber(int n, short improv)
	{
		clock_t tStart = clock();
		GenerateSolution();
		GenerateAllNeighbors();
		final_solution = current_solution;
		for (int i = 0; i < n; i++)
		{	
			SingleHillClimber(improv);
			if (final_solution.solution_value > current_solution.solution_value)
				final_solution = current_solution;
			if ((double)((clock() - tStart) / CLOCKS_PER_SEC) / 60 > 200)
				return;
			GenerateSolution();
			GenerateAllNeighbors();
			
		}
		
	}

	float ComputeInitialTemperature(float hi)
	{
		if (function_used == nullptr)
		{
			std::cout << "function not set!";
			exit(1);
		}
		GenerateSolution();
		GenerateAllNeighbors();
		float dE=0;
		for (int i = 0; i < 5; i++)
		{
			
			Solution best =FirstImprovement();
			dE += abs(best.solution_value - current_solution.solution_value);
			current_solution = best;
		}
		return -1* (dE / 5) / (log(hi));
	}


	void SimulatedAnnealing()
	{
		float T = ComputeInitialTemperature(0.8);
		GenerateSolution();
		GenerateAllNeighbors();
		std::random_device rd;
		std::uniform_int_distribution<int> dist(0, full_length - 1);
		std::uniform_real_distribution<float> dist2(0, 1.0);

		int k = 1;
		do {
			float s = 0;
			int i = 0;
			int y = 0;
			int worse = 0;
			do
			{	
				//choose at random
				int index = dist(rd);
				Solution rdm_neighbor = current_solution.all_neighbors[index];
				s += exp((float)(-1 * abs(rdm_neighbor.solution_value - current_solution.solution_value) / T));
				if (rdm_neighbor.solution_value < current_solution.solution_value)
				{
					current_solution = rdm_neighbor;
					GenerateAllNeighbors();
				}
				else
				{
					
					//chosen neighbor is worse, add it?
					//worse++;
					if (dist2(rd) < exp((float)(-1 * abs(rdm_neighbor.solution_value - current_solution.solution_value) / T)))
					{
						//y++;
						//std::cout << "   " << (exp((float)(-1 * abs(rdm_neighbor.solution_value - current_solution.solution_value) / T))) << std::endl;
						//std::cout << rdm_neighbor.solution_value << ";   " << current_solution.solution_value << std::endl;
						current_solution = rdm_neighbor; //UPDATES THE NEIGHBOR
						GenerateAllNeighbors();
						
					}
				}
				
				i++;
			} while (i<1000);
			//std::cout << " temp " << T <<" nr asgnari "<<y<<" din worse "<<worse<< std::endl;
			//y = 0;
			T =T * 0.8;
			//k++;
			//std::cout << "pentru temp " << T << " prob = " << s / 1000 << std::endl;
		} while (T > 0.0005);
		final_solution = current_solution;
	}

};

