#pragma once
#include <math.h>
#include <random>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "Functions.h"
#include "Solution.h"
#include <time.h>
#include <algorithm>
using namespace std;
#define precision 5
#define selection_pressure 2
#define P_c 0.6
#define P_m_0 0.01
#define population_size 100.00

random_device device;
uniform_int_distribution<int> distribution(0, 1);

random_device device3;
uniform_real_distribution<float> distribution3(0, 1.0);
typedef std::vector< boost::dynamic_bitset<> > bitset_vector;

//MIND THE DIFFERENCE BETWEEN BIT VECTOR AND BITSET VECTOR !!!! (QUADRATIC)

class Problem
{
public:
	float a;
	float b;
	Solution* global_minimum;
	boost::dynamic_bitset<> global_minimum_representation;

	int dimension;
	int nr_of_intervals;
	int one_argument_length;
	int full_length;
	float (*function)(const vector<float>& args);
	vector<Solution*> Population;

public:
	Problem(float (*function)(const vector<float>& args), float a, float b, int dimension)
		: a(a), b(b),
		dimension(dimension),
		function(function)
	{
		ComputeSpecs();
		Population.reserve(population_size);
		global_minimum = new Solution(full_length, dimension);
		for (int i = 0; i < population_size; i++)
		{
			Solution* chromosome = new Solution(full_length, dimension);
			Population.emplace_back(chromosome);
		}
	}
	boost::dynamic_bitset<> getRandomSequence()
	{
		boost::dynamic_bitset<> res(full_length);
		for (int i = 0; i < full_length; i++)
		{
			res[i] = distribution(device);
		}
		return res;
	}
	void GeneratePopulation()
	{
		for (auto& chromosome : Population)
		{
			chromosome = new Solution(full_length, dimension);
			chromosome->bits_representation = getRandomSequence();
			chromosome->solution_args = ConvertBitsetToVector(chromosome->bits_representation);
			chromosome->solution_value = function(chromosome->solution_args);
			//cout << chromosome->bits_representation << endl;
		}
	}
	void ComputeSpecs()
	{
		nr_of_intervals = (b - a) * pow(10, precision);
		one_argument_length = (int)abs(ceil(log2(nr_of_intervals)));
		full_length = one_argument_length * dimension;
		global_minimum = 0;
		cout << "Specificatiile situatiei:\n";
		cout << "numar intervale: " << nr_of_intervals << "\n";
		cout << "lungimea in biti a unei valori: " << one_argument_length << " " << pow(2, one_argument_length) << "\n";
		cout << "lungimea in biti a vectorului solutii: " << full_length << "\n";
	}
	unsigned long NumberByIndex(int index, const boost::dynamic_bitset<>& bit_vector)
	{
		if (index >= dimension || index < 0)
		{
			cout << "Invalid Index";
			exit(2);
		}
		boost::dynamic_bitset<> atomic_number(one_argument_length);
		for (int i = 0; i < one_argument_length; i++)
		{
			atomic_number[one_argument_length - 1 - i] = bit_vector[full_length - 1 - (index * one_argument_length) - i];
		}

		return atomic_number.to_ulong(); // convert to unsigned long
	}
	vector<float> ConvertBitsetToVector(const boost::dynamic_bitset<>& bit_vector)//from binary to float arguments
	{
		vector<float> float_vector;
		float argument;
		//boost::dynamic_bitset<> atomic_number(one_argument_length);
		for (int i = 0; i < dimension; i++)
		{
			argument = (float)NumberByIndex(i, bit_vector);
			argument = (argument / (pow(2, (one_argument_length)) - 1));
			argument *= (b - a);
			argument += a;
			float_vector.emplace_back(argument);
		}
		return float_vector;
	}
	void crossover_two_cromosomes(Solution*& chrm1, Solution*& chrm2)
	{
		float r, probability_choose_best = 0.7;
		boost::dynamic_bitset<> aux1(full_length);
		boost::dynamic_bitset<> aux2(full_length);

		if (chrm1->solution_value > chrm2->solution_value)
		{
			aux1 = chrm1->bits_representation;
			chrm1->bits_representation = chrm2->bits_representation;
			chrm2->bits_representation = aux1;
		}
		// full_length-1 limit for slower advance (faster convergence) - most significant bit
		for (int i = 0; i < full_length; i++)
		{
			r = distribution3(device3);
			if (r < probability_choose_best)
				aux1[i] = chrm1->bits_representation[i];
			else
				aux1[i] = chrm2->bits_representation[i];

			r = distribution3(device3);
			if (r < probability_choose_best)
				aux2[i] = chrm1->bits_representation[i];
			else
				aux2[i] = chrm2->bits_representation[i];
		}
		chrm1->bits_representation = aux1;
		chrm2->bits_representation = aux2;
		chrm1->solution_args = ConvertBitsetToVector(chrm1->bits_representation);
		chrm1->solution_value = function(chrm1->solution_args);
		chrm2->solution_args = ConvertBitsetToVector(chrm2->bits_representation);
		chrm2->solution_value = function(chrm2->solution_args);
	}
	void CrossoverPopulation()
	{
		vector<int> chosen_chromosomes;
		auto rng = std::default_random_engine{};
		int i = 0;

		for (i = 0; i < population_size; i++)
		{
			if (distribution3(device3) < P_c)
			{
				chosen_chromosomes.emplace_back(i);	// add indices of chromosomes from population for crossover
			}
		}
		for (i = 0; i <= chosen_chromosomes.size() / 2; i += 2)
		{
			//	cout << Population[i]->bits_representation << endl << Population[i + 1]->bits_representation << endl << endl;
			crossover_two_cromosomes(Population[i], Population[i + 1]);
			Population[i]->solution_args = ConvertBitsetToVector(Population[i]->bits_representation);
			Population[i]->solution_value = function(Population[i]->solution_args);
			Population[i + 1]->solution_args = ConvertBitsetToVector(Population[i + 1]->bits_representation);
			Population[i + 1]->solution_value = function(Population[i + 1]->solution_args);
			//	cout << Population[i]->bits_representation << endl << Population[i + 1]->bits_representation << endl << endl;
		}
	}
	struct AscendingOptimum
	{
		bool operator()(Solution* chrm1, Solution* chrm2)
		{
			return (chrm1->solution_value < chrm2->solution_value);
		}
	};
	void RankSelection()
	{
		sort(Population.begin(), Population.end(), AscendingOptimum());		//sort Population
		//cout << Population[0]->solution_value << endl;
		//index of chromosome i is i+1;
		// 0,..,n-1 have the indexes 1,..,n

		//compute probabilities.
		float cumulativeProbs[(int)population_size] = {};
		float p = 0;
		//cout << endl;

		for (int rank = 0; rank < population_size; rank++)
		{
			p = (1 / population_size);
			p *= (selection_pressure - (float)rank * 2 * (selection_pressure - 1) / (population_size - 1));
			//cout  << rank + 1 << "    " << p << endl;

			//cout << "for rank " << rank+1 << "with fitness: "
			//	<< Population[rank]->solution_value
			//	<< "the probability is: " << p << endl;

			if (rank == 0)
				cumulativeProbs[rank] = p;
			else
				cumulativeProbs[rank] = cumulativeProbs[rank - 1] + p;
		}
		//cout << "The cumulated probabilities look like this: \n";
		//for (float p : cumulativeProbs)
		//	cout << p << "; ";
		//let's choose some of the chrm's which we will later crossover

		int j = 0;
		float start = 0, end = 0;
		vector<Solution*> Population2;
		Population2.reserve(population_size);
		while (Population2.size() != (int)population_size)
		{
			j = 0;
			p = distribution3(device3);
			//start = cumulativeProbs[0];
			//end = cumulativeProbs[(int)population_size];

			while (cumulativeProbs[j] < p)
				j++;
			//select j chromosome of population
			if (j < population_size)
			{
				Solution* selected_chromosome = new Solution(full_length, dimension);
				selected_chromosome->bits_representation = Population[j]->bits_representation;
				selected_chromosome->solution_args = Population[j]->solution_args;
				selected_chromosome->solution_value = Population[j]->solution_value;
				Population2.emplace_back(selected_chromosome);
			}
		}
		for (auto chromosome : Population)
			delete chromosome;
		Population = Population2;
	}
	void Mutation()
	{
		sort(Population.begin(), Population.end(), AscendingOptimum());
		int N = population_size;
		//cout << "======================= " << endl;

		for (int i = 0; i < population_size; i++)
		{
			//r = (i+1) is the rank
			float p = 0;
			p = P_m_0 * (1 - (float)(100 - i - 1) / (N - 1));
			for (int j = 0; j < full_length; j++)
			{
				float r = distribution3(device3);
				if (r < p)
				{
					Population[i]->bits_representation[j] = (int)!Population[i]->bits_representation[j];
					Population[i]->solution_args = ConvertBitsetToVector(Population[i]->bits_representation);
					Population[i]->solution_value = function(Population[i]->solution_args);
				}
			}
		}
		//cout << "======================= ";

	}
	void Evaluate(int generation)
	{
		Solution* Population_minimum=new Solution(full_length, dimension);
		Population_minimum->bits_representation = Population[0]->bits_representation;
		Population_minimum->solution_args = Population[0]->solution_args;
		Population_minimum->solution_value = Population[0]->solution_value;

		if (generation == 0)
		{
			global_minimum = Population_minimum;
		}

		for (int i = 0; i < population_size; i++)
		{
			if (Population_minimum->solution_value > Population[i]->solution_value)
			{
				Population_minimum->bits_representation = Population[i]->bits_representation;
				Population_minimum->solution_args = Population[i]->solution_args;
				Population_minimum->solution_value = Population[i]->solution_value;
			}
		}
		if (global_minimum > Population_minimum)
		{
			//global_minimum = Population_minimum;
			global_minimum->bits_representation = Population_minimum->bits_representation;
			global_minimum->solution_args = Population_minimum->solution_args;
			global_minimum->solution_value = Population_minimum->solution_value;
		}
	}
	float GeneticAlgorithm()
	{
		float t = 0;
		GeneratePopulation();
		/*cout << "start chromosomes: " << endl;
		for (auto elem : Population)
			cout << elem->solution_value << endl;*/

		Evaluate(t);
		for (t; t < 1000; t++)
		{
			RankSelection();

			CrossoverPopulation();
			Mutation();
			Evaluate(t);
		}
		global_minimum->all_neighbors.reserve(full_length);

		float best = 0;

		best = global_minimum->solution_value;

		Solution* starting_point = new Solution(full_length, dimension);

		starting_point->bits_representation = global_minimum->bits_representation;
		starting_point->solution_args = global_minimum->solution_args;
		starting_point->solution_value = global_minimum->solution_value;

		for (int i = 0; i < 1000; i++)
		{
			GenerateAllNeighbors();
			SingleHillClimber();
			if (best < global_minimum->solution_value)
			{
				best = global_minimum->solution_value;

			}
			global_minimum->bits_representation = starting_point->bits_representation;
			global_minimum->solution_args = starting_point->solution_args;
			global_minimum->solution_value = starting_point->solution_value;
		}
		

		return best;
	}

	Solution* BestImprovement()
	{

		if (function == nullptr)
		{
			std::cout << "function not set!";
			exit(1);
		}
		Solution *min=new Solution(full_length, dimension);
		min = global_minimum;
		for (Solution* sol : global_minimum->all_neighbors)
		{
			if (sol->solution_value < min->solution_value)
			{
				global_minimum->bits_representation = min->bits_representation;
				global_minimum->solution_args = min->solution_args;
				global_minimum->solution_value = min->solution_value;
			}
				
		}
		return min;

	}

	void SingleHillClimber()
	{
		while (1)
		{
			Solution *sol=new Solution(full_length, dimension);
			sol = BestImprovement();
			if (sol->solution_value < global_minimum->solution_value)
			{
				global_minimum = sol;
				
				GenerateAllNeighbors(); //UPDATES THE NEIGHBORS
			}
			else
				return;
		}
	}

	void GenerateAllNeighbors()
	{
		bitset_vector neighbors;
		std::vector<Solution*> neighbors_vector;
		for (int i = 0; i < full_length; i++)
		{
			Solution *new_neighbor=new Solution(full_length, dimension);
			//negating a bit
			neighbors.push_back(global_minimum->bits_representation);
			neighbors[i][i] = ~neighbors[i][i];

			new_neighbor->bits_representation = neighbors[i];

			auto bitset_converted_to_normal_args = ConvertBitsetToVector(neighbors[i]);
			new_neighbor->solution_args = bitset_converted_to_normal_args;

			new_neighbor->solution_value = function(new_neighbor->solution_args);

			neighbors_vector.push_back(new_neighbor);
			//std::cout << " function value:" << Rastrigin(bitset_converted_to_normal_args) << std::endl; //print to debug
		}
		//std::cout << "IN TOTAL AVEM " << float_neighbors.size() << std::endl;
		global_minimum->all_neighbors = neighbors_vector;

	}


	
};

