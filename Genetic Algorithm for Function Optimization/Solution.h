#pragma once
#include "Problem.h"
using namespace std;
typedef vector< vector<float> > float_matrix;
class Solution
{
public:
    boost::dynamic_bitset<> bits_representation;
    vector<float> solution_args;
    float solution_value;
    std::vector<Solution*> all_neighbors;


    Solution(float full_length, float dimension)
    {
        all_neighbors.reserve(1);
        this->bits_representation.reserve(full_length);
        this->solution_args.reserve(dimension);
        this->solution_value = 0;
    }

     Solution operator=( const Solution* rvalue)
    {
        bits_representation = rvalue->bits_representation;
        solution_args = rvalue->solution_args;
        solution_value = rvalue->solution_value;
        return *this;
    }

};

