#include <iostream>
#include "Problem.h"
#include <Windows.h>
#include <cstdlib>
//debug function
void print(std::vector<float> x)
{
	for (auto el : x)
	{
		std::cout << el << "; ";
	}
	std::cout << std::endl;
}

void KeepMonitorActive() {
	// Enable away mode and prevent the sleep idle time-out.
	SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED | ES_AWAYMODE_REQUIRED);
}

void RestoreMonitorSettings() {
	// Clear EXECUTION_STATE flags to disable away mode and allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);
}


int main()

{	//0 - first
	//1 - best
	//2 worst
	KeepMonitorActive();
	atexit(RestoreMonitorSettings);




	std::ofstream myfile4;
	myfile4.open("D:\\uni2.0\\AlgoritmiGenetici\\proiect\\MIHA.txt", std::ios::app);
	myfile4 << std::endl;
	Problem* dj0 = new Problem(0, M_PI, 5, 5);
	dj0->SetFunction(Michalewicz);
	for (int i = 0; i < 20; i++)
	{
		clock_t tStart = clock();
		dj0->IterateHillClimber(10000, 1);
		myfile4 << " BEST IMPROVEMENT: iterated hill climber, finished with the minimum of " << (float)(dj0->final_solution.solution_value) << " after " << (double)((clock() - tStart) / CLOCKS_PER_SEC) / 60 << " min " <<
			"on dimension " << 5 << std::endl;
	}

	myfile4 << std::endl;
	myfile4 << "=============================";
	myfile4 << std::endl;




	return 0;
	 
}