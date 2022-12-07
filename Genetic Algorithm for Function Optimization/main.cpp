#include "Problem.h"
#include "Functions.h"
#include <Windows.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>



void KeepMonitorActive() {
	// Enable away mode and prevent the sleep idle time-out.
	SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED | ES_AWAYMODE_REQUIRED);
}

void RestoreMonitorSettings() {
	// Clear EXECUTION_STATE flags to disable away mode and allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);
}
int main()
{
	KeepMonitorActive();
	atexit(RestoreMonitorSettings);

	
	Problem* dj1 = new Problem(Rastrigin, -5.12, 5.12, 5);
	for(int i= 0; i<10; i++)
	std::cout << dj1->GeneticAlgorithm() << endl;

	
	return 0;
}