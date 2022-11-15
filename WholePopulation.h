#pragma once

#include "WholeIndividual.h"

class WholePopulation
{
public:
	WholePopulation();
	~WholePopulation();
	void newGeneration();
	void evaluation();
	void distance();

	WholeIndividual* pop[WPOP_SIZE];	// 個体集団

private:
	void changePC(int a, int b);
	void sort(int lb, int ub);
	void congestion_sort(int lb, int ub);

	WholeIndividual* childpop[CHILDPOP];	// 個体集団
	double childfit[CHILDPOP];				// 個体集団の適応度(スケーリング法により変換した後の値)
	int m;
};
