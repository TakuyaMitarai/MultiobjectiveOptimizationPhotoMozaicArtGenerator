#pragma once

#include "PartialIndividual.h"

class PartialPopulation
{
public:
	PartialPopulation();
	~PartialPopulation();
	void newGeneration();
	void evaluation();
	void evalinit();

	PartialIndividual* pop[PPOP_SIZE];
	PartialIndividual* distance[PPOP_SIZE];

private:
	void sort(int lb, int ub);
};
