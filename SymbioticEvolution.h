#pragma once

#include "WholePopulation.h"

class SymbioticEvolution
{
public:
	SymbioticEvolution();
	~SymbioticEvolution();
	void solve();

	double used_photo_cnt[PEOPLE_NUM];
	int best[MAX_NUM];
	double bestfit;

	WholePopulation* wpop;
	PartialPopulation* ppop[WCHROM_LEN];
};
