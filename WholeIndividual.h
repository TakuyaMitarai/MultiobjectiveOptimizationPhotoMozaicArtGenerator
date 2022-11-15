#pragma once

#include "PartialPopulation.h"

class WholeIndividual
{
public:
	static PartialPopulation** ppop;

	WholeIndividual();
	~WholeIndividual();
	void newGeneration();
	void newGeneration(WholeIndividual* p);
	void newGeneration(WholeIndividual* p1, WholeIndividual* p2, int index1, int index2);
	void evaluation();
	void evaluation1();
	void evaluation2();
	int used_photo_cnt[PEOPLE_NUM];
	int photo_sub[ArtPhotoNum];

	PartialIndividual* chrom[WCHROM_LEN];	// ���F��
	double fitness;							// �K���x
	double fitness1, fitness2;

private:
	void mutate();
};
