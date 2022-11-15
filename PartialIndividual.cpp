#include "PartialIndividual.h"

// �R���X�g���N�^
PartialIndividual::PartialIndividual()
{
	int i, ch;
	
	for(i = 0; i < PCHROM_LEN; i++) {
		chrom[i] = rand() % MATERIAL_PHOTO_NUM;
	}
	fitness = DBL_MAX;
}

// �f�X�g���N�^
PartialIndividual::~PartialIndividual()
{

}

// ����??���i??R?????j
void PartialIndividual::newGeneration()
{
	mutate();
	fitness = DBL_MAX;
}

// ����??���i�R�s�[��??R???j
// p: �R�s�[��?????|�C���^
void PartialIndividual::newGeneration(PartialIndividual* p)
{
	int i;
	
	for(i = 0; i < PCHROM_LEN; i++)
		chrom[i] = p->chrom[i];
	mutate();
	fitness = DBL_MAX;
}

// ����??���i�Q�_������??R???j
// p1: �e???��?P
// p2: �e???��?Q
// index1: �����_�i0�`PCHROM_LEN-1�j
// index2: �����_�i0�`PCHROM_LEN-1�j
void PartialIndividual::newGeneration(PartialIndividual* p1, PartialIndividual* p2, int index1, int index2)
{
	int i, min, max;

	if(index1 < index2) {
		min = index1;
		max = index2;
	} else {
		min = index2;
		max = index1;
	}
	for(i = 0; i < min; i++)
		chrom[i] = p1->chrom[i];
	for( ; i < max; i++){
		chrom[i] = p2->chrom[i];
	}
	for( ; i < PCHROM_LEN; i++)
		chrom[i] = p1->chrom[i];
	mutate();
	fitness = DBL_MAX;
}

// ??R???
void PartialIndividual::mutate()
{
	int i, c, r1, r2;

	for(i = 0; i < PCHROM_LEN; i++) {
		if(rand() < MUTATE_PROB)
		{
			chrom[i] = rand() % MATERIAL_PHOTO_NUM;
		}
	}
}
