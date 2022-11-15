#include "SymbioticEvolution.h"

// コンストラクタ
SymbioticEvolution::SymbioticEvolution()
{
	int i, j;
	// 初期集団生成
	for(i = 0; i < WCHROM_LEN; i++)
		ppop[i] = new PartialPopulation();
	WholeIndividual::ppop = ppop;
	wpop = new WholePopulation();
	wpop->evaluation();
	for(i = 0; i < WCHROM_LEN; i++)
		ppop[i]->evaluation();
	// 最適解初期化
	bestfit = wpop->pop[0]->fitness;
	for(i = 0; i < WCHROM_LEN; i++) {
		for(j = 0; j < PCHROM_LEN; j++)
			best[i*PCHROM_LEN+j] = wpop->pop[0]->chrom[i]->chrom[j];
	}
	printf("初期世代\t最良個体の適応度は%f\n", bestfit);
}

// デストラクタ
SymbioticEvolution::~SymbioticEvolution()
{
	int i;

	for(i = 0; i < WCHROM_LEN; i++)
		delete ppop[i];
	delete wpop;
}

// 進化
void SymbioticEvolution::solve(void)
{
	int gen, i, j;
	double ave = 0, diff = 0;

	for(gen = 1; gen <= GENERATION_MAX; gen++) {
		// 次世代生成
		
		for(i = 0; i < WCHROM_LEN; i++)
			ppop[i]->newGeneration();
		wpop->newGeneration();
		// 評価
		for(i = 0; i < WCHROM_LEN; i++)
			ppop[i]->evalinit();
		wpop->evaluation();
		for(i = 0; i < WCHROM_LEN; i++)
			ppop[i]->evaluation();
		wpop->distance();
		// 最適解更新
		if(wpop->pop[0]->fitness < bestfit) {
			bestfit = wpop->pop[0]->fitness;
		}
		cout << gen << "世代目" << endl;
		for(i = 0; i < 200; i++) {
			cout << wpop->pop[i]->fitness1 << "," << wpop->pop[i]->fitness2 << endl;
		}
		cout << endl;
	}
}
