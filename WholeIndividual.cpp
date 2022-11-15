#include "WholeIndividual.h"
#include <vector>
#include <iostream>
using namespace std;
PartialPopulation** WholeIndividual::ppop = NULL;

extern vector<int> B, G, R;     //構成画像のRGBそれぞれの中央値
extern vector<int> OB, OG, OR;  //元画像の各領域のRGBのそれぞれの中央値

// コンストラクタ
WholeIndividual::WholeIndividual()
{
	int i;

	for(i = 0; i < WCHROM_LEN; i++)
		chrom[i] = ppop[i]->pop[rand() % PPOP_SIZE];
	fitness = DBL_MAX;
}

// デストラクタ
WholeIndividual::~WholeIndividual()
{

}

// 次世代生成（ランダム）
void WholeIndividual::newGeneration()
{
	int i;

	for(i = 0; i < WCHROM_LEN; i++)
		chrom[i] = ppop[i]->pop[rand() % PPOP_SIZE];
	fitness = DBL_MAX;
}

// 次世代生成（コピー＆突然変異）
// p: コピー元個体へのポインタ
void WholeIndividual::newGeneration(WholeIndividual* p)
{
	int i;
	
	for(i = 0; i < WCHROM_LEN; i++)
		chrom[i] = p->chrom[i];
	mutate();
	fitness = DBL_MAX;
}

// 次世代生成（２点交叉＆突然変異）
// p1: 親個体その１
// p2: 親個体その２
// index1: 交叉点（0～WCHROM_LEN-1）
// index2: 交叉点（0～WCHROM_LEN-1）
void WholeIndividual::newGeneration(WholeIndividual* p1, WholeIndividual* p2, int index1, int index2)
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
	for( ; i < max; i++)
		chrom[i] = p2->chrom[i];
	for( ; i < PCHROM_LEN; i++)
		chrom[i] = p1->chrom[i];
	mutate();
	fitness = DBL_MAX;
}

// 突然変異
void WholeIndividual::mutate()
{
	int i;

	for(i = 0; i < WCHROM_LEN; i++) {
		if(rand() < MUTATE_PROBI)
			chrom[i] = ppop[i]->pop[rand() % PPOP_SIZE];
	}
}
// 評価1
void WholeIndividual::evaluation1()
{
	int i, j, k;

	fitness1 = 0.0;
	//適応度算出
	for(i = 0; i < ArtPhotoNum; i++) {
		double dis = (B[photo_sub[i]] - OB[i]) * (B[photo_sub[i]] - OB[i]) + (G[photo_sub[i]] - OG[i]) * (G[photo_sub[i]] - OG[i]) + (R[photo_sub[i]] - OR[i]) * (R[photo_sub[i]] - OR[i]);
		dis = sqrt(dis);
		fitness1 += dis;
	}
}
// 評価2
void WholeIndividual::evaluation2()
{
	int i, j, k;
	double ave = 0, diff = 0;

	fitness2 = 0.0;
	//初期化
	for(i = 0; i < PEOPLE_NUM; i++) {
		used_photo_cnt[i] = 0;
	}
	for(i = 0; i < WCHROM_LEN; i++) {
		for(j = 0; j < PCHROM_LEN; j++) {
			photo_sub[i * PCHROM_LEN + j] = chrom[i]->chrom[j];
			used_photo_cnt[chrom[i]->chrom[j] / (MATERIAL_PHOTO_NUM / PEOPLE_NUM)]++;
		}
	}
	//平均
	for(i = 0; i < PEOPLE_NUM; i++) {
		ave += used_photo_cnt[i];
	}
	ave /= PEOPLE_NUM;

	//分散
	for(i = 0; i < PEOPLE_NUM; i++) {
		diff += (used_photo_cnt[i] - ave) * (used_photo_cnt[i] - ave);
	}
	diff /= PEOPLE_NUM;
	fitness2 = diff * SCALE;
}

// 評価
void WholeIndividual::evaluation()
{
	int i, j, k;
	double ave = 0, diff = 0;

	fitness = 0.0;
	//初期化
	for(i = 0; i < PEOPLE_NUM; i++) {
		used_photo_cnt[i] = 0;
	}
	for(i = 0; i < WCHROM_LEN; i++) {
		for(j = 0; j < PCHROM_LEN; j++) {
			photo_sub[i * PCHROM_LEN + j] = chrom[i]->chrom[j];
			used_photo_cnt[chrom[i]->chrom[j] / (MATERIAL_PHOTO_NUM / PEOPLE_NUM)]++;
		}
	}
	//平均
	for(i = 0; i < PEOPLE_NUM; i++) {
		ave += used_photo_cnt[i];
	}
	ave /= PEOPLE_NUM;

	//分散
	for(i = 0; i < PEOPLE_NUM; i++) {
		diff += (used_photo_cnt[i] - ave) * (used_photo_cnt[i] - ave);
	}
	diff /= PEOPLE_NUM;

	//適応度算出
	for(i = 0; i < ArtPhotoNum; i++) {
		double dis = (B[photo_sub[i]] - OB[i]) * (B[photo_sub[i]] - OB[i]) + (G[photo_sub[i]] - OG[i]) * (G[photo_sub[i]] - OG[i]) + (R[photo_sub[i]] - OR[i]) * (R[photo_sub[i]] - OR[i]);
		dis = sqrt(dis);
		fitness += dis;
	}
	fitness += diff * SCALE;
	fitness = fabs(fitness);
}
