#pragma once

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <vector>
#include <iostream>
using namespace std;

// 定数の定義
#define GENERATION_MAX	40000					// 世代交代数
#define WPOP_SIZE		3000					// 全体解集団内の個体の数
#define PPOP_SIZE		3000					// 部分解集団内の個体の数
#define MUTATE_PROB		0.00001 * RAND_MAX		// 突然変異確率×乱数の最大値1
#define MUTATE_PROBI	0.0001 * RAND_MAX		// 突然変異確率×乱数の最大値2
#define MAX_NUM			WCHROM_LEN*PCHROM_LEN	// フロイド問題におけるN
#define MATERIAL_PHOTO_NUM 40000				//構成画像の枚数
#define PEOPLE_NUM		625						//人数 500 or 1000
#define WCHROM_LEN		75						// 全体解の染色体の長さ
#define PCHROM_LEN		75						// 部分解の染色体の長さ
#define CHILDPOP		6
#define ArtPhotoNum     5625                 	//モザイクアートの構成写真枚数
#define SCALE			7000					//分散の重み
#define MaxRankUsingCR  5						//混雑度ランキングが使用されるパレートランクの最大値
#define RANKCHOICE      300						//ランキング選択で選ばれる上位数

class PartialIndividual
{
public:
	PartialIndividual();
	~PartialIndividual();
	void newGeneration();
	void newGeneration(PartialIndividual* p);
	void newGeneration(PartialIndividual* p1, PartialIndividual* p2, int index1, int index2);

	int chrom[PCHROM_LEN];	// 染色体
	double fitness;			// 適応度

private:
	void mutate();
};
