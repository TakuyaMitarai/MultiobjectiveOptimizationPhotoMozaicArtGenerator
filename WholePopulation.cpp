#include "WholePopulation.h"

// コンストラクタ
WholePopulation::WholePopulation()
{
	int i;

	for(i = 0; i < WPOP_SIZE; i++)
		pop[i] = new WholeIndividual();
	for(i = 0; i < CHILDPOP; i++)
		childpop[i] = new WholeIndividual();
}

// デストラクタ
WholePopulation::~WholePopulation()
{
	int i;

	for(i = 0; i < WPOP_SIZE; i++)
		delete pop[i];
	for(i = 0; i < CHILDPOP; i++)
		delete childpop[i];
}

// 次世代生成
void WholePopulation::newGeneration()
{
	int i, a, b, index1, index2;
	
	/*
	a = rand() % WPOP_SIZE;
	b = (rand() % (WPOP_SIZE - 1) + a + 1) % WPOP_SIZE;
	index1 = rand() % WCHROM_LEN;
	index2 = ((rand() % (WCHROM_LEN - 1)) + index1) % WCHROM_LEN;
	childpop[0]->newGeneration(pop[a], pop[b], index1, index2);
	childpop[1]->newGeneration(pop[b], pop[a], index1, index2);
	childpop[2]->newGeneration(pop[a]);
	childpop[3]->newGeneration(pop[b]);
	childpop[4]->newGeneration();
	childpop[5]->newGeneration();
	changePC(a, b);
	for(i = 0; i < WPOP_SIZE; i++) {
		if((i != a) && (i != b))
			pop[i]->newGeneration();
	}
	*/
	for(i = WPOP_SIZE / 4 * 3; i < WPOP_SIZE; i++) {
		pop[i]->newGeneration();
	}
}

// 親子選択
void WholePopulation::changePC(int a, int b)
{
	int i, best, worst;
	double alpha, beta, fitsum, r;
	WholeIndividual* tmp;

	// 最良個体を残す
	childpop[0]->evaluation();
	for(best = 0, i = 1; i < CHILDPOP; i++) {
		childpop[i]->evaluation();
		if(childpop[best]->fitness > childpop[i]->fitness)
			best = i;
	}
	if(best) {
		tmp = childpop[0];
		childpop[0] = childpop[best];
		childpop[best] = tmp;
	}
	tmp = pop[a];
	pop[a] = childpop[0];
	childpop[0] = tmp;
	
	// スケーリング法で適応度を変換
	for(best = 1, worst = 1, i = 2; i < CHILDPOP; i++) {
		if(childpop[worst]->fitness < childpop[i]->fitness)
			worst = i;
		if(childpop[best]->fitness > childpop[i]->fitness)
			best = i;
	}
	alpha = (100000 / (double)(CHILDPOP - 1) - 100) / (childpop[best]->fitness - childpop[worst]->fitness);
	beta = 100 - alpha * childpop[worst]->fitness;
	fitsum = 0.0;
	for(i = 1; i < CHILDPOP; i++) {
		childfit[i] = alpha * childpop[i]->fitness + beta;
		fitsum += childfit[i];
	} 
	for(i = 1; i < CHILDPOP; i++)
		childfit[i] /= fitsum;

	// 残りからルーレットで１つ選んで残す
	r = (double)rand() / RAND_MAX;
	for(i = 1; i < CHILDPOP - 1; i++) {
		if(r < childfit[i])
			break;
		r -= childfit[i];
	}
	tmp = pop[b];
	pop[b] = childpop[i];
	childpop[i] = tmp;
}

// 評価
void WholePopulation::evaluation()
{
	int i, j, k;
	int cnt = 0;
	int current_cnt;
	int flag;
	int rank = 1;
	int fj;
	double fit1, fit2;
	double min1, min2;
	int num;
	vector<int> now_s;
	vector<int> next_s;

	//初期化
	for(i = 0; i < WPOP_SIZE; i++) {
		next_s.push_back(i);
	}

	for(i = 0; i < WPOP_SIZE; i++) {
		pop[i]->evaluation1();
		pop[i]->evaluation2();
	}
	
	fit1 = -1;
	fit2 = 10000000;

	while(cnt < WPOP_SIZE) {
		min1 = 10000000;
		flag = 1;
		for(i = 0; i < WPOP_SIZE - cnt; i++) {
			if(fit1 < pop[next_s[i]]->fitness1 && fit2 > pop[next_s[i]]->fitness2) {
				if(min1 > pop[next_s[i]]->fitness1) {
					min1 = pop[next_s[i]]->fitness1;
					num = next_s[i];
					flag = 0;
				}
			}
		}
		if(flag == 0) {
			pop[num]->fitness = rank;		
			fit1 = pop[num]->fitness1;
			fit2 = pop[num]->fitness2;
			for(i = 0; i < WPOP_SIZE - cnt; i++) {
				if(next_s[i] == num) {
					next_s.erase(next_s.begin() + i);
					break;
				}
			}
			cnt++;
		} else {
			rank++;
			fit1 = -1;
			fit2 = 10000000;
		}
	}
	sort(0, WPOP_SIZE-1);

	rank = 1;
	j = 0;
	for(i = 1; i <= MaxRankUsingCR; i++) {
		fj = j;
		while(rank == i && j <= WPOP_SIZE) {
			if(pop[j]->fitness != i) {
				rank++;
			}
			j++;
		}
		if(j <= WPOP_SIZE)
			congestion_sort(fj, j - 1);
	}
	m = 400;
	for(i = 0; i < WPOP_SIZE; i++) {
		// 参照先の個体の適応度を算出
		for(j = 0; j < WCHROM_LEN; j++) {
			if(pop[i]->chrom[j]->fitness > pop[i]->fitness) {
				pop[i]->chrom[j]->fitness = pop[i]->fitness;
			}
		}
	}
}

void WholePopulation::distance()
{
	int i, j, k;
	int dis1, dis2, dis;
	int min1;
	int num[m];
	for(i = 0; i < WCHROM_LEN; i++) {
		for(j = 0; j < PPOP_SIZE; j++) {
			WholeIndividual::ppop[i]->distance[j] = NULL;
		}
	}
	for(i = 0; i < m; i++) {
		min1 = 10000000;
		for(k = 0; k < m; k++) {
			dis1 = pop[i]->fitness1 - pop[k]->fitness1;
			dis2 = pop[i]->fitness2 - pop[k]->fitness2;
			if(dis1 < 0) dis1 = -dis1;
			if(dis2 < 0) dis2 = -dis2;
			dis = dis1 + dis2;
			if(dis <= min1 && i != k) {
				min1 = dis;
				num[i] = k;
			}
		}
	}
	for(i = 0; i < m; i++) {
		for(j = 0; j < WCHROM_LEN; j++) {
			for(k = 0; k <= PPOP_SIZE; k++) {
				if(pop[i]->chrom[j] == WholeIndividual::ppop[j]->pop[k]) break;
			}
			WholeIndividual::ppop[j]->distance[k] = pop[num[i]]->chrom[j];
		}
	}
}

// 昇順にクイックソート
// lb : 開始点添字
// ub : 終了点添字
void WholePopulation::sort(int lb, int ub)
{
	int i, j, k;
	double pivot;
	WholeIndividual* swap;

	if(lb < ub) {
		k = (lb + ub) / 2;
		pivot = pop[k]->fitness;
		i = lb;
		j = ub;
		do {
			while(pop[i]->fitness < pivot)
				i++;
			while(pop[j]->fitness > pivot)
				j--;
			if(i <= j) {
				swap = pop[i];
				pop[i] = pop[j];
				pop[j] = swap;
				i++;
				j--;
			}
		} while(i <= j);
		sort(lb, j);
		sort(i, ub);
	}
}

void WholePopulation::congestion_sort(int lb, int ub)
{
	int i, j, k;
	double pivot;
	WholeIndividual* swap;
	double min1;
	double min2, min3;
	double dis1, dis2, dis;

	for(i = lb; i <= ub; i++) {
		min1 = 10000000;
	 	min2 = 10000000;
		min3 = 10000000;
		for(j = lb; j<= ub; j++) {
			if(i != j) {
				dis1 = pop[i]->fitness1 - pop[j]->fitness1;
				dis2 = pop[i]->fitness2 - pop[j]->fitness2;
				if(dis1 < 0) dis1 = -dis1;
				if(dis2 < 0) dis2 = -dis2;
				dis = dis1 + dis2;
				if(dis < min1) {
					min1 = dis ;
				} else if(dis < min2) {
					min2 = dis;
				} else if(dis < min3) {
					min3 = dis;
				}
			}
		}
		if(min1 == 10000000) min1 = 0;
		if(min2 == 10000000) min2 = 0;
		if(min3 == 10000000) min3 = 0;
		pop[i]->fitness += -(min1 + min2 + min3) / 3 / 1000000;
	}

	if(lb < ub) {
		k = (lb + ub) / 2;
		pivot = pop[k]->fitness;
		i = lb;
		j = ub;
		do {
			while(pop[i]->fitness < pivot)
				i++;
			while(pop[j]->fitness > pivot)
				j--;
			if(i <= j) {
				swap = pop[i];
				pop[i] = pop[j];
				pop[j] = swap;
				i++;
				j--;
			}
		} while(i <= j);
		sort(lb, j);
		sort(i, ub);
	}
}
