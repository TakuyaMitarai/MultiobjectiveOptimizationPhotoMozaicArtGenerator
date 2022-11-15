// floydSE.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include "SymbioticEvolution.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
using namespace std;

vector<int> B, G, R;           //構成画像のRGBそれぞれの中央値
vector<int> OB, OG, OR;        //元画像の各領域のRGBそれぞれの中央値       
int pchrom_sub[PCHROM_LEN];    //部分解写真番号
int pchrom_num[PCHROM_LEN];    //部分解の順序
int photo_sub[ArtPhotoNum];    //全体解の写真番号

//カンマ区切りで読み込む関数
vector<string> split(string& input, char delimiter)
{
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

//部分解の順序のソート関数
void sortn(int lb, int ub)
{
	int i, j, k;
	double pivot;
	int tmp;

	if(lb < ub) {
		k = (lb + ub) / 2;
		pivot = pchrom_num[k];
		i = lb;
		j = ub;
		do {
			while(pchrom_num[i] < pivot) {
				i++;
			}
			while(pchrom_num[j] > pivot) {
				j--;
			}
			if(i <= j) {
				tmp = pchrom_num[i];
				pchrom_num[i] = pchrom_num[j];
				pchrom_num[j] = tmp;
				tmp = pchrom_sub[i];
				pchrom_sub[i] = pchrom_sub[j];
               	pchrom_sub[j] = tmp;
				i++;
				j--;
			}
		} while(i <= j);
		sortn(lb, j);
		sortn(i, ub);
	}
}

// メインルーチン
int main()
{
	int i, j;
	float comptime;
	clock_t start;
	double used_photo_cnt[PEOPLE_NUM];
	double ave = 0, diff = 0;
	SymbioticEvolution *se;

	// 時間計測スタート
	start = clock();

	// 乱数のタネの設定
	srand((unsigned int)time(NULL));

	ifstream ifs("../OriginalProcessing/OriginalResult.csv");

    //画素値読込
    string line;
    while (getline(ifs, line)) {

        vector<string> strvec = split(line, ',');
       
        OB.push_back(stoi(strvec.at(0)));
        OG.push_back(stoi(strvec.at(1)));
        OR.push_back(stoi(strvec.at(2)));
    }

    //素材画像の画素値の中央値を記述したファイルをopen
    ifstream ifs2("../MaterialProcessing/MaterialResult.csv");

    //画素値読込
    string line2;
    while (getline(ifs2, line2)) {

        vector<string> strvec = split(line2, ',');
       
        B.push_back(stoi(strvec.at(0)));
        G.push_back(stoi(strvec.at(1)));
        R.push_back(stoi(strvec.at(2)));
    }
	
	FILE *fp;
    fp = fopen("Result.csv", "w");

	// 解探索
	se = new SymbioticEvolution();
	
	se->solve();

	// 処理時間計算
	comptime = (float)(clock() - start) / CLOCKS_PER_SEC;

	// 結果表示
	//初期化
	for(i = 0; i < PEOPLE_NUM; i++) {
		used_photo_cnt[i] = 0;
	}

	for(i = 0; i < WCHROM_LEN; i++) {
		for(j = 0; j < PCHROM_LEN; j++) {
			used_photo_cnt[se->best[i * PCHROM_LEN + j] / (MATERIAL_PHOTO_NUM / PEOPLE_NUM)]++;
		}
	}
	//平均
	for(i = 0; i < PEOPLE_NUM; i++) {
		printf("%d人目　%f枚\n",i, used_photo_cnt[i]);
		ave += used_photo_cnt[i];
	}
	ave /= PEOPLE_NUM;

	//分散
	for(i = 0; i < PEOPLE_NUM; i++) {
		diff += (used_photo_cnt[i] - ave) * (used_photo_cnt[i] - ave);
	}
	diff /= PEOPLE_NUM;
	//出力
	for(i = 0; i < ArtPhotoNum; i++) {
		//cout << i << " " << se->best[i] << endl;
		fprintf(fp, "%d,", se->best[i]);
	}

    fclose(fp);

	printf("\n画像としての良さ：%f, 分散：%f\n", se->bestfit - diff * SCALE, diff);
	printf("最良個体の適応度：%f\n", se->bestfit);
	printf("処理時間：%f秒\n", comptime);

	// 後処理
	delete se;

	return 0;
}
