#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
//#include <random>
//#include <ctime>

using namespace std;

//assumes power > 0 and that map is bijective
vector<unsigned int> mapPower(const vector<unsigned int>& map, int power)
{
	if (power == 1) { return map; }
	if (power == 0) {
		vector<unsigned int> newMap = vector<unsigned int>(map.size());
		for (int i = 0; i < map.size(); i++) { newMap[i] = i; }
	}
	if (power < 0)
	{
		vector<unsigned int> newMap = vector<unsigned int>(map.size());
		for (int i = 0; i < map.size(); i++) { newMap[map[i]] = i; }
		return mapPower(newMap, -power);
	}


	vector<unsigned int> newMap = vector<unsigned int>(map.size());
	vector<bool> visited = vector<bool>(map.size(), false);
	for (unsigned int i = 0; i < map.size(); i++)
	{
		if (!visited[i])
		{
			visited[i] = true;
			vector<unsigned int> queue = vector<unsigned int>(power + 1);
			int j = 1;
			queue[0] = i;
			unsigned int k = map[i];
			while ((k != i) && (j < power))
			{
				visited[k] = true;
				queue[j] = k;
				j++;
				k = map[k];
			}
			if (j == power)
			{
				visited[k] = true;
				newMap[i] = k;
				queue[power] = k;
				unsigned int p1 = power + 1;
				unsigned int l = queue[1];
				int m = power + 2;
				while (l != i)
				{
					int z = (m - 1) % p1;
					queue[z] = map[queue[(m - 2) % p1]];
					newMap[l] = queue[z];
					m++;
					l = queue[(m) % p1];
					visited[queue[z]] = true;
				}
			}
			else
			{
				for (int l = 0; l < j; l++)
				{
					newMap[queue[l]] = queue[(l + power) % j];
				}
			}
		}
	}
	return newMap;
}

int main2(int argc, char* argv[])
{
	//srand(time(NULL));
	vector<unsigned int> powers2 = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };
	unsigned int M = 1;
	unsigned int N = 2;
	unsigned int L = 25;
	//cout << "Display maps? Enter 1 for yes: ";
	//bool displayMaps = false;
	//cin >> displayMaps;
	//The string abcdefg corresponds to the binary number whose ones digit is a-1, etc.
	//Step 1: generate the 25th map. Discard the earlier ones
	//The maps are slightly different than before. For digits x, y,
	//mapx[i] = C(x, i). lengthMapx[i] = (bool) Length(E[x,iy])%2 (independent of the digit y)
	while (N < powers2[L/2])
	{
		//cout << "Testing N = " << N << ".\n";
		vector<unsigned int> mapM = vector<unsigned int>{ 1, 0 };
		vector<unsigned int> mapN = vector<unsigned int>{ 1, 0 };
		unsigned int length = 2;
		while (length <= L)
		{
			unsigned int pow2l = powers2[length];
			unsigned int pow2l_1 = powers2[length - 1];
			vector<unsigned int> newMapM = vector<unsigned int>(pow2l);
			vector<unsigned int> newMapN = vector<unsigned int>(pow2l);
			vector<unsigned int> powMapMM = mapPower(mapM, M);
			vector<unsigned int> powMapMN = mapPower(mapM, N);
			vector<unsigned int> powMapNM = mapPower(mapN, M);
			vector<unsigned int> powMapNN = mapPower(mapN, N);
			for (unsigned int i = 0; i < pow2l_1; i++)
			{
				unsigned int j = 2 * i;
				newMapM[j] = 1 + 2 * powMapMM[i];
				newMapM[j + 1] = 2 * powMapNM[i];
				newMapN[j] = 1 + 2 * powMapMN[i];
				newMapN[j + 1] = 2 * powMapNN[i];
			}
			mapM = newMapM;
			mapN = newMapN;
			//bool orbit2 = false;
			/*if (displayMaps)
			{
				for (int i = 0; i < map1.size(); i++)
				{
					cout << i << "->" << map1[i] << " ";
				}
				cout << "\n";
				for (int i = 0; i < map1.size(); i++)
				{
					cout << i << "->" << map2[i] << " ";
					//if (map2[map2[i]] == i) { orbit2 = true; }
				}
			}/**/

			length++;
		}
		unsigned int i = mapM[0];
		unsigned int j = 1;
		while (i != 0) { i = mapM[i]; j++; }
		if (j != powers2[(L+1)/2])
		{
			cout << "N = " << N << ": 0's orbit " << j << ". ";
		}
		else
		{
			cout << N << " ";
		}

		/*unsigned int i = 0b101010101010101010101010;
		unsigned int j = 1;
		while (i != 0) { i = map1[i]; j++; }
		if (j != powers2[(L+1)/2])
		{
			cout << "N = " << N << ": orbit " << j << ". ";
		}
		else
		{
			cout << N << " ";
		}/**/



		N += 2;
	}
	int a; cin >> a;
	return 0;
}