#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <bitset>
//#include <random>
//#include <ctime>

using namespace std;

//assumes map is bijective
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
		if ( ! visited[i] )
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


int main(int argc, char* argv[])
{
	vector<unsigned int> powers2 = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };
	/*
	unsigned int M = 1;
	unsigned int N = 2;
	//cout << "Enter N: ";
	//cin >> N;
	//cout << "Enter M: ";
	//cin >> M;
	for (M=1; M < 10; M += 2)
	{
		for (N=2; N < 33; N += 2)
		{
			vector<unsigned int> mapM = vector<unsigned int>{ 1, 0 };
			vector<unsigned int> mapN = vector<unsigned int>{ 1, 0 };
			unsigned int length = 2;
			while (length < 25)
			{
				unsigned int pow2l = powers2[length];
				unsigned int pow2l_1 = powers2[length - 1];
				vector<unsigned int> newMapM = vector<unsigned int>(pow2l);
				vector<unsigned int> newMapN = vector<unsigned int>(pow2l);
				vector<unsigned int> powMapMM = power(mapM, M);
				vector<unsigned int> powMapMN = power(mapM, N);
				vector<unsigned int> powMapNM = power(mapN, M);
				vector<unsigned int> powMapNN = power(mapN, N);
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
				//cout << "Finished length " << length << "\n";
				length++;
			}
			///look at distribution of starting points for the first 'terms' terms
			unsigned int terms = 4000000000;
			vector<unsigned int> startsDist = vector<unsigned int>(powers2[length-1], 0);
			unsigned int currentStart = 0;
			vector<bool> KMN = vector<bool>(terms + N + M);
			unsigned int i;
			for (i = 0; i < M; i++)
			{
				KMN[i] = false;
			}
			KMN[M] = true;
			bool value = true;
			unsigned int block = 1; // 0-indexed
			while (i < terms)
			{
				if (KMN[block])
				{
					for (int j = 0; j < N; j++)
					{
						KMN[i] = value; i++;
					}
				}
				else
				{
					for (int j = 0; j < M; j++)
					{
						KMN[i] = value; i++;
					}
				}
				value = !value;
				block++;
				//cout << "i:" << i << " ";
				//cout << block << " ";
			}
			//for (i = 0; i < 1000; i++){cout << KMN[i];}
			for (i = 0; i < terms; i++)
			{
				startsDist[currentStart]++;
				if (KMN[i]) { currentStart = mapN[currentStart]; }
				else { currentStart = mapM[currentStart]; }
			}
			//cout << "First 1000 terms of Distribution: \n";
			//for (i = 0; i < 1000; i++){cout << startsDist[i] << " ";}
			double chi2 = 0;
			double exp = (double)terms / startsDist.size();
			//cout << exp << " ";
			int min = exp;
			int max = 0;
			for (i = 0; i < startsDist.size(); i++)
			{
				unsigned int j = startsDist[i];
				if (j > max) { max = j; }
				if (j < min) { min = j; }
				chi2 += (j - exp)*(j - exp) / exp;
			}
			cout << "\nFor (m, n) = (" << M << ", " << N << "), Chi2/dF: " << chi2/(startsDist.size() - 2) << ". min: " << min << ". max: " << max;

		}
	}/**/

	//===================================
	int M = 1;
	//cout << "Enter m: ";
	//cin >> M;
	//cout << " Entered " << M;
	constexpr unsigned int L = 23;
	for (M = 3; M < powers2[L/2+1]; M+=2)
	{
		unsigned int N = 2;
		unsigned int p2L = powers2[L];
		//cout << "Enter N: ";
		//cin >> N;
		//cout << "Enter M: ";
		//cin >> M;
		vector<unsigned int> hasMaxOrbit = vector<unsigned int>(p2L); for (int i = 0; i < p2L; i++) { hasMaxOrbit[i] = i; }
		//for (M = -1; M < 6; M += 2)
		//{
		for (N = 2; N < 33; N += 2)
		{
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
				if (length > 23) { cout << length << " "; }
				length++;
			}
			vector<unsigned int> map3 = mapPower(mapM, powers2[L / 2]);

			vector<unsigned int> newHMO = vector<unsigned int>();
			//cout << "HMO size " << hasMaxOrbit.size();
			for (int i = 0; i < hasMaxOrbit.size(); i++)
			{
				if (map3[hasMaxOrbit[i]] != hasMaxOrbit[i]) { newHMO.push_back(hasMaxOrbit[i]); }
			}
			hasMaxOrbit = newHMO;
			if (hasMaxOrbit.size() == 0) { break; }
			/**/
			//cout << N << " ";
		}
		if (hasMaxOrbit.size() > 0)
		{
			cout << "\nFor m = " << M << ", the following have max orbit for all n so far: ";
			for (int i = 0; i < __min(100, hasMaxOrbit.size()); i++)
			{
				cout << hasMaxOrbit[i] << " ";
			}
			cout << "\n in binary: ";
			for (int i = 0; i < __min(100, hasMaxOrbit.size()); i++)
			{
				std::bitset<L> t(hasMaxOrbit[i]);
				cout << t << " ";
			}
			cout << "\n";
			int a; cin >> a;
		} 
		else { cout << M << " "; }
	}
}