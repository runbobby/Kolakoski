#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>
#include <iomanip>

using namespace std; 
//change to unsigned long long for higher counts
typedef unsigned long long countt;

namespace
{
	vector<unsigned int> powers2 = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };

	unsigned int gcd(const unsigned int a, const unsigned int b) {
		return b == 0 ? a : gcd(b, a % b);
	}
	long long gcd(const long long a, const long long b) {
		return b == 0 ? a : gcd(b, a % b);
	}

	void fastTestRatio(const unsigned int length, const vector<unsigned int>& mapM, const vector<unsigned int>& mapN, const vector<countt>& mapMM, const vector<countt>& mapMN, const vector<countt>& mapNM, const vector<countt>& mapNN,
		vector<unsigned int>& cycle, vector<bool>& cycleDigits, const long long& num, const long long& den, unsigned long long& newNum, unsigned long long& newDen)
	{
		unsigned int pow2l = powers2[length];
		unsigned int pow2l_1 = powers2[length - 1];
		//We are testing the ratio num/den (should be >= 1/2) for a cycle whose 1-density > 1/2
		unsigned int iteration = 0;
		//vector<bool> i_visited = vector<bool>(powers2[l], false);
		//i_visited[0]=true;
		bool unstable = true;
		bool zeroUpdated = false;
		vector<long long> weightsM = vector<long long>(pow2l, 0); // M-edge weights times den
		vector<long long> weightsN = vector<long long>(pow2l, 0);
		for (unsigned int i = 0; i < pow2l; i++)
		{
			//unsigned int j = 3 * i;
			weightsM[i] = num*(mapMM[i] + mapMN[i]) - den*(mapMM[i]);
			weightsN[i] = num*(mapNM[i] + mapNN[i]) - den*(mapNM[i]);
		}
		vector<long long> potential = vector<long long>(pow2l, 0x3fffffffffffffff);
		potential[0] = 0;
		vector<unsigned int> parents = vector<unsigned int>(pow2l, 0); //may need to change 2nd arg
		while (iteration <= 10 * length && unstable && !zeroUpdated)
		{
			unstable = false;
			for (unsigned int i = 0; i < pow2l; i++)
			{
				//unsigned int j = 3 * i;
				long long m = potential[i] + weightsM[i];
				unsigned int k = mapM[i];
				if (potential[k] > m) //may need to test for virtual infinity
				{
					potential[k] = m;
					parents[k] = i;
					unstable = true;
					//if (k == 0) { zeroUpdated = true; }
				}
				m = potential[i] + weightsN[i];
				k = mapN[i];
				if (potential[k] > m) //may need to test for virtual infinity
				{
					potential[k] = m;
					parents[k] = i;
					unstable = true;
					//if (k == 0) { zeroUpdated = true; }
				}
			}
			iteration++;
		}
		if (!unstable) { cout << "No negative cycles (Stable?). Previous ratio is probably the answer. \n"; newDen = 0; return; }
		//form cycle
		unsigned int cycleStart = UINT32_MAX;
		long long difference = 0;
		for (unsigned int i = 0; i < pow2l; i++)
		{
			//unsigned int j = 3 * i;
			long long m = potential[i] + weightsM[i];
			unsigned int k = mapM[i];
			if (potential[k] > m + difference) //may need to test for virtual infinity
			{
				difference = potential[k] - m;
				cycleStart = i;
			}
			m = potential[i] + weightsN[i];
			k = mapN[i];
			if (potential[k] > m) //may need to test for virtual infinity
			{
				difference = potential[k] - m;
				cycleStart = i;
			}
		}
		if (cycleStart == UINT32_MAX) { cout << "No negative cycles found so far."; newDen = 0; return; }

		//m, n are unrelated to M, N
		vector<unsigned int> tempCycle = vector<unsigned int>();
		vector<bool> tempCycleDigits = vector<bool>(); //false: M true: N
		vector<long long> cycleCountM = vector<long long>();
		vector<long long> cycleCountN = vector<long long>();
		long long countM = 0;
		long long countN = 0;
		unsigned int m = parents[cycleStart];
		unsigned int n = cycleStart;
		vector<bool> not_visited = vector<bool>(pow2l, true);
		while (not_visited[m])
		{
			not_visited[m] = false;
			//unsigned int m3 = 3 * m;
			if (mapM[m] == n)
			{
				tempCycle.push_back(m);
				tempCycleDigits.push_back(false);
				countM += mapMM[m];
				countN += mapMN[m];
				cycleCountM.push_back(countM);
				cycleCountN.push_back(countN);
			}
			else
			{
				tempCycle.push_back(m);
				tempCycleDigits.push_back(true);
				countM += mapNM[m];
				countN += mapNN[m];
				cycleCountM.push_back(countM);
				cycleCountN.push_back(countN);
			}
			n = m;
			m = parents[n];
		}
		//another one
		//unsigned int m3 = 3 * m;
		if (mapM[m] == n)
		{
			tempCycle.push_back(m);
			tempCycleDigits.push_back(false);
			countM += mapMM[m];
			countN += mapMN[m];
			cycleCountM.push_back(countM);
			cycleCountN.push_back(countN);
		}
		else
		{
			tempCycle.push_back(m);
			tempCycleDigits.push_back(true);
			countM += mapNM[m];
			countN += mapNN[m];
			cycleCountM.push_back(countM);
			cycleCountN.push_back(countN);
		}
		n = m;
		m = parents[n];
		unsigned int i = 0;
		while (tempCycle[i] != m) { i++; }
		unsigned int s = tempCycle.size();
		cycle = vector<unsigned int>(s - i);
		cycleDigits = vector<bool>(s - i);
		for (unsigned int j = 0; j < s - i; j++)
		{
			cycle[j] = tempCycle[s - j - 1];
			cycleDigits[j] = tempCycleDigits[s - j - 1];
		}
		countM = 0;
		countN = 0;
		if (i > 0)
		{
			countM = cycleCountM[i - 1];
			countN = cycleCountN[i - 1];
		}
		long long n1 = cycleCountM[s - 1] - countM;
		long long d1 = n1 + cycleCountN[s - 1] - countN;
		unsigned int g = gcd(n1, d1);
		newNum = n1 / g;
		newDen = d1 / g;
		for (unsigned int i = 0; i < cycle.size(); i++)
		{
			cout << cycle[i] << " ";
		}
		cout << " End cycle. \nRatio: " << newNum << "/" << newDen << "\n";
	}
	/**/
}

int main(int argc, char* argv[])
{
	// Example program
	/*
#include <iostream>
#include <string>
#include<random>
#include<ctime>
#include<vector>
	while (true)
	{
		srand(time(NULL));
		int samples = 0;
		int score = 0;
		std::vector<int> scores(101, 0);
		std::vector<int> count(101, 0);
		std::vector<int> exp_index(100, 0);
		while (samples < 1000000)
		{
			samples++;
			std::vector<int> flips(100);
			int heads = 0;
			for (int i = 0; i < 100; i++) { flips[i] = (rand()) % 2; heads += flips[i]; }
			int heads_seen = 0;
			int guess = 0;
			int this_score = 0;
			for (int i = 0; i < 100; i++)
			{
				if (2 * (heads - heads_seen)>100 - i) { guess = 1; }
				else { guess = 0; }
				if (flips[i] == guess) { this_score++; exp_index[i]++; }
				heads_seen += flips[i];
			}
			count[heads]++;
			scores[heads] += this_score;
			score += this_score;
		}
		std::cout << "Expected score: " << (float)score / samples << ". Score for each # heads\n";
		for (int i = 0; i <= 100; i++)
		{
			if (count[i] != 0)
			{
				std::cout << i << ": " << (float)scores[i] / count[i] << "\t";
			}
		}
		cout << "\nScore by index\n";
		for (int i = 0; i < 100; i++)
		{
			std::cout << i << ": " << exp_index[i] << "\t";
		}
		int a; cin >> a;
	}
	return 0;
	/**/

	/////////////
	unsigned int M = 1;
	unsigned int N = 2;
	cout << "This program computes, for K(m, n), the maps t->C(m, t), t-> #m in E(m, t), t-> # n in E(m, t) etc.\n Enter n: ";
	cin >> N;
	cout << "Enter m: ";
	cin >> M;
	//cout << "Display maps? Enter 1 for yes: ";
	//bool displayMaps = false;
	//cin >> displayMaps;

	//The string abcdefg corresponds to the binary number whose ones digit is (0 if a=m and 1 if a=n), etc.
	//mapM is the map t->C(m, t).
	//mapMM is the map t-> #m in E(m, t).
	//mapMN is the map t-> #n in E(m, t), etc.
	//change the type in mapMM, mapMN etc. to unsigned long long if necessary. 
	//For convenience, the sizes of these arrays has type countt
	//These are initialized for length 1

	vector<unsigned int> mapM = vector<unsigned int>{ 1, 0 };
	vector<unsigned int> mapN = vector<unsigned int>{ 1, 0 };
	vector<countt> mapMM = vector<countt>{ M, 0 };
	vector<countt> mapMN = vector<countt>{ 0, M };
	vector<countt> mapNM = vector<countt>{ N, 0 };
	vector<countt> mapNN = vector<countt>{ 0, N };

	unsigned int length = 2;
	//change this clause to look at higher lengths
	while (length < 20)
	{
		unsigned int pow2l = powers2[length];
		unsigned int pow2l_1 = powers2[length - 1];
		//currently, this algorithm needs to allocate the vectors for both length-1 and length at the same time
		//It may be possible to only have one at a time

		vector<unsigned int> newMapM = vector<unsigned int>(pow2l);
		vector<unsigned int> newMapN = vector<unsigned int>(pow2l);
		vector<countt> newMapMM = vector<countt>(pow2l);
		vector<countt> newMapMN = vector<countt>(pow2l);
		vector<countt> newMapNM = vector<countt>(pow2l);
		vector<countt> newMapNN = vector<countt>(pow2l);

		for (unsigned int i = 0; i < pow2l_1; i++)
		{
			unsigned int j = 2 * i;
			unsigned int k = 0, i1 = i, i2 = i;
			countt c1 = 0, c2 = 0, c3 = 0, c4 = 0;
			while (k < M) { c1 += mapMM[i1]; c2 += mapNM[i2]; c3 += mapMN[i1]; c4 += mapNN[i2]; i1 = mapM[i1]; i2 = mapN[i2]; k++; }
			newMapM[j] = 1 + 2 * i1;
			newMapM[j + 1] = 2 * i2;
			newMapMM[j] = c1;		//#m in E(m, mi)
			newMapMM[j + 1] = c2;	//#m in E(m, ni)
			newMapMN[j] = c3;		//#n in E(m, mi)
			newMapMN[j + 1] = c4;	//#n in E(m, ni)

			k = 0, i1 = i, i2 = i;
			c1 = 0, c2 = 0, c3 = 0, c4 = 0;
			while (k < N) { c1 += mapMM[i1]; c2 += mapNM[i2]; c3 += mapMN[i1]; c4 += mapNN[i2]; i1 = mapM[i1]; i2 = mapN[i2]; k++; }
			newMapN[j] = 1 + 2 * i1;
			newMapN[j + 1] = 2 * i2;
			newMapNM[j] = c1;		//#m in E(n, mi)
			newMapNM[j + 1] = c2;	//#m in E(n, ni)
			newMapNN[j] = c3;		//#n in E(n, mi)
			newMapNN[j + 1] = c4;	//#n in E(n, ni)
		}
		mapM = newMapM;
		mapN = newMapN;
		mapMM = newMapMM;
		mapMN = newMapMN;
		mapNM = newMapNM;
		mapNN = newMapNN;

		//process the arrays here. example:
		cout << "i: mapM, mapMM, mapMN, mapN, mapNM, mapNN\n";
		for (int i = 0; i < mapM.size(); i++)
		{
			cout << std::setw(10) << i << " " << std::setw(10) << mapM[i] << " " << std::setw(20) << mapMM[i] << " " << std::setw(20) << mapMN[i] << " " << std::setw(10) << mapN[i] << " " << std::setw(20) << mapNM[i] << " " << std::setw(20) << mapNN[i] << "\n";
		}/**/
		int b; cin >> b;
		/*
		unsigned long long currentNum = 1; unsigned long long currentDen = 2;
		vector<unsigned int> cycle{};
		vector<bool> cycleDigits{};
		while (currentDen > 0)
		{
			fastTestRatio(length, mapM, mapN, mapMM, mapMN, mapNM, mapNN, cycle, cycleDigits, currentNum, currentDen, currentNum, currentDen);
			cout << (double)(currentNum) / currentDen << "\t";
		}
		/**/

		cout << "\n Finished length " << length << ".\n";

		int a; cin >> a; //(pause)

		length++;
	}
	return 0;
}