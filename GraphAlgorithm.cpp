#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <random>

using namespace std;

vector<unsigned int> powers2 = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };

namespace
{

	unsigned int gcd(const unsigned int a, const unsigned int b) {
		return b == 0 ? a : gcd(b, a % b);
	}
	long long gcd(const long long a, const long long b) {
		return b == 0 ? a : gcd(b, a % b);
	}
	//probably need unsigned long long gcd

	void fastTestRatio(const unsigned int length, const vector<unsigned int>& map1, const vector<unsigned int>& map2,
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
		vector<long long> weights1 = vector<long long>(pow2l, 0); // 1-edge weights times den
		vector<long long> weights2 = vector<long long>(pow2l, 0);
		for (unsigned int i = 0; i < pow2l; i++)
		{
			unsigned int j = 3 * i;
			weights1[i] = num*(map1[j + 1] + map1[j + 2]) - den*(map1[j + 1]);
			weights2[i] = num*(map2[j + 1] + map2[j + 2]) - den*(map2[j + 1]);
		}
		vector<long long> potential = vector<long long>(pow2l, 0xffffffffffffff);
		potential[0] = 0;
		vector<unsigned int> parents = vector<unsigned int>(pow2l, 0); //may need to change 2nd arg
		while (iteration <= 10 * length && unstable && !zeroUpdated)
		{
			unstable = false;
			for (unsigned int i = 0; i < pow2l; i++)
			{
				unsigned int j = 3 * i;
				long long m = potential[i] + weights1[i];
				unsigned int k = map1[j];
				if (potential[k] > m) //may need to test for virtual infinity
				{
					potential[k] = m;
					parents[k] = i;
					unstable = true;
					//if (k == 0) { zeroUpdated = true; }
				}
				m = potential[i] + weights2[i];
				k = map2[j];
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
			unsigned int j = 3 * i;
			long long m = potential[i] + weights1[i];
			unsigned int k = map1[j];
			if (potential[k] > m + difference) //may need to test for virtual infinity
			{
				difference = potential[k] - m;
				cycleStart = i;
			}
			m = potential[i] + weights2[i];
			k = map2[j];
			if (potential[k] > m) //may need to test for virtual infinity
			{
				difference = potential[k] - m;
				cycleStart = i;
			}
		}
		if (cycleStart == UINT32_MAX) { cout << "No negative cycles found so far."; newDen = 0; return; }

		vector<unsigned int> tempCycle = vector<unsigned int>();
		vector<bool> tempCycleDigits = vector<bool>();
		vector<long long> cycleCount1 = vector<long long>();
		vector<long long> cycleCount2 = vector<long long>();
		long long count1 = 0;
		long long count2 = 0;
		unsigned int m = parents[cycleStart];
		unsigned int n = cycleStart;
		vector<bool> not_visited = vector<bool>(pow2l, true);
		while (not_visited[m])
		{
			not_visited[m] = false;
			unsigned int m3 = 3 * m;
			if (map1[m3] == n)
			{
				tempCycle.push_back(m);
				tempCycleDigits.push_back(false);
				count1 += map1[m3 + 1];
				count2 += map1[m3 + 2];
				cycleCount1.push_back(count1);
				cycleCount2.push_back(count2);
			}
			else
			{
				tempCycle.push_back(m);
				tempCycleDigits.push_back(true);
				count1 += map2[m3 + 1];
				count2 += map2[m3 + 2];
				cycleCount1.push_back(count1);
				cycleCount2.push_back(count2);
			}
			n = m;
			m = parents[n];
		}
		//another one
		unsigned int m3 = 3 * m;
		if (map1[m3] == n)
		{
			tempCycle.push_back(m);
			tempCycleDigits.push_back(false);
			count1 += map1[m3 + 1];
			count2 += map1[m3 + 2];
			cycleCount1.push_back(count1);
			cycleCount2.push_back(count2);
		}
		else
		{
			tempCycle.push_back(m);
			tempCycleDigits.push_back(true);
			count1 += map2[m3 + 1];
			count2 += map2[m3 + 2];
			cycleCount1.push_back(count1);
			cycleCount2.push_back(count2);
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
		count1 = 0;
		count2 = 0;
		if (i > 0)
		{
			count1 = cycleCount1[i - 1];
			count2 = cycleCount2[i - 1];
		}
		long long n1 = cycleCount1[s - 1] - count1;
		long long d1 = n1 + cycleCount2[s - 1] - count2;
		unsigned int g = gcd(n1, d1);
		newNum = n1 / g;
		newDen = d1 / g;
		for (unsigned int i = 0; i < cycle.size(); i++)
		{
			cout << cycle[i] << " ";
		}
		cout << " End cycle. \nRatio: " << newNum << "/" << newDen << "\n";
	}/**/

	//fastTestRatio MN version?
	
	/*
	typedef unsigned long long countt;

	
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
		vector<long long> potential = vector<long long>(pow2l, 0xfffffffffffffff);
		potential[0] = 0;
		vector<unsigned int> parents = vector<unsigned int>(pow2l, 0); //may need to change 2nd arg
		while (iteration <= 10 * length && unstable && !zeroUpdated)
		{
			unstable = false;
			for (unsigned int i = 0; i < pow2l; i++)
			{
				unsigned int j = 3 * i;
				long long m = potential[i] + weightsM[i];
				unsigned int k = mapM[j];
				if (potential[k] > m) //may need to test for virtual infinity
				{
					potential[k] = m;
					parents[k] = i;
					unstable = true;
					//if (k == 0) { zeroUpdated = true; }
				}
				m = potential[i] + weightsN[i];
				k = mapN[j];
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
			unsigned int j = 3 * i;
			long long m = potential[i] + weightsM[i];
			unsigned int k = mapM[j];
			if (potential[k] > m + difference) //may need to test for virtual infinity
			{
				difference = potential[k] - m;
				cycleStart = i;
			}
			m = potential[i] + weightsN[i];
			k = mapN[j];
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
	

	///The string abcdefg corresponds to the binary number whose ones digit is a-1, etc.
	vector<vector<unsigned int>> maps1{};
	vector<vector<unsigned int>> maps2{};
	maps1.push_back(vector<unsigned int>{ 1,1,0,0,0,1 });
	maps2.push_back(vector<unsigned int>{ 1,2,0,0,0,2 });
	unsigned int length = 2;	
	while (length < 29)
	{
		unsigned int pow2l = powers2[length];
		unsigned int pow2l_1 = powers2[length - 1];
		vector<unsigned int> newmap1 = vector<unsigned int>(3 * pow2l);
		vector<unsigned int> newmap2 = vector<unsigned int>(3 * pow2l);
		vector<unsigned int> lastmap1 = maps1[length - 2];
		vector<unsigned int> lastmap2 = maps2[length - 2];
		for (unsigned int i = 0; i < pow2l_1; i++)
		{
			unsigned int j = 6 * i;
			unsigned int k = 3 * i;
			newmap1[j] = 1 + 2 * lastmap1[k];
			newmap1[j + 1] = lastmap1[k + 1];
			newmap1[j + 2] = lastmap1[k + 2];
			newmap1[j + 3] = 2 * lastmap2[k];
			newmap1[j + 4] = lastmap2[k + 1];
			newmap1[j + 5] = lastmap2[k + 2];

			newmap2[j] = 1 + 2 * lastmap1[3 * lastmap1[k]];
			newmap2[j + 1] = lastmap1[k + 1] + lastmap1[3 * (lastmap1[k]) + 1];
			newmap2[j + 2] = lastmap1[k + 2] + lastmap1[3 * (lastmap1[k]) + 2];
			newmap2[j + 3] = 2 * lastmap2[3 * lastmap2[k]];
			newmap2[j + 4] = lastmap2[k + 1] + lastmap2[3 * (lastmap2[k]) + 1];
			newmap2[j + 5] = lastmap2[k + 2] + lastmap2[3 * (lastmap2[k]) + 2];
		}
		maps1.push_back(newmap1);
		maps2.push_back(newmap2);

		//for (unsigned int i = 0; i < newmap1.size(); i++) { cout << newmap1[i] << " "; }
		//cout << "\n";
		//for (unsigned int i = 0; i < newmap1.size(); i++) { cout << newmap2[i] << " "; }
		//cout << "\n";
		/*if (length > 20)
		{
			//unsigned long long currentNum = 1001; unsigned long long currentDen = 2000;
			unsigned long long currentNum = 1; unsigned long long currentDen = 2;
			vector<unsigned int> cycle{};
			vector<bool> cycleDigits{};
			while (currentDen > 0)
			{
				fastTestRatio(length, newmap1, newmap2, cycle, cycleDigits, currentNum, currentDen, currentNum, currentDen);
				//fastTestRatio(length, mapM, mapN, mapMM, mapMN, mapNM, mapNN, cycle, cycleDigits, currentNum, currentDen, currentNum, currentDen);
			}
		}/**/

		vector<bool> visited = vector<bool>(pow2l, false);
		visited[0] = true;
		unsigned int numLeft = pow2l - 1;
		while (numLeft > 0)
		{
			for (int i = 0; i < visited.size(); i++)
			{
				if (visited[i])
				{
					if (!visited[newmap1[3 * i]])
					{
						visited[newmap1[3 * i]] = true; numLeft--;
					}
					if (!visited[newmap2[3 * i]])
					{
						visited[newmap2[3 * i]] = true; numLeft--;
					}
				}
			}
			cout << "\n number Left: " << numLeft;
		}

		cout << "Finished length " << length << ".\n ";

		length++;
	}
	///possibly print to file.

	cout << "Enter anything to end.";
	cin >> length;
	return 0;	
}


