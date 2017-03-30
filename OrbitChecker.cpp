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
	unsigned int maxLength = 26;
	unsigned int numIterations = 4096;
}

int main(int argc, char* argv[])
{
	//The string abcdefg corresponds to the binary number whose ones digit is a-1, etc.
	//Step 1: generate the 25th map. Discard the earlier ones
	//The maps are slightly different than before. For digits x, y,
	//mapx[i] = C(x, i). lengthMapx[i] = (bool) Length(E[x,iy])%2 (independent of the digit y)
	vector<unsigned int> map1 = vector<unsigned int>{ 1, 0 };
	vector<unsigned int> map2 = vector<unsigned int>{ 1, 0 };
	vector<bool> lengthMap1 = vector<bool>{ false, true };
	vector<bool> lengthMap2 = vector<bool>{ true, true };
	unsigned int length = 2;
	while (length < maxLength)
	{
		unsigned int pow2l = powers2[length];
		unsigned int pow2l_1 = powers2[length - 1];
		vector<unsigned int> newMap1 = vector<unsigned int>(pow2l);
		vector<unsigned int> newMap2 = vector<unsigned int>(pow2l);
		vector<bool> newLMap1 = vector<bool>(pow2l);
		vector<bool> newLMap2 = vector<bool>(pow2l);
		for (unsigned int i = 0; i < pow2l_1; i++)
		{
			unsigned int j = 2 * i;
			newMap1[j] = 1 + 2 * map1[i];
			newMap1[j+1] = 2 * map2[i];
			newLMap1[j] = lengthMap1[i];
			newLMap1[j + 1] = lengthMap2[i];

			newMap2[j] = 1 + 2 * map1[map1[i]];
			newMap2[j + 1] = 2 * map2[map2[i]];
			newLMap2[j] = lengthMap1[i] != lengthMap1[map1[i]];
			newLMap2[j+1] = lengthMap2[i] != lengthMap2[map2[i]];
		}
		map1 = newMap1;
		map2 = newMap2;
		lengthMap1 = newLMap1;
		lengthMap2 = newLMap2;
		cout << "Finished length " << length << ".\n";

		length++;
	}
	//possibly print to file.
	//Step 2: Let s^1 be the string 12211. s^{i+1} is the string E(s^i s^i*, 11)
	//where s* is the complement of s (i.e. s^1* is 21122)
	//For i from 1 to 12, we compute Length(E(s^i s^i* s^i s^i* ... s^i s^i*, 111....111))%2
	//where the first argument has 2^13 terms and the second argument has length 26.
	vector<bool> initialString = vector<bool>{ false, true, true, false, false };
	for (int i = 1; i <= 12; i++)
	{
		unsigned int startingPoints = 0;
		bool length = false;
		for (int j = 0; j < numIterations; j++)
		{
			for (int k = 0; k < initialString.size(); k++)
			{
				if (initialString[k])
				{
					length = length != lengthMap2[startingPoints];
					startingPoints = map2[startingPoints];
				}
				else
				{
					length = length != lengthMap1[startingPoints];
					startingPoints = map1[startingPoints];
				}
			}
			for (int k = 0; k < initialString.size(); k++)
			{
				if (initialString[k])
				{
					length = length != lengthMap1[startingPoints];
					startingPoints = map1[startingPoints];
				}
				else
				{
					length = length != lengthMap2[startingPoints];
					startingPoints = map2[startingPoints];
				}
			}
		}
		cout << "i = " << i << ": length is ";
		if (length) { cout << "odd. \n"; }
		else { cout << "even. \n"; }

		
		//generate string s^{i+1}; The following is perhaps too optimized. Not sorry
		int size = initialString.size();
		for (int i = 0; i < size; i++) { initialString.push_back(!initialString[i]); }
		
		vector<bool> newString{};
		newString.reserve(initialString.size()*1.6);
		bool nextDigit = false;
		for (int i = 0; i < initialString.size(); i++)
		{
			newString.push_back(nextDigit);
			if (initialString[i]) { newString.push_back(nextDigit); }
			nextDigit = !nextDigit;
		}
		initialString = newString;
		newString = vector<bool>{};
		newString.reserve(initialString.size()*1.6);
		nextDigit = false;
		for (int i = 0; i < initialString.size(); i++)
		{
			newString.push_back(nextDigit);
			if (initialString[i]) { newString.push_back(nextDigit); }
			nextDigit = !nextDigit;
		}
		initialString = newString;
		//for (int j = 0; j < initialString.size(); j++) { cout << initialString[j]; }
		//int a; cin >> a;
	}
	cout << "Enter anything to end.";
	cin >> length;
	return 0;
}


