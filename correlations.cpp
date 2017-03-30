#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>


using namespace std;

//assumes power > 0 and that map is bijective

struct seqs_with_lens
{
	vector<vector<bool>> seqs;
	vector<int> lengths;
};

seqs_with_lens expand(vector<bool> seq)
{
	vector<vector<bool>> results;
	vector<int> lengths;
	vector<bool> seq1;
	vector<bool> seq2;
	int length1 = 0;
	bool b = true;
	for (unsigned int k = 0; k < seq.size(); k++)
	{
		seq1.push_back(b);
		seq2.push_back(!(b));
		length1 = length1 + 1 + b;
		if (seq[k])
		{
			seq1.push_back(b);
			seq2.push_back(!(b));
			length1 = length1 + 1 + b;
		}
		b = !b;
	}
	int length2 = 3 * seq1.size() - length1;

	vector<bool> seq1_01 = seq1; seq1_01.push_back(false);
	vector<bool> seq1_10 = seq1; seq1_10.insert(seq1_10.begin(), false);
	vector<bool> seq1_11 = seq1_10; seq1_11.push_back(false);
	vector<bool> seq2_01 = seq2; seq2_01.push_back(false);
	vector<bool> seq2_10 = seq2; seq2_10.insert(seq2_10.begin(), false);
	vector<bool> seq2_11 = seq2_10; seq2_11.push_back(false);
	if (seq1[0])
	{
		if (seq1[seq1.size() - 1])
		{
			results.push_back(seq1);
			results.push_back(seq1_01);
			results.push_back(seq1_10);
			results.push_back(seq1_11);
			lengths.push_back(length1);
			lengths.push_back(length1 + 1);
			lengths.push_back(length1 + 1);
			lengths.push_back(length1 + 2);

			results.push_back(seq2_11);
			lengths.push_back(length2 + 2);
		}
		else
		{
			results.push_back(seq1_01);
			results.push_back(seq1_11);
			lengths.push_back(length1 + 1);
			lengths.push_back(length1 + 2);

			results.push_back(seq2_10);
			results.push_back(seq2_11);
			lengths.push_back(length2 + 1);
			lengths.push_back(length2 + 2);
		}
	}
	else
	{
		if (seq1[seq1.size() - 1])
		{
			results.push_back(seq1_10);
			results.push_back(seq1_11);
			lengths.push_back(length1 + 1);
			lengths.push_back(length1 + 2);

			results.push_back(seq2_01);
			results.push_back(seq2_11);
			lengths.push_back(length2 + 1);
			lengths.push_back(length2 + 2);
		}
		else
		{
			results.push_back(seq1_11);
			lengths.push_back(length1 + 2);

			results.push_back(seq2);
			results.push_back(seq2_01);
			results.push_back(seq2_10);
			results.push_back(seq2_11);
			lengths.push_back(length2);
			lengths.push_back(length2 + 1);
			lengths.push_back(length2 + 1);
			lengths.push_back(length2 + 2);
		}
	}
	seqs_with_lens data;
	data.seqs = results;
	data.lengths = lengths;
	return data;
}


seqs_with_lens expand2(vector<bool> seq) //expand twice
{
	seqs_with_lens data1 = expand(seq);

	vector<vector<bool>> seqs2;
	vector<int> lengths2;
	for (int i = 0; i < data1.seqs.size(); i++)
	{
		seqs_with_lens data_i = expand(data1.seqs[i]);
		vector<vector<bool>> seqs_i = data_i.seqs;
		vector<int> lengths_i = data_i.lengths;
		for (int j = 0; j < seqs_i.size(); j++)
		{
			seqs2.push_back(seqs_i[j]);
			lengths2.push_back(lengths_i[j]);
		}
	}
	seqs_with_lens data2;
	data2.seqs = seqs2;
	data2.lengths = lengths2;
	return data2;
}




// exact for K(1,2)
int main1(int argc, char* argv[])
{
	//srand(time(NULL));
	vector<unsigned int> powers2 = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };
	vector<unsigned int> powers3 = { 1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049,177147, 531441, 1594323, 4782969, 14348907, 43046721, 129140163, 387420489, 1162261467, 3486784401 };

	unsigned int M = 1;
	unsigned int N = 2;

	int size_1 = 501;

	vector<vector<vector<bool>>> sequences(size_1); // don't store full sequences for length > 221 to save space
	// how much space even is a bool vector? lol
	sequences[0] = { {} };
	sequences[1] = { {false} };
	sequences[2] = { {false, false} };
	//sequences[2] = { { false, false },{ true } };
	//sequences[3] = { {false, false, false} }; //incomplete; to be completed in recursive construction

	vector<vector<bool>> numerators(size_1); // ints are definitely bigger than necessary, but we keep it anyways
	numerators[0] = { false };
	numerators[1] = { false };
	numerators[2] = {  true };
	//numerators[3] = { 1 };

	vector<vector<char>> denominators(size_1); // actual denominators are 3 ^ (int stored here)
	denominators[0] = { 0 };
	denominators[1] = { 0 };
	denominators[2] = { 1 };
	//denominators[3] = { 1 };

	vector<vector<bool>> start_equals_end(80);

	for (unsigned int i = 1; i <= 320; i++)
	{

		cout << i << "\n";
		for (unsigned int j = 0; j < sequences[i].size(); j++)
		{
			//if (i == 2 && j == 0) continue; //small case
			vector<bool> seq = sequences[i][j];
			bool new_num = numerators[i][j];
			char new_den = char(denominators[i][j] + 1);
			//expand seq both ways
			vector<bool> seq1;
			vector<bool> seq2;
			int length1 = 0;
			int length2 = 0;
			bool b = true;
			for (unsigned int k = 0; k < seq.size(); k++)
			{
				seq1.push_back(b);
				seq2.push_back(!(b));
				length1 = length1 + 1 + b;
				if (seq[k])
				{
					seq1.push_back(b);
					seq2.push_back(!(b));
					length1 = length1 + 1 + b;
				}
				b = !b;
			}
			length2 = 3 * i - length1;
			vector<bool> seq1_01 = seq1; seq1_01.push_back(false);
			vector<bool> seq1_10 = seq1; seq1_10.insert(seq1_10.begin(), false);
			vector<bool> seq1_11 = seq1_10; seq1_11.push_back(false);
			vector<bool> seq2_01 = seq2; seq2_01.push_back(false);
			vector<bool> seq2_10 = seq2; seq2_10.insert(seq2_10.begin(), false);
			vector<bool> seq2_11 = seq2_10; seq2_11.push_back(false);


			// now determine what augments of seq1, seq2 should be added to the data
			if (seq1[0])
			{
				if (seq1[i - 1])
				{
					sequences[length1].push_back(seq1);
					numerators[length1].push_back(new_num);
					denominators[length1].push_back(new_den);
					sequences[length1 + 1].push_back(seq1_01);
					numerators[length1 + 1].push_back(new_num);
					denominators[length1 + 1].push_back(new_den);
					sequences[length1 + 1].push_back(seq1_10);
					numerators[length1 + 1].push_back(new_num);
					denominators[length1 + 1].push_back(new_den);
					sequences[length1 + 2].push_back(seq1_11);
					numerators[length1 + 2].push_back(new_num);
					denominators[length1 + 2].push_back(new_den);

					sequences[length2 + 2].push_back(seq2_11);
					numerators[length2 + 2].push_back(new_num);
					denominators[length2 + 2].push_back(new_den);
				}
				else
				{
					sequences[length1 + 1].push_back(seq1_01);
					numerators[length1 + 1].push_back(new_num);
					denominators[length1 + 1].push_back(new_den);
					sequences[length1 + 2].push_back(seq1_11);
					numerators[length2 + 2].push_back(new_num);
					denominators[length1 + 2].push_back(new_den);

					sequences[length2 + 1].push_back(seq2_10);
					numerators[length2 + 1].push_back(new_num);
					denominators[length2 + 1].push_back(new_den);
					sequences[length2 + 2].push_back(seq2_11);
					numerators[length2 + 2].push_back(new_num);
					denominators[length2 + 2].push_back(new_den);
				}
			}
			else
			{
				if (seq1[i - 1])
				{
					sequences[length1 + 1].push_back(seq1_10);
					numerators[length1 + 1].push_back(new_num);
					denominators[length1 + 1].push_back(new_den);
					sequences[length1 + 2].push_back(seq1_11);
					numerators[length1 + 2].push_back(new_num);
					denominators[length1 + 2].push_back(new_den);

					sequences[length2 + 1].push_back(seq2_01);
					numerators[length2 + 1].push_back(new_num);
					denominators[length2 + 1].push_back(new_den);
					sequences[length2 + 2].push_back(seq2_11);
					numerators[length2 + 2].push_back(new_num);
					denominators[length2 + 2].push_back(new_den);
				}
				else
				{
					sequences[length1 + 2].push_back(seq1_11);
					numerators[length1 + 2].push_back(new_num);
					denominators[length1 + 2].push_back(new_den);

					sequences[length2].push_back(seq2);
					numerators[length2].push_back(new_num);
					denominators[length2].push_back(new_den);
					sequences[length2 + 1].push_back(seq2_01);
					numerators[length2 + 1].push_back(new_num);
					denominators[length2 + 1].push_back(new_den);
					sequences[length2 + 1].push_back(seq2_10);
					numerators[length2 + 1].push_back(new_num);
					denominators[length2 + 1].push_back(new_den);
					sequences[length2 + 2].push_back(seq2_11);
					numerators[length2 + 2].push_back(new_num);
					denominators[length2 + 2].push_back(new_den);
				}
			}
		}

		//cout << i;
	}
	/*cout << "number of possible subsequences by length << \n";
	for (int i = 0; i < sequences.size(); i++)
	{
		cout << sequences[i].size() << " ";
	}
	cout << "\n";
	int a; cin >> a;*/

	vector<int> num_list;
	vector<int> den_list;
	for (unsigned int d = 2; d <= 501; d++)
	{

		int num = 0;
		int den = 1;
		int num2 = 0;
		int den2 = 0;
		for (unsigned int j = 0; j < sequences[d].size(); j++)
		{
			//cout << numerators[d][j] << " " << denominators[d][j] << " " << sequences[d][j].size() << "     ";
			if (sequences[d][j].size() % 2 == 1)
			{
				int num_j = 1 + int(numerators[d][j]);
				int den_j = denominators[d][j];
				int max_den = max(den, den_j);
				num = num * powers3[max_den - den] + num_j * powers3[max_den - den_j];
				den = max_den;
			}
			if (sequences[d][j].size() % 2 == 0)
			{
				int num_j = 1 + int(numerators[d][j]);
				int den_j = denominators[d][j];
				int max_den = max(den2, den_j);
				num2 = num2 * powers3[max_den - den2] + num_j * powers3[max_den - den_j];
				den2 = max_den;
			}
		}
		cout << "correlation frequency for distance " << d - 1 << " is " << num << " / " << powers3[den] << " = " << num * 1.0 / powers3[den] << "\n";
		cout << "complement check: " << num2 << " / " << den2 << "\n";
		num_list.push_back(num);
		den_list.push_back(powers3[den]);

	}
	int a2; cin >> a2;
	for (int i = 0; i < num_list.size(); i++)
	{
		cout << num_list[i] << " / " << den_list[i] << " , ";

	}
	int a3; cin >> a3;


	//ad-hoc way to compute correlation frequency exactly 782 apart
	// remember a factor of 9
	/*
	int num = 0;
	int den = 1;
	int num2 = 0;
	int den2 = 1;


    for (int h = 0; h <= 51; h++)
	{
		int g = h / 2;
		int i = 0;
		if (h % 2 == 1)
		{
			i = 343 - g;
		}
		else
		{
			i = 344 + g;
		}/*
	/*for (int i = 40; i <= 60; i++)
	{/*
		cout << i << "\n";
		int count_i = 0;
		for (int j = 0; j < sequences[i].size(); j++)
		{
			vector<bool> seq = sequences[i][j];
			int num_j = 1 + int(numerators[i][j]);
			//if (rand() < 21474836) cout << num_j << " ";
			int den_j = denominators[i][j];

			//seqs_with_lens data = expand2(seq);
			seqs_with_lens data = expand2(seq);


			for (int k = 0; k < data.lengths.size(); k++)
			{
				if (data.lengths[k] == 783 && (data.seqs[k].size() % 2 == 1))
				{
					int max_den = max(den, den_j);
					num = num * powers3[max_den - den] + num_j * powers3[max_den - den_j];
					den = max_den;


				}
				if (data.lengths[k] == 783 && (data.seqs[k].size() % 2 == 0))
				{
					int max_den = max(den2, den_j);
					num2 = num2 * powers3[max_den - den2] + num_j * powers3[max_den - den_j];
					den2 = max_den;
				}
			}
			
		}
		//cout << "count_i equals " << count_i << " ";
		cout << "current estimate for distance " << 782 << " is " << num << " / " << powers3[den + 2] << " = " << num * 1.0 / powers3[den + 2] << "\n";
		cout << num2 * 1.0 / powers3[den2 + 2] << "\n";
	}
	cout << "correlation frequency for distance " << 782 << " is " << num << " / " << powers3[den + 2] << " = " << num * 1.0 / powers3[den + 2] << "\n";
	cout << num2 * 1.0 / powers3[den2 + 2] << "\n";
	int a; cin >> a; /**/

	return 0;
}


// approximate for all (M,N,d)
// obsoleted tbh
int main2(int argc, char* argv[])
{
	int M = 1;
	int N = 2;
	long long length = 10000000;
	vector<bool> kolakoski(length, false);
	if (M == 1) kolakoski[1] = true;
	//for (int i = 0; i < M; i++) kolakoski[i] = false;
	long long index = 1;
	long long write_index = M;
	bool b = true;
	while (write_index < length - max(M, N))
	{
		bool c = kolakoski[index];

		// write the bool b, N times if c and M times if !c
		if (c)
		{
			for (int i = 0; i < N; i++)
			{
				kolakoski[write_index] = b;
				write_index++;
			}
		}
		else
		{
			for (int i = 0; i < M; i++)
			{
				kolakoski[write_index] = b;
				write_index++;
			}
		}
		b = !b;
		index++;
	}

	/*for (int i = 0; i < kolakoski.size(); i++)
	{
		cout << kolakoski[i];
	}
	cout << "\n";*/

	int max_dist = 5000;
	vector<double> frequencies(max_dist, 0);
	for (int d = 1; d < max_dist + 1; d++)
		//for (int d = 781; d <= 783; d++)
	{
		long long match = 0;
		long long total = 0;
		for (long long i = 0; i < length - d; i++)
		{
			if (kolakoski[i] == kolakoski[i + d])
			{
				match++;
			}
		}
		total = length - d;
		cout << "\ncorrelation at distance " << d << " is estimated as " << match * 1.0 / total;
		frequencies[d] = match*1.0 / total;
	}
	int a1; cin >> a1;
	cout << "\nIndicator for being greater than 0.5: ";
	for (int d = 1; d < max_dist; d++)
	{
		cout << (frequencies[d] > 0.5);
	}
	cout << "\n estimated decimal values: ";
	cin >> a1;
	cout.precision(17);
	for (int d = 1; d < max_dist; d++)
	{
		cout << frequencies[d] << ", ";
	}
	int a; cin >> a;

}

// exact for M = (1,2), small d, but the bound on d is not really easy to calculate; w
// we regard its output as starting exact, then approximate for an unknown cutoff

// ad-hoc method to compute Eulerian cycle on 2-regular digraph. assumes connected
// input: two vector<int>of same length; 
// vertex set is common index set; each vertex has out-degree 2, and those 2 edges are given by v0, v1
// output is a sequence of 0s and 1s only s.t. the path starting from vertex 0 and following either
// edge from v1 or edge from v0 according to the output is an Eulerian cycle
vector<int> EulerianCycle(vector<int> v0, vector<int> v1)
{
	int n = v0.size();
	vector<int> path(0);
	//vector<int> temp_path(0);
	vector<bool> used0(n, false);
	vector<bool> used1(n, false);
	vector<bool> iccc(n, false); // in current connected component
	vector<int> ccc(1, 0); // vector which contains the vertices in current connected component
	//int next_v = 0; // least vertex (index) which doesn't have both out-edges in current path
	vector<int> v_to_p(n, 0);
	vector<int> next_edge(2 * n, -1);
	vector<int> prev_edge(2 * n, -1);

	iccc[0] = true;
	int ni = 0; // next index to test in current connected component; all earlier 

	// just making a single cycle for god's sake
	int v = 0;
	int e = rand() % 2;
	int edge0 = e;
	if (e == 0) v = v0[v]; 
	if (e == 1) v = v1[v];
	int edge_prev = e;
	while (v != 0)
	{
		if (!iccc[v]) { ccc.push_back(v); iccc[v] = true; }
		e = rand() % 2;
		if (next_edge[2 * v + 0] >= 0) e = 1;
		if (next_edge[2 * v + 1] >= 0) e = 0;
		int edge_next = 2 * v + e;
		next_edge[edge_prev] = edge_next;
		prev_edge[edge_next] = edge_prev;
		// shift variables
		if (e == 0) v = v0[v];
		if (e == 1) v = v1[v];
		edge_prev = edge_next;
	}
	// complete the loop; last value of edge_prev is the last sege
	next_edge[edge_prev] = edge0;
	prev_edge[edge0] = edge_prev;

	/*for (int i = 0; i < 2 * n; i++)
	{
		cout << next_edge[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < 2 * n; i++)
	{
		cout << prev_edge[i] << " ";
	}*/
	cout << endl;



	while (ni < ccc.size() )
	{
		//cout << "path size (1) = " << path.size() << " ";
		// find vtx (index) to extend cycle by splicing
		while (ni < ccc.size() && (next_edge[2*ccc[ni]+0]>=0 && next_edge[2*ccc[ni]+1]>=0)) ni++; 

		
		if (ni >= ccc.size()) {  
			if (ccc.size() == n) continue;
			cout << "\ndata did not make a connected graph\n"; return path;
		}
		int v = ccc[ni];
		//cout << v << endl;
		int vcopy = v;
		//find the edges to displace. say our new cycle is edge3, edge4, ... edge10. we set ne[e2] = e3; pe[e3] = e2; ne[e10] = e1; pe[e1] = e2
		int edge1 = 0;
		if (next_edge[2 * v + 0] >= 0) edge1 = 2 * v;
		if (next_edge[2 * v + 1] >= 0) edge1 = 2 * v + 1;
		int edge2 = prev_edge[edge1];

		int edge_prev = edge2;
		bool are_we_there_yet = false;
		
		while (v != vcopy || !are_we_there_yet)
		{
			if (!iccc[v]) { ccc.push_back(v); iccc[v] = true; }
			are_we_there_yet = true;
			int e = rand() % 2;
			if (next_edge[2 * v + e] >= 0) e = 1 - e;
			next_edge[edge_prev] = 2 * v + e;
			prev_edge[2 * v + e] = edge_prev;
			// shift variables
			edge_prev = 2 * v + e;
			if (e == 0) v = v0[v];
			if (e == 1) v = v1[v];
		}
		// last value of edge_prev is the last edge in the cycle; now rejoin
		next_edge[edge_prev] = edge1;
		prev_edge[edge1] = edge_prev;




	}
	if (ccc.size() != n) cout << "\ndata did not make a connected graph\n";// return vector<int>();



	int edge = 0;
	do
	{
		path.push_back(edge % 2);
		edge = next_edge[edge];
	} while (edge != 0);

	return path;
}

//input: M < N
//returns E_{ m,n }(s, t)
vector<bool> expand_run(vector<bool> s, int t, int M, int N)
{
	t = t % 2;
	vector<bool> r(0);
	bool b = bool(t);
	for (int i = 0; i < s.size(); i++)
	{
		for (int j = 0; j < M; j++)
			r.push_back(b);
		if (s[i])
		{
			for (int j = 0; j < N - M; j++)
				r.push_back(b);
		}
		b = !b;
	}
	return r;
}


// implements FFT mod 10*3^18+1; assumes that argument length is 2*power of 3
// the pair w ^ (|v| / 3), w^(2 |v| / 3) is always equal to 2094712797, 1779492093 mod 10*3^18+1 in one of the two orders.
// note 32^ (2*3^17) == 2094712797
// option flag:
// If w ^ (|v| / 3) ==  2094712797, then flag must be set to true
// if == 1779492093, flag must be set to false
// the recursive calls always have the same flag


// 32 is a root of unith of order (2*3^18), and that's the maximum length we will try FFT on
// 3874204891 = 10*3^18+1


vector<unsigned int> FFT(vector<unsigned int> v, unsigned long long w, bool flag)
{
	unsigned long long a = 0;
	unsigned long long b = 0;
	unsigned long long c = 3874204891;
	if (flag) { a = 2094712797; b = 1779492093; }
	if (!flag) { b = 2094712797; a = 1779492093; }
	if (v.size() == 2)
	{
		return vector<unsigned int>({ unsigned int((unsigned long long(v[0]) + unsigned long long(v[1])) % c), unsigned int((unsigned long long(v[0]) + unsigned long long(v[1]) * w) % c) });
		//vector<unsigned int> u(0);
		// this is what modular arithmetic modulo a large non-power of 2 looks like on a pleb laptop
	}
	unsigned int n = v.size();
	//otherwise, n is a multiple of 3


	vector<unsigned int> v1(n / 3, 0);
	vector<unsigned int> v2(n / 3, 0);
	vector<unsigned int> v3(n / 3, 0);
	unsigned long long pow_w_i = 1;
	unsigned long long pow_w_2i = 1;
	for (unsigned int i = 0; i < n / 3; i++)
	{
		v1[i] = int(((long long(v[i]) + long long(v[i + n / 3]))%c + long long(v[i + 2 * n / 3])) % c);
		v2[i] = int(((((long long(v[i]) + (long long(v[i + n / 3]) * a) % c + (long long(v[i + 2 * n / 3]) * b) % c))%c) * pow_w_i ) % c);
		v3[i] = int(((((long long(v[i]) + (long long(v[i + n / 3]) * b) % c + (long long(v[i + 2 * n / 3]) * a) % c))%c) * pow_w_2i) % c);
		pow_w_i = (pow_w_i * w) % c;
		pow_w_2i = (pow_w_i * pow_w_i) % c;
	}
	unsigned long long w3 = (w * ((w * w) % c)) % c;

	v1 = FFT(v1, w3, flag);
	//v1 = {}; //free up space
	v2 = FFT(v2, w3, flag);
	//v2 = {};
	v3 = FFT(v3, w3, flag);
	//v3 = {};

	// overwrite v. lol does this save memory?
	for (unsigned int i = 0; i < n / 3; i++)
	{
		v[3 * i + 0] = v1[i];
		v[3 * i + 1] = v2[i];
		v[3 * i + 2] = v3[i];
	}
	return v; 
}

// We want to compute the convolution of universal_sequence with itself.
// Ideally, we use DFT for vectors of exactly the correct length
// unfortunately, FFT is annoying to implement for many different lengths
// therefore, we will use FFT for length 2*3^18 or 2*3^k
// we pad with zeros to get this length
// When we compute the convolution of the padded vectors, we get a small error
// we will correct this error for the first 10^6 terms (and ignore the rest; don't even output them)
// fortunately, this error term also has the form of a convolution
// this convolution involves taking slices of universal_sequence of length 1000000

int main(int argc, char* argv[])
{

	///////////// main options
	int M = 2;
	int N = 5;
	char filepath[] = "cfp_25_run2.txt";
	// options implemented: (m,n) in { (1,2), (1,4), (1,6), (2,3), (2,5), (3,4) }
	int L;
	if (M + N == 3) L = 18;
	if (M + N == 5) L = 12;
	if (M + N == 7) L = 10;
	
	int D = 1000000; // max autocorrelation distance calculated. should not be more than 3^13 = 1594323

	unsigned long long C = 3874204891;
	srand(time(NULL));

	// precomputation
	vector<unsigned int> powers2 = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 };
	vector<unsigned int> powers3 = { 1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049,177147, 531441, 1594323, 4782969, 14348907, 43046721, 129140163, 387420489, 1162261467, 3486784401 };



	
	vector<vector<int>> maps1(L+1); // maps1[k][s] for seqs s of length k equals C(1,s), except sequences are encoded as integer. 1,2,1,1,1 -> 00010 binary -> 2
	vector<vector<int>> maps2(L+1);
	maps1[1] = { 1,0 };
	maps2[1] = { 1,0 };
	
	for (int k = 1; k < L; k++)
	{
		maps1[k + 1] = vector<int>(powers2[k + 1]);
		maps2[k + 1] = vector<int>(powers2[k + 1]);
		// cout << k;
		for (int i = 0; i < powers2[k]; i++)
		{
			int j = i;
			for (int l = 0; l < M; l++) { j = maps1[k][j]; }
			maps1[k + 1][0 + 2 * i] = 1 + 2 * j;
			j = i;
			for (int l = 0; l < M; l++) { j = maps2[k][j]; }
			maps1[k + 1][1 + 2 * i] = 2 * j;
			j = i;
			for (int l = 0; l < N; l++) { j = maps1[k][j]; }
			maps2[k + 1][0 + 2 * i] = 1 + 2 * j;
			j = i;
			for (int l = 0; l < N; l++) { j = maps2[k][j]; }
			maps2[k + 1][1 + 2 * i] = 2 * j;
		}
	}
	cout << "\n finished computing C_{m,n} maps";

	vector<int> path = EulerianCycle(maps1[L], maps2[L]);
	cout << "\n finished finding Eulerian Cycle";

	// universal_sequence is a finite sequence of {0,1} such that replacing with {m,n}, then repeating infinitely 
	// gives a periodic approximation to K(m,n)
	vector<bool> universal_sequence(0);
	int v = 0;
	for (int i = 0; i < path.size(); i++)
	{
		vector<bool> block = { bool(path[i]) };
		int u = v;
		for (int k = 0; k < L; k++)
		{
			block = expand_run(block, u % 2, M, N);
			u = u / 2;
		}
		if (path[i] == 0) v = maps1[L][v];
		if (path[i] == 1) v = maps2[L][v];
		universal_sequence.insert(universal_sequence.end(), block.begin(), block.end());
	}
	cout << "\n finished generating universal sequence";
	

	vector<unsigned int> universal_sequence_i(774840978, 0);
	// universal_sequence_i is universal sequence, cast to integers, with zeros replaced with 10*3^18, THEN padded with zeros to length 2*3^18


	for (int i = 0; i < universal_sequence.size(); i++)
	{
		unsigned int a = int(universal_sequence[i]);
		if (a == 0) a = 3874204890;
		//
		universal_sequence_i[i] = a;
	}
	//for (int i = 0; i < universal_sequence.size(); i++) cout << universal_sequence[i];
	//free up space
	maps1 = { {} };
	maps2 = { {} };

	/*unsigned int M = 3;
	unsigned int N = 4;
	vector<unsigned int> universal_sequence_i(774840978 + M + N, 1);
	if (M == 1) universal_sequence_i[2] = 3874204890;
	int index = 0;
	int r_index = 0;
	unsigned int value = 1;
	while (index < 774840978)
	{
		int a = universal_sequence_i[r_index];
		if (a == 1)
		{
			for (int i = 0; i < M; i++)
			{
				universal_sequence_i[index] = value;
				index++;
			}
		}
		else
		{
			for (int i = 0; i < N; i++)
			{
				universal_sequence_i[index] = value;
				index++;
			}
		}
		r_index++;
		value = C - value;
	}
	universal_sequence_i.resize(774840978);
	cout << "finished generating sequence \n";
	/**/


	// computing the residue 
	// this involves taking slices and colvolving them
	// for efficiency, we convolve using FFT for vectors of length 2*3^13 
	// 2*3^13 = 3188646

	vector<unsigned int> slice1(3188646, 0);
	vector<unsigned int> slice2(3188646, 0);
	int s = 564950497;
	if (M + N == 5) s = 488281249;
	for (int i = 0; i < D; i++)
	{
		slice1[i] = universal_sequence_i[i];
		// 2 * 5^13 = 488281250
		// 2 * 7^10 = 564950498

		slice2[i] = universal_sequence_i[s - i];
	}

	// L = 18 -> omega = 32
	// L = 13 -> omega = 1058122295
	// L = 11 -> omega = 2745970270

	slice1 = FFT(slice1, 1058122295, true);
	slice2 = FFT(slice2, 1058122295, true);
	for (int i = 0; i < 3188646; i++)
	{
		slice1[i] = unsigned int((unsigned long long(slice1[i]) * unsigned long long(slice2[i])) % C);
	}
	slice1 = FFT(slice1, 2091729343, false);
	for (int i = 0; i < slice1.size(); i++)
	{
		slice1[i] = unsigned int((unsigned long long(slice1[i]) * 3874203676) % C);
	}

	// slice1 is the residue; slice1[d-1] is the residue when computing the correlation at distance d

	vector<unsigned int> fft = FFT(universal_sequence_i, 32, true);
	fft[0] = unsigned int((unsigned long long(fft[0]) * unsigned long long(fft[0])) % C);
	// pointwise multiplication of fft with reverse(fft)
	for (int i = 1; i <= (fft.size() / 2); i++)
	{
		fft[i] = unsigned int((unsigned long long(fft[i]) * unsigned long long(fft[fft.size() - i])) % C);
		fft[fft.size() - i] = fft[i];
	}
	// L = 18 -> omega^(-1) = 1573895737
	// L = 13 -> omega^(-1) = 2091729343
	// L = 11 -> omega^(-1) = 3380164432
	//vector<unsigned int> fft2 = FFT(fft, 1573895737, false);

	fft = FFT(fft, 1573895737, false);
	// multiply by N^(-1)
	// L = 18 -> N^(-1) = 3874204886
	// L = 13 -> N^(-1) = 3874203676
	// L = 11 -> N^(-1) = 3874193956

	for (int i = 0; i < fft.size(); i++)
	{
		fft[i] = unsigned int((unsigned long long(fft[i]) * 3874204886) % C);
	}
	// if M = 1, N = 2, there is no residue, so don't add it. Computing the fake residue is so fast that I decided to compute it unconditionally
	if (M + N > 3)
	{
		for (int i = 1; i <= D; i++)
		{
			fft[i] = unsigned int(unsigned long long(fft[i]) + unsigned long long(slice1[i - 1])) % C;
		}
	}



	ofstream myfile;
	myfile.open(filepath);
	for (int i = 0; i <= D; i++)
	{
		signed long long a = fft[i];
		if (a > (C / 2)) { a = a - C; }

		if (a < -300000000) { a = a + 420762405; } // I hate that I have to include this line


		myfile << a << ",";

		if (i < 1000)
			cout << a << ", "; // print first 1000 terms to cout for fun

		//currently, this loop turns some intended negative values x into x - 2^32 + C, which is bad
		// I just manually fix this in mathematica, where the data analysis is performed
	}
	
	myfile.close();
	cout << "\nfinished printing to file";


	/*vector<unsigned int> fft = FFT(universal_sequence_i, 2696225585, true);
	cout << " " << fft.size() << " ";
	// multiply by its reverse, not itself
	fft[0] = unsigned int((unsigned long long(fft[0]) * unsigned long long(fft[0])) % 3874204891);
	for (int i = 1; i <= (fft.size() / 2); i++)
	{
		fft[i] = unsigned int((unsigned long long(fft[i]) * unsigned long long(fft[fft.size() - i])) % 3874204891);
		fft[fft.size() - i] = fft[i];
	}
	//for (int i = 0; i < fft.size(); i++) fft[i] = unsigned int((unsigned long long(fft[i]) * unsigned long long(fft[i])) % 3874204891);

	vector<unsigned int> fft3 = FFT(fft, 58785636, false);
	// divide by (2*3^L) mod (10^3*18+1); L = 5 in this test
	for (int i = 0; i < fft3.size(); i++)
	{
		fft3[i] = unsigned int((unsigned long long(fft3[i]) * unsigned long long(3866233276)) % 3874204891);
	}
	cout << "results of fft test: \n";
	for (int i = 0; i < fft3.size(); i++)
	{
		int a = fft3[i]; if (a >(3874204891 / 2)) a = a - 3874204891;
		cout << a << ", ";
	}/**/

	cin.get();

}