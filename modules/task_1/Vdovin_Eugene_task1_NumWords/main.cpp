#include <iostream>
#include <mpi.h>
#include <string.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>
using namespace std;

int converter_in_number(const string &s)
{
	int len = s.length();
	int a = 0;
	int i = 0;
	for (i = 0; (i < len); i++)
		a = a * 10 + (s[i] - '0');
	return a;
}

int main(int argc, char *argv[])
{
	srand((int)time(0));

	string size = argv[1];
	string l = "qwertyuiop asdfghjkl zxcvbnm";
	int lenl = l.length();
	int len = converter_in_number(size);
	len = len * 3;
	char* s;

	int ProcNum, ProcRank;
	double stime = 0.0;
	double etime = 0.0;
	double stime_ = 0.0;
	double etime_ = 0.0;
	int nword = 0;
	int nword_ = 0;
	int nresword = 0;
	int residue = 0;

	MPI_Init(&argc, &argv);
	MPI_Status ss;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		s = new char[len + 1];
		int k = 0;
		while (k < len - 2)
		{
			s[k] = 'a';
			s[k + 1] = l[(int)(rand() % lenl)];
			s[k + 2] = 'a';
			k = k + 3;
		}
		stime = MPI_Wtime();
		for (int i = 0; i < len; i++)
			if (s[i] == ' ')
				nword++;
		nword++;
		etime = MPI_Wtime();
		cout << "1 process" << endl;
		cout << "Number words: " << nword << endl;
		cout << "Time: " << etime - stime << " sec" << endl;
		stime_ = MPI_Wtime();
		for (int i = 1; i < ProcNum; i++)
			MPI_Send(s + len / ProcNum * i, len / ProcNum, MPI_CHAR, i, 0, MPI_COMM_WORLD);
	}
	else
	{
		s = new char[len / ProcNum + 1];
		MPI_Recv(s, len / ProcNum, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &ss);
	}

	for (int i = 0; i < len / ProcNum; i++)
		if (s[i] == ' ')
			nword_++;

	MPI_Reduce(&nword_, &nresword, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		residue = len % ProcNum;
		if (residue)
			for (int i = len - residue; i < len; i++)
				if (s[i] == ' ')
					nresword++;
		nresword++;
		etime_ = MPI_Wtime();
		cout << ProcNum << " process" << endl;
		cout << "Number words: " << nresword << endl;
		cout << "Time: " << etime_ - stime_ << " sec" << endl;
	}

	MPI_Finalize();
	delete[] s;
	return 0;
}