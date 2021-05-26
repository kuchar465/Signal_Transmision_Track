#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <cstddef>
#include <bitset>
#include <vector> 
#include <numeric>
#include <list>
#include <complex>
#ifndef Encoder_h
#define Encoder_h

class Encoder {
	class hammingMatrixes {
	public:
		std::vector<std::vector<int>> G;
		std::vector<std::vector<int>> H;
		hammingMatrixes();
	};
private:

	std::vector<int> encoding(const std::vector<int>& data, int sec = 0);
	std::vector<std::vector<int>> CreatePackets(const std::vector<int>& data);
	void bitN(std::vector<int>& mt, int n);
	std::vector<char> generateCLK(int bitamount);
	double signalgen(double t, double A = 1.0, double fm = 8.0, double fi = 0);
public:
	std::complex<double>* dft(std::vector<double> inputs, int n);
	std::complex<double>* dft(std::vector<int> inputs, int n);
	std::vector<int> NRZI(const std::vector<int> signal);
	std::vector<int> TTL(const std::vector<int> signal);
	std::vector<int> manchester(const std::vector<int> signal, int bitamount);
	std::vector<int> BAMI(const std::vector<int> signal);
	std::byte* S2BS(std::string table, std::string sw);
	std::vector<int> encoderFullH(const std::vector<int>& mt);
	std::vector<char> S2BS_CHAR_CUT(int bitamount, const std::byte *ob);
	std::vector<int> S2BS_INT(const std::vector<char>& data);
	void print(const std::vector<int>& data);
	std::vector<double> ASK(const std::vector<int> signal, int cz, double freq, double A1, double A2, double fm);
	std::vector<double> PSK(const std::vector<int> signal, int cz, double freq, double A, double fm, double fi0, double fi1);
	std::vector<double> FSK(const std::vector<int> signal, int cz, double freq, double A, double f0, double f1);
	void zapisdft(std::complex<double>* inputs, int n, std::string nazwa, double freq);
};


#endif