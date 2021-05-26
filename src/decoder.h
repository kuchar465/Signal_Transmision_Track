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
#ifndef Decoder_h
#define Decoder_h

class Decoder {
	class hammingMatrixes {
	public:
		std::vector<std::vector<int>> G;
		std::vector<std::vector<int>> H;
		hammingMatrixes();
	};
private:
	std::vector<char> generateCLK(int bitamount);
	
	double signalgen(double t, double A = 1.0, double fm = 8.0, double fi = 0);
public:
	double findH(std::vector<double> tab);
	std::vector<double> prep1(std::vector<double> sig, double A, double fm, double fi, std::string mod, double freq, int bitamount);
	std::vector<std::vector<double>> prep2(std::vector<double> sig, double A, double freq, int bitamount);
	std::vector<int> decodeTTL(const std::vector<int> TTL);
	std::vector<int> decodeNRZI(const std::vector<int> nrzi);
	std::vector<int> decodeBAMI(const std::vector<int> bami);
	std::vector<int> decodeManchester(const std::vector<int> manchester, int bitamount);
	std::vector<int> decoding(const std::vector<int>& data);
	std::vector<std::vector<int>> CreatePackets7(const std::vector<int>& data);
	std::vector<int> decoderFullH7(const std::vector<int>& data);
	std::vector<char> bitToWord(const std::vector<int>& data);
	std::vector<double> calka(int nrOfBits, int TB, std::vector<double> mul);
	std::vector<int> calkaToBit(std::vector<double> tab);
	std::vector<int> calkaToBit(std::vector<double> tab, double h);
	void bitN(std::vector<int>& mt, int n);
	void print(const std::vector<int>& data);
};

#endif