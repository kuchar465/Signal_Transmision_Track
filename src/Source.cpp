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
#include <cstdlib>
#include <ctime>
#include "encoder.h"
#include "decoder.h"
#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)

const double BER_rate(const std::vector<int>& sygnal1, const std::vector<int>& sygnal2) {
	double BER = 0;
	for (size_t i = 0; i < sygnal1.size(); i++) {
		if (sygnal1[i] ^ sygnal2[i])
			BER++;
	}
	return BER / sygnal1.size();
}

template <typename T>
void save(std::string filename, double TS, T data) {
	std::ofstream plik;
	plik.open(filename);
	double tmp=0;
	plik << "x" << " " << "y" << std::endl;
	for (int i = 0; i < data.size(); i++) {
		plik << tmp << " " << data[i] << std::endl;
		tmp += TS;
	}
	plik.close();

}

std::vector<double> gen_noise(const int& numSamples) {
	const static int q = 15;
	const static float c1 = (1 << q) - 1;
	const static float c2 = ((int)(c1 / 3)) + 1;
	const static float c3 = 1.f / c1;

	float random = 0.f;
	float noise = 0.f;

	std::vector<double> out;

	for (int i = 0; i < numSamples; i++)
	{
		random = ((float)rand() / (float)(RAND_MAX + 1));
		noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
		out.push_back(noise);
	}
	return out;
}

void add_noise(std::vector<double>& sygnal, double alfa, const int& wielkoscSzumu = 100) {
	std::vector<double> szum = gen_noise(wielkoscSzumu);
	int nn = 0, NN = szum.size();

	for (int n = 0; n < sygnal.size(); n++)
	{
		sygnal[n] = (sygnal[n] * alfa) + (szum[nn] * (1.0 - alfa));
		nn++;

		if (nn >= NN)
			nn = 0;
	}
}


double signalgen(double t, double A, double fm, double fi)
{
	return A * sin(2 * M_PI * fm * t + fi);
}

int main() {
	{


		Encoder koder;
		Decoder dekoder;
		std::string word = "ALA";
		int bitamount = word.size() * 8;

		std::byte* ob = koder.S2BS(word, "littleEndian");
		std::vector<char> mt = koder.S2BS_CHAR_CUT(bitamount, ob);
		delete ob;
		std::vector<int> mtint = koder.S2BS_INT(mt);

		std::cout << std::dec << std::endl;
		std::vector<int> wynik = koder.encoderFullH(mtint);

		double tb = 0.1;
		int N = wynik.size() * (1 / tb);
		double f = N * pow(tb, -1);
		double f_fsk[] = { (N + 1.0) / tb, (N + 2.0) / tb };
		double A = 1.0;
		double A_ask[] = { 0.0, 1.0 };
		double fi = 1.0;
		double fi_psk[] = { 0.0, M_PI };
		double freq = float(1) / float(1000);

		save("signal.txt", 0.1, mtint);
		save("HEncoded.txt", 0.1, wynik);
		koder.print(wynik);

		std::vector<double> ask, fsk, psk, probki;

		for (size_t i = 0; i < wynik.size(); i++) {
			for (double j = i; j < i + 1.0; j += freq) {
				double sample = j;
				if (wynik[i] == 0) {
					ask.push_back(signalgen(sample, A_ask[0], f, fi));
					fsk.push_back(signalgen(sample, A, f_fsk[0], fi));
					psk.push_back(signalgen(sample, A, f, fi_psk[0]));
				}
				else {
					ask.push_back(signalgen(sample, A_ask[1], f, fi));
					fsk.push_back(signalgen(sample, A, f_fsk[1], fi));
					psk.push_back(signalgen(sample, A, f, fi_psk[1]));
				}
			}
		}

		save("ask.txt", freq, ask);
		save("fsk.txt", freq, fsk);
		save("psk.txt", freq, psk);

		bool widmapodstawa = false;

		if (widmapodstawa) {
			std::complex<double>* askwidmo = koder.dft(ask, ask.size());
			std::complex<double>* pskwidmo = koder.dft(psk, psk.size());
			std::complex<double>* fskwidmo = koder.dft(fsk, fsk.size());

			koder.zapisdft(askwidmo, ask.size(), "askWidmoPod.txt", 0.1);
			koder.zapisdft(pskwidmo, psk.size(), "pskWidmoPod.txt", 0.1);
			koder.zapisdft(fskwidmo, fsk.size(), "fskWidmoPod.txt", 0.1);
		}

		double alfa = 0.01;

		while (alfa < 1) {

		std::cout << "alfa : " << alfa << std::endl;
		


		std::vector<double>	asktmp = ask;
		std::vector<double>	psktmp = psk;
		std::vector<double>	fsktmp = fsk;

		add_noise(asktmp, alfa, 100);
		add_noise(psktmp, alfa, 100);
		add_noise(fsktmp, alfa, 100);

		
		/*if (alfa == 0.01) {
			std::complex<double>* askwidmo = koder.dft(asktmp, ask.size());
			std::complex<double>* pskwidmo = koder.dft(psktmp, psk.size());
			std::complex<double>* fskwidmo = koder.dft(fsk, fsk.size());

			koder.zapisdft(askwidmo, ask.size(), "askWidmo1.txt", 0.1);
			koder.zapisdft(pskwidmo, psk.size(), "pskWidmo1.txt", 0.1);
			koder.zapisdft(fskwidmo, fsk.size(), "fskWidmo1.txt", 0.1);

			delete askwidmo;
			delete pskwidmo;
			delete fskwidmo;
		}*/

		/*if (alfa == 0.39) {
			std::complex<double>* askwidmo = koder.dft(asktmp, ask.size());
			std::complex<double>* pskwidmo = koder.dft(psktmp, psk.size());
			std::complex<double>* fskwidmo = koder.dft(fsktmp, fsk.size());

			koder.zapisdft(askwidmo, ask.size(), "askWidmo2.txt", 0.1);
			koder.zapisdft(pskwidmo, psk.size(), "pskWidmo2.txt", 0.1);
			koder.zapisdft(fskwidmo, fsk.size(), "fskWidmo2.txt", 0.1);

			delete askwidmo;
			delete pskwidmo;
			delete fskwidmo;
		}*/

		/*if (alfa == 0.5) {
			std::complex<double>* askwidmo = koder.dft(asktmp, ask.size());
			std::complex<double>* pskwidmo = koder.dft(psktmp, psk.size());
			std::complex<double>* fskwidmo = koder.dft(fsktmp, fsk.size());

			koder.zapisdft(askwidmo, ask.size(), "askWidmo3.txt", 0.1);
			koder.zapisdft(pskwidmo, psk.size(), "pskWidmo3.txt", 0.1);
			koder.zapisdft(fskwidmo, fsk.size(), "fskWidmo3.txt", 0.1);

			delete askwidmo;
			delete pskwidmo;
			delete fskwidmo;
		}*/
		// ber non zero for alfa =0.01, 0.39, 0.5 

		std::vector<double> s_ask, s_fsk1, s_fsk2, s_psk, x_ask, x_fsk1, x_fsk2, x_psk;

		for (size_t i = 0, k = 0; i < wynik.size(); i++) {
			for (double j = i; j < i + 1.0; j += freq) {
				double sample = j;
				s_ask.push_back(signalgen(sample, A_ask[1], f, fi));
				x_ask.push_back(s_ask[k] * asktmp[k]);

				s_fsk1.push_back(signalgen(sample, A, f_fsk[0], fi));
				x_fsk1.push_back(s_fsk1[k] * fsktmp[k]);
				s_fsk2.push_back(signalgen(sample, A, f_fsk[1], fi));
				x_fsk2.push_back(s_fsk2[k] * fsktmp[k]);

				s_psk.push_back(signalgen(sample, A, f, fi_psk[1]));
				x_psk.push_back(s_psk[k] * psktmp[k]);

				k++;
			}
		}



		std::vector<double> ccalka1 = dekoder.calka(wynik.size(), 1000, x_ask);
		std::vector<double> ccalka2 = dekoder.calka(wynik.size(), 1000, x_fsk1);
		std::vector<double> ccalka3 = dekoder.calka(wynik.size(), 1000, x_fsk2);
		std::vector<double> ccalka4 = dekoder.calka(wynik.size(), 1000, x_psk);

		std::vector<double> sc2;
		for (int i = 0; i < wynik.size(); i++) {
			sc2.push_back(ccalka3[i] - ccalka2[i]);
		}

		std::vector<int> wASK2 = dekoder.calkaToBit(ccalka1);
		std::vector<int> wPSK2 = dekoder.calkaToBit(ccalka4);
		std::vector<int> wFSK2 = dekoder.calkaToBit(sc2);

		std::vector<int> wwASK2 = dekoder.decoderFullH7(wASK2);
		std::vector<int> wwPSK2 = dekoder.decoderFullH7(wPSK2);
		std::vector<int> wwFSK2 = dekoder.decoderFullH7(wFSK2);

		std::vector<char> wordASK2 = dekoder.bitToWord(wwASK2);
		std::vector<char> wordPSK2 = dekoder.bitToWord(wwPSK2);
		std::vector<char> wordFSK2 = dekoder.bitToWord(wwFSK2);

		for (int i = 0; i < wordASK2.size(); i++) {
			std::cout << wordASK2[i];
		}
		std::cout << std::endl;
		for (int i = 0; i < wordPSK2.size(); i++) {
			std::cout << wordPSK2[i];
		}
		std::cout << std::endl;
		for (int i = 0; i < wordFSK2.size(); i++) {
			std::cout << wordFSK2[i];
		}
		std::cout << std::endl;
		for (int i = 0; i < mtint.size(); i++) {
			std::cout << mtint[i] << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < wwASK2.size(); i++) {
			std::cout << wwASK2[i] << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < wwPSK2.size(); i++) {
			std::cout << wwPSK2[i] << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < wwFSK2.size(); i++) {
			std::cout << wwFSK2[i] << " ";
		}
		std::cout << std::endl;

		std::cout << std::endl << "BER ask: " << BER_rate(mtint, wwASK2) << std::endl;
		std::cout << std::endl << "BER psk: " << BER_rate(mtint, wwPSK2) << std::endl;
		std::cout << std::endl << "BER fsk: " << BER_rate(mtint, wwFSK2) << std::endl;
		alfa += 0.01;
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		}
		
	}

	int tmp = 12>>2;
	std::cout << "---------" << std::endl;
	std::cout << tmp << std::endl;
	system("PAUSE");
	_CrtDumpMemoryLeaks();
	return 0;
}

// signal mt: 0 1 0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 1 0 0 1 1 1 1 0 1 0 0 1 0 1 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 1 0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 1 0 0 1 1 0 0 0 1 0 0 0 0 0 1
// encoded mt: 1 0 0 1 1 0 0 1 1 0 1 0 0 1 0 1 0 0 1 0 1 1 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 1 1 1 0 0 1 1 0 0 0 1 1 0 0 1 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 1 1 0 0 1 1 0 1 0 0 1 1 0 0 1 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 1 1 0 0 1 1 0 1 0 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0 0 1 0 0 1 1 0 0 1 1 0 1 0 0 1