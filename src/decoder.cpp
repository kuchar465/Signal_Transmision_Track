#include "decoder.h"

std::vector<char> Decoder::generateCLK(int bitamount)
{
	std::vector<char> CLK;
	for (int i = 0; i < bitamount * 2 + 1; i++) {
		if (i % 2 == 0) {
			CLK.push_back('0');
		}
		else {
			CLK.push_back('1');
		}
	}
	return CLK;
}

double Decoder::findH(std::vector<double> tab)
{
	double suma = 0;
	for (int i = 0; i < tab.size(); i++) {
		suma += tab[i];
	}
	return suma / double(tab.size());
}

std::vector<int> Decoder::decodeTTL(const std::vector<int> TTL)
{
	std::vector<int> decodedTTL;
	for (int i = 0; i < TTL.size(); i++) {
		if (TTL[i] == 1) {
			decodedTTL.push_back(1);
		}
		else {
			decodedTTL.push_back(0);
		}
	}
	return decodedTTL;
}

std::vector<int> Decoder::decodeNRZI(const std::vector<int> nrzi)
{
	int tmp, tmp2, tmp3;
	std::vector<int> decodedNRZI;
	for (int i = 0; i < nrzi.size(); i++) {
		if (i == 0) {
			tmp2 = nrzi[0];
			tmp3 = 1;
			if (tmp2 == -1) {
				tmp2 = 0;
			}
			tmp = tmp2 ^ tmp3;
			decodedNRZI.push_back(tmp);
		}
		else {
			tmp2 = nrzi[i];
			tmp3 = nrzi[i - 1];
			if (tmp2 == -1) {
				tmp2 = 0;
			}
			if (tmp3 == -1) {
				tmp3 = 0;
			}
			tmp = tmp2 ^ tmp3;
			decodedNRZI.push_back(tmp);
		}
	}
	return decodedNRZI;
}

std::vector<int> Decoder::decodeBAMI(const std::vector<int> bami)
{
	std::vector<int> decodedBAMI;
	for (int i = 0; i < bami.size(); i++) {
		if (bami[i] == -1 || bami[i] == 1) {
			decodedBAMI.push_back(1);
		}
		else {
			decodedBAMI.push_back(0);
		}
	}
	return decodedBAMI;
}

std::vector<int> Decoder::decodeManchester(const std::vector<int> manchester, int bitamount)
{
	int tmp = 0, tmp2, tmp3;
	std::vector<char> CLK = generateCLK(bitamount);
	std::vector<int> decodedManchester;
	for (int i = 0; i < CLK.size() - 1; i += 2) {
		tmp2 = manchester[i];
		tmp3 = manchester[i + 1];
		if (tmp2 == 1 && tmp3 == -1) {
			decodedManchester.push_back(1);
		}
		else if (tmp2 == -1 && tmp3 == 1) {
			decodedManchester.push_back(0);
		}
	}
	return decodedManchester;
}

std::vector<int> Decoder::decoding(const std::vector<int>& data)
{
	hammingMatrixes mat;
	std::vector<int> h, d;

	for (int i = 0; i < mat.H.size(); i++) {
		int wynik = 0;
		for (int j = 0; j < mat.H[i].size(); j++)
			wynik += data[j] * mat.H[i][j];
		wynik %= 2;
		h.push_back(wynik);
	}

	d.push_back(data[2]);
	d.push_back(data[4]);
	d.push_back(data[5]);
	d.push_back(data[6]);

	int index_korekty = 0;
	for (int i = 0; i < h.size(); i++)
		index_korekty += h[i] * std::pow(2, i);
	if (index_korekty == 0)
		return d;

	bitN(d, index_korekty - 1);

	return d;
}

std::vector<std::vector<int>> Decoder::CreatePackets7(const std::vector<int>& data)
{
	std::vector<std::vector<int>> matrix7;
	std::vector<int> vectors7;
	for (size_t i = 1; i <= data.size(); i++) {
		vectors7.push_back(data[i - 1]);
		if (i % 7 == 0) {
			matrix7.push_back(vectors7);
			vectors7.clear();
		}
	}
	return matrix7;
}

std::vector<int> Decoder::decoderFullH7(const std::vector<int>& data)
{
	std::vector<std::vector<int>> matrix=CreatePackets7(data);
	std::vector<int> decoded;
	for (size_t i = 0; i < matrix.size(); i++) {
		std::vector<int> multiplied = decoding(matrix[i]);
		for (int i = 0; i < multiplied.size(); i++) {
			decoded.push_back(multiplied[i]);
		}
	}
	return decoded;
}

std::vector<char> Decoder::bitToWord(const std::vector<int>& data)
{
	std::string tmpword;
	std::vector<char> word(data.size()/8);
	int tmp2=0;
	int tmp=0;
	int tmp3=0;
	std::bitset<8> binary=0;
	std::cout << std::endl;
	for (int i = 0; i < data.size(); i++) {
		tmp++;
		if (data[i] == 1) {
			binary[7-tmp2] = 1;
			tmp2++;
		}
		else {
			binary[7-tmp2] = 0;
			tmp2++;
		}
		if (tmp % 8 == 0 && tmp > 0) {
			word[(data.size() / 8)-1-tmp3]=(int(binary.to_ullong()));
			tmp2 = 0;
			tmp = 0;
			tmp3++;
		}

	}
	return word;
}

std::vector<double> Decoder::calka(int nrOfBits, int TB, std::vector<double> mul)
{
	std::vector<double> calka;
	for (int i = 0; i < nrOfBits; i++) {
		calka.push_back(0);
	}
	double a;
	double b;
	double H;
	for (int i = 1; i <= nrOfBits; i++) {
		a = double(i) - 1.0;
		b = double(i);
		H = abs(b - a) / double(TB);
		int ogrd = 0 + ((i - 1) * TB);
		int ogrg = i * TB;
		for (int k = ogrd; k < ogrg; k++) {
			calka[i - 1] = calka[i - 1] + mul[k];
		}
		calka[i - 1] *= H;
	}
	return calka;
}

double Decoder::signalgen(double t, double A, double fm, double fi)
{
	return A * sin(2 * M_PI * fm * t + fi);
}

std::vector<double> Decoder::prep1(std::vector<double> sig, double A, double fm, double fi, std::string mod, double freq, int bitamount)
{
	int tmp = 0;
	if (mod == "ASK") {
		std::vector<double> San;
		std::vector<double> xa;
		for (double k = 0; k < bitamount; k += freq, tmp++) {
			San.push_back(signalgen(k, A, fm, 1));
			xa.push_back(San[tmp] * sig[tmp]);
		}
		return xa;
	}
	else {
		std::vector<double> Sfin;
		std::vector<double> xp;
		for (double k = 0; k < bitamount; k += freq, tmp++) {
			Sfin.push_back(signalgen(k, A, fm, fi));
			xp.push_back(Sfin[tmp] * sig[tmp]);
		}
		return xp;
	}
}

std::vector<std::vector<double>> Decoder::prep2(std::vector<double> sig, double A, double freq, int bitamount)
{
	int tmp=0;
	std::vector<std::vector<double>> out;
	std::vector<double> Sf1;
	std::vector<double> Sf2;
	std::vector<double> xf1;
	std::vector<double> xf2;
	for (double k = 0; k < bitamount; k += freq, tmp++) {
		Sf1.push_back(signalgen(k, A, 1, 1));
		Sf2.push_back(signalgen(k, A, 2, 1));
		xf1.push_back(Sf1[tmp] * sig[tmp]);
		xf2.push_back(Sf2[tmp] * sig[tmp]);
	}
	out.push_back(xf1);
	out.push_back(xf2);
	return out;
}

std::vector<int> Decoder::calkaToBit(std::vector<double> tab)
{
	std::vector<int> out;
	double H = findH(tab)/4;
	for (int i = 0; i < tab.size(); i++) {
		if (tab[i] >= H) {
			out.push_back(1);
		}
		else {
			out.push_back(0);
		}
	}
	return out;
}

std::vector<int> Decoder::calkaToBit(std::vector<double> tab, double h)
{
	std::vector<int> out;
	double H = h;
	for (int i = 0; i < tab.size(); i++) {
		if (tab[i] >= H) {
			out.push_back(1);
		}
		else {
			out.push_back(0);
		}
	}
	return out;
}

void Decoder::bitN(std::vector<int>& mt, int n)
{
	if (mt[n] == 1) {
		mt[n] = 0;
	}
	else {
		mt[n] = 1;
	}
}

void Decoder::print(const std::vector<int>& data)
{
	for (int i = 0; i < data.size(); i++) {
		std::cout << data[i] << " ";
	}
	std::cout << std::endl;
}

Decoder::hammingMatrixes::hammingMatrixes()
{
	G.push_back(std::vector<int> {1, 1, 0, 1});
	G.push_back(std::vector<int> {1, 0, 1, 1});
	G.push_back(std::vector<int> {1, 0, 0, 0});
	G.push_back(std::vector<int> {0, 1, 1, 1});
	G.push_back(std::vector<int> {0, 1, 0, 0});
	G.push_back(std::vector<int> {0, 0, 1, 0});
	G.push_back(std::vector<int> {0, 0, 0, 1});

	H.push_back(std::vector<int> {1, 0, 1, 0, 1, 0, 1});
	H.push_back(std::vector<int> {0, 1, 1, 0, 0, 1, 1});
	H.push_back(std::vector<int> {0, 0, 0, 1, 1, 1, 1});
}
