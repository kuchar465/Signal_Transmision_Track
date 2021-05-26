#include "encoder.h"

std::complex<double>* Encoder::dft(std::vector<double> inputs, int n)
{
	std::complex<double>* outputs = new std::complex<double>[n];
	for (int k = 0; k < n; ++k) { //czas
		outputs[k] = std::complex<double>(0.0, 0.0);
		for (int t = 0; t < n; ++t) {
			std::complex<double> dft_sino = std::polar(1.0, -(t * k * (2 * M_PI) / double(n))); //wzor
			if (fabs(real(dft_sino)) < std::numeric_limits<float>::epsilon())
				dft_sino = std::complex<double>(0.0, imag(dft_sino));
			if (fabs(imag(dft_sino)) < std::numeric_limits<float>::epsilon())
				dft_sino = std::complex<double>(real(dft_sino), 0.0);
			outputs[k] += inputs[t] * dft_sino;
		}
	}
	return outputs;
}

std::complex<double>* Encoder::dft(std::vector<int> inputs, int n)
{
	std::complex<double>* outputs = new std::complex<double>[n];
	for (int k = 0; k < n; ++k) { //czas
		outputs[k] = std::complex<double>(0.0, 0.0);
		for (int t = 0; t < n; ++t) {
			std::complex<double> dft_sino = std::polar(1.0, -(t * k * (2 * M_PI) / double(n))); //wzor
			if (fabs(real(dft_sino)) < std::numeric_limits<float>::epsilon())
				dft_sino = std::complex<double>(0.0, imag(dft_sino));
			if (fabs(imag(dft_sino)) < std::numeric_limits<float>::epsilon())
				dft_sino = std::complex<double>(real(dft_sino), 0.0);
			outputs[k] += double(inputs[t]) * dft_sino;
		}
	}
	return outputs;
}

std::vector<int> Encoder::NRZI(const std::vector<int> signal)
{
	std::vector<int> nzri;
	int voltage;
	if (signal[0] == 1) {
		voltage = 1;
	}
	else
		voltage = -1;
	for (int i = 0; i < signal.size(); i++) {
		if (signal[i] == 1) {
			voltage = -voltage;
			nzri.push_back(voltage);
		}
		else {
			nzri.push_back(voltage);

		}
	}
	return nzri;
}

std::vector<int> Encoder::TTL(const std::vector<int> signal)
{
	std::vector<int> ttl;
	for (int i = 0; i < signal.size(); i++) {
		if (signal[i] == 1) {
			ttl.push_back(1);
		}
		else {
			ttl.push_back(-1);
		}
	}
	return ttl;
}

std::vector<int> Encoder::manchester(const std::vector<int> signal, int bitamount)
{
	std::vector<char> CLK = generateCLK(bitamount);
	std::vector<int> Manchester;
	int tmp = 0;
	for (int i = 0; i < CLK.size() - 1; i++) {
		if (i % 2 == 0 && i > 0 && i < CLK.size() - 1) {
			tmp++;
		}
		if (CLK[i] == '0') {
			if (signal[tmp] == 0) {
				Manchester.push_back(-1);
			}
			else {
				Manchester.push_back(1);
			}
		}
		else {
			if (signal[tmp] == 1) {
				Manchester.push_back(-1);
			}
			else {
				Manchester.push_back(1);
			}
		}

	}

	return Manchester;
}

std::vector<int> Encoder::BAMI(const std::vector<int> signal)
{
	std::vector<int> bami;
	int voltage = 1;
	for (int i = 0; i < signal.size(); i++) {
		if (signal[i] == 1) {
			bami.push_back(voltage);
			voltage = -voltage;
		}
		else {
			bami.push_back(0);
		}
	}
	return bami;
}

std::byte* Encoder::S2BS(std::string table, std::string sw)
{
	int tmp = table.length();
	std::byte* tab = new std::byte[tmp];
	if (sw == "littleEndian")
	{
		for (int i = 0; i < tmp; i++) {
				tab[i] = std::byte{ table.at(tmp - i - 1) };
		}
		/*for (int i = 0; i < tmp; i++) {
			std::cout << std::hex << std::to_integer<int>(tab[i]);
			std::cout << " ";
		}*/
		std::cout << std::endl;
		return tab;
	}
	else if (sw == "bigEndian") {
		for (int i = 0; i < tmp; i++) {
			tab[i] = std::byte{ table.at(i) };
		}
		/*for (int i = 0; i < tmp; i++) {
			std::cout << std::hex << std::to_integer<int>(tab[i]);
			std::cout << " ";
		}*/
		std::cout << std::endl;
		return tab;
	}
	else {
		std::cout << "wrong argument" << std::endl;
	}
}

void Encoder::bitN(std::vector<int>& mt, int n)
{
	if (mt[n] == 1) {
		mt[n] = 0;
	}
	else {
		mt[n] = 1;
	}
}

std::vector<char> Encoder::generateCLK(int bitamount)
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

double Encoder::signalgen(double t, double A, double fm, double fi)
{
	return A * sin(2 * M_PI * fm * t + fi);
}

std::vector<int> Encoder::encoding(const std::vector<int>& data, int sec)
{
	hammingMatrixes mat;
	std::vector<int> h;
	for (int i = 0; i < 7; i++) {
		int wynik = 0;
		for (int j = 0; j < 4; j++)
			wynik += data[j] * mat.G[i][j];
		wynik %= 2;
		h.push_back(wynik);
	}

	if (sec) {
		int p4 = std::accumulate(h.begin(), h.end(), 0) % 2;
		h.push_back(p4);
	}

	return h;
}

std::vector<int> Encoder::encoderFullH(const std::vector<int> &mt)
{
	std::vector<std::vector<int>> matrix = CreatePackets(mt);
	std::vector<int> encoded;
	for (size_t i = 0; i < matrix.size(); i++) {
		std::vector<int> multiplied = encoding(matrix[i]);
		for (int i = 0; i < multiplied.size(); i++) {
			encoded.push_back(multiplied[i]);
		}
	}
	return encoded;
}

std::vector<char> Encoder::S2BS_CHAR_CUT(int bitamount, const std::byte *ob)
{
	std::vector<char> mt;
	int tmp = 0;
	int tmp2 = 0;
	std::string tmpword;
	for (int i = 0; i < bitamount; i++) {
		if (i % 8 == 0) {
			std::bitset<8> binary(std::to_integer<int>(ob[tmp]));
			tmpword = binary.to_string();
			tmp2 = 0;
			tmp++;
		}
		mt.push_back(tmpword.at(tmp2));
		tmp2++;
	}
	return mt;
}

std::vector<int> Encoder::S2BS_INT(const std::vector<char>& data)
{
	std::vector<int> mtint;
	for (int i = 0; i < data.size(); i++) {
		if (data[i] == '1') {
			mtint.push_back(1);
		}
		else {
			mtint.push_back(0);
		}
	}
	return mtint;
}

std::vector<std::vector<int>> Encoder::CreatePackets(const std::vector<int>& data)
{
	std::vector<std::vector<int>> matrix;
	std::vector<int> vectors;

	for (int i = 1; i <= data.size(); i++) {
		vectors.push_back(data[i - 1]);
		if (i % 4 == 0) {
			matrix.push_back(vectors);
			vectors.clear();
		}
	}
	return matrix;
}

void Encoder::print(const std::vector<int>& data)
{
	for (int i = 0; i < data.size(); i++) {
		std::cout << data[i] << " ";
	}
	std::cout << std::endl;
}

std::vector<double> Encoder::ASK(const std::vector<int> signal, int cz, double freq, double A1, double A2, double fm)
{
	std::vector<double> out;
	for (int i = 0; i < cz; i++) {
		if (signal[i] == '0') {
			for (double k = i; k < (i + 1); k += freq) {
				out.push_back(signalgen(k, A1, fm, 1));
			}
		}
		else {
			for (double k = i; k < (i + 1); k += freq) {
				out.push_back(signalgen(k, A2, fm, 1));
			}
		}
	}
	return out;
}

std::vector<double> Encoder::PSK(const std::vector<int> signal, int cz, double freq, double A, double fm, double fi0, double fi1)
{
	std::vector<double> out;
	for (int i = 0; i < cz; i++) {
		if (signal[i] == '0') {
			for (double k = i; k < (i + 1); k += freq) {
				out.push_back(signalgen(k, A, fm, fi0));
			}
		}
		else {
			for (double k = i; k < (i + 1); k += freq) {
				out.push_back(signalgen(k, A, fm, fi1));
			}
		}
	}
	return out;
}

std::vector<double> Encoder::FSK(const std::vector<int> signal, int cz, double freq, double A, double f0, double f1)
{
	std::vector<double> out;
	for (int i = 0; i < cz; i++) {
		if (signal[i] == '0') {
			for (double k = i; k < (i + 1); k += freq) {
				out.push_back(signalgen(k, A, f0, 1));
			}
		}
		else {
			for (double k = i; k < (i + 1); k += freq) {
				out.push_back(signalgen(k, A, f1, 1));
			}
		}
	}
	return out;
}

void Encoder::zapisdft(std::complex<double>* inputs, int n, std::string nazwa, double freq)
{
	std::ofstream plik;
	plik.open(nazwa);
	double* M = new double[n];
	double* M2 = new double[n];
	double fk = 0.0;
	plik << "x" << " " << "y" << std::endl;
	for (int i = 0; i < n; i++) {
		M[i] = sqrt(pow(inputs[i].real(), 2) + pow(inputs[i].imag(), 2));
		M2[i] = 10 * log10(M[i]);
		fk = i * freq / double(n);
		plik << fk << " " << M2[i] << std::endl;
	}
	plik.close();
	delete M;
	delete M2;
}

Encoder::hammingMatrixes::hammingMatrixes()
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
