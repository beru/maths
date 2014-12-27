
#include <stdio.h>
#include <stdint.h>
#include <intrin.h>
#include <math.h>

double newton_recipro(double a)
{
	size_t nlz = __lzcnt((int)a);
	double initial = 1.0 / (1 << (32 - nlz));
	double x = initial;
	for (size_t i=0; i<2; ++i) {
		x = 2 * x - a * x * x;
	}
	return x;
}

static uint32_t g_reciproTable[256];

uint32_t newton_recipro_FixedPoint(uint32_t a, size_t* q)
{
	static const size_t nIte = 1;
	size_t nlz = __lzcnt((int)a);
	size_t nBits = 32 - nlz;
	uint32_t x;
	if (nBits > 15) {
		size_t offsetBits = nBits - 15;
		uint32_t idx = a >> (nBits - 8);
		uint32_t initial2 = g_reciproTable[idx];
		x = initial2;
		for (size_t i=0; i<nIte; ++i) {
			uint32_t xmulx = ((uint64_t)x * x) >> (31 + offsetBits);
			uint64_t right = (uint64_t)a * xmulx;
			uint64_t xmul2 = (uint64_t)x << 16;
			x = (xmul2 - right) >> 15;
		}
		*q = 46 + offsetBits;
	}else {
		uint32_t initial = 1u << (nBits << 1);
		x = initial;
		for (size_t i=0; i<nIte; ++i) {
			uint32_t xmulx = (((uint64_t)x * x) >> (nBits*2)) - 1;
			uint64_t right = (uint64_t)a * xmulx + xmulx;
			uint64_t xmul2 = (uint64_t)x << (1 + nBits);
			x = (xmul2 - right) >> nBits;
		}
		*q = nBits * 3;
	}
	return x;
}

int main(int argc, char* argv[])
{
	for (size_t i=1; i<=256; ++i) {
		uint64_t tmp = (1ULL << 62) / i;
		size_t lnz = __lzcnt64(tmp);
		tmp >>= 32 - lnz;
		g_reciproTable[i-1] = tmp;
	}

	size_t start = (1 << 24);
	size_t end = start + (1 << 16);
	size_t step = (end-start) / 256;
	for (size_t i=start; i<end; i+=step) {
//		double recipro = newton_recipro(i);
//		printf("%d %.9f\n", i, recipro);
		size_t q;
		uint32_t recipro2 = newton_recipro_FixedPoint(i, &q);
		double recipro2f = recipro2 / (double)(1LL << q);
		double recipro = 1/(double)i;
		uint64_t recipro1 = recipro * (1LL << q);
		int64_t diff = (int64_t)recipro1 - (int64_t)recipro2;
		printf("%u %llu %u %lld %f\n", i, recipro1, recipro2, diff, (diff*100.0)/recipro1);
	}

	return 0;
}