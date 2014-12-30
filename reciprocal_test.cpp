/*
The MIT License (MIT)

Copyright (c) 2014 berupon@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

// Newtonñ@Ç≈ãtêîÇãÅÇﬂÇÈÉvÉçÉOÉâÉÄ

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

static uint32_t g_zero = 0;
static uint32_t g_reciproTable[256];

uint32_t newton_recipro_FixedPoint(uint32_t a, size_t* q)
{
	if (a > 256) {
		size_t nBits = 32 - __lzcnt(a);
		uint32_t idx = a >> (nBits - 8);
		uint32_t initial = g_reciproTable[idx];
		uint32_t x = initial;
		int32_t intOffsetBits = (int)nBits - 15;
		*q = 46 + intOffsetBits;
		size_t absOffsetBits = abs(intOffsetBits);
		a <<= (absOffsetBits << 1) * (nBits < 16);
		static const size_t nIte = 2;
		for (size_t i=0; i<nIte; ++i) {
			uint32_t xmulx = ((uint64_t)x * x) >> 32;
			uint64_t right = (uint64_t)a * xmulx;
			uint64_t xmul2 = ((uint64_t)x << (15 + absOffsetBits)) - 1;
			x = (xmul2 - right) >> (14 + absOffsetBits);
		}
		return x;
	}else {
		if (a == 0) {
			return 0;
		}
		*q = 63 - __lzcnt(a - 1);
		return g_reciproTable[a - 1];
	}
}

int main(int argc, char* argv[])
{
	for (size_t i=1; i<=256; ++i) {
		uint64_t tmp = (1ULL << 62) / i;
		size_t lnz = __lzcnt64(tmp);
		tmp >>= 32 - lnz;
		g_reciproTable[i-1] = tmp;
	}

	size_t start = 1U << 13;
	size_t end = start + (1 << 4);
	size_t step = 1;
	for (size_t i=start; i<end; i+=step) {
//		double recipro = newton_recipro(i);
//		printf("%d %.9f\n", i, recipro);
		size_t q;
		uint32_t recipro2 = newton_recipro_FixedPoint(i, &q);
		double recipro2f = recipro2 / (double)(1LL << q);
		double recipro = (i == 0) ? 0 : (1/(double)i);
		uint64_t recipro1 = recipro * (1LLU << q);
		int64_t diff = (int64_t)recipro1 - (int64_t)recipro2;
		printf("%u %llu %u %lld %f\n", i, recipro1, recipro2, diff, (diff*100.0)/recipro1);
	}

	return 0;
}