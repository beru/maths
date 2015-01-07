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
#if 0
		// newton method
		x = 2 * x - a * x * x;
#elif 1
		// 1 / A ÇÃÇSéüé˚ë©ÇÃëQâªéÆ
		double h = 1 - a * x;
		x = x * (1 + h) * (1 + h * h);
#endif
	}
	return x;
}

void testFloat()
{
	size_t start = 1U << 1;
	size_t end = start + (1 << 8);
	size_t step = 1;
	for (size_t i=start; i<end; i+=step) {
		double recipro0 = 1.0 / i;
		double recipro1 = newton_recipro(i);
		printf("%d %.9f %.9f %.9f\n", i, recipro0, recipro1, recipro1 - recipro0);
	}
}

static uint32_t g_zero = 0;
static uint32_t g_reciproTable[256];

uint32_t newton_recipro_FixedPoint(uint32_t a, size_t* q)
{
	if (a < 256) {
		if (a == 0) {
			return 0;
		}
		*q = 63 - __lzcnt(a - 1);
		return g_reciproTable[a - 1];
	}else {
		// newton method
		size_t clz = __lzcnt(a);
		size_t nBits = 32 - clz;
#if 0
		uint32_t initial = (1llu << 31) - 1;
#else
		uint32_t idx = a >> (nBits - 8);
		uint32_t initial = g_reciproTable[idx];
#endif
		uint32_t x = initial;
		static const size_t nIte = 2;
		for (size_t i=0; i<nIte; ++i) {
			uint32_t xmulx = ((uint64_t)x * x) >> 32;
			uint64_t right = (uint64_t)a * xmulx;
			uint64_t xmul2 = (uint64_t)x << nBits;
			x = ((xmul2 - right) >> (nBits - 1)) + 0;
		}
		*q = 31 + nBits;
		return x;
	}
}

void testFixedPoint()
{
	for (size_t i=1; i<=256; ++i) {
		uint64_t tmp = (1ULL << 62) / i;
		size_t lnz = __lzcnt64(tmp);
		tmp >>= 32 - lnz;
		g_reciproTable[i-1] = tmp;
	}

	size_t start = 1U << 11;
	size_t end = start + (1 << 10);
	size_t step = 1;
	for (size_t i=start; i<end; i+=step) {
		size_t q;
		uint32_t recipro2 = newton_recipro_FixedPoint(i, &q);
		double recipro2f = recipro2 / (double)(1LL << q);
		double recipro = (i == 0) ? 0 : (1/(double)i);
		uint64_t recipro1 = recipro * (1LLU << q);
		int64_t diff = (int64_t)recipro1 - (int64_t)recipro2;
		printf("%u %llu %u %lld %f\n", i, recipro1, recipro2, diff, (diff*100.0)/recipro1);
	}
}

int main(int argc, char* argv[])
{
//	testFloat();

	testFixedPoint();
	return 0;
}