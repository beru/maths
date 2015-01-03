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

// 平方根（近似値）を求めるプログラム

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <intrin.h>
#include <math.h>
#include <float.h>
#include <assert.h>

double g_isqrt_initial_values[256+1];

double calc_recipro_sqrt(double s)
{
	int exp;
	double frac = frexp(s, &exp);
	frac -= 0.5;
	double initial;
	if (exp >= 3) {
		exp -= 3;
		frac += 0.5 * (exp & 1);
		size_t idx = (_countof(g_isqrt_initial_values) - 1) * frac;
		initial = g_isqrt_initial_values[idx];
		int scale = exp >> 1;
		initial /= (2 << scale);
	}else {
		if (s == 0) {
			return 0;
		}else {
			exp = 3 - exp;
			frac += 0.5 * (exp & 1);
			size_t idx = (_countof(g_isqrt_initial_values) - 1) * frac;
			initial = g_isqrt_initial_values[idx];
			int scale = exp >> 1;
			initial *= (1 << (scale + (exp & 1) - 1));
		}
	}

	double x = initial;
#if 0
	// quartic algorithm
	for (size_t i=0; i<2; ++i) {
		double h = 1 - s * x * x;
		double h2 = h * h;
		double h3 = h * h2;
		x = x + x * (8 * h + 6 * h2 + 5 * h3) / 16;
//		x = x + x * h * (8 + h * (6 + 5 * h)) / 16;
	}
#elif 0
	// cubic iteration
	for (size_t i=0; i<2; ++i) {
		double h = 1 - s * x * x;
		double h2 = h * h;
		x = x + x * (4 * h + 3 * h2) / 8;
	}
#else
	for (size_t i=0; i<4; ++i) {
		double h = 1 - s * x * x;
		x = x + x * h / 2;
	}
#endif
	return x;
}

void testFloat()
{
	// setup table
	{
		double begin = 1;
		double end = begin + 3;
		double step = (end - begin) / (_countof(g_isqrt_initial_values) - 1);
		size_t cnt = 0;
		for (double x = begin; x <= end; x += step) {
			double isqrt = 1.0 / sqrt(x);
			g_isqrt_initial_values[cnt++] = isqrt;
		}
		int hoge = 0;
	}

	double begin = 0;
	double end = 64;
	double step = (end - begin) / 32;
	for (double x = begin; x < end; x += step) {
		double recipro_sqrt = calc_recipro_sqrt(x);
		double square_root0 = sqrt(x);
		double square_root1 = x * recipro_sqrt;
		double square_root2 = 1.0/recipro_sqrt;
		printf("%f %.13f %.13f\n", x, square_root0, square_root2);
	}
}

// 引数 s  : 負でない整数
// 引数 oq : 戻り値の固定小数点のQ値（小数部ビット数）
// 戻り値  : 平方根の逆数 固定小数点形式
uint32_t calc_recipro_sqrt_integer(uint32_t s, size_t* oq)
{
	if (s == 0) {
		return 0;
	}
	size_t clz = __lzcnt(s);
	size_t nBits = (32 - clz) & (~1);
	s <<= clz;
	size_t offset = 31 - (clz >> 1);
	uint32_t x = 1u << 31;
	for (size_t i=0; i<4; ++i) {
		uint32_t s_mul_x = ((uint64_t)s * x) >> 31;								// (Q.clz * Q.31) >> 31 = Q.clz
		uint32_t x_mul_x = ((uint64_t)x * x) >> 31;								// (Q.31 * Q.31) >> 31 = Q.31
		uint64_t s_mul_x_mul_x_mul_x = ((uint64_t)s_mul_x * x_mul_x) >> (clz + nBits);	// (Q.clz * Q.31) >> (clz + nBits) = Q.(31 - nBits)
		x = x + ((x - s_mul_x_mul_x_mul_x) >> 1);
	}
	*oq = 16 + offset - (clz & 1);
	return x;
}

void testFixedPoint()
{

//	calc_recipro_sqrt_integer(125348);

	uint32_t begin = 1 << 0;
	uint32_t end = begin + 64;
	size_t q = 0;
	for (uint32_t i = begin; i < end; ++i) {
		size_t q;
		uint32_t x = calc_recipro_sqrt_integer(i, &q);
		double square_root0 = sqrt((double)i);
		double square_root1 = ((uint64_t)i * x) / (double)(1llu << q);
		double square_root2 = 1.0 / (x / (double)(1llu << q));
		printf("%u %u %.13f %.13f\n", i, x, square_root0, square_root2);
	}
}

int main(int argc, char* argv[])
{
	testFloat();
	testFixedPoint();
	return 0;
}

