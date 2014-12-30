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
	double s10 = s * 10;
	double ss3 = s * s * 3;
	for (size_t i=0; i<2; ++i) {
		double xs = x * x;
//		double y = s * xs;
//		x = x / 8 * (15 - y * (10 - 3 * y));
		x *= (15 - xs * (s10 - ss3 * xs));
		x /= 8;
	}
	return x;
}



int main(int argc, char* argv[])
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
	double end = 1024.0;
	double step = (end - begin) / 10240;
	for (double x = begin; x < end; x += step) {
		double recipro_sqrt = calc_recipro_sqrt(x);
		printf("%f %f %.10f\n", x, x * recipro_sqrt, x*x*recipro_sqrt*recipro_sqrt);
	}
	
	return 0;
}

