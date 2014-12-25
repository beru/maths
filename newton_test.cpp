
#include <stdio.h>
#include <intrin.h>


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

int main(int argc, char* argv[])
{
	for (size_t i=0; i<1024; ++i) {
		double recipro = newton_recipro(i);
		printf("%d %.9f\n", i, recipro);
	}

	return 0;
}