
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <math.h>

#define PI M_PI
#define DEG2RAD(deg) (deg * M_PI / 180.0)

float sine(float x)
{
    // Convert the input value to a range of -1 to 1
    x = x * (1.0f / PI);

    // Wrap around
    volatile float z = (x + 25165824.0f);
    x = x - (z - 25165824.0f);

    #if LOW_SINE_PRECISION
        return 4.0f * (x - x * abs(x));
    #else
        float y = x - x * abs(x);

        const float Q = 3.1f;
        const float P = 3.6f;

        return y * (Q + P * abs(y));
    #endif
}

void testSine()
{
	printf("angle\tsin\tline\tdiff\n");
	for (double i=0; i<90.0; i+=1.0) {
		double angle = DEG2RAD(i);
		double resultA = sin(angle);
		double resultB = sine(angle);
		double resultC = 0.5 * i / 30.0;
		printf("%f %f %f %f %f %f\n",
			i,
			resultA,
			resultB,
			resultC,
			resultB - resultA,
			resultC - resultB
		);

	}

}


int main(int argc, char* argv[])
{
	testSine();
	return 0;
}