
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <math.h>

#define PI M_PI
#define DEG2RAD(deg) (deg * M_PI / 180.0)

// http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648/
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

double sine2(double x)
{
	bool bMinus = false;
	if (x < 0) {
		x = -x;
		bMinus = true;
	}
	x = fmod(x, 2 * PI);
	if (x > PI) {
		x -= PI;
		bMinus = !bMinus;
	}
	if (x > PI/2) {
		x = PI - x;
	}
#if 0
	double x2 = x * x;
	double x3 = x * x2;
	double x5 = x2 * x3;
	double x7 = x2 * x5;
	double x9 = x2 * x7;
	return
		x
		- x3/6
		+ x5/120
		- x7/5040
		+ x9/362880
	;
#else
	double x2 = x * x;
	double y = x * (1 + x2 * (-1.0/6.0 + x2 * (+1.0/120.0 + x2 * (-1.0/5040.0 + x2 * 1.0/362880))));
#endif
	return bMinus ? -y : y;
}

static double minimaxsin(double x)
{
    static const
    double a0 =  1.0,
           a1 = -1.666666666640169148537065260055e-1,
           a2 =  8.333333316490113523036717102793e-3,
           a3 = -1.984126600659171392655484413285e-4,
           a4 =  2.755690114917374804474016589137e-6,
           a5 = -2.502845227292692953118686710787e-8,
           a6 =  1.538730635926417598443354215485e-10;
    double x2 = x * x;
    return x * (a0 + x2 * (a1 + x2 * (a2 + x2
             * (a3 + x2 * (a4 + x2 * (a5 + x2 * a6))))));
}

void testSine()
{
	printf("angle\tsin\tline\tdiff\n");
	double start = 1330;
	for (double i=start - 360; i<start + 360.0; i+=1.0) {
		double angle = DEG2RAD(i);
		double resultA = sin(angle);
		double resultC = sine2(angle);
		printf("%f %f %f %e\n",
			i,
			resultA,
			resultC,
			resultC - resultA
		);

	}

}


int main(int argc, char* argv[])
{	
	testSine();
	return 0;
}