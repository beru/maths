
// 三角関数（Sine, Cosine）の近似法プログラム

#include <stdio.h>
#include <stdint.h>
#include <limits.h>

#define _USE_MATH_DEFINES
#include <math.h>

#define PI M_PI
#define DEG2RAD (PI / 180.0)

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

double approxSineByRadian(double x)
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
	double x2 = x * x;
#if 1
	double x3 = x * x2;
	double x5 = x2 * x3;
	double x7 = x2 * x5;
	double x9 = x2 * x7;
	double y =
		x
		- x3/6
		+ x5/120
		- x7/5040
		+ x9/362880
	;
#else
	double y = x * (1 + x2 * (-1.0/6.0 + x2 * (+1.0/120.0 + x2 * (-1.0/5040.0 + x2 * 1.0/362880))));
#endif
	return bMinus ? -y : y;
}

double approxSineByDegree(double x)
{
	bool bMinus = false;
	if (x < 0) {
		x = -x;
		bMinus = true;
	}
	x = fmod(x, 360);
	if (x > 180) {
		x -= 180;
		bMinus = !bMinus;
	}
	if (x > 90) {
		x = 180 - x;
	}
	double c0 = DEG2RAD;						// +(pi/180)
	double c1 = -8.86096155701298015989E-7;		// -(pi/180)^3 / 3!
	double c2 = 1.349601623163255010593E-11;	// +(pi/180)^5 / 5!
	double c3 = -9.78838486161772760954E-17;	// -(pi/180)^7 / 7!
	double c4 = 4.14126741725732068529E-22;		// +(pi/180)^9 / 9!

	double x2 = x * x;
#if 0
	double x3 = x * x2;
	double x5 = x2 * x3;
	double x7 = x2 * x5;
	double x9 = x2 * x7;
	double y =
		x * c0
		+ x3 * c1
		+ x5 * c2
		+ x7 * c3
		+ x9 * c4
	;
#else
	double y = x * (c0 + x2 * (c1 + x2 * (c2 + x2 * (c3 + x2 * c4))));
#endif
	return bMinus ? -y : y;
}

double approxSineByNormalized(double x)
{
	bool bMinus = false;
	if (x < 0) {
		x = -x;
		bMinus = true;
	}
	x = fmod(x, 4.0);
	if (x > 2.0) {
		x = 4 - x;
		bMinus = !bMinus;
	}
	if (x > 1.0) {
		x = 2.0 - x;
	}
	double c0 = 1.570796326794896619231;		// +(pi/2)
	double c1 = -0.6459640975062462536558;		// -(pi/2)^3 / 3!
	double c2 = 0.07969262624616704512051;		// +(pi/2)^5 / 5!
	double c3 = -0.004681754135318688100685;	// -(pi/2)^7 / 7!
	double c4 = 1.604411847873598218727E-4;		// +(pi/2)^9 / 9!

	double x2 = x * x;
#if 0
	double x3 = x * x2;
	double x5 = x2 * x3;
	double x7 = x2 * x5;
	double x9 = x2 * x7;
	double y =
		x * c0
		+ x3 * c1
		+ x5 * c2
		+ x7 * c3
		+ x9 * c4
	;
#else
	double tmp;
	tmp = c3 + x2 * c4;
	tmp = c2 + x2 * tmp;
	tmp = c1 + x2 * tmp;
	tmp = c0 + x2 * tmp;
	tmp = x * tmp;
	double y = tmp;
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

typedef int fp16;
fp16 ne3d_sintp_02pi(fp16 angle)
{
	fp16 Xi2,Xi3,Yi,Yi2,Yi3,Yi4;

	angle = angle >> 1;
	Xi2 = ((angle >> 1) * (angle >> 2)) >> 13;
	Xi3 = (6640 * Xi2) >> 12;
	Yi =  (((20860 * angle) >> 13) - Xi3) >> 1;
	Yi2 = (Yi * Yi) >> 14;
	Yi3 = (14746 * Yi2) >> 16;
	Yi4 = (25395 * Yi ) >> 14;
	Yi = Yi4 + Yi3;

	return Yi;
}

// 分岐処理無しで絶対値を求める
// Compute the integer absolute value (abs) without branching
// https://graphics.stanford.edu/~seander/bithacks.html#IntegerAbs
static inline
unsigned int getAbsolute(
	int v	// we want to find the absolute value of v
	)
{
	int const mask = v >> sizeof(int) * CHAR_BIT - 1;
	unsigned int r;  // the result goes here 
	r = (v + mask) ^ mask;
	return r;
}

// 分岐処理無しで符号反転
// Conditionally negate a value without branching
// https://graphics.stanford.edu/~seander/bithacks.html#ConditionalNegate
static inline
int negate(
	int v,			// Input value to negate if fNegate is true.
	bool fNegate	// Flag indicating if we should negate v.
	)
{
	int r;         // result = fNegate ? -v : v;
	r = (v ^ -fNegate) + fNegate;
	return r;
}

// 整数の符号ビット取得
// 入力が負の整数だったら 1 が返る
// 正の整数だったら 0 が返る
// Compute the sign of an integer 
// https://graphics.stanford.edu/~seander/bithacks.html#CopyIntegerSign
template <typename IntegerT>
static inline
int getSignBit(IntegerT v)
{
	return v >> (sizeof(IntegerT) * CHAR_BIT - 1);
}

// input : fixed point radian angle. 1rad == 180deg == (PI << InputFracBitsCount) == (3.14.... << InputFracBitsCount)
// output : fixed point normalized angle. from -1.0(-2rad,-360deg) to +1.0(+2rad,+360deg)
template <size_t InputFixedPointQ, size_t OutputFixedPointQ>
static inline
int32_t normalizeRadianAngleBy2PI(int32_t radianFixedPoint)
{
	return 0;
}

template <size_t FixedPointQ, typename IntegerT>
IntegerT getReciprocal(int64_t dividiend, IntegerT divisor)
{
	int64_t result = ((divisor/2) + (dividiend << FixedPointQ)) / divisor;
	return (IntegerT) result;
}

template <typename IntegerT>
IntegerT minval(IntegerT x, IntegerT y)
{
	return y ^ ((x ^ y) & -(x < y));
}

// Degree単位の角度値の正規化
// 入力：度数法の値、固定小数点形式（Q.InputFixedPointQ）
// 出力：正規化された角度、固定小数点形式（Q.InputFixedPointQ）
// input : fixed point degree angle. ex 180deg == (180 << InputFracBitsCount)
// output : fixed point normalized angle. from -1.0(-2rad,-360deg) to +1.0(+2rad,+360deg)
template <size_t InputFixedPointQ, size_t OutputFixedPointQ>
static inline
int32_t normalizeDegreeAngleBy2PI(int32_t degreeFixedPoint)
{
	static const int invPI2_FixedPointQ = 39;
	static const int invPI2 = ((1ULL << invPI2_FixedPointQ) + 180) / 360;
	// multiply input value with invPI2 and right shift to achieve input / PI2
	int32_t result = ((int64_t)degreeFixedPoint * invPI2) >> ((InputFixedPointQ + invPI2_FixedPointQ) - OutputFixedPointQ);

	// extract fractional bits only
	static const int outputFracBitsMask = (1 << OutputFixedPointQ) - 1;
	return result & outputFracBitsMask;
}

// 三角関数 sine の近似関数
// 浮動小数点演算器を持たないプロセッサ用の固定小数点演算
// 入力：度数法の角度、固定小数点形式（Q.InputFixedPointQ）
// 出力：sine の結果、固定小数点形式 Q.31
// input : fixed point decimal angle.
// output : fixed point sine value. from -1.0 to 1.0
template <
	size_t InputFixedPointQ,	// 入力の固定小数点数の小数部ビット数
	size_t WorkFixedPointQ = 24	// 処理変数の固定小数点数の小数部ビット数
>
static inline
int32_t approxSineByDegreeFixedPoint(int32_t x)
{
	// 角度が負の値の場合、処理結果が負になる。（絶対値は同じ）
	// 角度を絶対値にして処理する。
	// 処理結果の値の符号を反転して上下反転する。
	bool bNegate = getSignBit(x);
	uint32_t ux = normalizeDegreeAngleBy2PI<InputFixedPointQ, WorkFixedPointQ>(getAbsolute(x));

	// 絶対値が 0.5(+π)以上の場合は入力角度を折り返し、処理結果の符号を反転させる
	bool isGreaterThanOrEqHalf = ux >> (WorkFixedPointQ - 1);
	bNegate ^= isGreaterThanOrEqHalf;
	static const uint32_t halfBitMask = 1 << (WorkFixedPointQ - 1);
	static const uint32_t halfFracBitsMask = halfBitMask - 1;
	ux &= halfFracBitsMask;

	// 0.25（π/2）を超えていたら、左右に反転する
	// 0.5（π）から引く事により左右反転が可能。
	// 正負符号は反転する必要はない。
	static const uint32_t quarterBitMask = halfBitMask >> 1;
	static const uint32_t quarterFracBitsMask = quarterBitMask - 1;
	int border = ux & quarterBitMask;
	int frac = ux & quarterFracBitsMask;
	ux = getAbsolute(border - frac);

	// 0～0.25(π/2)までに限定したのでここから近似計算

	// 1.0 = 2rad = 360deg = (1 << InputFixedPointQ)
	// 0.25 = rad/2 = 90deg = (1 << (InputFixedPointQ - 2))

	int64_t x2 = (uint64_t)ux * ux;	// Q.(InputFixedPointQ + InputFixedPointQ - 4)
	x2 >>= (WorkFixedPointQ + WorkFixedPointQ - 4) - 31;	// Q.35
	// 近似計算

	// double y = x * (c0 + x2 * (c1 + x2 * (c2 + x2 * (c3 + x2 * c4))));
	// sin x(deg) = deg2rad * x - (deg2rad^3 * x^3/3!) + (deg2rad^5 * x^5/5!) - (deg2rad^7 * x^7/7!)
	// = x * (nrm2rad + x^2 * (-(nrm2rad^3)/6 + x^2 * (+nrm2rad^5/120 + x^2 * (-nrm2rad^7/5040 + x^2 * nrm2rad^9/362880))))
	static const int64_t c0 = +3373259426LL;	// +(2^31) * (pi/2)
	static const int32_t c1 = -1387197337;		// -(2^31) * ((pi/2)^3) / 3!
	static const int32_t c2 = +1369108894;		// +(2^34) * ((pi/2)^5) / 5!
	static const int32_t c3 = -1286910778;		// -(2^38) * ((pi/2)^7) / 7!
	static const int32_t c4 = +1411255586;		// +(2^43) * ((pi/2)^9) / 9!

	int64_t tmp;
	tmp = c3 + ((x2 * c4) >> 36);	// Q.38 + ((Q.35 * Q.43) >> 40) = Q.38
	tmp = c2 + ((x2 * tmp) >> 35);	// Q.34 + ((Q.35 * Q.38) >> 39) = Q.34
	tmp = c1 + ((x2 * tmp) >> 34);	// Q.31 + ((Q.35 * Q.34) >> 38) = Q.31
	tmp = c0 + ((x2 * tmp) >> 31);	// Q.31 + ((Q.35 * Q.31) >> 35) = Q.31
	tmp = (ux * tmp) >> (WorkFixedPointQ - 2);	// Q.InputFixedPointQ * Q.31
		// シフト数の - 2 は、2pi = 320deg = 1.0 に正規化された入力の値を、90deg = 1.0 の値として扱う為
	int result = minval(tmp, (int64_t)INT32_MAX);

	// 対称性を考慮した符号反転
	return negate(result, bNegate);
}

void testSine()
{
	printf("degree\n");
	double start = 0;
	for (double i=start - 360; i<start + 360.0; i+=1) {
		double radAngle = DEG2RAD * i;
		double resultA = sin(radAngle);
		double resultB = approxSineByRadian(radAngle);
		double resultC = approxSineByDegree(i);
		double resultD = approxSineByNormalized(i / 90);
//		int32_t iResultE = approxSineByDegreeFixedPoint<16>(i * (1 << 16));
		int32_t iResultE = ne3d_sintp_02pi(radAngle * (1 << 16));
		double resultE = iResultE / (double)(1LL << 16);
		printf("%f %f %f %f %f %f %e %e %e %e\n",
			i,
			resultA,
			resultB,
			resultC,
			resultD,
			resultE,
			resultB - resultA,
			resultC - resultA,
			resultD - resultA,
			resultE - resultA
		);

	}

}

int main(int argc, char* argv[])
{
	testSine();
	return 0;
}

