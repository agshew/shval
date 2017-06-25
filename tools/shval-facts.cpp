/*
 * Copyright (c) 2018, Andrew Shewmaker.
 * All rights reserved.
 *
 * This file is part of SHVAL. For details, see https://github.com/lam2mo/shval
 *
 * Please also see the LICENSE file for our notice and the LGPL.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*
 * Shadow value analysis of the distribution of numbers
 */

#include <cassert>

#include <sys/types.h>	// getpid
#include <unistd.h>	// getpid

#include <math.h>

#include "pin.H"

#include <map>
using namespace std;

#include <cmath>
#include <cfloat>

KNOB<UINT32> KnobOpsPerStep(KNOB_MODE_WRITEONCE, "pintool",
	"step", "100000", "record stats every X operations (default=100000)");

static UINT32 _ops_per_step;

UINT64 ops = 0;
int exponent;

FILE	*normalsOut,
	*subnormalsOut,
	*fractionsOut,
	*exponentsOut,
	*poiOut;

map<double, size_t> normals, subnormals, fractions;
map<double, size_t>::iterator dbl_bucket;

map<int, size_t> exponents;
map<int, size_t>::iterator exp_bucket;

// The histogram of normals lacks bins for subnormals
#define NORMAL_BINS 97
double normal_breaks[NORMAL_BINS] =
{-DBL_MAX,
 -10000000000000000000000000.0,
 -1000000000000000000000000.0,
 -100000000000000000000000.0,
 -10000000000000000000000.0,
 -1000000000000000000000.0,
 -100000000000000000000.0,
 -10000000000000000000.0,
 -1000000000000000000.0,
 -100000000000000000.0,
 -10000000000000000.0,
 -1000000000000000.0,
 -100000000000000.0,
 -10000000000000.0,
 -1000000000000.0,
 -100000000000.0,
 -10000000000.0
 -1000000000.0,
 -100000000.0,
 -10000000.0,
 -1000000.0,
 -100000.0,
 -10000.0,
 -1000.0,
 -100.0,
 -10.0,
 -8.0,
 -4.0,
 -2.0,
 -1.0,
 -0.5,
 -0.25,
 -0.125,
 -0.0625,
 -0.03125,
 -0.015625,
 -0.0078125,
 -0.00390625,
 -0.001953125,
 -0.0009765625,
 -1.9073486328125e-06,
 -9.5367431640625e-07,
 -4.76837158203125e-07,
 -2.384185791015625e-07,
 -1.1920928955078125e-07,
 -5.960464477539063e-08,
 -2.9802322387695312e-08,
 -DBL_MIN,
 0.0,
 DBL_MIN,		/* smallest normal double, 2.22507e-308*/
 2.9802322387695312e-08,/* 1.0 / pow(2,25) */
 5.960464477539063e-08,	/* 1.0 / pow(2,24) */
 1.1920928955078125e-07,/* 1.0 / pow(2,23) */
 2.384185791015625e-07,	/* 1.0 / pow(2,22) */
 4.76837158203125e-07,	/* 1.0 / pow(2,21) */
 9.5367431640625e-07,	/* 1.0 / pow(2,20) */
 1.9073486328125e-06,	/* 1.0 / pow(2,19) */
 0.0009765625,		/* 1.0 / pow(2,10) */
 0.001953125,
 0.00390625,
 0.0078125,
 0.015625,
 0.03125,
 0.0625,
 0.125,
 0.25,
 0.5,
 1.0,
 2.0,
 4.0,
 8.0,
 10.0,
 100.0,
 1000.0,
 10000.0,
 100000.0,
 1000000.0,
 10000000.0,
 100000000.0,
 1000000000.0,
 10000000000.0,			/* 1e+10 */
 100000000000.0,
 1000000000000.0,
 10000000000000.0,
 100000000000000.0,
 1000000000000000.0,
 10000000000000000.0,
 100000000000000000.0,
 1000000000000000000.0,
 10000000000000000000.0,
 100000000000000000000.0,	/* 1e+20 */
 1000000000000000000000.0,
 10000000000000000000000.0,
 100000000000000000000000.0,
 1000000000000000000000000.0,
 10000000000000000000000000.0,	/* 1e+25 */
 DBL_MAX};			/* largest normal double, 1.79769e+308 */

// The histogram of subnormals lacks bins for normals
#define SUBNORMAL_BINS 35
double subnormal_breaks[SUBNORMAL_BINS] =
{-DBL_MIN,
 -2e-309,
 -2e-310,
 -2e-311,
 -2e-312,
 -2e-313,
 -2e-314,
 -2e-315,
 -2e-316,
 -2e-317,
 -2e-318,
 -2e-319,
 -2e-320,
 -2e-321,
 -2e-322,
 -2e-323,
 -5e-324,
 0.0,
 5e-324,		/* smallest subnormal double */
 2e-323,
 2e-322,
 2e-321,
 2e-320,
 2e-319,
 2e-318,
 2e-317,
 2e-316,
 2e-315,
 2e-314,
 2e-313,
 2e-312,
 2e-311,
 2e-310,
 2e-309,
 DBL_MIN};		/* smallest normal double, 2.22507e-308*/

#define FRACTION_BINS 202
double fraction_breaks[FRACTION_BINS] =
{-1.0,
 -0.995,
 -0.99,
 -0.985,
 -0.98,
 -0.975,
 -0.97,
 -0.965,
 -0.96,
 -0.955,
 -0.95,
 -0.945,
 -0.94,
 -0.935,
 -0.93,
 -0.925,
 -0.92,
 -0.915,
 -0.91,
 -0.905,
 -0.9,
 -0.895,
 -0.89,
 -0.885,
 -0.88,
 -0.875,
 -0.87,
 -0.865,
 -0.86,
 -0.855,
 -0.85,
 -0.845,
 -0.84,
 -0.835,
 -0.83,
 -0.825,
 -0.82,
 -0.815,
 -0.81,
 -0.805,
 -0.8,
 -0.795,
 -0.79,
 -0.785,
 -0.78,
 -0.775,
 -0.77,
 -0.765,
 -0.76,
 -0.755,
 -0.75,
 -0.745,
 -0.74,
 -0.735,
 -0.73,
 -0.725,
 -0.72,
 -0.715,
 -0.71,
 -0.705,
 -0.7,
 -0.695,
 -0.69,
 -0.685,
 -0.68,
 -0.675,
 -0.67,
 -0.665,
 -0.66,
 -0.655,
 -0.65,
 -0.645,
 -0.64,
 -0.635,
 -0.63,
 -0.625,
 -0.62,
 -0.615,
 -0.61,
 -0.605,
 -0.6,
 -0.595,
 -0.59,
 -0.585,
 -0.58,
 -0.575,
 -0.57,
 -0.565,
 -0.56,
 -0.555,
 -0.55,
 -0.545,
 -0.54,
 -0.535,
 -0.53,
 -0.525,
 -0.52,
 -0.515,
 -0.51,
 -0.505,
 -0.5,
  0.5,
  0.505,
  0.51,
  0.515,
  0.52,
  0.525,
  0.53,
  0.535,
  0.54,
  0.545,
  0.55,
  0.555,
  0.56,
  0.565,
  0.57,
  0.575,
  0.58,
  0.585,
  0.59,
  0.595,
  0.6,
  0.605,
  0.61,
  0.615,
  0.62,
  0.625,
  0.63,
  0.635,
  0.64,
  0.645,
  0.65,
  0.655,
  0.66,
  0.665,
  0.67,
  0.675,
  0.68,
  0.685,
  0.69,
  0.695,
  0.7,
  0.705,
  0.71,
  0.715,
  0.72,
  0.725,
  0.73,
  0.735,
  0.74,
  0.745,
  0.75,
  0.755,
  0.76,
  0.765,
  0.77,
  0.775,
  0.78,
  0.785,
  0.79,
  0.795,
  0.8,
  0.805,
  0.81,
  0.815,
  0.82,
  0.825,
  0.83,
  0.835,
  0.84,
  0.845,
  0.85,
  0.855,
  0.86,
  0.865,
  0.87,
  0.875,
  0.88,
  0.885,
  0.89,
  0.895,
  0.9,
  0.905,
  0.91,
  0.915,
  0.92,
  0.925,
  0.93,
  0.935,
  0.94,
  0.945,
  0.95,
  0.955,
  0.96,
  0.965,
  0.97,
  0.975,
  0.98,
  0.985,
  0.99,
  0.995,
  1.0};

// Points of Interest (POI)
static map<double, const char*> poi_names;
map<double, size_t> poi;

// POI not defined elsewhere
#define ONE_THIRD 1.0 / 3.0
#define GAMMA 0.57721566490153286060651209008240243
#define LAPLACE_LIMIT 0.66274341934918158097474209710925290
#define TWO_THIRD 2.0 / 3.0
#define ZETA3 1.2020569
#define PHI 1.61803398874989484820458683436563811
#define SQRT3 1.73205080756887729352744634150587236
#define SQRT5 2.23606797749978969640917366873127623
#define TWO_PI 2 * M_PI

#define POI_BINS 41
double poi_values[POI_BINS] =
{0.0,
 0.1,
 0.2,
 0.3,
 M_1_PI,
 ONE_THIRD,
 0.4,
 M_LOG10E,
 GAMMA,
 0.6,
 M_2_PI,
 LAPLACE_LIMIT,
 TWO_THIRD,
 M_LN2,
 0.7,
 M_SQRT1_2,
 M_PI_4,
 0.8,
 0.9,
 1.0,
 1.1,
 M_2_SQRTPI,
 1.2,
 ZETA3,
 1.3,
 1.4,
 M_SQRT2,
 M_PI_2,
 1.6,
 PHI,
 1.7,
 SQRT3,
 1.8,
 1.9,
 SQRT5,
 M_LN10,
 M_E,
 M_PI,
 TWO_PI,
 HUGE_VAL,
 INFINITY};

void clear_histograms()
{
    for (int i = 0; i < NORMAL_BINS; i++)
        normals[normal_breaks[i]] = 0;

    for (int i = 0; i < SUBNORMAL_BINS; i++)
        subnormals[subnormal_breaks[i]] = 0;

    for (int i = 0; i < FRACTION_BINS; i++)
        fractions[fraction_breaks[i]] = 0;

    for (int i = DBL_MIN_EXP; i < DBL_MAX_EXP; i++)
        exponents[i] = 0;

    for (int i = 0; i < POI_BINS; i++) {
        poi[poi_values[i]] = 0;
        poi[-poi_values[i]] = 0;
    }
}

FILE* initialize_histogram_file(int pid, const char* name, FILE* histogramOut)
{
    char fileName[128];

    sprintf(fileName, "%s-%d.dat", name, pid);
    histogramOut = fopen(fileName, "a");
    fprintf(histogramOut, "Ops Value Count\n");

    return histogramOut;
}

FILE* initialize_poi_file(int pid, FILE* histogramOut)
{
    char fileName[128];

    sprintf(fileName, "poi-%d.dat", pid);
    histogramOut = fopen(fileName, "a");
    fprintf(histogramOut, "Ops Value Name Count\n");

    return histogramOut;
}

void initialize()
{
    int pid = getpid();

    normalsOut = initialize_histogram_file(pid, "normals", normalsOut);
    subnormalsOut = initialize_histogram_file(pid, "subnormals", subnormalsOut);
    fractionsOut = initialize_histogram_file(pid, "fractions", fractionsOut);
    exponentsOut = initialize_histogram_file(pid, "exponents", exponentsOut);
    poiOut = initialize_poi_file(pid, poiOut);

    clear_histograms();

    poi_names[0.0] = "0";
    poi_names[0.1] = "0.1";
    poi_names[0.2] = "0.2";
    poi_names[0.3] = "0.3";
    poi_names[M_1_PI] = "1/pi";
    poi_names[ONE_THIRD] = "1/3";
    poi_names[0.4] = "0.4";
    poi_names[M_LOG10E] = "log10(E)";
    poi_names[GAMMA] = "gamma";
    poi_names[0.6] = "0.6";
    poi_names[M_2_PI] = "2/pi";
    poi_names[LAPLACE_LIMIT] = "Laplace limit";
    poi_names[TWO_THIRD] = "2/3";
    poi_names[M_LN2] = "ln(2)";
    poi_names[0.7] = "0.7";
    poi_names[M_SQRT1_2] = "1/sqrt(2)";
    poi_names[M_PI_4] = "pi/4";
    poi_names[0.8] = "0.8";
    poi_names[0.9] = "0.9";
    poi_names[1.0] = "1";
    poi_names[1.1] = "1.1";
    poi_names[M_2_SQRTPI] = "2/sqrt(pi)";
    poi_names[1.2] = "1.2";
    poi_names[ZETA3] = "zeta(3)";
    poi_names[1.3] = "1.3";
    poi_names[1.4] = "1.4";
    poi_names[M_SQRT2] = "sqrt(2)";
    poi_names[M_PI_2] = "pi/2";
    poi_names[1.6] = "1.6";
    poi_names[PHI] = "phi";
    poi_names[1.7] = "1.7";
    poi_names[SQRT3] = "sqrt(3)";
    poi_names[1.8] = "1.8";
    poi_names[1.9] = "1.9";
    poi_names[SQRT5] = "sqrt(5)";
    poi_names[M_LN10] = "ln(10)";
    poi_names[M_E] = "e";
    poi_names[M_PI] = "pi";
    poi_names[TWO_PI] = "2*pi";
    poi_names[HUGE_VAL] = "HUGE_VAL";
    poi_names[INFINITY] = "infinity";
}

bool isPOI(double num)
{
    for (int i = 0; i < POI_BINS; i++)
        if (fabs(num) == poi_values[i])
            return true;

    if (isinf(num))
            return true;

    return false;
}

void print_float_histogram(UINT64 ops, FILE* histogramOut, map<double, size_t> &histogram)
{
    for(map<double, size_t>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
        fprintf(histogramOut, "%lu %g %lu\n", ops, it->first, it->second);
}

void print_int_histogram(UINT64 ops, FILE* histogramOut, map<int, size_t> &histogram)
{
    for(map<int, size_t>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
        fprintf(histogramOut, "%lu %d %lu\n", ops, it->first, it->second);
}

void print_poi(UINT64 ops, FILE* histogramOut, map<double, size_t> &histogram)
{
    for(map<double, size_t>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
        if (signbit(it->first) == 0)
            fprintf(histogramOut, "%lu %g \"%s\" %lu\n", ops, it->first, poi_names[it->first], it->second);
        else
            fprintf(histogramOut, "%lu %g \"-%s\" %lu\n", ops, it->first, poi_names[-it->first], it->second);
}

void print_step(UINT64 ops)
{
    print_float_histogram(ops, normalsOut, normals);
    print_float_histogram(ops, subnormalsOut, subnormals);
    print_float_histogram(ops, fractionsOut, fractions);
    print_int_histogram(ops, exponentsOut, exponents);
    print_poi(ops, poiOut, poi);
}

void finalize(UINT64 ops)
{
    print_step(ops);
    fclose(normalsOut);
    fclose(subnormalsOut);
    fclose(fractionsOut);
    fclose(exponentsOut);
    fclose(poiOut);
}

#define UPDATE(X)   ops++; \
                    switch(fpclassify(X)) { \
                        case FP_NORMAL:    dbl_bucket = normals.lower_bound(X); \
                                           if (dbl_bucket != normals.end()) \
                                               dbl_bucket->second++; \
                                           break; \
                        case FP_SUBNORMAL: dbl_bucket = subnormals.lower_bound(X); \
                                           if (dbl_bucket != subnormals.end()) \
                                               dbl_bucket->second++; \
                                           break; \
                    } \
                    dbl_bucket = fractions.lower_bound(frexp(X, &exponent)); \
                    if (dbl_bucket != fractions.end()) \
                        dbl_bucket->second++; \
                    exponents[exponent]++; \
                    if (isPOI(X)) poi[X]++; \
                    if (ops % _ops_per_step == 0) {\
                        print_step(ops); \
                        clear_histograms(); \
                    }

#define SH_TYPE         double
#define SH_INFO         "statistics of values"
#define SH_PARAMS       (decstr(_ops_per_step))

#define SH_INIT         _ops_per_step = KnobOpsPerStep.Value(); \
                        initialize();

#define SH_ALLOC(V)     ;
#define SH_FREE(V)      ;
#define SH_SET(V,X)     (V)=(X); UPDATE(X)
#define SH_COPY(V,S)    (V)=(S); UPDATE(S)
#define SH_OUTPUT(O,V)  O << fltstr(V,15); UPDATE(V)
#define SH_FINI         finalize(ops);

#define SH_PACKED_TYPE  double
#define SH_PACK(P,V)    (P)=(V)
#define SH_UNPACK(P,V)  (V)=(P)

#define SH_ADD(V,S)     (V)=((V)+(S)); UPDATE(V); UPDATE(S)
#define SH_SUB(V,S)     (V)=((V)-(S)); UPDATE(V); UPDATE(S)
#define SH_MUL(V,S)     (V)=((V)*(S)); UPDATE(V); UPDATE(S)
#define SH_DIV(V,S)     (V)=((V)/(S)); UPDATE(V); UPDATE(S)
#define SH_MIN(V,S)     (V)=fmin((V),(S)); UPDATE(V); UPDATE(S)
#define SH_MAX(V,S)     (V)=fmax((V),(S)); UPDATE(V); UPDATE(S)
#define SH_SQRT(V,S)    (V)=sqrt(S); UPDATE(S)
#define SH_ABS(V,S)     (V)=fabs(S); UPDATE(S)

// Excluding from statistics since they will be included later
// if actually used.
#define SH_NEG(V,S)     (V)=(-(S))
#define SH_RND(V,S)     (V)=round(S)

#define SH_DBL(V)       (V)

inline double _relerr(double v, double x)
{
    return ((x == 0.0) ? (v == 0.0 ? 0.0 : 1.0) : (fabs((double)x-v) / fabs(x)));
}

// native64 only: must match exactly or both be NaNs
//
inline bool _iserr(double v, double x)
{
    bool vnan = isnan(v);
    bool xnan = isnan(x);
    return ( vnan && !xnan) ||
           (!vnan &&  xnan) ||
           (!vnan && !xnan && v != x);
}

#define SH_RELERR(V,X)  (_relerr((V),(X)))
#define SH_ISERR(V,X)   (_iserr((V),(X)))

#include "shval.cpp"

