/* NICK CAIN */
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
double gasdev();
void init_genrand(unsigned long s);

void rush_rinzel_1step(double V_start, double *V_new, double b_start, double *b_new, double n_start, double *n_new, 
		double shared, unsigned long *i_seed, double sharedInput, double dt, double tMax, double g_a, double mu, double  sigma) {
	//init state
	double V, B, N;
	double C = 1;
	double c_m = 0.01;
	double g_l = 0.003;
	double g_na = 1.2;
	double g_k = 0.20;
	double E_L = -17;
	double E_na = 55;
	double E_k = -72;
	double E_a = -75;
	double a_m = 0.38*(V_start +29.7)/(1-exp(-0.1*(V_start +29.7)));
	double a_n = 0.02*(V_start +45.7)/(1-exp(-0.1*(V_start +45.7)));
	double b_m = 15.2*exp(-0.0556*(V_start+54.7));
	double b_n = 0.25*exp(-0.0125*(V_start+55.7));
	double tau_n = 1/(a_n + b_n);;
	double tau_b = .124 + 2.678/(1+exp(0.0624*(V_start+50)/16.027));
	double m_inf = a_m / (a_m + b_m);
	double n_inf = a_n / (a_n + b_n);
	double a_inf = pow(0.0761*exp(0.0314*(V_start+94.22))/(1+exp(0.0346*(V_start+1.17))),1.0/3.0);
	double b_inf = pow(1/(1+exp(0.0688*(V_start+53.3))),4);
	double dndt = (n_inf - n_start)/tau_n;
	double dbdt = (b_inf - b_start)/tau_b;
	double I = g_l*(V_start-E_L)  + g_na*pow(m_inf,3)*(0.9 - 1.2*n_start)*(V_start-E_na) + g_k*pow(n_start,4)*(V_start-E_k) + g_a*pow(a_inf,3)*b_start*(V_start-E_a);
	//init rand
	double unique_stoch = gasdev();
	//printf("%f %f\n", shared, unique_stoch);
//output: k = [1]
	V = V_start + (dt*((-I + mu/C)/c_m) + sqrt(dt)*sigma*(sqrt(sharedInput)*shared + sqrt(1-sharedInput)*unique_stoch));
	*V_new = V;
//output: k = [3]
	B = b_start + (dt*(dbdt));
	*b_new = B;
//output: k = [2]
	N = n_start + (dt*(dndt));
	*n_new = N;
}

int main(int argc, char** argv){
	double V_start[2], V_new[2];
	double b_start[2], b_new[2];
	double n_start[2], n_new[2];

        //double dt=0.01,tMax = 500;
        //double mu=.47; double sigma=0.10;double sharedInput = 0.50;
        int i,j;
	unsigned long i_seed= ((unsigned long)time(NULL));
	double shared;
	if (argc < 5) { 
		printf("oops, not enough args\n");
		printf("usage: %s sharedInput  dt  tMax  g_a i_seed\n",argv[0]);
		printf("stdin vars: mu  sigma\n");
		printf("configured actions: printspiketimes\n");
		return 1;
	}
	char *endp;
	double sharedInput = strtod(argv[1],&endp);
	double  dt = strtod(argv[2],&endp);
	double  tMax = strtod(argv[3],&endp);
	double  g_a = strtod(argv[4],&endp);
	//unsigned long i_seed = atoi(argv[5]);
	init_genrand(i_seed);
	double mu;
	double sigma;
	scanf("%lf %lf",&mu, &sigma);
/* initial conditions */
	 for (j=0; j<2; j++) {
		V_start[j]=-65;
		b_start[j]=.5;
		n_start[j]=.5;
		V_new[j]=-65;
		b_new[j]=.5;
		n_new[j]=.5;
	}
	for (i=0; i*dt <= tMax; i++) {
		for (j=0; j<2; j++) {
			 /* Optional Spike Detection*/
//settings = {'static_g_k': '0.20', 'static_E_a': '-75', 'static_g_l': '0.003', 'static_a_inf': 'pow(0.0761*exp(0.0314*(V_start+94.22))/(1+exp(0.0346*(V_start+1.17))),1.0/3.0)', 'actions': 'printspiketimes', 'staticparams': 'C, c_m, g_l, g_na, g_k, E_L, E_na, E_k, E_a, a_m, a_n, b_m, b_n, tau_n, tau_b, m_inf, n_inf, a_inf, b_inf, dndt, dbdt, I', 'static_a_m': '0.38*(V_start +29.7)/(1-exp(-0.1*(V_start +29.7)))', 'static_a_n': '0.02*(V_start +45.7)/(1-exp(-0.1*(V_start +45.7)))', 'static_b_n': '0.25*exp(-0.0125*(V_start+55.7))', 'static_b_inf': 'pow(1/(1+exp(0.0688*(V_start+53.3))),4)', 'static_b_m': '15.2*exp(-0.0556*(V_start+54.7))', 'static_tau_n': '1/(a_n + b_n);', 'static_E_L': '-17', 'static_tau_b': '.124 + 2.678/(1+exp(0.0624*(V_start+50)/16.027))', 'static_g_na': '1.2', 'static_C': '1', 'stdinparams': 'mu, sigma', 'static_dbdt': '(b_inf - b_start)/tau_b', 'static_n_inf': 'a_n / (a_n + b_n)', 'static_I': 'g_l*(V_start-E_L)  + g_na*pow(m_inf,3)*(0.9 - 1.2*n_start)*(V_start-E_na) + g_k*pow(n_start,4)*(V_start-E_k) + g_a*pow(a_inf,3)*b_start*(V_start-E_a)', 'scale_eq1': '', 'scale_eq2': '', 'scale_eq3': '', 'static_dndt': '(n_inf - n_start)/tau_n', 'spike_var': '@V', 'static_E_k': '-72', 'static_m_inf': 'a_m / (a_m + b_m)', 'static_c_m': '0.01', 'spike_threshold': '0', 'static_E_na': '55', 'ic_var1': '-65', 'ic_var3': '.5', 'ic_var2': '.5', 'spike_trigger': 'threshold'}
			if ((V_new[j]>0)&&(V_start[j]<=0)) {
				printf("%d %f\n",j+1,dt*i);
			}
			V_start[j] = V_new[j];
			b_start[j] = b_new[j];
			n_start[j] = n_new[j];
		}

	shared = gasdev();
		for (j=0; j<2; j++) {
			rush_rinzel_1step(V_start[j], &V_new[j], b_start[j], &b_new[j], n_start[j], &n_new[j],
				shared, &i_seed, sharedInput, dt, tMax, g_a, mu, sigma);

                }
        }
        return 0;
}

/*
A C-program for MT19937, with initialization improved 2002/1/26.
Coded by Takuji Nishimura and Makoto Matsumoto.

Before using, initialize the state by using init_genrand(seed)
or init_by_array(init_key, key_length).

Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. The names of its contributors may not be used to endorse or promote
products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Any feedback is very welcome.
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
mt[0]= s & 0xffffffffUL;
for (mti=1; mti<N; mti++) {
mt[mti] =
(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
/* In the previous versions, MSBs of the seed affect   */
/* only MSBs of the array mt[].                        */
/* 2002/01/09 modified by Makoto Matsumoto             */
mt[mti] &= 0xffffffffUL;
/* for >32 bit machines */
}
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
int i, j, k;
init_genrand(19650218UL);
i=1; j=0;
k = (N>key_length ? N : key_length);
for (; k; k--) {
mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
+ init_key[j] + j; /* non linear */
mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
i++; j++;
if (i>=N) { mt[0] = mt[N-1]; i=1; }
if (j>=key_length) j=0;
}
for (k=N-1; k; k--) {
mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
- i; /* non linear */
mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
i++;
if (i>=N) { mt[0] = mt[N-1]; i=1; }
}

mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
unsigned long y;
static unsigned long mag01[2]={0x0UL, MATRIX_A};
/* mag01[x] = x * MATRIX_A  for x=0,1 */

if (mti >= N) { /* generate N words at one time */
int kk;

if (mti == N+1)   /* if init_genrand() has not been called, */
init_genrand(5489UL); /* a default initial seed is used */

for (kk=0;kk<N-M;kk++) {
y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
}
for (;kk<N-1;kk++) {
y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
}
y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

mti = 0;
}

y = mt[mti++];

/* Tempering */
y ^= (y >> 11);
y ^= (y << 7) & 0x9d2c5680UL;
y ^= (y << 15) & 0xefc60000UL;
y ^= (y >> 18);

return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
return genrand_int32()*(1.0/4294967295.0);
/* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
return genrand_int32()*(1.0/4294967296.0);
/* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
/* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
double gasdev()
{
double genrand_res53();
static int iset=0;
static double gset;
double fac,rsq,v1,v2;

if (iset == 0) {
do {
v1=2.0*genrand_res53()-1.0;
v2=2.0*genrand_res53()-1.0;
rsq=v1*v1+v2*v2;
} while (rsq >= 1.0 || rsq == 0.0);
fac=sqrt(-2.0*log(rsq)/rsq);
gset=v1*fac;
iset=1;
return (double)(v2*fac);
} else {
iset=0;
return (double)gset;
}
}

/* These real versions are due to Isaku Wada, 2002/01/09 added */
/*
int main(void)
{
int i;
unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
/*init_by_array(init, length);* /
long i_seed= -1*((unsigned int)time(NULL));
init_genrand(i_seed);

printf("1000 outputs of genrand_int32()\n");
for (i=0; i<1000; i++) {
printf("%10lu ", genrand_int32());
if (i%5==4) printf("\n");
}
printf("\n1000 outputs of genrand_real2()\n");
for (i=0; i<1000; i++) {
printf("%10.8f ", genrand_real2());
if (i%5==4) printf("\n");
}
return 0;
}
*/
