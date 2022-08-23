#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pgm.h>
#include <time.h>

using namespace std;


float ricker (float t, float fpeak)
{

float x,xx;

x = M_PI*fpeak*(t-(1.0/fpeak));
xx = x*x;
/* return (-6.0+24.0*xx-8.0*xx*xx)*exp(-xx); */
/* return PI*fpeak*(4.0*xx*x-6.0*x)*exp(-xx); */
return exp(-xx)*(1.0-2.0*xx);
}

// rdo w czasie i przestrzeni

void ptsrc (float xs, float zs,
	        int nx  , float dx,
            int nz  , float dz,
            float t , float fpeak, float **s)
{

int ix,iz,ixs,izs;
float ts,xn,zn,xsn,zsn;


dtr2 = dtr*dtr;
ds2 = ds*ds;



/* zero source array */
for (ix=0; ix<=nx; ++ix)
for (iz=0; iz<=nz; ++iz)
s[ix][iz] = 0.0;


ts = ricker(t,fpeak);


xsn = (xs)/dx;
zsn = (zs)/dz;

//zamiana na inty
ixs = (int)(xsn);
izs = (int)(zsn);



for (ix=ixs-3; ix<=ixs+3; ++ix)
{
 for (iz=izs-3; iz<=izs+3; ++iz)
 {
 xn = ix-xsn;
 zn = iz-zsn;
 s[ix][iz] = ts*exp(-xn*xn-zn*zn);
 }
}


}//koniec funkcji

void fdmod1 (int nx, int nz, float ds,
             float xs, float zs,
             float fpeak, float dt, float et,
             float **V,int surf,int prec,
             float **seis)
{
int i,j,k,kk,kkk;
float dx=ds,dz=ds,**s,**pm,**p,**pp,dtr,vmax=-9999,t,w1,w2;
FILE *d;
char nazwa[25];
struct timeval tv;
float start,stop;

p=(float  **)malloc ((nx+1)*sizeof(float *));
    for (i=0;i<=nx;i++)
    p[i]=(float *)malloc ((nz+1)*sizeof(float));

pm=(float  **)malloc ((nx+1)*sizeof(float *));
    for (i=0;i<=nx;i++)
    pm[i]=(float *)malloc ((nz+1)*sizeof(float));

pp=(float  **)malloc ((nx+1)*sizeof(float *));
    for (i=0;i<=nx;i++)
    pp[i]=(float *)malloc ((nz+1)*sizeof(float));

s = (float  **)malloc ((nx+1)*sizeof(float *));
for (i=0;i<=nx;i++)
s[i] = (float *)malloc ((nz+1)*sizeof(float));


for (i=0;i<=nx;i++) //szukanie vmax
for (j=0;j<=nz;j++)
if (vmax<V[i][j]) vmax=V[i][j];


for (i=0;i<=nx;i++) //zerowanie sejsmogramu
for (j=0;j<=et/dt+1;j++)
seis[i][j] = 0.0;


dtr = ds/(2.0*vmax);   // dt do obliczen dla stabilnosci
w2 = 0;

while(1)// TU ZMIENIC NA WARUNEK Z IFA
{
w2 = w2+1;
w1 = dt/w2;

 if (w1<=dtr)
 {
 dtr = w1;
 break;
 }

}


ptsrc (xs,zs,nx,dx,nz,dz,0,fpeak,s);  // rozwizanie dola czasu 0
for (i=0;i<=nx;i++)
for (j=0;j<=nz;j++)
pm[i][j] = s[i][j];

ptsrc (xs,zs,nx,dx,nz,dz,dtr,fpeak,s); // rozwiazanie dla czasu dtr
for (i=0;i<=nx;i++)
for (j=0;j<=nz;j++)
p[i][j] = s[i][j];

gettimeofday(&tv, NULL);
start = tv.tv_sec+tv.tv_usec/1000000.0;
//printf ("%f \n",start);
kk = 1; // ilosc dtr na jedna dt
kkk = 0; // ilosc dt
k = 1;   // ilosc dtr

while (1)
{

k++;
kk++;
t = (float)(k)*dtr;
ptsrc (xs,zs,nx,dx,nz,dz,t,fpeak,s);


//gwna p?la do modelowa

if (prec==0)//2th order
{
 for (i=2; i<=nx-2; i++)
 {
  for (j=2; j<=nz-2; j++)
  {


	pp[i][j] = 2.0*p[i][j]-pm[i][j] + ((dtr2)/(ds))*V[i][j]*V[i][j]* ( p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]-4.0*p[i][j]) + s[i][j];// TU PEWNIE ZOPTYMALIZOWAC

  }
 }
}
else  //4th order
	{

    #pragma omp parallel for private(i,j)
 	for (i=2; i<=nx-2; i++)
 		{
  			for (j=2; j<=nz-2; j++)
  		{
 			 pp[i][j] = ( 2.0 - 5.0*pow( (V[i][j]*dtr)/ds , 2 ) ) * p[i][j] - pm[i][j] +  (4.0/3.0) * pow( (V[i][j]*dtr)/ds , 2 ) *(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]) - (1.0/12.0)* pow( (V[i][j]*dtr)/ds , 2 ) * (p[i+2][j]+p[i-2][j]+p[i][j+2]+p[i][j-2]) + s[i][j];// TU POPATRZEC
  }
 }
}

//pomicncze p?le
for (j=1; j<=nz-1; j++)
{
	i=1;
//2th order
		pp[i][j] = 2.0*p[i][j]-pm[i][j] + ((dtr2)/(ds2))*V[i][j]*V[i][j]* ( p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]-4.0*p[i][j]) + s[i][j];
	i=nx-1;
		pp[i][j] = 2.0*p[i][j]-pm[i][j] + ((dtr2)/(ds2))*V[i][j]*V[i][j]*( p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]-4.0*p[i][j]) + s[i][j];
}

for (i=1; i<=nx-1; i++)
{
	j=1;
//2th order
		pp[i][j] = 2.0*p[i][j]-pm[i][j] + ((dtr2)/(ds2))*V[i][j]*V[i][j]*( p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]-4.0*p[i][j]) + s[i][j];
        j = nz-1;
		pp[i][j] = 2.0*p[i][j]-pm[i][j] + ((dtr2)/(ds2))*V[i][j]*V[i][j]*( p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]-4.0*p[i][j]) +  s[i][j];
}

//trasparent lewa
for (j=1;j<=nz;j++)
{
	pp[0][j] = p[0][j] + p[1][j] - pm[1][j] +  V[0][j]*(dtr/ds)*(p[1][j]-p[0][j]-(pm[2][j]-pm[1][j]));
}

//trasparent prawa
for (j=1;j<=nz;j++)
{
	pp[nx][j]= p[nx][j] + p[nx-1][j] - pm[nx-1][j] +  V[nx][j]*(dtr/ds)*(p[nx-1][j]-p[nx][j]-(pm[nx-2][j]-pm[nx-1][j]));
}

//trasparent gora
if (surf==0)
{
 for (i=0;i<=nx;i++)
 {
 	pp[i][0] = p[i][0] + p[i][1] - pm[i][1] +  V[i][0]*(dtr/ds)*(p[i][1]-p[i][0]-(pm[i][2]-pm[i][1]));
 }
}
//trasparent dol

for (i=0;i<=nx;i++)
{
	pp[i][nz] = p[i][nz] + p[i][nz-1] - pm[i][nz-1] + V[i][nz]*(dtr/ds)*(p[i][nz-1]-p[i][nz]-(pm[i][nz-2]-pm[i][nz-1]));
}


for (i=0; i<=nx; i++)
{
 for (j=0; j<=nz; j++)
 {
	pm[i][j] = p[i][j];
 	p[i][j] = pp[i][j];
 }
}

if ( (float)kk*dtr + dtr/10.0 >= dt )
{
kk = 0;
kkk++;

  for (i=1;i<=nx-1;i++) //dodanie probek sejsmogramu
 	 seis[i][kkk] = pp[i][1];

if ((float)kkk*dt > et) break;
}

}//czas

gettimeofday(&tv, NULL);
stop = tv.tv_sec+tv.tv_usec/1000000.0;
printf ("%f \n",stop);
printf ("CPU TIME: %f \n",stop-start);

 for (i=0;i<=nx;i++)
 	free(p[i]);
 	free(p);

 for (i=0;i<=nx;i++)
	free(pm[i]);
 	free(pm);

 for (i=0;i<=nx;i++)
 	free(pp[i]);
 	free(pp);

 for (i=0;i<=nx;i++)
	free(s[i]);
 	free(s);


}//od funkcji


int main (int argc, char *argv[])
{


double start,stop;
struct timeval tv;

int nx,nz; //wymiary modelu
float ds,xs,zs,fpeak,dt,et;
//ds - krok, xs,zs - lokalizacja zrodla,
//fpeak czestotliwosc dt - krok czasowy rejestracji et - koniec rejestracji


float **V;
//V- model predkosciowy osrodka


int surf, prec;
//surf - warunki brzegowe na powierzchnim, prec- rozwazanie drugiego lub czwartego zedu


float **seis;
//seis macierz dla wyników

int i, j, nt; //nt - ilosc probek czasowych



//dane poczatkowe
nz = atoi(argv[1]);   //malutki model testowy -> zwykle nx nz to K000
nx = atoi(argv[1]);
fpeak = 30.0;
dt = 0.002;
et = 0.5;
surf = 0;
prec = 0;
ds = 1.0;
xs = nx/2.0;
zs = 25;

nt = et/dt + 1;

V=(float  **)malloc ((nx+1)*sizeof(float *));
	 for (i=0;i<=nx;i++)
		 V[i] = (float *)malloc ((nz+1)*sizeof(float));

seis=(float  **)malloc ((nx+1)*sizeof(float *));

 for (i=0;i<=nx;i++)
	 seis[i]=(float *)malloc ((nt+1)*sizeof(float));


for (i=0;i<=nx;i++)			//model z dwoma warstwami
{
	for (j=0;j<=nz/2;j++)
        V[i][j]=1000.0;
	for (j=nz/2+1;j<=nz;j++)
        V[i][j]=2000.0;
}


gettimeofday(&tv, NULL);
start = tv.tv_sec+tv.tv_usec/1000000.0;
fdmod1 (nx,nz,ds, xs,zs, fpeak, dt, et, V,surf,prec, seis);
gettimeofday(&tv, NULL);
stop = tv.tv_sec+tv.tv_usec/1000000.0;
printf ("%lf\n",stop-start);

//-----Co nam wyszlo-------
gray **p;
float min,max;
FILE *d;

p = pgm_allocarray(nx+1,nt+1);

min = seis[1][1];
max = seis[1][1];

 for (i=0;i<=nx;i++)
 {
  for (j=0;j<=nt;j++)
  {
  if (seis[i][j]<min) min = seis[i][j];
  if (seis[i][j]>max) max = seis[i][j];
  }
 }

printf ("%f %f\n",min,max);

 for (i=0;i<=nx;i++)
 {
  for (j=0;j<=nt;j++)
  {
  p[j][i]=255.0-(seis[i][j]-min)/(max-min)  * 255.0;
  }
 }

d = popen("display","w");
pgm_writepgm(d,p,nx,nt,255,0);
pclose (d);

//-------------------------------

pgm_freearray(p,nt+1);

for (i=0;i<=nx;i++)
free(V[i]);
free(V);

for (i=0;i<=nx;i++)
free(seis[i]);
free(seis);

return 0;
}






