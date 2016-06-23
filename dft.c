#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <malloc.h>

#define M_PI 3.141592653589793

typedef complex double cplx;

void print_complex_array(cplx input[],int length);
double *bitrevorder(double *i,int n);
cplx *recursive_fft(cplx *a,int n);
void print_complex_array_index_value(cplx input[],int length,int index);

int main()
{
    cplx *dft=(cplx *)malloc(512*sizeof(cplx));
    cplx *c=(cplx *)malloc(512*sizeof(cplx));
    int i;
    for(i=0;i<8;i++)
        *(c+i)=1;
//    print_array(a,8);
    dft=recursive_fft(c,8);
    print_complex_array_index_value(c,8,2);
    print_complex_array_index_value(dft,8,2);
    return 0;
}

void print_complex_array(cplx input[],int length)
{
    int i;
    for(i=0;i<length;i++)
        printf("%f+j%f\n",creal(*(input+i)),cimag(*(input+i)));
}

void print_complex_array_index_value(cplx input[],int length,int index)
{
    if(index<length && index>=0)
        printf("%f+j%f\n",creal(*(input+index)),cimag(*(input+index)));
}

cplx *recursive_fft(cplx *a,int n)
{
    cplx *y=(cplx *)malloc(n*sizeof(cplx));
    cplx wn=cexp(2*M_PI*I/n);
    cplx w=1;
    cplx p,q,*real;
    cplx *y0=(cplx *)malloc((n/2)*sizeof(cplx));
    cplx *y1=(cplx *)malloc((n/2)*sizeof(cplx));
    cplx *a0=(cplx *)malloc((n/2)*sizeof(cplx));
    cplx *a1=(cplx *)malloc((n/2)*sizeof(cplx));
    int i,j,k;

    if(n==1)
    {
        return a;
    }

    j=0;
    for(i=0;i<n/2;i++)
    {
        *(a0+i)=*(a+j);
        j=j+2;
    }

    j=1;
    for(i=0;i<n/2;i++)
    {
        *(a1+i)=*(a+j);
        j=j+2;
    }

    y0=recursive_fft(a0,n/2);
    y1=recursive_fft(a1,n/2);

    free(a0);
    free(a1);

    for(k=0;k<(n/2);k++)
    {
        p=*(y0+k);
        q=*(y1+k);
        *(y+k)=p+w*q;
        *(y+k+(n/2))=p-w*q;
        w=w*wn;
    }
//    free(a0);
//    free(a1);
    free(y0);
    free(y1);
    return y;
}
