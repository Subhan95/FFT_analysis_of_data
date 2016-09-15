#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

#define M_PI 3.14159265358979323846

typedef complex double cplx;

void print_array(double *input,int length);
void print_array_index_value(double *input,double length,int index);
void print_complex_array(cplx input[],int length);
void print_complex_array_index_value(cplx input[],int length,int index);
void file_write(char s[],double a[],int n);

double *linspace(double start,double stop,double points);
double *generate_carrier(double amplitude,double fc1,double fc2,double fc3,double t[],double len_of_t);
double *dot_multiply(double a[],double b[],double n,int increment);
double *fftshift(double *a,int n);

double drand48(void);
double randn(double mean, double var);

cplx *recursive_fft(cplx *a,int n);



int main()
{
	int h,i,j,k,m;
    int radix,increment;
    int start_index,stop_index;
    double val;
	double snrdb;
	double signal_power,noise_power;
	double amplitude;
	double fs=2*pow(10,9);
    double fc1=50*pow(10,6),fc2=100*pow(10,6),fc3=200*pow(10,6);
    double len_of_t=6000;
    double threshold=0.52;
    double *t=linspace(0,3*pow(10,-6),6*pow(10,3));
    double *carrier=(double *)malloc((int)len_of_t*sizeof(double));
    double data[(int)len_of_t];

    printf("Enter radix of fft (16,32,64,128,256 or 512): ");
    scanf("%d",&radix);

    if((int)len_of_t % radix == 0)
        increment=0;
    else
    {
        increment=radix-((int)len_of_t % radix);
    }


    double *signal=(double *)malloc(((int)len_of_t+increment)*sizeof(double));
    double noise[(int)len_of_t+increment];
    double *actual_signal=(double *)malloc(((int)len_of_t+increment)*sizeof(double));
    double *quant_signal=(double *)malloc(((int)len_of_t+increment)*sizeof(double));
    double *shifted_signal=(double *)malloc(((int)len_of_t+increment)*sizeof(double));
    double *quant=linspace(-0.25,0.25,256);
    double neg_freq[((int)len_of_t+increment)/2];
    double pulse_width;

    cplx quantized_data[radix];
    cplx *fft=(cplx *)malloc(radix*sizeof(cplx));

    FILE *fp,*fd;

	noise_power=1.907*pow(10,-10);

	printf("Enter snr in dB: ");
    scanf("%lf",&snrdb);

    signal_power=noise_power*pow(10,(0.1*snrdb));
    amplitude=sqrt((2.0/3)*signal_power)*200;

    file_write("/home/subhan/t.txt",t,(int)len_of_t);

    carrier=generate_carrier(amplitude,fc1,fc2,fc3,t,len_of_t);

    file_write("/home/subhan/carrier.txt",carrier,(int)len_of_t);

    for(i=0;i<(int)len_of_t;i=i+1)
    {
        if(i<(int)len_of_t/3)
            data[i]=0;
        else if(i<2.0*(int)len_of_t/3)
            data[i]=1;
        else if(i<(int)len_of_t)
            data[i]=0;
    }

    file_write("/home/subhan/data.txt",data,(int)len_of_t);

    signal=dot_multiply(data,carrier,(int)len_of_t,increment);

    free(carrier);

    file_write("/home/subhan/signal.txt",signal,(int)len_of_t+increment);

    srand( (unsigned)time(NULL));
    for(i=0;i<len_of_t+increment;i++)
        noise[i]=sqrt(noise_power)*randn(0,1);

    file_write("/home/subhan/noise.txt",noise,(int)len_of_t+increment);

    for(i=0;i<(int)len_of_t+increment;i++)
    {
        actual_signal[i]=signal[i]+noise[i];
    }

    file_write("/home/subhan/actual_signal.txt",actual_signal,(int)len_of_t+increment);

    file_write("/home/subhan/quant.txt",quant,256);

    for(i=0;i<(int)len_of_t+increment;i++)
    {
    	for(j=0;j<256;j++)
    	{
    		if(j==255)
            {
                quant_signal[i]=quant[255];
            }
    		else if(actual_signal[i]>=quant[j] && actual_signal[i]<quant[j+1])
    		{
    			quant_signal[i]=quant[j+1];
    			break;    		
            }
    	}
    }

    file_write("/home/subhan/quantized.txt",quant_signal,(int)len_of_t+increment);

    fp=fopen("/home/subhan/quantized.txt","r");
    for(i=0;i<12;i++)
    {
        start_index=i*radix;
        stop_index=start_index+radix-1;
        k=0;
        for(j=start_index;j<=stop_index;j++)
        {
            fscanf(fp,"%16le\n",&val);
            quantized_data[k++]=val+0*I;
        }
        fft=recursive_fft(quantized_data,radix);
        if(i==0)
        {
            fd=fopen("/home/subhan/fftdata.txt","w");
            for(h=0;h<radix;h++)
                fprintf(fd,"%e\n",cabs(*(fft+h)));
            fclose(fd);
        }
        else
        {
            fd=fopen("/home/subhan/fftdata.txt","a");
            for(h=0;h<radix;h++)
                fprintf(fd,"%e\n",cabs(*(fft+h)));
            fclose(fd);
        }
    }
    fclose(fp);

    k=0;
    fd=fopen("/home/subhan/fftdata.txt","r");
    for(i=0;i<(int)len_of_t+increment;i++)
    {
        fscanf(fd,"%le\n",&val);
        actual_signal[k++]=val;
    }
    fclose(fd);
    shifted_signal=fftshift(actual_signal,(int)len_of_t+increment);

    file_write("/home/subhan/fftshift.txt",shifted_signal,(int)len_of_t+increment);

    for(i=0;i<12;i++)
    {
        start_index=i*radix;
        stop_index=start_index+radix-1;
        for(j=start_index;j<=stop_index;j++)
        {
            if(shifted_signal[j]>threshold)
            {
                m=m+1;
                break;   
            }
        }
    }

    pulse_width=(m*radix)/fs;
    printf("Pulse width is %e\n",pulse_width);

	return 0;
}

double *generate_carrier(double amplitude,double fc1,double fc2,double fc3,double t[],double len_of_t)
{
    int i;
    double *ans = (double *)malloc((len_of_t+1)*sizeof(double));
    for(i=0;i<len_of_t;i++)
        *(ans+i)=amplitude*(sin(2*M_PI*fc1*t[i])+sin(2*M_PI*fc2*t[i])+sin(2*M_PI*fc3*t[i]));
    return ans;
}

double *linspace(double start,double stop,double points)
{
    int i;
    double *ans = malloc(points*sizeof(double));
    double increment=(stop-start)/(points-1);

    for(i=0;i<points;i++)
        *(ans+i)=start+i*increment;
    return ans;
}

void print_array(double *input,int length)
{
    int i;
    for(i=0;i<length;i++)
    {
        printf("%e\n",*(input+i));
    }
}

void print_array_index_value(double *input,double length,int index)
{
    if(index<length && index>=0)
        printf("%e\n",*(input+index));
}

double *dot_multiply(double a[],double b[],double n,int increment)
{
    int i;
    double *ans =(double *)malloc(((int)n+increment)*sizeof(double));
    for(i=0;i<n;i++)
        *(ans+i)=a[i]*b[i];
    for(i=n;i<(int)n+increment;i++)
    	*(ans+i)=0;
    return ans;
}

double drand48(void)
{
    return (float)rand()/(float)RAND_MAX;
}

double randn(double mean, double var)
{
    return sqrt(var)*cos(2*M_PI*drand48())*sqrt(abs(2.0*log(drand48())))+mean;
}

void print_complex_array(cplx input[],int length)
{
   int i;
   for(i=0;i<length;i++)
       printf("%e+j%e\n",creal(*(input+i)),cimag(*(input+i)));
}

void print_complex_array_index_value(cplx input[],int length,int index)
{
   if(index<length && index>=0)
   {
        printf("%e+j%e\n",creal(*(input+index)),cimag(*(input+index)));
        printf("%e\n",cabs(*(input+index)));
   }
       
}

cplx *recursive_fft(cplx *a,int n)
{
    cplx *y=(cplx *)malloc(n*sizeof(cplx));
    cplx wn=cexp(2*M_PI*I/n);
    cplx w=1;
    cplx p,q;
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

    for(k=0;k<(n/2);k++)
    {

        p=*(y0+k);
        q=*(y1+k);
        *(y+k)=p+w*q;
        *(y+k+(n/2))=p-w*q;
        w=w*wn;
    }
    return y;
}

void file_write(char s[],double a[],int n)
{
    FILE *fp;
    int i;
    fp=fopen(s,"w");
    for(i=0;i<n;i++)
        fprintf(fp,"%e\n",a[i]);
    fclose(fp);
}

double *fftshift(double *a,int n)
{
    int i;
    double key;
    double *ans=(double *)malloc(n*sizeof(double));
    int p=n/2;
    if(n%2==0)
    {
        for(i=0;i<p;i++)
        {
            key=*(a+i);
            *(ans+i)=*(a+i+p);
            *(ans+i+p)=key;
        }
    }
    
    else
    {
        for(i=0;i<p;i++)
        {
            key=*(a+i);
            *(ans+i)=*(a+i+p+1);
            *(ans+i+p+1)=key;
        }
        *(ans+i)=*(a+p);

        key=*(a+p);
        for(i=p;i<n-1;i++)
        {
            *(ans+i)=*(ans+i+1);
        }
        *(ans+i)=key;
    }
    return ans;
}
