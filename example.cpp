#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "fft.h"
#define MAX_N          20
FFT1::complex dat[MAX_N];
void genTestWave(FFT1::complex *dat, int sz);
void printComplex(const FFT0::complex *dat, int sz);
void printComplex(const FFT1::complex *dat, int w, int h);
void demo_FFT0(int sz_index)
{
    FFT0 fft;                   // Use FFT0 class to run fft of 2^sz_index.
    assert(1 << sz_index <= MAX_N);
    // Our test data
    genTestWave(dat, 1 << sz_index);
    puts("FFT0 : Source");
    printComplex(dat, 1 << sz_index);
    // Init. After that, you can use "Resize" to change size.
    fft.InitFFT(sz_index);
    // For FFT0:
    // Step1. BRO
    fft.BRO(dat);
    // Step2. FFT
    fft.FFT(dat);
    puts("FFT0 : FFT");
    printComplex(dat, 1 << sz_index);
    // Then, Let's do IFFT.
    // Step1. BRO
    fft.BRO(dat);
    // Step2. IFFT
    fft.IFFT(dat);
    puts("FFT0 : IFFT");
    printComplex(dat, 1 << sz_index);
    // fft.~FFT0();             // You can explicitly deconstruct it to release resource.
}
void demo_FFT1(int sz)
{
    FFT1 fft;                   // Use FFT0 class to run fft of any size.
    assert(sz <= MAX_N);
    // Our test data
    genTestWave(dat, sz);
    puts("FFT1 : Source");
    printComplex(dat, sz);
    // Init. After that, you can use "Resize" to change size.
    fft.InitFFT(sz);
    // FFT
    fft.FFT(dat);
    puts("FFT1 : FFT");
    printComplex(dat, sz);
    // Let's do IFFT.
    fft.IFFT(dat);
    puts("FFT1 : IFFT");
    printComplex(dat, sz);
    // fft.~FFT1();
}
void demo_FFT2(int w, int h)
{
    FFT2 fft;                   // Use FFT0 class to run 2d-fft of any size.
    assert(w * h <= MAX_N);
    // Our test data
    // We use 1d data on 2d fft. Please note that we want to show how to use it.
    // Instead of being a complete application
    genTestWave(dat, w * h);
    puts("FFT2 : Source");
    printComplex(dat, w, h);
    // Init. After that, you can use "Resize" to change size.
    fft.InitFFT(w, h);
    // FFT
    // Do you want to use buffer? 
    // Step1. You can use "fft.BufferLoadReal" to push data.
    // Step2. Then, use "fft.FFT()" or "fft.FFT(fft.GetBufferPtr())" to run FFT
    // Step3. In addition, use "fft.BufferGetReal" or "fft.GetBufferPtr" to get data.
    fft.FFT(dat);
    // Using fft.Shift() to shift data
    fft.Shift();
    puts("FFT2 : FFT");
    printComplex(dat, w, h);
    // Due to we use fft.Shift, so, shift data firstly.
    fft.Shift();
    // Then, we can use IFFT
    fft.IFFT(dat);
    puts("FFT2 : IFFT");
    printComplex(dat, w, h);
    // fft.~FFT2();
}
// Generate test data
// Result like these:
// Real:
// |\ Value
// -
// -   |             |
// -   |             |
// -   | |         | |
// -+-+-+-+-+-+-+-+-+-+--> freq
// Imag:
// |\ Value
// -
// -   |             |
// -   |             |
// -   |             |
// --+-+-+-+-+-+-+-+-+--> freq
void genTestWave(FFT1::complex *dat, int sz)
{
    int i;
    double t;
    for (i = 0; i < sz; i++)
    {
        t = i * FFT1::PI2 / sz;
        dat->real = sin(t) + sin(3 * t) / 3 + cos(t);
        dat->imag = 0;
        dat++;
    }
}
// Print output
void printComplex(const FFT1::complex *dat, int sz)
{
    while(sz-- > 0)
    {
        if (dat->imag < 0)
            printf("%+lf - %lfi\n", dat->real, -dat->imag);
        else
            printf("%+lf + %lfi\n", dat->real, +dat->imag);
        dat++;
    }
}
void printComplex(const FFT1::complex *dat, int w, int h)
{
    int k;
    while(h -- > 0)
    {
        k = w;
        while (k-- > 0)
        {
            if (dat->imag < 0)
                printf("%+lf - %lfi, ", dat->real, -dat->imag);
            else
                printf("%+lf + %lfi, ", dat->real, +dat->imag);
            dat++;
        }
        putchar('\n');
    }
}
// Main proc
int main()
{
    demo_FFT0(3);
    demo_FFT1(10);
    demo_FFT2(3, 4);

    puts("Press 'Enter' to exit...");
    getchar();
    return 0;
}
