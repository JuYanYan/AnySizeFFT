# AnySizeFFT

Any size FFT! It's not just an integer power of two.
 1. class FFT0 : Radix-2 DIT FFT.
 2. class FFT1 : Any size FFT, using Bluestein.
 3. class FFT2 : Two dimensional FFT with a buffer.

I wanted to use it in my image process library. This library was completed at my birthday. However, after that, I have decided to use DCT. So, I public it.

## Performance

Use MS VC++ Compiler with O2 Optimization:

![Performance](https://github.com/JuYanYan/AnySizeFFT/blob/master/classFFT1.png)
