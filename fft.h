/*
 | Mark
 | 文件名称: fft.h
 | 文件作用: 快速傅里叶变换
 | 创建日期: 2020-03-03
 | 更新日期: 2020-03-05
 | 开发人员: JuYan
 +----------------------------
 MIT License

 Copyright (C) JuYan

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
*/
#ifndef _INCLUDE_FFT_H_
#define _INCLUDE_FFT_H_
// define FFT2_USE_BUFF to use buffer in FFT2 class.
// #define FFT2_USE_BUFF
// define yourself memory mananger.
#define memalloc(typ, count)     (typ*)malloc((count) * sizeof(typ))
#define memfree(p)               free(p)
// -------------------------------
class FFT0                                          // Radix-2 lifting FFT
{
public:
    FFT0();
    ~FFT0();
    using data = double;                            // FFT data type
    using iodata = float;                           // Buffer data type
    struct complex                                  // 复数类型
    {
        data   real;
        data   imag;
        complex & operator+=(const  complex & b)
        {
            real += b.real;
            imag += b.imag;
            return *this;
        }
        complex & operator-=(const  complex & b)
        {
            real -= b.real;
            imag -= b.imag;
            return *this;
        }
        complex  operator*(const  complex & b)
        {
            complex t;
            t.real = real * b.real - imag * b.imag;
            t.imag = imag * b.real + real * b.imag;
            return t;
        }
    };
    void   InitFFT(int sz_index);
    void   Resize(int new_sz);
    void   BRO(complex *dat);
    void   FFT(complex *dat);
    void   IFFT(complex *dat);
    int    GetBitRevPos(int pos);
    static constexpr auto PI2 = 6.28318530717958647692;
private:
    int   n;
    int   wnkSz;
    int   nIndex;
    data *wnkTable;
    struct liftingDat                               // 提升运算参数
    {
        data     p;
        data     q;
    };
    void  liftingWnkCalc(liftingDat *fl, int k);
};
class FFT1 : public FFT0                            // chirp Z变换算法, 任意大小的一维FFT
{
public:
    FFT1();
    ~FFT1();
    void   InitFFT(int sz);
    void   Resize(int new_sz);
    void   FFT(complex *dat);
    void   IFFT(complex *dat);
    inline int GetFFTSize() const
    {
        return fftSz;
    }
    static constexpr auto PI = 3.1415926535897932384626;
private:
    int      fftSz, upnIndex, covSz;
    bool     runRadix2FFT;
    complex *wnkTable, *tmp1, *tmp2;
    int      ilog2UpwordRound(int n);
};
class FFT2 : public FFT1                            // 基于FFT1的二维FFT, 继承的类作为宽度方向FFT
{
public:
    FFT2();
    ~FFT2();
    enum class GetType
    {
        Real, Norm, NormWithScale, NormWithLog,
    };
    void   InitFFT(int w, int h);
    void   Resize(int new_w, int new_h);
#ifdef FFT2_USE_BUFF
    void   FFT(complex *dat = NULL);
    void   IFFT(complex *dat = NULL);
    void   Shift(complex *dat = NULL);
    void   BufferLoadReal(const iodata *rdat);
    void   BufferGetReal(iodata *rdat, GetType flag);
    inline complex * GetBufferPtr()
    {
        return dataBuffer;
    }
#else
    void   FFT(complex *dat);
    void   IFFT(complex *dat);
    void   Shift(complex *dat);
#endif
private:
    FFT1     hFFT;                                  // 高度方向的FFT
    int      matW, matH;
    complex *tmp, *dataBuffer;
};
#endif // _INCLUDE_FFT_H_
