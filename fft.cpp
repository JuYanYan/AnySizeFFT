/*
 | Mark
 | 文件名称: fft.cpp
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
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include "fft.h"
FFT0::FFT0()
{
    n = 0;
    nIndex = 0;
}
FFT0::~FFT0()
{
    if (n > 0)
    {
        memfree(wnkTable);
    }
    n = 0;
    nIndex = 0;
}
FFT1::FFT1()
{
    fftSz = 0;
    upnIndex = 0;
}
FFT1::~FFT1()
{
    if (fftSz > 0 && runRadix2FFT == false)
    {
        memfree(tmp1);
        memfree(tmp2);
        memfree(wnkTable);
    }
    fftSz = 0;
    upnIndex = 0;
}
FFT2::FFT2()
{
    matW = 0;
    matH = 0;
}
FFT2::~FFT2()
{
    if (matW > 0)
    {
        hFFT.~FFT1();                                               // 显示析构掉用于高度方向计算的FFT
        memfree(tmp);
#ifdef FFT2_USE_BUFF
        memfree(dataBuffer);
#endif
    }
    matW = 0;
    matH = 0;
}
// ------------------------
// class FFT0
// Radix-2 DIT Fast Fourier Transform
// ------------------------
// 初始化, sz_index表达2的n次方
void FFT0::InitFFT(int sz_index)
{
    int i;
    liftingDat d;
    assert(sz_index > 0 && sz_index < 31);
    n = 1 << sz_index;
    wnkSz = n;
    nIndex = sz_index;
    wnkTable = memalloc(data, n);                                   // 这个表只保存一半的旋转因子, 每个需要2个数字, 因此是:(N / 2) * 2 = N
    for (i = 0; i < n / 2; i++)
    {
        liftingWnkCalc(&d, i);
        wnkTable[i * 2 + 0] = d.p;
        wnkTable[i * 2 + 1] = d.q;
    }
}
// 重设FFT的大小
void FFT0::Resize(int new_szindex)
{
    n = 1 << new_szindex;
    nIndex = new_szindex;
    if (n > wnkSz)                                                  // 需要更改旋转因子表的大小
    {
        int i;
        liftingDat d;
        wnkSz = n;
        memfree(wnkTable);
        wnkTable = memalloc(data, n);
        for (i = 0; i < n / 2; i++)
        {
            liftingWnkCalc(&d, i);
            wnkTable[i * 2 + 0] = d.p;
            wnkTable[i * 2 + 1] = d.q;
        }
    }
}
// 取得某个数的反比特顺序
int FFT0::GetBitRevPos(int pos)
{
    int j, r = 0;
    for (j = 0; j < nIndex; j++)
    {
        int value = 1 << j;
        if ((pos & value) == value)                                 // 发现第j+1位是1(因为j是从0开始的...)
            r |= 1 << (nIndex - j - 1);                             // 对应输出的颠倒过来的那个位设置为1
    }
    return r;
}
// 计算Wnk的提升运算参数, k < n/2
void FFT0::liftingWnkCalc(liftingDat *fl, int k)
{
    data d = (data)((PI2 / wnkSz) * k);                             // k(2Pi/N)
    data c, s;
    c = (data)+cos(d);                                              // real = cos(k(2Pi/N))
    s = (data)-sin(d);                                              // imag = -sin(k(2Pi/N)) ??? ??? ??? ??? ??? ???
    assert(k < (wnkSz / 2));
    if (k == 0)                                                     // k可能等于0, 此时无效, 而另一种情况是k==n/2,被禁止
    {
        fl->p = 0;
        fl->q = 0;
        return;
    }
    fl->p = (c - 1) / s;
    fl->q = s;
}
// BRO
void FFT0::BRO(complex *dat)
{
    int i, j;
    for (i = 1, j = 0; i < n - 1; i++)
    {
        j = GetBitRevPos(i);
        if (i < j)
        {
            std::swap(dat[i], dat[j]);
        }
    }
}
// FFT
void FFT0::FFT(complex *dat)
{
    int        lay, curN, mergeN, j, k;
    int        midN, pos;
    complex    t;
    liftingDat lt;
    assert(n > 0);
    assert(nIndex > 0);
    for (lay = 1; lay <= nIndex; lay++)                             // 从最底层开始计算DFT, 一共存在着log(2,N)层
    {
        curN = 1 << lay;                                            // 当前合并的DFT大小
        midN = curN >> 1;                                           // 被合并的DFT大小, 它是合并结果大小curN的一半
        mergeN = 1 << (nIndex - lay);                               // 需要合并的次数
        for (k = 0; k < midN; k++)                                  // 处理这些小DFT中的数据
        {
            pos = (wnkSz / curN) * k;
            lt.p = wnkTable[0 + (pos << 1)];                        // 查表找到提升运算的数据
            lt.q = wnkTable[1 + (pos << 1)];
            for (j = k; j < n; j += curN)                           // 处理每一个DFT相同的数据点
            {
                t = dat[j + midN];
                t.real += lt.p * t.imag;                            // 完成提升运算, 以取代复数乘法
                t.imag += lt.q * t.real;
                t.real += lt.p * t.imag;
                dat[j + midN].real = dat[j].real - t.real;          // a = (-1) * t.real,  dat[..] = dat[j] + a
                dat[j + midN].imag = dat[j].imag - t.imag;          // 展开后+变成-

                dat[j] += t;
            }
        }
    }
}
// IFFT
void FFT0::IFFT(complex *dat)
{
    int        lay, curN, mergeN, j, k;
    int        midN, pos;
    data       tmpVal;
    complex    t;
    liftingDat lt;
    assert(n > 0);
    assert(nIndex > 0);
    for (lay = 1; lay <= nIndex; lay++)
    {
        curN = 1 << lay;
        midN = curN >> 1;
        mergeN = 1 << (nIndex - lay);
        for (k = 0; k < midN; k++)
        {
            pos = (wnkSz / curN) * k;
            lt.p = wnkTable[0 + (pos << 1)];
            lt.q = wnkTable[1 + (pos << 1)];
            for (j = k; j < n; j += curN)
            {
                t = dat[j + midN];
                t.real -= lt.p * t.imag;                            // 完成提升运算, 以取代复数乘法
                t.imag -= lt.q * t.real;
                t.real -= lt.p * t.imag;
                dat[j + midN].real = dat[j].real - t.real;
                dat[j + midN].imag = dat[j].imag - t.imag;

                dat[j] += t;
            }
        }
    }
    tmpVal = (data)(1.0 / n);
    for (k = 0; k < n; k++)                                         // * 1/N 来完成IFFT
    {
        (dat + k)->real *= tmpVal;
        (dat + k)->imag *= tmpVal;
    }
}
// ------------------------
// class FFT1
// Bluestein (chirp Z Transform) FFT
// ------------------------
// 初始化
void FFT1::InitFFT(int sz)
{
    fftSz = sz;
    upnIndex = ilog2UpwordRound(sz);
    if (sz == (1 << upnIndex))                                      // 只有非2的整数幂FFT才会被执行chirp Z
    {
        runRadix2FFT = true;
    }
    else {
        int  i, k;
        data d;
        upnIndex += 1;                                              // 从1-N...N+1, 刚好是对称的2N个点
        tmp1 = memalloc(complex, (1 << upnIndex));
        tmp2 = memalloc(complex, (1 << upnIndex));
        wnkTable = memalloc(complex, sz * 2);
        for (i = 0; i < sz * 2; i++)
        {
            k = (i * i) % (sz * 2);
            d = (data)((PI / sz) * k);
            (wnkTable + i)->real = +(data)cos(d);                   // o.real = +cos(k * (2Pi / N))
            (wnkTable + i)->imag = -(data)sin(d);                   // o.imag = -sin(k * (2Pi / N))
        }
        runRadix2FFT = false;

    }
    covSz = 1 << upnIndex;
    FFT0::InitFFT(upnIndex);
}
// 变更大小
void FFT1::Resize(int new_sz)
{
    if (new_sz != fftSz)                                            // 只有大小不同才重新初始化
    {
        this->~FFT1();
        InitFFT(new_sz);
    }
}
// FFT
void FFT1::FFT(complex *dat)
{
    int i;
    if (runRadix2FFT)
    {
        FFT0::BRO(dat);
        FFT0::FFT(dat);
        return;
    }                                                               // 下面我们将用Bluestein算法进行FFT
    memset(tmp1, 0, sizeof(complex) * covSz);
    memset(tmp2, 0, sizeof(complex) * covSz);
    tmp1[0].real = 1;
    tmp1[0].imag = 0;
    for (i = 1; i < fftSz; i++)                                     // w(n), n \in [1 - N, N - 1]
    {
        tmp1[i] = wnkTable[i];
        tmp1[covSz - i] = wnkTable[i];
    }
    FFT0::BRO(tmp1);
    FFT0::FFT(tmp1);

    for (i = 0; i < fftSz; i++)
    {
        tmp2[i].real = wnkTable[i].real * dat[i].real + wnkTable[i].imag * dat[i].imag;
        tmp2[i].imag = wnkTable[i].real * dat[i].imag - wnkTable[i].imag * dat[i].real;
    }
    FFT0::BRO(tmp2);
    FFT0::FFT(tmp2);

    for (i = 0; i < covSz; i++)
    {
        tmp2[i] = tmp2[i] * tmp1[i];
    }
    FFT0::BRO(tmp2);
    FFT0::IFFT(tmp2);

    for (i = 0; i < fftSz; i++)                                     // 最后, 我们拷贝数据
    {
        dat[i].real = wnkTable[i].real * tmp2[i].real + wnkTable[i].imag * tmp2[i].imag;
        dat[i].imag = wnkTable[i].real * tmp2[i].imag - wnkTable[i].imag * tmp2[i].real;
    }
}
// IFFT
void FFT1::IFFT(complex *dat)
{
    int  i;
    data v;
    complex t;
    if (runRadix2FFT)
    {
        FFT0::BRO(dat);
        FFT0::IFFT(dat);
        return;
    }
    memset(tmp1, 0, sizeof(complex) * covSz);                       // 注意往下的所有的Wnk全部都是虚部取负数的
    memset(tmp2, 0, sizeof(complex) * covSz);
    tmp1[0].real = 1;
    tmp1[0].imag = 0;
    for (i = 1; i < fftSz; i++)                                     // w(n), n \in [1 - N, N - 1]
    {
        t = wnkTable[i];
        t.imag = -t.imag;
        tmp1[i] = t;
        tmp1[covSz - i] = t;
    }
    FFT0::BRO(tmp1);
    FFT0::FFT(tmp1);

    for (i = 0; i < fftSz; i++)
    {
        tmp2[i].real = wnkTable[i].real * dat[i].real - wnkTable[i].imag * dat[i].imag;
        tmp2[i].imag = wnkTable[i].real * dat[i].imag + wnkTable[i].imag * dat[i].real;
    }
    FFT0::BRO(tmp2);
    FFT0::FFT(tmp2);

    for (i = 0; i < covSz; i++)
    {
        tmp2[i] = tmp2[i] * tmp1[i];
    }
    FFT0::BRO(tmp2);
    FFT0::IFFT(tmp2);

    v = (data)(1.0 / fftSz);
    for (i = 0; i < fftSz; i++)
    {
        dat[i].real = v * (wnkTable[i].real * tmp2[i].real - wnkTable[i].imag * tmp2[i].imag);
        dat[i].imag = v * (wnkTable[i].real * tmp2[i].imag + wnkTable[i].imag * tmp2[i].real);
    }
}
// 取得离它最近的log2的值, 比如3返回2(2^2 = 4), 5返回3(2^3 = 8)
int FFT1::ilog2UpwordRound(int n)
{
    int count = 0;
    int t = n;
    assert(n > 0);
    do
    {
        ++count;
        n >>= 1;
    } while (n > 1);
    if ((1 << count) != t)                                          // 为了保证n为2的整数幂的时候返回正确的值
    {                                                               // 不然向上舍入
        count++;
    }
    return count;
}
// ------------------------
// class FFT2
// FFT 2D
// ------------------------
// 初始化
void FFT2::InitFFT(int w, int h)
{
    matW = w;
    matH = h;
    hFFT.InitFFT(h);
    FFT1::InitFFT(w);
    tmp = memalloc(complex, h);
#ifdef FFT2_USE_BUFF
    dataBuffer = memalloc(complex, w * h);
#else
    dataBuffer = NULL;
#endif
}
// 改变大小
void FFT2::Resize(int new_w, int new_h)
{
    if (new_w != matW || new_h != matH)                             // 改变尺寸不一致时才重新初始化
    {
        this->~FFT2();
        InitFFT(new_w, new_h);
    }
}
// FFT2
void FFT2::FFT(complex *dat)
{
    int  i, k;
#ifdef FFT2_USE_BUFF
    if (dat == NULL)
    {
        dat = dataBuffer;
    }
#endif
    complex *pdat = dat;
    for (i = 0; i < matH; i++)                                      // 对行进行FFT
    {
        FFT1::FFT(pdat);                                            // 因为数组是行主序的, 因此直接在dat上进行
        pdat += matW;
    }
    for (i = 0; i < matW; i++)                                      // 对每一列进行FFT
    {
        pdat = dat + i;
        for (k = 0; k < matH; k++)
        {
            tmp[k] = *pdat;
            pdat += matW;                                           // 这样来得到每一列的元素
        }
        hFFT.FFT(tmp);                                              // 对当前列进行FFT
        pdat = dat + i;
        for (k = 0; k < matH; k++)                                  // 然后把数据拷贝回去
        {
            *pdat = tmp[k];
            pdat += matW;
        }
    }
}
// IFFT2
void FFT2::IFFT(complex *dat)
{
    int  i, k;
#ifdef FFT2_USE_BUFF
    if (dat == NULL)
    {
        dat = dataBuffer;
    }
#endif
    complex *pdat = dat;
    for (i = 0; i < matH; i++)                                      // 对行进行IFFT
    {
        FFT1::IFFT(pdat);
        pdat += matW;
    }
    for (i = 0; i < matW; i++)                                      // 对每一列进行IFFT
    {
        pdat = dat + i;
        for (k = 0; k < matH; k++)
        {
            tmp[k] = *pdat;
            pdat += matW;
        }
        hFFT.IFFT(tmp);
        pdat = dat + i;
        for (k = 0; k < matH; k++)                                  // 数据拷贝
        {
            *pdat = tmp[k];
            pdat += matW;
        }
    }
}
// 把低频集中到中央的移位
void FFT2::Shift(complex *dat)
{
    int i, j;
    int swapW = matW / 4;
    int swapH = matH / 4;
#ifdef FFT2_USE_BUFF
    if (dat == NULL)
    {
        dat = dataBuffer;
    }
#endif
    complex *pdat = dat;
    for (j = 0; j < matH; j++)
    {
        for (i = 0; i < swapW; i++)
        {
            std::swap(pdat[i], pdat[matW / 2 - i - 1]);
            std::swap(pdat[matW / 2 + i], pdat[matW - i - 1]);
        }
        pdat += matW;                                               // 跨行
    }
    for (j = 0; j < matW; j++)
    {
        pdat = dat + j;                                             // 跨列
        for (i = 0; i < swapH; i++)
        {
            std::swap(pdat[matW * i], pdat[matW * (matH / 2 - i - 1)]);
            std::swap(pdat[matW * (matH / 2 + i)], pdat[matW * (matH - i - 1)]);
        }
    }
}
#ifdef FFT2_USE_BUFF
// 向缓冲区加载实数
void FFT2::BufferLoadReal(const iodata *rdat)
{
    int sz = matW * matH;
    complex *pTmp = dataBuffer;
    while (sz-- > 0)
    {
        pTmp->real = *rdat;
        pTmp->imag = 0;
        pTmp++;
        rdat++;
    }
}
// 从缓冲区取得实数
void FFT2::BufferGetReal(iodata *rdat, GetType flag)
{
    int sz = matW * matH;
    data r;
    complex *pTmp = dataBuffer;
    while (sz-- > 0)
    {
        switch (flag)
        {
        case GetType::Real:
            r = pTmp->real;
            break;
        case GetType::Norm:
            r = sqrt(pTmp->real * pTmp->real + pTmp->imag * pTmp->imag);
            break;
        case GetType::NormWithScale:
            r = (sqrt(pTmp->real * pTmp->real + pTmp->imag * pTmp->imag) / (matW * matH * 0.25));
            break;
        case GetType::NormWithLog:
            r = (sqrt(pTmp->real * pTmp->real + pTmp->imag * pTmp->imag) / (matW * matH * 0.25));
            if (fabs(r) < FLT_EPSILON)
            {
                r = 0;
            }
            else {
                r = -20 * log10(r);
                if (r > 80)
                {
                    r = 80;
                }
                r = 1.0f - r / 80.0f;
            }
            break;
        }
        *rdat = (iodata)r;
        pTmp++;
        rdat++;
    }
}
#endif // FFT2_USE_BUFF

