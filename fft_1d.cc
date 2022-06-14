#include "complex.h"
#include <vector>
#include <complex>
#include <stddef.h>
#include <utility> // std::swap in c++11
#include <assert.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <string.h>
#include <functional>

#define USING_DOUBLE_PRECISION 1

#if USING_DOUBLE_PRECISION
	typedef double d_type;
#else
	typedef float d_type;
#endif

inline complex_t<d_type> TwiddleFactor(size_t N, size_t k, int sign)
{
    d_type r = (d_type)1;
    d_type theta = sign*C_2PI*k/N;
    return std::polar2<d_type>(r, theta);
}

template<typename T>
inline void ButterFly(T &a, T &b, T w)
{
    T tmp = a;
    a = a + b*w;
    b = tmp - b*w;
}

// below function produce  https://oeis.org/A030109
void bit_reverse_permute(size_t radix2_num, std::vector<size_t> &arr)
{
    size_t k;
    arr.resize(std::pow(2,radix2_num));
    arr[0] = 0;
    for(k=0;k<radix2_num;k++){
       size_t last_k_len = std::pow(2, k);
       size_t last_k;
       for(last_k = 0; last_k < last_k_len; last_k++){
           arr[last_k] = 2*arr[last_k];
           arr[last_k_len+last_k] = arr[last_k]+1;
       }
    }
}

template<typename T>
void bit_reserse_radix2(T *seq, size_t length)
{
    std::vector<size_t> index_r;
    int nbits = std::log2(length);

    // generate the sequence of bit-reversed index
    // https://en.wikipedia.org/wiki/Bit-reversal_permutation
    for(size_t i = 0; i < length; i++)
    {
        size_t ir = 0;
        for(int b = 0; b < nbits; b++)
        {
            if((i >> b) & 1)
            {
                ir |= 1 << (nbits - 1 - b);
            }
        }
        index_r.push_back(ir);
    }

    for(size_t i = 0; i < length; i++)
    {
        size_t ir = index_r[i];
        if(i < ir)
            std::swap(seq[i], seq[ir]);
    }
}

template<typename T>
void rand_vec(std::vector<complex_t<T>> & seq){
    static std::random_device rd;   // seed

    static std::mt19937 mt(rd());
    static std::uniform_real_distribution<T> dist(-2.0, 2.0);

    size_t i;
    for(i = 0; i < seq.size(); i++)
    {
        seq[i] = complex_t<T>(dist(mt), dist(mt));
    }
}

template<typename T>
void rand_vec(std::vector<T> & seq){
    static std::random_device rd;   // seed

    static std::mt19937 mt(rd());
    static std::uniform_real_distribution<T> dist(-2.0, 2.0);

    size_t i;
    for(i = 0; i < seq.size(); i++)
    {
        seq[i] = dist(mt);
    }
}

#ifndef ABS
#define ABS(x) ((x)>0?(x):-1*(x))
#endif

template<typename T>
int valid_vector(const std::vector<complex_t<T>> & lhs, const std::vector<complex_t<T>> & rhs, T delta = (T)0.001){
    assert(lhs.size() == rhs.size());
    size_t i;
    int err_cnt = 0;
    for(i = 0;i < lhs.size(); i++){
        T d_re = std::real(lhs[i]) - std::real(rhs[i]);
        T d_im = std::imag(lhs[i]) - std::imag(rhs[i]);
        d_re = ABS(d_re);
        d_im = ABS(d_im);
        if(d_re > delta || d_im > delta){
            std::cout<<" diff at "<<i<<", lhs:"<<lhs[i]<<", rhs:"<<rhs[i]<<std::endl;
            err_cnt++;
        }
    }
    return err_cnt;
}

template<typename T>
void DFT_naive(std::vector<complex_t<T>> &t_seq, std::vector<complex_t<T>> &f_seq, size_t length, bool is_inverse)
{
    // https://en.wikipedia.org/wiki/Discrete_Fourier_transform

    int sign = is_inverse? 1:-1;
    for(size_t k = 0; k < length; k++)
    {
        complex_t<T> omega_k = TwiddleFactor(length, k, sign);
        complex_t<T> X_k;

        for(size_t n = 0; n < length; n++)
        {
            X_k += t_seq[n]*std::pow(omega_k, (T)n);
        }
        
        if(is_inverse) X_k /= (T)length;
        f_seq.push_back(X_k);
    }
}

/*
template<typename T>
void iDFT_naive(std::vector<complex_t<T>> &f_seq, std::vector<complex_t<T>> &t_seq, size_t length)
{
    for(size_t n = 0; n < length; n++)
    {
        complex_t<T> omega_n = TwiddleFactor(length, n, 1);
        complex_t<T> x_n;

        for(size_t k = 0; k < length; k++)
        {
            x_n += f_seq[k]*std::pow(omega_n, (T)k);
        }
        x_n /= (T)length;
        t_seq.push_back(x_n);
    }
}
*/

template<typename T>
void fft_cooley_tukey(complex_t<T> *seq, size_t length, bool is_inverse)
{
    if(length == 1) return; // f_seq = t_seq

    // arrange seq to a bit-reversed order
    bit_reserse_radix2(seq, length);

    // pre-computed half list of omege_k
    std::vector<complex_t<T>> omega_k;
    int sign = is_inverse? 1:-1;
    for(size_t k = 0; k < length/2; k++)
    {
        omega_k.push_back(TwiddleFactor(length, k, sign));
    }

    // pt is the number of points in each fft group
    for(size_t pt = 2; pt <= length; pt *= 2)
    {
        size_t N_group = length/pt;
        size_t stride = pt/2;
        for(size_t g = 0; g < N_group; g++)
        {
            for(size_t s = 0; s < stride; s++)
            {
                ButterFly(seq[g*N_group + s], seq[g*N_group + s + stride], omega_k[N_group*s]);
            }
        }
    }

    if(is_inverse)
    {
        for(size_t k = 0; k < length; k++)
        {
            seq[k] /= (T)length;
        }
    }
}

template<typename T>
void FFT_R2C(T *t_seq, complex_t<T> *f_seq, size_t length)
{
    if(length == 1) return;
    assert( ( (length & (length - 1)) == 0 ) && "the length mush be a power of 2");

    // TI "Implementing FFT Algorithms of Real-Valued Sequence" 
    std::vector<complex_t<T>> seq;
    for(size_t n = 0; n < length/2; n++)
    {
        seq.push_back(complex_t<T>(t_seq[2*n], t_seq[2*n + 1]));
    }
    fft_cooley_tukey(seq.data(), length/2, 0);

    seq.push_back(seq[0]); // seq[length/2] = seq[0]
    complex_t<T> A, B;
    for(size_t k = 0; k < length/2; k++)
    {
        complex_t<T> tmp1((T)0.5, (T)0);
        complex_t<T> tmp2((T)0, (T)0.5);
        complex_t<T> omega_k = TwiddleFactor(length, k, -1);
        A = tmp1 - tmp2*omega_k;
        B = tmp1 + tmp2*omega_k;

        f_seq[k] = A*seq[k] + std::conj(seq[length/2 - k])*B;
    }
    f_seq[length/2] = complex_t<T>( std::real(seq[0])-std::imag(seq[0]), (T)0); // why?
    
    // conjugate symmetric
    for(size_t k = 1; k < length/2; k++)
    {
        f_seq[length - k] = std::conj(f_seq[k]);
    }
}

template<typename T>
void iFFT_C2R(complex_t<T> *f_seq, T *t_seq, size_t length)
{
    if(length == 1) return;
    assert( ( (length & (length - 1)) == 0 ) && "the length mush be a power of 2");

    // TI "Implementing FFT Algorithms of Real-Valued Sequence" 
    std::vector<complex_t<T>> seq;
    complex_t<T> A, B;
    for(size_t k = 0; k < length/2; k++)
    {
        complex_t<T> tmp1((T)0.5, (T)0);
        complex_t<T> tmp2((T)0, (T)0.5);
        complex_t<T> omega_k = TwiddleFactor(length, k, -1);
        A = tmp1 - tmp2*omega_k;
        B = tmp1 + tmp2*omega_k;

        seq.push_back(f_seq[k]*std::conj(A) + std::conj(f_seq[length/2 - k])*std::conj(B));
    }
    fft_cooley_tukey(seq.data(), length/2, 1);

    for(size_t n = 0; n < length/2; n++)
    {
        t_seq[2*n] = std::real(seq[n]);
        t_seq[2*n+1] = std::imag(seq[n]);
    }
}

int main()
{
    const size_t max_size = 8;
    for(size_t size = 2; size<=max_size; size *= 2)
    {
        std::vector<complex_t<d_type>> t_seq;
        std::vector<complex_t<d_type>> f_seq;
        std::vector<complex_t<d_type>> t_seq_r;

        std::vector<d_type> seq_fwd_real;
        std::vector<d_type> seq_bwd_real;
        seq_fwd_real.resize(size);
        seq_bwd_real.resize(size);

        rand_vec(seq_fwd_real);

        for(size_t i = 0; i < size; i++)
            t_seq.push_back(complex_t<d_type>(seq_fwd_real[i], (d_type)0));

        DFT_naive(t_seq, f_seq, size, 0);
        DFT_naive(f_seq, t_seq_r, size, 1); //inverse

        std::vector<complex_t<d_type>> seq_fwd;
        std::vector<complex_t<d_type>> seq_bwd;
        seq_fwd.resize(size);

        FFT_R2C(seq_fwd_real.data(), seq_fwd.data(), size);
        iFFT_C2R(f_seq.data(), seq_bwd_real.data(), size);
        for(size_t i = 0; i < size; i++)
        {
            seq_bwd.push_back(complex_t<d_type>(seq_bwd_real[i], (d_type)0));
        }
        
        int err_cnt = valid_vector(f_seq, seq_fwd);
        int ierr_cnt = valid_vector(t_seq_r, seq_bwd);
        std::cout<<"length:"<<size<<", r2c fwd valid:"<< ( (err_cnt==0)?"y":"n" ) <<
            ", c2r bwd valid:"<<( (ierr_cnt==0)?"y":"n" ) <<std::endl;
    }

    return 0;
}