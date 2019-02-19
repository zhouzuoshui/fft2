#ifndef FFT_HPP
#define FFT_HPP

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <iostream>

#define EPSILON (1e-5)

namespace myFFT{
    template <typename V>
    struct Complex{
        V real;
        V imag;

        Complex(): real(0), imag(0){}
        Complex(V r): real(r), imag(0){}
        Complex(V r, V i):
            real(r), imag(i){}

        Complex conj(){
            Complex T;
            T.real = real; T.imag = -imag;
            return T;
        }

        friend Complex operator+(Complex oprd1, Complex oprd2){
            Complex T;
            T.real = oprd1.real + oprd2.real;
            T.imag = oprd1.imag + oprd2.imag;
            return T;
        }
        friend Complex operator-(Complex oprd1, Complex oprd2){
            Complex T;
            T.real = oprd1.real - oprd2.real;
            T.imag = oprd1.imag - oprd2.imag;
            return T;
        }
        friend Complex operator*(Complex oprd1, Complex oprd2){
            Complex T;
            T.real = oprd1.real * oprd2.real - oprd1.imag * oprd2.imag;
            T.imag = oprd1.real * oprd2.imag + oprd1.imag * oprd2.real;
            return T;
        }
        friend Complex operator/(Complex oprd1, Complex oprd2){
            Complex T;
            T = oprd1 * oprd2.conj();
            V mod = oprd2.real * oprd2.real + oprd2.imag * oprd2.imag;
            assert(mod > EPSILON);
            T.real = T.real / mod;
            T.imag = T.imag / mod;
            return T;
        }
        
    };


    // 参数是当前下标以及最长下标加1，也就是原序列的长度
    int rev(int x, int num){
        int n = std::log2(num); // max. length of binary
        std::vector<int> ans;
        while(x != 0){
            ans.push_back(x%2);
            x /= 2;
        }
        while(ans.size() < n){
            ans.push_back(0); // 0 padding
        }
        int reval = 0;
        // convert binary to decimal 
        int power = std::pow(2, n-1);
        
        for(int i: ans){
            reval += i * power;
            power /= 2;
        }

        return reval;
    }

    template <typename V>
    void bitReverseCopy(Complex<V>* in, int size){
        Complex<V>* data = new Complex<V> [size];
        for(int i=0; i < size; ++i){
            auto revc = rev(i, size);
            data[revc] = in[i];
        }
        for(int i=0; i < size; ++i){
            in[i] = data[i];
        }
        delete [] data;
    }

    template <typename V>
    void _CooleyTukey(Complex<V>* data, int size, bool inverse=false){
        bitReverseCopy(data, size);
        
        // for each layer
        for(int s=1; s < std::log2(size)+1; ++s){
            int m = std::pow(2, s);
            float theta = (-2*M_PI)/m;
            Complex<V> wm{std::cos(theta), std::sin(theta)};

            for(int k=0; k < size; k += m){
                Complex<V> w = 1;
                for(int j=0; j < m/2; ++j){
                    // butterfly
                    Complex<V> u    = data[k+j];
                    Complex<V> t    = w * data[k+j+m/2];
                    data[k+j]       = u + t;
                    data[k+j+m/2]   = u - t;

                    w = inverse? w/wm: w*wm;
                }
            }
        }
        // for IFFT
        if(inverse)
            for(int i=0; i < size; ++i){
                data[i].real /= size;
                data[i].imag /= size;
            }
    }

    template <typename V>
    void _Bluestein(Complex<V>* data, int size, bool inverse=false){
        // to avoid data winding, length should be twice
        // to use Cooley-Tukey, padding 0 to power of 2
        int Npow2 = std::exp2(std::ceil(std::log2(size*2)));
        Complex<V>* w = new Complex<V> [size*2];
        Complex<V>* y = new Complex<V> [Npow2];
        Complex<V>* b = new Complex<V> [Npow2];
        // positive real point of unit circle, the start of rotation
        w[0] = {1,0};
        y[0] = {1,0};
        // rotation angle, in complex format
        for(int i=1; i < size; ++i){
            double theta = M_PI / size *i*i;
            w[i].real = std::cos(theta);
            w[i].imag = inverse? std::sin(theta): -std::sin(theta);
            y[i] = w[i].conj();
            y[Npow2 - i] = y[i];
        }

        // preparation in frequence domain
        _CooleyTukey(y, Npow2);

        // multiply in time domain
        for(int i=0; i < size; ++i)
            b[i]= data[i] * w[i];

        // tranform into frequency domain
        _CooleyTukey(b, Npow2);

        // Convolution in time domain is product in F domain
        for(int i=0; i < Npow2; ++i)
            b[i] = b[i] * y[i];

        // back to time domain
        _CooleyTukey(b, Npow2, true);

        // multiply
        for(int i=0; i < size; ++i){
            data[i] = w[i] * b[i];
        }

        delete [] w;
        delete [] y;
        delete [] b;

        if(inverse)
            for(int i=0; i < size; ++i){
                data[i].real /= size;
                data[i].imag /= size;
            }
    }

    // for internal call
    template <typename V>
    void fft(Complex<V>* in, int size){
        if(std::ceil(std::log2(size)) - std::floor(std::log2(size)) > EPSILON) // depth is not an integer
            return _Bluestein(in, size);
        else 
            return _CooleyTukey(in, size);
    }

    // for external call, same for the rest
    template <typename V>
    Complex<V>* fft(V* in, int size){
        Complex<V> *newData = new Complex<V> [size];
        std::memset(newData, 0, size * sizeof newData[0]);
        for(int i=0; i < size; ++i)
            newData[i].real = in[i];
        fft(newData, size);

        return newData;    
    }

    template <typename V>
    void ifft(Complex<V>* data, int size){
        if(std::ceil(std::log2(size)) - std::floor(std::log2(size)) > EPSILON) // depth is not an integer
            return _Bluestein(data, size, true);
        else 
            return _CooleyTukey(data, size, true);
    }

    template <typename V>
    Complex<V>* ifft(V* in, int size){
        Complex<V> *newData = new Complex<V> [size];
        std::memset(newData, 0, size * sizeof newData[0]);
        for(int i=0; i < size; ++i)
            newData[i].real = in[i];
        ifft(newData, size, true);

        return newData;
    }

    template <typename V>
    void tranpose(V* data, int width, int height){
        V* dataNew = new V[width * height];

        for(int j=0; j < height; ++j)
            for(int i=0; i < width; ++i)
                dataNew[i* height + j] = data[j * width + i];

        for(int j=0; j < height; ++j)
            for(int i=0; i < width; ++i)
                data[j * width + i] = dataNew[j * width + i];

        delete []dataNew;
    }

    // FFT2: seperable, first 1d FFT in each row, and then each columns.
    template <typename V>
    void _fft2(Complex<V>* data, int width, int height, bool inverse=false){
        for(int j=0; j < height; ++j){
            Complex<V>* tar = &data[j*width];
            inverse? ifft(tar, width): fft(tar, width);
            for(int i=0; i < width; ++i){
                data[j*width+i] = tar[i];
            }
        }
       
        tranpose(data, width, height);
        
        for(int j=0; j < width; ++j){
            Complex<V>* tar = &data[j*height];
            inverse? ifft(tar, height): fft(tar, height);
            for(int i=0; i < height; ++i){
                data[j*height+i] = tar[i];
            }
        }
        
        tranpose(data, height, width);
    };

    // for potential internal call
    template <typename V>
    void fft2(Complex<V>* data, int width, int height){
        return _fft2(data, width, height, false);
    }

    template <typename V>
    Complex<V>* fft2(V* data, int width, int height){
        int size = width * height;
        Complex<V> *newData = new Complex<V> [size];
        std::memset(newData, 0, size * sizeof newData[0]);
        for(int j=0; j < height; ++j){
            for(int i=0; i < width; ++i){
                newData[j * width + i].real = data[j * width + i];
            }
        }
        _fft2(newData, width, height, false);

        return newData;
    }

    template <typename V>
    void ifft2(Complex<V>* data, int width, int height){
        return _fft2(data, width, height, true);
    }

    template <typename V>
    Complex<V>* ifft2(V* data, int width, int height){
        int size = width * height;
        Complex<V> *newData = new Complex<V> [size];
        std::memset(newData, 0, size * sizeof newData[0]);
        for(int j=0; j < height; ++j){
            for(int i=0; i < width; ++i){
                newData[j * width + i].real = data[j * width + i];
            }
        }
        _fft2(newData, width, height, true);

        return newData;
    }
}

#undef EPSILON
#endif
