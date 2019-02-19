#include "fft.hpp"
#include <iomanip>
using namespace std;

int main(){
    int width = 5, height = 5;
    int size = width * height;
    float *data = new float[size];

    std::memset(data, 0, size * sizeof data[0]);

    for(int i=0; i < size; ++i)
        data[i] = i;

    for(int i=0; i < height; ++i){
        for(int j=0; j < width; ++j){
            data[i*width+j] = i;
        }
    }

    cout.setf(ios::fixed | ios::showpos);

    cout << "original elements of vector are:" << endl;
    for(int i=0; i < height; ++i){
        for(int j=0; j < width; ++j){
            cout << setprecision(2) << data[i*width+j] << ", ";
        }
        cout << endl;
    }

    auto nw = myFFT::fft2(data, width, height);

    cout << "after FFT2:" << endl;
    for(int i=0; i < height; ++i){
        for(int j=0; j < width; ++j){
            cout << setprecision(2) << nw[i*width+j].real << nw[i*width+j].imag << "j, ";
        }
        cout << endl;
    }

    myFFT::ifft2(nw, width, height);

    cout << "after IFFT2:" << endl;
    for(int i=0; i < height; ++i){
        for(int j=0; j < width; ++j){
            cout << setprecision(2) << nw[i*width+j].real << nw[i*width+j].imag << "j, ";
        }
        cout << endl;
    }

    return 0;
}
