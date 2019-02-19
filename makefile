CC = clang++
CXXFLAGS = -g -stdlib=libc++ -std=c++14

fft: fft.hpp fft.cpp
	$(CC) $(CXXFLAGS) -o fft fft.cpp
