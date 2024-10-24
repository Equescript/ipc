#pragma once

template <typename T>
class Array {
public:
    T* data;
    int length;
    Array(): data(0) {}
    Array(T* data, int length): data(data), length(length) {}
    T& operator[](int i) const {
        return data[i];
    };
};