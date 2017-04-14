#ifndef RINGBUFFER_HPP
#define RINGBUFFER_HPP

#include <stdint.h>
#include <vector>
#include <stdexcept>

template <class T>
class RingBuffer {
    public:
        RingBuffer(int size, T defaultValue);
        void add(T value);
        unsigned int getPosition() {
            return bufferPosition;
        }
        void setPosition(unsigned int position) {
            if (position >= (unsigned int)buffer.size()) return;
            bufferPosition = position;
        }
        unsigned int size() {
            return (unsigned int)buffer.size();
        }
        std::vector<T> &getBuffer() {
            return buffer;
        }
        // only implemented for float & double
        T getMedian() {
            throw std::domain_error("RingBuffer::getMedian called on unsupported type");
            return defaultValue;
        }
    private:
        std::vector<T> buffer;
        unsigned int bufferPosition;
        T defaultValue;
};

template <class T>
RingBuffer<T>::RingBuffer(int size, T defaultValue) {
    this->defaultValue = defaultValue;
    for (int i=0; i<size; i++) {
        buffer.push_back(defaultValue);
    }
    bufferPosition = 0;
}

template <class T>
void RingBuffer<T>::add(T value) {
    buffer[bufferPosition] = value;
    bufferPosition++;
    if (bufferPosition >= buffer.size()) bufferPosition = 0;
}

template <>
float RingBuffer<float>::getMedian() {
    float median = 0.0;
    for (int i=0; i<(int)buffer.size(); i++) {
        median += buffer[i];
    }
    return median / (float)buffer.size();
}

template <>
double RingBuffer<double>::getMedian() {
    double median = 0.0;
    for (int i=0; i<(int)buffer.size(); i++) {
        median += buffer[i];
    }
    return median / (double)buffer.size();
}

template <>
long double RingBuffer<long double>::getMedian() {
    long double median = 0.0;
    for (int i=0; i<(int)buffer.size(); i++) {
        median += buffer[i];
    }
    return median / (long double)buffer.size();
}

#endif
