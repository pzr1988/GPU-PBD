#ifndef FIXED_VECTOR_H
#define FIXED_VECTOR_H

#include "Pragma.h"

template<typename T, size_t Capacity>
class FixedVector {
 private:
  T _data[Capacity];
  size_t _size = 0;

 public:
  DEVICE_HOST void push_back(const T& value) {
    assert(_size < Capacity);
    _data[_size++] = value;
  }
  DEVICE_HOST const T& operator[](size_t index) const {
    assert(index < _size);
    return _data[index];
  }
  DEVICE_HOST T& operator[](size_t index) {
    assert(index < _size);
    return _data[index];
  }
  DEVICE_HOST size_t size() const {
    return _size;
  }
  DEVICE_HOST bool empty() const {
    return _size == 0;
  }
  DEVICE_HOST size_t capacity() const {
    return Capacity;
  }
};

#endif