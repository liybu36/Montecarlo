#ifndef SIZED_ARRAY_H
#define SIZED_ARRAY_H
#include <algorithm>

template <typename V> class sized_array {
  public:
    sized_array (int size): _size(size) { _ptr = new V[_size]; }
    sized_array (): _size(0),_ptr(0){}
    ~sized_array () { delete [] _ptr; }
    sized_array (const sized_array<V>& r): _size(r._size) { _ptr = new V[size]; std::copy (r._ptr, r._ptr + _size, _ptr); }
    void allocate (int size) {if(!_ptr) {_size = size ; _ptr = new V[_size];}} 
    void reallocate (int size) { V* new_ptr = new V[size]; if(_ptr) { std::copy (_ptr, _ptr + std::min(_size, size), new_ptr);} _size = size; _ptr = new_ptr; } 
    const sized_array<V>& operator= (const sized_array<V>& r) { delete [] _ptr; _size = r._size;  _ptr = new V[size]; std::copy (r._ptr, r._ptr + _size, _ptr); }
    operator V* () { return _ptr; } 
    operator const V* () const { return _ptr; } 
    V& operator[] (int i) { return _ptr[i]; }
    const V& operator[] (int i) const { return _ptr[i]; }
    int size () const { return _size; }
  private:
    V* _ptr;
    int _size;
};

#endif

/*
 * $Log: DSSizedArray.hh,v $
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
