/*
MIT License

Copyright (c) 2018 Jérémy Anger, Gabriele Facciolo, Mauricio Delbracio

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
#pragma once
#include <limits>
#include <fftw3.h>

template <class T>
class fftw_alloc {
    public:
        typedef T        value_type;
        typedef T*       pointer;
        typedef const T* const_pointer;
        typedef T&       reference;
        typedef const T& const_reference;
        typedef std::size_t    size_type;
        typedef std::ptrdiff_t difference_type;

        template <class U>
        struct rebind {
            typedef fftw_alloc<U> other;
        };

        pointer address (reference value) const {
            return &value;
        }
        const_pointer address (const_reference value) const {
            return &value;
        }

        fftw_alloc() throw() {
        }
        fftw_alloc(const fftw_alloc&) throw() {
        }
        template <class U>
        fftw_alloc (const fftw_alloc<U>&) throw() {
        }
        ~fftw_alloc() throw() {
        }

        size_type max_size () const throw() {
            return std::numeric_limits<std::size_t>::max() / sizeof(T);
        }

        pointer allocate (size_type num, const void* = 0) {
            void* ptr;
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
// Replace the fftwf line with the fftw line if type double is required
// Also requires linking against libfftw_3 as well as libfftw_3f
//            ptr = fftw_malloc(num*sizeof(T));
            ptr = fftwf_malloc(num*sizeof(T));
            return (pointer) ptr;
        }

        void construct (pointer p, const T& value) {
            new((void*)p)T(value);
        }

        void destroy (pointer p) {
            p->~T();
        }

        void deallocate (pointer p, size_type num) {
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
// Replace the fftwf line with the fftw line if type double is required
// Also requires linking against libfftw_3 as well as libfftw_3f
//            fftw_free(p);
            fftwf_free(p);
            num = num;
        }
};

template <class T1, class T2>
bool operator== (const fftw_alloc<T1>&,
                 const fftw_alloc<T2>&) throw() {
    return true;
}

template <class T1, class T2>
bool operator!= (const fftw_alloc<T1>&,
                 const fftw_alloc<T2>&) throw() {
    return false;
}
