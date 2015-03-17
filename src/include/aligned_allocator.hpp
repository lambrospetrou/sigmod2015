
#include "LPUtils.hpp"
#include <cstdint>
#include <vector>
#include <iostream>
#include <cstddef>
#include <stdexcept>
#include <xmmintrin.h>

/**
 * Allocator for aligned data.
 *
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */
template <typename T, std::size_t Alignment>
class aligned_allocator
{
	public:
 
		// The following will be the same for virtually all allocators.
		typedef T * pointer;
		typedef const T * const_pointer;
		typedef T& reference;
		typedef const T& const_reference;
		typedef T value_type;
		typedef std::size_t size_type;
		typedef ptrdiff_t difference_type;
 
		T * address(T& r) const { return &r; }
 
		const T * address(const T& s) const { return &s; }
 
		std::size_t max_size() const {
			// The following has been carefully written to be independent of
			// the definition of size_t and to avoid signed/unsigned warnings.
			return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
		}
 
 
		// The following must be the same for all allocators.
		template <typename U>
		struct rebind {
			typedef aligned_allocator<U, Alignment> other;
		};
 
		bool operator!=(const aligned_allocator& other) const {
			return !(*this == other);
		}
 
		ALWAYS_INLINE void construct(T * const p, const T& t) const {
			void * const pv = static_cast<void *>(p);
			new (pv) T(t);
		}
 
		void destroy(T * const p) const	{
			p->~T();
		}
 
		// Returns true if and only if storage allocated from *this
		// can be deallocated from other, and vice versa.
		// Always returns true for stateless allocators.
		bool operator==(const aligned_allocator& other) const { (void)other;
			return true;
		}
 
 
		// Default constructor, copy constructor, rebinding constructor, and destructor.
		// Empty for stateless allocators.
		aligned_allocator() { }
 
		aligned_allocator(const aligned_allocator&) { }
 
		template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }
 
		~aligned_allocator() { }
 
 
		// The following will be different for each allocator.
		T * allocate(const std::size_t n) const	{
			// The return value of allocate(0) is unspecified.
			// Mallocator returns NULL in order to avoid depending
			// on malloc(0)'s implementation-defined behavior
			// (the implementation can define malloc(0) to return NULL,
			// in which case the bad_alloc check below would fire).
			// All allocators can return NULL in this case.
			if (!n) { return NULL; }
 
			// All allocators should contain an integer overflow check.
			// The Standardization Committee recommends that std::length_error
			// be thrown in the case of integer overflow.
			if (unlikely(n > max_size())) {
				throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
			}
 
			// Mallocator wraps malloc().
			void * const pv = _mm_malloc(n * sizeof(T), Alignment);
 
			// Allocators should throw std::bad_alloc in the case of memory allocation failure.
			if (unlikely(pv == nullptr)) {
				throw std::bad_alloc();
			}
 
			return static_cast<T *>(pv);
		}
 
		void deallocate(T * const p, const std::size_t n) const { (void)n;
			_mm_free(p);
		}
 
 
		// The following will be the same for all allocators that ignore hints.
		template <typename U>
		T * allocate(const std::size_t n, const U * /* const hint */) const	{
			return allocate(n);
		}
 
 
		// Allocators are not required to be assignable, so
		// all allocators should have a private unimplemented
		// assignment operator. Note that this will trigger the
		// off-by-default (enabled under /Wall) warning C4626
		// "assignment operator could not be generated because a
		// base class assignment operator is inaccessible" within
		// the STL headers, but that warning is useless.
	private:
		aligned_allocator& operator=(const aligned_allocator&);
};
/*
int main()
{
	typedef std::vector<__m128, aligned_allocator<__m128, sizeof(__m128)> > aligned_vector;
	aligned_vector lhs;
	aligned_vector rhs;

	float a = 1.0f;
	float b = 2.0f;
	float c = 3.0f;
	float d = 4.0f;

	float e = 5.0f;
	float f = 6.0f;
	float g = 7.0f;
	float h = 8.0f;

	for (std::size_t i = 0; i < 1000; ++i)
	{
		lhs.push_back(_mm_set_ps(a, b, c, d));
		rhs.push_back(_mm_set_ps(e, f, g, h));

		a += 1.0f; b += 1.0f; c += 1.0f; d += 1.0f;
		e += 1.0f; f += 1.0f; g += 1.0f; h += 1.0f;
	}

	__m128 mul = _mm_mul_ps(lhs[10], rhs[10]);
}
*/
