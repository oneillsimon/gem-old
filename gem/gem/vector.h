#ifndef __GEM_VECTOR_H__
#define __GEM_VECTOR_H__

#include "matrix.h"

namespace gem
{
	template<typename T, unsigned int N>
	class vector : public gem::matrix<T, 1, N>
	{
	public:
		vector() : gem::matrix<T, 1, N>() {}
		vector(const std::initializer_list<T>& l) : gem::matrix<T, 1, N>(l) {}
		vector(const std::vector<T>& v) : gem::matrix<T, 1, N>(v) {}

		template<typename U>
		vector(const gem::matrix<U, 1, N>& m)
		{
			for(unsigned i = 0; i < N; ++i)
			{
				(*this)(i) = m(0, i);
			}
		}
		
		inline T x() const { return (*this)(0); }
		inline T y() const { return (*this)(1); }
		inline T z() const { return (*this)(2); }
		inline T w() const { return (*this)(3); }

		inline unsigned int size() const { return cols() };

		inline void set_x(T t) { (*this)(0) = t; }
		inline void set_y(T t) { (*this)(0) = t; }
		inline void set_z(T t) { (*this)(0) = t; }
		inline void set_w(T t) { (*this)(0) = t; }

		template<typename U>
		vector<T, N> operator *(const vector<U, N>& v) const
		{
			auto r = vector<T, N>();

			for(unsigned i = 0; i < N; ++i)
			{
				r(i) = (*this)(i) * v(i);
			}

			return r;
		}

		template<typename U>
		vector<T, N> operator /(const vector<U, N>& v) const
		{
			auto r = vector<T, N>();

			for(unsigned i = 0; i < N; ++i)
			{
				r(i) = (*this)(i) / v(i);
			}

			return r;
		}

		template<typename U>
		void operator *=(const vector<U, N>& v)
		{
			*this = (*this * v);
		}
		
		template<typename U>
		void operator /=(const vector<U, N>& v)
		{
			*this = *this / v;
		}

		T operator ()(unsigned int i) const { return m_values[0][i]; }
		T& operator ()(unsigned int i) { return m_values[0][i]; }
	};

	template<typename T>
	using vector2 = vector<T, 2>;

	typedef vector2<float> vector2f;
	typedef vector2<int> vector2i;

	template<typename T>
	using vector3 = vector<T, 3>;

	typedef vector3<float> vector3f;
	typedef vector3<int> vector3i;

	template<typename T>
	using vector4 = vector<T, 4>;

	typedef vector4<float> vector4f;
	typedef vector4<int> vector4i;
};
#endif