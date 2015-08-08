#ifndef __GEM_MATRIX_H__
#define __GEM_MATRIX_H__

#include <initializer_list>
#include <iostream>
#include <sstream>
#include <vector>

#include "gem.h"

namespace gem
{
	template<class T, unsigned int M, unsigned int N>
	class matrix
	{
	public:
		matrix() { gem::fill<T, M, N>(*this, 0); }

		matrix(const std::initializer_list<T>& l) : matrix(std::vector<T>(l)) {}

		matrix(const std::vector<T>& v)
		{
			unsigned index = 0;
			for(unsigned i = 0; i < M; ++i)
			{
				for(unsigned j = 0; j < N; ++j)
				{
					m_values[i][j] = v[index++];
				}
			}
		}

		template<typename U>
		matrix(const matrix<U, M, N>& m)
		{
			*this = m;
		}

		inline unsigned int cols() const { return N; }
		inline unsigned int rows() const { return M; }
		inline unsigned int size() const { return cols() * rows(); }

		template<typename U>
		matrix<T, M, N>& operator =(const matrix<U, M, N>& m)
		{
			for(unsigned i = 0; i < M; ++i)
			{
				for(unsigned j = 0; j < N; ++j)
				{
					(*this)(i, j) = (T)m(i, j);
				}
			}

			return *this;
		}

		template<class U, unsigned int O, unsigned int P>
		bool operator ==(const matrix<U, O, P>& m) const
		{
			if(!gem::same_dimensions(*this, m))
			{
				return false;
			}

			for(unsigned i = 0; i < rows(); ++i)
			{
				for(unsigned j = 0; j < cols(); ++j)
				{
					if((*this)(i, j) != m(i, j))
					{
						return false;
					}
				}
			}

			return true;
		}

		template<class U, unsigned int O, unsigned int P>
		bool operator !=(const matrix<U, O, P>& m) const
		{
			return !(*this == m);
		}

		matrix<T, M, N> operator -() const
		{
			auto r = matrix<T, M, N>();

			for(unsigned i = 0; i < rows(); ++i)
			{
				for(unsigned j = 0; j < cols(); ++j)
				{
					r(i, j) = -(*this)(i, j);
				}
			}

			return r;
		}

		template<class U, unsigned int O, unsigned int P>
		matrix<T, M, N> operator +(const matrix<U, O, P>& m) const
		{
			auto r = matrix<T, M, N>();

			for(unsigned i = 0; i < rows(); ++i)
			{
				for(unsigned j = 0; j < cols(); ++j)
				{
					r(i, j) = (*this)(i, j) + m(i, j);
				}
			}

			return r;
		}

		template<class U, unsigned int O, unsigned int P>
		matrix<T, M, N> operator -(const matrix<U, O, P>& m) const
		{
			return *this + -m;
		}

		template<class U, unsigned int O, unsigned int P>
		matrix<T, M, P> operator *(const matrix<U, O, P>& m) const
		{
			auto r = matrix<T, M, P>();

			for(unsigned i = 0; i < rows(); ++i)
			{
				for(unsigned j = 0; j < m.cols(); ++j)
				{
					for(unsigned k = 0; k < m.rows(); ++k)
					{
						r(i, j) += (*this)(i, k) * m(k, j);
					}
				}
			}

			return r;
		}

		template<class U, unsigned int O, unsigned int P>
		matrix<T, M, N> operator /(const matrix<U, O, P>& m) const
		{
			return *this * gem::inverse(m);
		}

		template<class U, unsigned int O, unsigned int P>
		void operator +=(const matrix<U, O, P>& m)
		{
			*this = *this + m;
		}

		template<class U, unsigned int O, unsigned int P>
		void operator -=(const matrix<U, O, P>& m)
		{
			*this = *this - m;
		}

		template<class U, unsigned int O, unsigned int P>
		void operator *=(const matrix<U, O, P>& m)
		{
			*this = *this * m;
		}

		template<class U, unsigned int O, unsigned int P>
		void operator /=(const matrix<U, O, P>& m)
		{
			*this = *this / m;
		}

		matrix<T, M, N> operator +(T t) const
		{
			return *this + matrix<T, M, N>(gem::filled_vector<T, M * N>(t));
		}

		matrix<T, M, N> operator -(T t) const
		{
			return *this + -t;
		}

		matrix<T, M, N> operator *(T t) const
		{
			auto r = *this;

			for(unsigned i = 0; i < rows(); ++i)
			{
				for(unsigned j = 0; j < cols(); ++j)
				{
					r(i, j) *= t;
				}
			}

			return r;
		}

		matrix<T, M, N> operator /(T t) const
		{
			auto a = (T)1 / t;
			return *this * ((T)1 / t);
		}

		void operator +=(T t)
		{
			*this = *this + t;
		}

		void operator -=(T t)
		{
			*this = *this - t;
		}

		void operator *=(T t)
		{
			*this = *this * t;
		}

		void operator /=(T t)
		{
			*this = *this / t;
		}

		T operator ()(unsigned int i, unsigned int j) const { return m_values[i][j]; }
		T& operator ()(unsigned int i, unsigned int j) { return m_values[i][j]; }

		std::string to_string() const
		{
			std::ostringstream s;

			for(unsigned i = 0; i < M; ++i)
			{
				s << "[ ";
				for(unsigned j = 0; j < N - 1; ++j)
				{
					s << (*this)(i, j) << ", ";
				}

				s << (*this)(i, N - 1) << " ]";

				if(i != M - 1) { s << std::endl; }
			}

			return s.str();
		}

		friend std::ostream&  operator <<(std::ostream& o, const matrix<T, M, N>& m)
		{
			return o << m.to_string();
		}

	protected:
		T m_values[M][N];
	};

	template<typename T>
	using matrix2 = matrix<T, 2, 2>;

	typedef matrix2<float> matrix2f;
	typedef matrix2<int> matrix2i;

	template<typename T>
	using matrix3 = matrix<T, 3, 3>;

	typedef matrix3<float> matrix3f;
	typedef matrix3<int> matrix3i;

	template<typename T>
	using matrix4 = matrix<T, 4, 4>;

	typedef matrix4<float> matrix4f;
	typedef matrix4<int> matrix4i;
};

#endif