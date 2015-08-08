#ifndef __GEM_H__
#define __GEM_H__

#include <assert.h>

#include "matrix.h"
#include "vector.h"

namespace gem
{
	namespace internal
	{
		template<typename T, unsigned int M, unsigned int N>
		T calc_determinant(const gem::matrix<T, M, N>& o, const gem::matrix<T, M, N>& m, int n)
		{
			T t = (T)0;
			gem::matrix<T, M, N> r = m;

			if(n == 1)
			{
				return r(0, 0);
			}
			else if(n == 2)
			{
				return (r(0, 0) * r(1, 1) - r(0, 1) * r(1, 0));
			}
			else
			{
				for(int p = 0; p < n; p++)
				{
					int h = 0;
					int k = 0;

					for(int i = 1; i < n; i++)
					{
						for(int j = 0; j < n; j++)
						{
							if(j == p)
							{
								continue;
							}

							r(h, k) = o(i, j);
							k++;

							if(k == n - 1)
							{
								h++;
								k = 0;
							}
						}
					}

					t += m(0, p) * (T)pow(-1, p) * calc_determinant(o, r, n - 1);
				}
			}

			return t;
		}
	};

	template<typename T, unsigned int M, unsigned int N>
	constexpr int is_square(const gem::matrix<T, M, N>& m)
	{
		return M == N;
	}

	template<typename T, unsigned int M, unsigned int N,
			 typename T1, unsigned int M1, unsigned int N1>
	constexpr int same_col_count(const gem::matrix<T, M, N>& m, const gem::matrix<T1, M1, N1>& n)
	{
		return N == N1;
	}

	template<typename T, unsigned int M, unsigned int N,
		     typename T1, unsigned int M1, unsigned int N1>
	constexpr int same_row_count(const gem::matrix<T, M, N>& m, const gem::matrix<T1, M1, N1>& n)
	{
		return M == M1;
	}

	template<typename T, unsigned int M, unsigned int N,
		     typename T1, unsigned int M1, unsigned int N1>
	constexpr int same_dimensions(const gem::matrix<T, M, N>& m, const gem::matrix<T1, M1, N1>& n)
	{
		return same_col_count(m, n) && same_row_count(m, n);
	}

	template<typename T, unsigned int M, unsigned int N>
	T determinant(const gem::matrix<T, M, N>& m)
	{
		return internal::calc_determinant(m, m, m.rows());
	}

	template<typename T, unsigned int M>
	gem::matrix<T, M, M> cofactor(const gem::matrix<T, M, M>& m)
	{
		assert(M > 1);
		auto r = gem::matrix<T, M, M>();

		if(r.rows() == 2)
		{
			r(0, 0) =  m(1, 1);
			r(0, 1) = -m(1, 0);
			r(1, 0) = -m(0, 1);
			r(1, 1) =  m(0, 0);
		}
		else
		{
			int r_count = r.rows();
			gem::matrix<T, M - 1, M - 1>*** temp = new gem::matrix<T, M - 1, M - 1>**[M];

			for(int i = 0; i < r_count; i++)
			{
				temp[i] = new gem::matrix<T, M - 1, M - 1>*[M];
			}

			for(int i = 0; i < r_count; i++)
			{
				for(int j = 0; j < r_count; j++)
				{
					temp[i][j] = new gem::matrix<T, M - 1, M - 1>();
				}
			}

			for(int k1 = 0; k1 < r_count; k1++)
			{
				for(int k2 = 0; k2 < r_count; k2++)
				{
					int i1 = 0;
					for(int i = 0; i < r_count; i++)
					{
						int j1 = 0;
						for(int j = 0; j < r_count; j++)
						{
							if(k1 == i || k2 == j)
							{
								continue;
							}

							temp[k1][k2][0](i1, j1++) = m(i, j);
						}

						if(k1 != i)
						{
							i1++;
						}
					}
				}
			}

			bool flag_positive = true;
			for(int k1 = 0; k1 < r_count; k1++)
			{
				flag_positive = ((k1 % 2) == 0);
				for(int k2 = 0; k2 < r_count; k2++)
				{
					if(flag_positive)
					{
						r(k1, k2) = determinant(*temp[k1][k2]);
						flag_positive = false;
					}
					else
					{
						r(k1, k2) = -determinant(*temp[k1][k2]);
						flag_positive = true;
					}
				}
			}

			for(int i = 0; i < r_count; i++)
			{
				for(int j = 0; j < r_count; j++)
				{
					delete temp[i][j];
				}
			}

			for(int i = 0; i < r_count; i++)
			{
				delete[] temp[i];
			}

			delete[] temp;
		}

		return r;
	}

	template<typename T, unsigned int M>
	gem::matrix<T, M, M> identity(const gem::matrix<T, M, M>& m)
	{
		auto r = gem::matrix<T, M, M>();

		for(unsigned i = 0; i < M; ++i)
		{
			r(i, i) = 1;
		}

		return r;
	}

	template<typename T, unsigned int M>
	gem::matrix<T, M, M> inverse(const gem::matrix<T, M, M>& m)
	{
		T d = determinant(m);
		assert(d != 0);

		return cofactor(transpose(m)) * ((T)1 / d);
	}

	template<typename T, unsigned int M>
	gem::matrix<T, M, M> matrix_of_minors(const gem::matrix<T, M, M>& m)
	{
		auto r = gem::matrix<T, M, M>();

		for(unsigned i = 0; i < r.rows(); ++i)
		{
			for(unsigned j = 0; j < r.cols(); ++j)
			{
				auto a = minor_matrix(m, i, j);
				r(i, j) = determinant(minor_matrix(m, i, j));
			}
		}

		return r;
	}

	template<typename T, int M>
	gem::matrix<T, M - 1, M - 1> minor_matrix(const gem::matrix<T, M, M>& m, int row, int col)
	{
		auto r = gem::matrix<T, M - 1, M - 1>();
		int row_count = 0;
		int col_count = 0;

		for(unsigned i = 0; i < M; ++i)
		{
			if(i != row)
			{
				for(unsigned j = 0; j < M; ++j)
				{
					if(j != col)
					{
						r(row_count, col_count++) = m(i, j);
					}
				}

				row_count++;
				col_count = 0;
			}
		}

		return r;
	}

	template<typename T, unsigned int M, unsigned int N>
	gem::matrix<T, N, M> transpose(const gem::matrix<T, M, N>& m)
	{
		auto r = gem::matrix<T, N, M>();

		for(unsigned i = 0; i < r.rows(); ++i)
		{
			for(unsigned j = 0; j < r.cols(); ++j)
			{
				r(i, j) = m(j, i);
			}
		}

		return r;
	}

	template<typename T, unsigned int M>
	std::vector<T> filled_vector(T t)
	{
		std::vector<T> v;

		for(unsigned i = 0; i < M; ++i)
		{
			v.push_back(t);
		}

		return v;
	}

	template<typename T, unsigned int M, unsigned int N>
	void fill(gem::matrix<T, M, N>& m, T t)
	{
		for(unsigned i = 0; i < m.rows(); ++i)
		{
			for(unsigned j = 0; j < m.cols(); ++j)
			{
				m(i, j) = t;
			}
		}
	}

	template<typename T, typename U>
	gem::matrix<T, 1, 3>cross(const gem::matrix<T, 1, 3>& a, const gem::matrix<U, 1, 3>& b)
	{
		return gem::matrix<T, 1, 3>{ a(0, 1) * b(0, 2) - a(0, 2) * b(0, 1),
			                         a(0, 2) * b(0, 0) - a(0, 0) * b(0, 2),
			                         a(0, 0) * b(0, 1) - a(0, 1) * b(0, 0) };
	}

	template<typename T, typename U, unsigned int N>
	double distance(const gem::matrix<T, 1, N>& a, const gem::matrix<U, 1, N>& b)
	{
		double coord_sum = 0;

		for(unsigned i = 0; i < N; ++i)
		{
			coord_sum += (a(0, i) - b(0, i)) * (a(0, i) - b(0, i));
		}

		return sqrt(coord_sum);
	}

	template<typename T, typename U, unsigned int N>
	double dot(const gem::matrix<T, 1, N>&a, const gem::matrix<U, 1, N>& b)
	{
		double d = 0;

		for(unsigned i = 0; i < N; ++i)
		{
			d += a(0, i) * b(0, i);
		}

		return d;
	}

	template<typename T, unsigned int N>
	double length(const gem::matrix<T, 1, N>& v)
	{
		double length = 0;

		for(unsigned i = 0; i < N; ++i)
		{
			length += v(0, i) * v(0, i);
		}

		return sqrt(length);
	}

	template<typename T, unsigned int N>
	gem::matrix<T, 1, N> normalise(const gem::matrix<T, 1, N>& v)
	{
		return v / length(v);
	}

	template<typename T, unsigned int M, unsigned int N>
	T max(const gem::matrix<T, M, N>& m)
	{
		T max = m(0, 0);

		for(unsigned i = 0; i < M; ++i)
		{
			for(unsigned j = 0; j < N; ++j)
			{
				max = m(i, j) > max ? m(i, j) : max;
			}
		}

		return max;
	}

	template<typename T, unsigned int M, unsigned int N>
	T min(const gem::matrix<T, M, N>& m)
	{
		T min = m(0, 0);

		for(unsigned i = 0; i < M; ++i)
		{
			for(unsigned j = 0; j < N; ++j)
			{
				min = m(i, j) < min ? m(i, j) : min;
			}
		}

		return min;
	}

	template<typename T, unsigned int M>
	T trace(const gem::matrix<T, M, M>& m)
	{
		T trace = 0;

		for(unsigned i = 0; i < M; ++i)
		{
			trace += m(i, i);
		}

		return trace;
	}
};

#endif