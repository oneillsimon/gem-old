#ifndef __GEM_QUATERNION_H__
#define __GEM_QUATERNION_H__

#include "vector.h"

namespace gem
{
	class quaternion : public gem::vector<float, 4>
	{
	public:
		quaternion() : gem::vector<float, 4>() {}
		quaternion(const std::initializer_list<float>& l) : gem::vector<float, 4>(l) {}
		quaternion(const std::vector<float>& v) : gem::vector<float, 4>(v) {}
		quaternion(const gem::matrix<float, 1, 4>& m) : gem::vector<float, 4>(m) {}

		quaternion(const gem::vector<float, 3>& axis, float angle)
		{
			float sinHalfAngle = sin(angle / 2.0f);
			float cosHalfAngle = cos(angle / 2.0f);

			set_x(axis.x() * sinHalfAngle);
			set_y(axis.y() * sinHalfAngle);
			set_z(axis.z() * sinHalfAngle);
			set_w(cosHalfAngle);
		}
	};
};

#endif