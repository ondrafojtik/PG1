#ifndef CAMERA_H_
#define CAMERA_H_

#include "vector3.h"
#include "matrix3x3.h"

/*! \class Camera
\brief A simple pin-hole camera.

\author Tom� Fabi�n
\version 1.0
\date 2018
*/
class Camera
{
public:
	Camera() { }

	Camera( const int width, const int height, const float fov_y,
		const Vector3 view_from, const Vector3 view_at );

	Camera(const int width, const int height, const float fov_y,
		const Vector3 view_from, const Vector3 view_at, float focal_point);

	/* generate primary ray, top-left pixel image coordinates (xi, yi) are in the range <0, 1) x <0, 1) */
	RTCRay GenerateRay( const float xi, const float yi ) const;
	Vector3 view_direction{ 0, 0, 0 };
	Matrix3x3 M_c_w_; // transformation matrix from CS -> WS	
	Vector3 view_from_; // ray origin or eye or O
	float fov_y_{ 0.785f }; // vertical field of view (rad)
	float f_y_{ 1.0f }; // focal lenght (px)
	int width_{ 640 }; // image width (px)
	int height_{ 480 };  // image height (px)

private:
	
	Vector3 view_at_; // target T
	Vector3 up_{ Vector3( 0.0f, 0.0f, 1.0f ) }; // up vector

	
	
};

#endif
