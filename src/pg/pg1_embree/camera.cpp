#include "stdafx.h"
#include "camera.h"

#include "vector3.h"
// stredove promitani 
// Z jde zpet (do oka)
// Y je nahoru

// CS - camera system
// WS - world system

Camera::Camera( const int width, const int height, const float fov_y,
	const Vector3 view_from, const Vector3 view_at )
{
	width_ = width;
	height_ = height;
	fov_y_ = fov_y;

	view_from_ = view_from;
	view_at_ = view_at;


	Vector3 _forward = view_from - view_at;
	_forward.Normalize();
	Vector3 x_cs = up_.CrossProduct(_forward);
	Vector3 y_cs = _forward.CrossProduct(x_cs);

	// TODO compute focal lenght based on the vertical field of view and the camera resolution
	f_y_ = height_ / (2.0f * tan(fov_y_ / 2.0f));
	// TODO build M_c_w_ matrix	
	M_c_w_ = Matrix3x3(x_cs, y_cs, _forward);
}

RTCRay Camera::GenerateRay( const float x_i, const float y_i ) const
{
	RTCRay ray = RTCRay();

	// TODO fill in ray structure and compute ray direction
	// ray.org_x = ...

	// ray setup
	ray.org_x = view_from_.x;
	ray.org_y = view_from_.y;
	ray.org_z = view_from_.z;
	ray.tnear = FLT_MIN; // start of ray segment
	Vector3 target{x_i - (width_ / 2),(height_ / 2) - y_i, -f_y_};
	target.Normalize();
	Vector3 direction = M_c_w_ * target;
	ray.dir_x = direction.x;
	ray.dir_y =	direction.y;
	ray.dir_z =	direction.z;
	ray.time = 0.0f; // time of this ray for motion blur
	ray.tfar = FLT_MAX; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}
