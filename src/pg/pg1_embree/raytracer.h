#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/
class Raytracer : public SimpleGuiDX11
{
public:
	Raytracer( const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3" );
	~Raytracer();

	int InitDeviceAndScene( const char * config );

	int ReleaseDeviceAndScene();

	void LoadScene( const std::string file_name );

	Color4f get_pixel( const int x, const int y, const float t = 0.0f ) override;

	int Ui();

private:
	Vector3 calc_diffuse(RTCRayHit ray_hit);
	Vector3 calc_normal(RTCRayHit ray_hit);
	Vector3 get_fragment_position(RTCRayHit ray_hit);
	Vector3 calc_light_dir(RTCRayHit ray_hit);
	Color4f calc_blinn_phong(RTCRayHit ray_hit);
	RTCRayHit generate_ray_hit(RTCRay ray);
	RTCRay generate_ray(Vector3 position, Vector3 direction, float tfar = FLT_MAX);
	bool generate_shadow_ray(Vector3 position, Vector3 direction);

private:
	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;

	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;

	// my decl
	Material* current_material = nullptr;
	//Vector3 light_position{ 50, 20, 50 };
	Vector3 light_position{ 175, 80, 100 };

};
