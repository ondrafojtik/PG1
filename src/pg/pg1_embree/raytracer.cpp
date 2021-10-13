#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include <math.h>

#define SPECULAR_STRENGTH 0.5

Raytracer::Raytracer( const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config ) : SimpleGuiDX11( width, height )
{
	InitDeviceAndScene( config );

	camera_ = Camera( width, height, fov_y, view_from, view_at );
	//background = new Texture("C:\\dev\\pg1_template_embree_vs2019\\data\\snowy_cemetery.jpg");
	//background = new Texture("C:\\dev\\pg1_template_embree_vs2019\\data\\large_corridor.jpg");
	background = new Texture("C:\\dev\\pg1_template_embree_vs2019\\data\\photo_studio_loft_hall.jpg");
	
}

Raytracer::~Raytracer()
{
	ReleaseDeviceAndScene();
}

int Raytracer::InitDeviceAndScene( const char * config )
{
	device_ = rtcNewDevice( config );
	error_handler( nullptr, rtcGetDeviceError( device_ ), "Unable to create a new device.\n" );
	rtcSetDeviceErrorFunction( device_, error_handler, nullptr );

	ssize_t triangle_supported = rtcGetDeviceProperty( device_, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED );

	// create a new scene bound to the specified device
	scene_ = rtcNewScene( device_ );

	return S_OK;
}

int Raytracer::ReleaseDeviceAndScene()
{
	rtcReleaseScene( scene_ );
	rtcReleaseDevice( device_ );

	return S_OK;
}

void Raytracer::LoadScene( const std::string file_name )
{
	const int no_surfaces = LoadOBJ( file_name.c_str(), surfaces_, materials_ );

	// surfaces loop
	for ( auto surface : surfaces_ )
	{
		RTCGeometry mesh = rtcNewGeometry( device_, RTC_GEOMETRY_TYPE_TRIANGLE );

		Vertex3f * vertices = ( Vertex3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof( Vertex3f ), 3 * surface->no_triangles() );

		Triangle3ui * triangles = ( Triangle3ui * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof( Triangle3ui ), surface->no_triangles() );

		rtcSetGeometryUserData( mesh, ( void* )( surface->get_material() ) );

		rtcSetGeometryVertexAttributeCount( mesh, 2 );

		Normal3f * normals = ( Normal3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof( Normal3f ), 3 * surface->no_triangles() );

		Coord2f * tex_coords = ( Coord2f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof( Coord2f ), 3 * surface->no_triangles() );		

		// triangles loop
		for ( int i = 0, k = 0; i < surface->no_triangles(); ++i )
		{
Triangle& triangle = surface->get_triangle(i);

// vertices loop
for (int j = 0; j < 3; ++j, ++k)
{
	const Vertex& vertex = triangle.vertex(j);

	vertices[k].x = vertex.position.x;
	vertices[k].y = vertex.position.y;
	vertices[k].z = vertex.position.z;

	normals[k].x = vertex.normal.x;
	normals[k].y = vertex.normal.y;
	normals[k].z = vertex.normal.z;

	tex_coords[k].u = vertex.texture_coords[0].u;
	tex_coords[k].v = vertex.texture_coords[0].v;
} // end of vertices loop

triangles[i].v0 = k - 3;
triangles[i].v1 = k - 2;
triangles[i].v2 = k - 1;
		} // end of triangles loop

		rtcCommitGeometry(mesh);
		unsigned int geom_id = rtcAttachGeometry(scene_, mesh);
		rtcReleaseGeometry(mesh);
	} // end of surfaces loop

	rtcCommitScene(scene_);
}

Color4f add_color(Color4f c1, Color4f c2)
{
	Color4f result{};
	result.r = c1.r + c2.r;
	result.g = c1.g + c2.g;
	result.b = c1.b + c2.b;
	result.a = 1.0f;

	return result;
}

Color4f multiply_color(Color4f c1, Color4f c2)
{
	Color4f result{};
	result.r = c1.r * c2.r;
	result.g = c1.g * c2.g;
	result.b = c1.b * c2.b;
	result.a = c1.a * c2.a;

	return result;
}

Color4f multiply_color(Color4f color, float value)
{
	Color4f result{};
	result.r = color.r * value;
	result.g = color.g * value;
	result.b = color.b * value;
	result.a = 1.0f;

	return result;
}

Color4f multiply_color(Color4f color, Vector3 value)
{
	Color4f result{};
	result.r = color.r * value.x;
	result.g = color.g * value.y;
	result.b = color.b * value.z;
	result.a = 1.0f;

	return result;
}

Vector3 reflect(Vector3 data, Vector3 normal)
{
	Vector3 reflect_dir = data - 2 * (data.DotProduct(normal) * normal);
	reflect_dir.Normalize();

	return reflect_dir;
}

Vector3 Raytracer::calc_diffuse(RTCRayHit ray_hit)
{
	current_material = surfaces_[ray_hit.hit.geomID]->get_material();

	Color3f diff{ 0, 0, 0 };
	float u = ray_hit.hit.u;
	float v = ray_hit.hit.v;

	RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);

	if (current_material->get_texture(current_material->kDiffuseMapSlot) != nullptr)
	{
		Coord2f tex_coord__;
		rtcInterpolate0(geometry, ray_hit.hit.primID, u, v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord__.u, 2);

		// TODO(ondra): facing wrong direction
		tex_coord__.v = abs(1 - tex_coord__.v);

		diff = current_material->get_texture(current_material->kDiffuseMapSlot)->get_texel(tex_coord__.u, tex_coord__.v);
	}
	else
		diff = { current_material->diffuse.x, current_material->diffuse.y, current_material->diffuse.z };


	return Vector3{ diff.r, diff.g, diff.b };
}

// THIS DOESNT WORK
Vector3 Raytracer::calc_normal(RTCRayHit ray_hit)
{
	RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);

	Vector3 normal_{};
	rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
		RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal_.x, 3);
	normal_.Normalize();

	if (normal_.x * ray_hit.ray.dir_x + normal_.y * ray_hit.ray.dir_y + normal_.z * ray_hit.ray.dir_z > 0.0f)
	{
		normal_.x *= -1;
		normal_.y *= -1;
		normal_.z *= -1;
	}

	//if (normal_.DotProduct({ ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z }) >= 0)
	//	normal_ = -normal_;

	return normal_;
}

// THIS DOESNT WORK
Vector3 Raytracer::get_fragment_position(RTCRayHit ray_hit)
{
#if 0
	Vertex A = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(0);
	Vertex B = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(1);
	Vertex C = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(2);
	
	float u = ray_hit.hit.v;
	float v = ray_hit.hit.u;
	float w = 1 - u - v;

	// https://spec.oneapi.io/oneart/latest/embree-spec.html

	Vector3 point_hit
	{
		(1.0f - u - v) * A.position.x + u * B.position.x + v * C.position.x,
		(1.0f - u - v) * A.position.y + u * B.position.y + v * C.position.y,
		(1.0f - u - v) * A.position.z + u * B.position.z + v * C.position.z
	};
	return point_hit;

#else
	Vector3 frag_pos
	{
		ray_hit.ray.org_x + (ray_hit.ray.dir_x * ray_hit.ray.tfar),
		ray_hit.ray.org_y + (ray_hit.ray.dir_y * ray_hit.ray.tfar),
		ray_hit.ray.org_z + (ray_hit.ray.dir_z * ray_hit.ray.tfar)
	};
	// calc normal
	
	Vector3 normal = calc_normal(ray_hit);
	Vector3 light_dir = light_position - frag_pos;
	light_dir.Normalize();

	if ((normal.x * light_dir.x + normal.y * light_dir.y + normal.z * light_dir.z) > 0.0f) {
		normal.x *= -1;
		normal.y *= -1;
		normal.z *= -1;
	}
	if (normal.DotProduct(light_dir) < 0)
		normal = -normal;

	frag_pos += normal * 0.001f;
	
	return frag_pos;

#endif
}

Vector3 Raytracer::calc_light_dir(RTCRayHit ray_hit)
{
	//Vector3 light_position{ -200, 150, -120 };
	//Vector3 light_position{ -50, 0, -40 };
#if 0
	Vector3 frag_pos = get_fragment_position(ray_hit);
	Vector3 light_direction
	{
		frag_pos.x - light_position.x,
		frag_pos.y - light_position.y,
		frag_pos.z - light_position.z
	};
	light_direction.Normalize();

	return light_direction;
#else
	// this should be the correct one..
	Vector3 frag_pos = get_fragment_position(ray_hit);
	Vector3 light_direction
	{
		light_position.x - frag_pos.x,
		light_position.y - frag_pos.y,
		light_position.z - frag_pos.z
	};
	light_direction.Normalize();

	return light_direction;
#endif

}

Color4f Raytracer::calc_blinn_phong(RTCRayHit ray_hit)
{
	Vector3 diff = calc_diffuse(ray_hit);
	Vector3 normal = calc_normal(ray_hit);
	Vector3 light_dir = calc_light_dir(ray_hit);

	// diffuse
	float diffuse_strength = light_dir.DotProduct(normal);
	diffuse_strength = abs(diffuse_strength);

	diff.x = diff.x * diffuse_strength;
	diff.y = diff.y * diffuse_strength;
	diff.z = diff.z * diffuse_strength;

	// reflecton
	if (normal.DotProduct(camera_.view_direction) < 0)
	{
		normal.x *= -1;
		normal.y *= -1;
		normal.z *= -1;
	}
	Vector3 reflect_dir = reflect(camera_.view_direction, normal);
	
	Vector3 halfway_dir = camera_.view_direction + light_dir;
	halfway_dir.Normalize();
	//float spec_ = pow(max(normal.DotProduct(halfway_dir), 0.0), current_material->shininess);
	float spec_ = pow(max(normal.DotProduct(reflect_dir), 0.0), current_material->shininess);
	Vector3 spec{};
	spec.x = SPECULAR_STRENGTH * spec_;
	spec.y = SPECULAR_STRENGTH * spec_;
	spec.z = SPECULAR_STRENGTH * spec_;

	//relection = direction - 2 (dot(direction, normal))normal;

	Color4f final_color{ 0, 0, 0, 1 };
	final_color.r = diff.x + spec.x;
	final_color.g = diff.y + spec.y;
	final_color.b = diff.z + spec.z;

	return final_color;
}

RTCRay Raytracer::generate_ray(Vector3 position, Vector3 direction, float tfar)
{
	RTCRay ray = RTCRay();

	ray.org_x = position.x;
	ray.org_y = position.y;
	ray.org_z = position.z;
	ray.tnear = 0.01f; // start of ray segment
	
	ray.dir_x = direction.x;
	ray.dir_y = direction.y;
	ray.dir_z = direction.z;
	ray.time = 0.0f; // time of this ray for motion blur
	ray.tfar = tfar; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}

bool Raytracer::generate_shadow_ray(Vector3 position, Vector3 light_position)
{
	Vector3 x1 = position;
	Vector3 x2 = light_position;

	Vector3 direction = x2 - x1;
	direction.Normalize();

	float distance = sqrt(pow((x2.x - x1.x), 2) + pow((x2.y - x1.y), 2) + pow((x2.z - x1.z), 2));
	RTCRay ray = generate_ray(x1, direction, distance);
	RTCRayHit ray_hit = generate_ray_hit(ray);

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		return 1;
	}
	return 0;
}

RTCRayHit Raytracer::generate_ray_hit(RTCRay ray)
{
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f; // geometry normal
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;

	// merge ray and hit structures
	RTCRayHit ray_hit;
	ray_hit.ray = ray;
	ray_hit.hit = hit;
	
	// intersect ray with the scene
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(scene_, &context, &ray_hit);
	return ray_hit;
}

Color4f Raytracer::shader(RTCRayHit ray_hit, float ior)
{
	// return diffuse?
	//if (depth == 5)
	//{
	//	depth = 0;
	//	return calc_blinn_phong(ray_hit);
	//}


	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		// what did we hit 
		if (surfaces_[ray_hit.hit.geomID]->get_material()->get_name() == "green_plastic_transparent" ||
			surfaces_[ray_hit.hit.geomID]->get_material()->get_name() == "wire_214229166")
		{
			if (depth == 10)
			{
				depth = 0;
				return { 1, 1, 1, 1 };//shader(ray_hit);
			}
			depth += 1;

			////////////////
			float n1 = ior;
			float n2;

			if (n1 == 1.0f)
			{
				n2 = 1.5f;
			}
			else
				n2 = 1.0f;

			Vector3 d = { ray_hit.ray.dir_x,ray_hit.ray.dir_y,ray_hit.ray.dir_z };
			Vector3 v = -d;
			RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
			Vector3 normal{};
			rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
				RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
			normal.Normalize();

			if ((normal.x * v.x + normal.y * v.y + normal.z * v.z) > 0.0f) {
				normal.x *= -1;
				normal.y *= -1;
				normal.z *= -1;
			}
			if (normal.DotProduct(v) < 0)
				normal = -normal;


			Vector3 n = normal;
			float cos_01 = n.DotProduct(v);
			if (cos_01 < 0)
				cos_01 = -normal.DotProduct(v);

			float cos_02 = sqrt(1 - pow((n1 / n2), 2) * (1 - pow(cos_01, 2)));
			Vector3 refracted_dir{};
			refracted_dir.x = (n1 / n2) * d.x + ((n1 / n2) * cos_01 - cos_02) * n.x;
			refracted_dir.y = (n1 / n2) * d.y + ((n1 / n2) * cos_01 - cos_02) * n.y;
			refracted_dir.z = (n1 / n2) * d.z + ((n1 / n2) * cos_01 - cos_02) * n.z;

			Vector3 reflected_dir{};
			reflected_dir.x = (2 * v.DotProduct(n)) * n.x - v.x;
			reflected_dir.y = (2 * v.DotProduct(n)) * n.y - v.y;
			reflected_dir.z = (2 * v.DotProduct(n)) * n.z - v.z;

			float Rs = pow((n2 * cos_02 - n1 * cos_02) / (n2 * cos_02 + n1 * cos_01), 2);
			float Rp = pow((n2 * cos_01 - n1 * cos_02) / (n2 * cos_01 + n1 * cos_02), 2);
			float R = (Rs + Rp) / 2;	// refl
			float T = 1 - R;			// refr

			Vector3 p{};
			p.x = ray_hit.ray.org_x + (ray_hit.ray.dir_x * ray_hit.ray.tfar);
			p.y = ray_hit.ray.org_y + (ray_hit.ray.dir_y * ray_hit.ray.tfar);
			p.z = ray_hit.ray.org_z + (ray_hit.ray.dir_z * ray_hit.ray.tfar);

			RTCRay refr_ray = generate_ray(p, refracted_dir);
			RTCRay refl_ray = generate_ray(p, reflected_dir);

			RTCRayHit refr_ray_hit = generate_ray_hit(refr_ray);
			RTCRayHit refl_ray_hit = generate_ray_hit(refl_ray);

			//Color4f r = shader(refl_ray_hit);
			Color4f t = shader(refr_ray_hit, n2);

			Vector3 att_coef = surfaces_[ray_hit.hit.geomID]->get_material()->emission;
			att_coef = { 1.0f, 0.1f, 1.0f };
			att_coef = { 0.2f, 0.2f, 0.2f };
			att_coef = { 0.0f, 0.0f, 0.0f };


			Vector3 length{};
			length.x = ray_hit.ray.org_x - p.x;
			length.y = ray_hit.ray.org_y - p.y;
			length.z = ray_hit.ray.org_z - p.z;
			Vector3 attenuation{};
			attenuation.x = exp(-1.0f * att_coef.x * length.L2Norm());
			attenuation.y = exp(-1.0f * att_coef.y * length.L2Norm());
			attenuation.z = exp(-1.0f * att_coef.z * length.L2Norm());
			//Color4f color = add_color(multiply_color(t, T), multiply_color(r, R));
			Color4f color = multiply_color(t, T);
			return multiply_color(color, attenuation);

			////////////////

			//
			//Vector3 position{};
			//position.x = ray_hit.ray.org_x * ray_hit.ray.tfar;
			//position.y = ray_hit.ray.org_y * ray_hit.ray.tfar;
			//position.z = ray_hit.ray.org_z * ray_hit.ray.tfar;
			//
			//Vector3 direction{};
			//direction.x = -ray_hit.ray.dir_x;
			//direction.y = -ray_hit.ray.dir_y;
			//direction.z = -ray_hit.ray.dir_z;
			//direction.Normalize();
			//
			//RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
			//Vector3 normal{};
			//rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
			//	RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
			//normal.Normalize();
			//
			//if ((normal.x * ray_hit.ray.dir_x + normal.y * ray_hit.ray.dir_y + normal.z * ray_hit.ray.dir_z) > 0.0f) {
			//	normal.x *= -1;
			//	normal.y *= -1;
			//	normal.z *= -1;
			//}
			//if (normal.DotProduct(direction) < 0)
			//	normal = -normal;
			//
			//
			//float n1 = ior;
			//float n2;
			//
			//if (n1 == 1.0f)
			//{
			//	n2 = 1.5f;
			//}
			//else
			//	n2 = 1.0f;
			//
			//// --TEST--
			///*
			//n1 = 1.5f;
			//n2 = 1.0f;
			//direction.x = -0.429;
			//direction.y = -0.903;
			//direction.z = 0;
			//
			//normal.x = 0;
			//normal.y = 1;
			//normal.z = 0;
			//*/
			//
			//// refract
			////350
			//
			//float div_n1n2 = (n1 / n2);
			//float dot_dn = -direction.DotProduct(normal);
			//
			//float cos_01 = (normal).DotProduct(direction);
			//if (cos_01 < 0) {
			//	cos_01 = (-normal).DotProduct(direction);
			//}
			//
			//float cos_02 = sqrt(1 - pow(div_n1n2, 2) * (1 - pow(cos_01, 2)));
			//float Rs = powf((n2 * cos_02 - n1 * cos_01) / (n2 * cos_02 + n1 * cos_01), 2);
			//float Rp = powf((n2 * cos_01 - n1 * cos_02) / (n2 * cos_01 + n1 * cos_02), 2);
			//float reflectivity__ = (Rs + Rp) * 0.5f; // R
			//float refractivity__ = 1.0f - reflectivity__; // T
			//
			//float phi = (div_n1n2 * dot_dn + sqrt(1 - pow(div_n1n2, 2) * (1 - pow(dot_dn, 2))));
			//Vector3 direction_refracted{};
			//// check the direction i normal  
			//direction_refracted.x = div_n1n2 * ray_hit.ray.dir_x + (div_n1n2 * cos_01 - cos_02) * normal.x;
			//direction_refracted.y = div_n1n2 * ray_hit.ray.dir_y + (div_n1n2 * cos_01 - cos_02) * normal.y;
			//direction_refracted.z = div_n1n2 * ray_hit.ray.dir_z + (div_n1n2 * cos_01 - cos_02) * normal.z;
			//
			////direction_refracted.x = div_n1n2 * -direction.x - phi * normal.x;
			////direction_refracted.y = div_n1n2 * -direction.y - phi * normal.y;
			////direction_refracted.z = div_n1n2 * -direction.z - phi * normal.z;
			//direction_refracted.Normalize();
			//
			//RTCRay ray_refracted = generate_ray(position, direction_refracted);
			//RTCRayHit refracted_ray_hit = generate_ray_hit(ray_refracted);
			//Color4f refracted_color = shader(refracted_ray_hit, n2);
			//
			//// reflect
			//Vector3 reflected{};
			//reflected.x = 2 * normal.DotProduct(direction) * normal.x - direction.x;
			//reflected.y = 2 * normal.DotProduct(direction) * normal.y - direction.y;
			//reflected.z = 2 * normal.DotProduct(direction) * normal.z - direction.z;
			//Vector3 d{ ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z };
			//reflected.x = d.x - 2 * (d.DotProduct(normal) * normal.x);
			//reflected.y = d.y - 2 * (d.DotProduct(normal) * normal.y);
			//reflected.z = d.z - 2 * (d.DotProduct(normal) * normal.z);
			//
			//reflected.Normalize();
			//
			//RTCRay ray = generate_ray(position, reflected);
			//RTCRayHit reflected_ray_hit = generate_ray_hit(ray);
			//Color4f blinn_phong = calc_blinn_phong(ray_hit);
			//Color4f reflected_color = shader(reflected_ray_hit);
			//float reflectivity = surfaces_[ray_hit.hit.geomID]->get_material()->reflectivity;
			////return multiply_color(blinn_phong, shader(reflected_ray_hit));
			////multiply_color(reflected_color, reflectivity__);
			//
			//Vector3 emission = surfaces_[ray_hit.hit.geomID]->get_material()->emission;
			//Vector3 a_pos{};
			//a_pos.x = ray_hit.ray.org_x - (ray_hit.ray.org_x * ray_hit.ray.tfar);
			//a_pos.y = ray_hit.ray.org_y - (ray_hit.ray.org_y * ray_hit.ray.tfar);
			//a_pos.z = ray_hit.ray.org_z - (ray_hit.ray.org_z * ray_hit.ray.tfar);
			//Color4f attentuation{};
			//attentuation.r = exp(-1 * emission.x * a_pos.L2Norm());
			//attentuation.g = exp(-1 * emission.y * a_pos.L2Norm());
			//attentuation.b = exp(-1 * emission.z * a_pos.L2Norm());
			//attentuation.a = 1.0f;
			//
			//return add_color(multiply_color(shader(refracted_ray_hit), refractivity__), 
			//				 multiply_color(shader(reflected_ray_hit), 1 - refractivity__));
			//return multiply_color(shader(refracted_ray_hit), attentuation);
			//
			////return refracted_color;
			//
			////Color4f amp = { 0.8, 1.0, 0.1, 1 };
			//Color4f amp = { 1, 1, 1, 1 };
			//
			//Color4f final__{};
			//final__.r = reflected_color.r * reflectivity__ + refracted_color.r * refractivity__ * amp.r;
			//final__.g = reflected_color.g * reflectivity__ + refracted_color.g * refractivity__ * amp.g;
			//final__.b = reflected_color.b * reflectivity__ + refracted_color.b * refractivity__ * amp.b;
			//final__.a = 1.0f;
			//return final__;
			//
			////return multiply_color(refracted_color, refractivity__);
			//return add_color(multiply_color(reflected_color, reflectivity__), multiply_color(refracted_color, refractivity__));
			//
			//return add_color(multiply_color(blinn_phong, reflectivity), multiply_color(reflected_color, 1 - reflectivity));
			//
			//
			//
			//return Color4f{ 1, 1, 1, 1 };
		}
		else if (surfaces_[ray_hit.hit.geomID]->get_material()->get_name() == "black_plastic")
		{
			if (depth == 5)
			{
				depth = 0;
				return Color4f{ 1, 1, 1, 1 };// shader(ray_hit);
			}
			depth += 1;

			Vector3 position{};
			position.x = ray_hit.ray.org_x * ray_hit.ray.tfar;
			position.y = ray_hit.ray.org_y * ray_hit.ray.tfar;
			position.z = ray_hit.ray.org_z * ray_hit.ray.tfar;

			Vector3 direction{};
			direction.x = -ray_hit.ray.dir_x;
			direction.y = -ray_hit.ray.dir_y;
			direction.z = -ray_hit.ray.dir_z;
			direction.Normalize();

			RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
			Vector3 normal{};
			rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
				RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
			normal.Normalize();

			/*
			if (normal.DotProduct(direction) < 0)
			{
				normal.x *= -1;
				normal.y *= -1;
				normal.z *= -1;
			}*/
			if ((normal.x * ray_hit.ray.dir_x + normal.y * ray_hit.ray.dir_y + normal.z * ray_hit.ray.dir_z) > 0.0f) {
				normal.x *= -1;
				normal.y *= -1;
				normal.z *= -1;
			}
			if (normal.DotProduct(direction) < 0)
				normal = -normal;

			Vector3 reflected{};
			reflected.x = 2 * normal.DotProduct(direction) * normal.x - direction.x;
			reflected.y = 2 * normal.DotProduct(direction) * normal.y - direction.y;
			reflected.z = 2 * normal.DotProduct(direction) * normal.z - direction.z;
			reflected.Normalize();

			float n2 = 1.0f;
			float n1 = 1.0f;

			float cos_01 = (normal).DotProduct(direction);
			if (cos_01 < 0) {
				cos_01 = (-normal).DotProduct(direction);
			}

			float n_d = n1 / n2;
			float sqrt_value = 1 - powf(n_d, 2) * (1 - powf(cos_01, 2));

			float cos_02 = sqrt(sqrt_value);

			float Rs = powf((n2 * cos_02 - n1 * cos_01) / (n2 * cos_02 + n1 * cos_01), 2);
			float Rp = powf((n2 * cos_01 - n1 * cos_02) / (n2 * cos_01 + n1 * cos_02), 2);
			float reflectivity__ = 0.5f * (Rs + Rp);
			float refractivity__ = 1.0f - reflectivity__;


			RTCRay ray = generate_ray(position, reflected);
			RTCRayHit reflected_ray_hit = generate_ray_hit(ray);
			Color4f blinn_phong = calc_blinn_phong(ray_hit);
			Color4f reflected_color = shader(reflected_ray_hit);
			float reflectivity = surfaces_[ray_hit.hit.geomID]->get_material()->reflectivity;
			reflectivity = 0.8f;
			reflectivity = 1.0f - reflectivity__;
			reflectivity = 0.85f;
			//return multiply_color(blinn_phong, shader(reflected_ray_hit));
			//multiply_color(reflected_color, reflectivity__);
			Vector3 specular = surfaces_[ray_hit.hit.geomID]->get_material()->specular;

			//return multiply_color(shader(reflected_ray_hit), reflectivity);

			return add_color(multiply_color(blinn_phong, reflectivity), multiply_color(reflected_color, 1-reflectivity));

			//Color4f blinn_phong = calc_blinn_phong(ray_hit);
			//Vector3 direction{ ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z };
			//direction.Normalize();
			//RTCRayHit reflect_ray_hit = generate_ray_hit(generate_ray(get_fragment_position(ray_hit), reflect(direction, //calc_normal(ray_hit))));
			//Color4f reflected = shader(reflect_ray_hit);
			//return multiply_color(blinn_phong, multiply_color(reflected, 0.9f));

		}
		else // white plastic
		{
			if (depth == 5)
			{
				depth = 0;
				return {1, 1, 1, 1}; shader(ray_hit);
			}
			depth += 1;


			Vector3 position{};
			position.x = ray_hit.ray.org_x * ray_hit.ray.tfar;
			position.y = ray_hit.ray.org_y * ray_hit.ray.tfar;
			position.z = ray_hit.ray.org_z * ray_hit.ray.tfar;

			Vector3 direction{};
			direction.x = -ray_hit.ray.dir_x;
			direction.y = -ray_hit.ray.dir_y;
			direction.z = -ray_hit.ray.dir_z;
			direction.Normalize();

			RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
			Vector3 normal{};
			rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
				RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
			normal.Normalize();

			/*if (normal.DotProduct(direction) < 0)
			{
				normal.x *= -1;
				normal.y *= -1;
				normal.z *= -1;
			}*/

			if ((normal.x * ray_hit.ray.dir_x + normal.y * ray_hit.ray.dir_y + normal.z * ray_hit.ray.dir_z) > 0.0f) {
				normal.x *= -1;
				normal.y *= -1;
				normal.z *= -1;
			}
			if (normal.DotProduct(direction) < 0)
				normal = -normal;

			Vector3 reflected{}; 
			reflected.x = 2 * normal.DotProduct(direction) * normal.x - direction.x;
			reflected.y = 2 * normal.DotProduct(direction) * normal.y - direction.y;
			reflected.z = 2 * normal.DotProduct(direction) * normal.z - direction.z;
			reflected.Normalize();

			RTCRay ray = generate_ray(position, reflected);
			RTCRayHit reflected_ray_hit = generate_ray_hit(ray);
			Color4f blinn_phong = calc_blinn_phong(ray_hit);
			Color4f reflected_color = shader(reflected_ray_hit);
			float reflectivity = surfaces_[ray_hit.hit.geomID]->get_material()->reflectivity;
			reflectivity = 0.8f;
			return add_color(multiply_color(blinn_phong, reflectivity), multiply_color(reflected_color, 1-reflectivity));


			//Color4f blinn_phong = calc_blinn_phong(ray_hit);
			//Vector3 direction{ ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z };
			//direction.Normalize();
			//RTCRayHit reflect_ray_hit = generate_ray_hit(generate_ray(get_fragment_position(ray_hit), reflect(direction, //calc_normal(ray_hit))));
			//Color4f reflected = shader(reflect_ray_hit);
			//return multiply_color(blinn_phong, multiply_color(reflected, 0.9f));
		}
	}
	depth = 0;
	// background
	Vector3 dir{ -ray_hit.ray.dir_x, -ray_hit.ray.dir_y, ray_hit.ray.dir_z };
	const float theta = acos(dir.z);
	const float phi = atan2f(dir.y, dir.x) + float(3.14159265358979323846);
	const float u = 1.0f - phi * 0.5f * float(0.318309886183790671538);
	const float v = theta * float(0.318309886183790671538);
	Color3f bg = background->get_texel(u, v);

	return Color4f{ bg.r, bg.g, bg.b, 1.0f };

	return Color4f{ 1, 1, 1, 1 };
	
}

// 640, 480
Color4f Raytracer::get_pixel( const int x, const int y, const float t )
{
	// TODO generate primary ray and perform ray cast on the scene
	RTCRay ray = camera_.GenerateRay(x, y);
	
	RTCRayHit ray_hit = generate_ray_hit(ray);
	return shader(ray_hit);

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{

		Color4f blinn_phong = calc_blinn_phong(ray_hit);
		Vector3 frag_pos = get_fragment_position(ray_hit);
		bool test = generate_shadow_ray(frag_pos, light_position);
	
		Color4f final_color{};
		final_color.a = 1.0f;

		current_material = surfaces_[ray_hit.hit.geomID]->get_material();
		Vector3 ambient = current_material->ambient;
		Vector3 diffuse = current_material->diffuse;
		Vector3 specular = current_material->specular;

		final_color.r = ambient.x + diffuse.x + (calc_normal(ray_hit).DotProduct(calc_light_dir(ray_hit)) + specular.x);
		final_color.g = ambient.y + diffuse.y + (calc_normal(ray_hit).DotProduct(calc_light_dir(ray_hit)) + specular.y);
		final_color.b = ambient.z + diffuse.z + (calc_normal(ray_hit).DotProduct(calc_light_dir(ray_hit)) + specular.z);
		return final_color;
		//final_color.r = blinn_phong.r;
		//final_color.g = blinn_phong.g;
		//final_color.b = blinn_phong.b;


		RTCRay test_ = generate_ray(get_fragment_position(ray_hit), reflect({ ray_hit.ray.dir_x,ray_hit.ray.dir_y,ray_hit.ray.dir_z }, calc_normal(ray_hit)));

		Color4f bonus{ 0,0,0,1 };

		float reflectivity = 0;
		RTCRayHit test__ = generate_ray_hit(test_);
		if (test__.hit.geomID != RTC_INVALID_GEOMETRY_ID)
		{
			bonus = calc_blinn_phong(test__);
			reflectivity = surfaces_[ray_hit.hit.geomID]->get_material()->reflectivity;
		
			// reflectivity
			final_color.r = (1 - reflectivity) * bonus.r + reflectivity * final_color.r;
			final_color.g = (1 - reflectivity) * bonus.g + reflectivity * final_color.g;
			final_color.b = (1 - reflectivity) * bonus.b + reflectivity * final_color.b;

		}
		
		float shadow_factor = 0.35f;
		if (test)
		{
			final_color.r = final_color.r * shadow_factor;
			final_color.g = final_color.g * shadow_factor;
			final_color.b = final_color.b * shadow_factor;
		}

		return final_color;
	}

	return Color4f{ (float)x/ 640.0f, 0.0f, (float)y / 480.0f, 1.0f };
}

int Raytracer::Ui()
{
	static float f = 0.0f;
	static int counter = 0;

	// Use a Begin/End pair to created a named window
	ImGui::Begin( "Ray Tracer Params" );
	
	ImGui::Text( "Surfaces = %d", surfaces_.size() );
	ImGui::Text( "Materials = %d", materials_.size() );
	ImGui::Separator();
	ImGui::Checkbox( "Vsync", &vsync_ );
	
	//ImGui::Checkbox( "Demo Window", &show_demo_window ); // Edit bools storing our window open/close state
	//ImGui::Checkbox( "Another Window", &show_another_window );

	ImGui::SliderFloat( "float", &f, 0.0f, 1.0f ); // Edit 1 float using a slider from 0.0f to 1.0f    
	//ImGui::ColorEdit3( "clear color", ( float* )&clear_color ); // Edit 3 floats representing a color

	// Buttons return true when clicked (most widgets return true when edited/activated)
	if ( ImGui::Button( "Button" ) )
		counter++;
	ImGui::SameLine();
	ImGui::Text( "counter = %d", counter );

	ImGui::Text( "Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate );
	ImGui::End();

	// 3. Show another simple window.
	/*if ( show_another_window )
	{
	ImGui::Begin( "Another Window", &show_another_window ); // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
	ImGui::Text( "Hello from another window!" );
	if ( ImGui::Button( "Close Me" ) )
	show_another_window = false;
	ImGui::End();
	}*/

	return 0;
}



/*
		DIFFERENT WAY OF CALCULATING THE HIT POINT

		Vertex A = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(0);
		Vertex B = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(1);
		Vertex C = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(2);

		float u = ray_hit.hit.v;
		float v = ray_hit.hit.u;
		float w = 1 - u - v;

		// https://spec.oneapi.io/oneart/latest/embree-spec.html


		Vector3 point_hit
		{
			(1.0f - u - v) * A.position.x + u * B.position.x + v * C.position.x,
			(1.0f - u - v) * A.position.y + u * B.position.y + v * C.position.y,
			(1.0f - u - v) * A.position.z + u * B.position.z + v * C.position.z
		}

		*/



		// Vratí texel o relativních souradnicich \a u a \a v, kde \f$(u,v)\in\left<0,1\right>^2\f$
/*
Vertex A = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(0);
Vertex B = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(1);
Vertex C = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(2);

float w = 1 - u - v;
float u_ = (A.texture_coords[0].u * w) + (B.texture_coords[0].u * u) + (C.texture_coords[0].u * v);
float v_ = (A.texture_coords[0].v * w) + (B.texture_coords[0].v * u) + (C.texture_coords[0].v * v);

//t_uv = (1-u-v)*t0 + u*t1 + v*t2
float t_u = w * A.texture_coords[0].u + u * B.texture_coords[0].u + v * C.texture_coords[0].u;
float t_v = w * A.texture_coords[0].v + u * B.texture_coords[0].v + v * C.texture_coords[0].v;

//t_u = (t_u + 1) / 2;
//t_v = (t_v + 1) / 2;
//t_u = sqrt(u_)-0.1;
//t_v = sqrt(v_ / 2);
//u_ = sqrt(u_);
//v_ = sqrt(v_ / 2);
//u_ = (u_ + 1) / 2;
//v_ = (v_ + 1) / 2;

// interpolation

*/