#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include <math.h>

#define SPECULAR_STRENGTH 0.5


namespace gamma_util
{
	float c_linear(float c_srgb, float gamma = 2.4f)
	{
		if (c_srgb <= 0.0f)
			return 0.0f;
		else if (c_srgb >= 1.0f)
			return 1.0f;
		assert((c_srgb >= 0.0f) && (c_srgb <= 1.0f));

		if (c_srgb <= 0.04045f)
			return c_srgb / 12.92f;
		else
		{
			const float a = 0.055f;
			return powf((c_srgb + a) / (1.0f + a), gamma);
		}
	}

	float c_srgb(float c_linear, float gamma = 2.4f)
	{
		if (c_linear <= 0.0f)
			return 0.0f;
		else if (c_linear >= 1.0f)
			return 1.0f;

		assert((c_linear >= 0.0f) && (c_linear <= 1.0f));

		if (c_linear <= 0.0031308f)
			return 12.92f * c_linear;
		else
		{
			const float a = 0.055f;
			return (1.0f + a) * powf(c_linear, 1.0f / gamma) - a;
		}
	}

}

Raytracer::Raytracer(const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char* config) : SimpleGuiDX11(width, height)
{
	InitDeviceAndScene(config);

	camera_ = Camera(width, height, fov_y, view_from, view_at);
	//background = new Texture("C:\\dev\\pg1_template_embree_vs2019\\data\\snowy_cemetery.jpg");
	//background = new Texture("C:\\dev\\pg1_template_embree_vs2019\\data\\large_corridor.jpg");
	background = new Texture("../../../data/large_corridor.jpg");
	//background = new Texture("C:\\dev\\pg1_template_embree_vs2019\\data\\photo_studio_loft_hall.jpg");

}

Raytracer::~Raytracer()
{
	ReleaseDeviceAndScene();
}

int Raytracer::InitDeviceAndScene(const char* config)
{
	device_ = rtcNewDevice(config);
	error_handler(nullptr, rtcGetDeviceError(device_), "Unable to create a new device.\n");
	rtcSetDeviceErrorFunction(device_, error_handler, nullptr);

	ssize_t triangle_supported = rtcGetDeviceProperty(device_, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED);

	// create a new scene bound to the specified device
	scene_ = rtcNewScene(device_);

	return S_OK;
}

int Raytracer::ReleaseDeviceAndScene()
{
	rtcReleaseScene(scene_);
	rtcReleaseDevice(device_);

	return S_OK;
}

void Raytracer::LoadScene(const std::string file_name)
{
	const int no_surfaces = LoadOBJ(file_name.c_str(), surfaces_, materials_);

	// surfaces loop
	for (auto surface : surfaces_)
	{
		RTCGeometry mesh = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_TRIANGLE);

		Vertex3f* vertices = (Vertex3f*)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof(Vertex3f), 3 * surface->no_triangles());

		Triangle3ui* triangles = (Triangle3ui*)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof(Triangle3ui), surface->no_triangles());

		rtcSetGeometryUserData(mesh, (void*)(surface->get_material()));

		rtcSetGeometryVertexAttributeCount(mesh, 2);

		Normal3f* normals = (Normal3f*)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof(Normal3f), 3 * surface->no_triangles());

		Coord2f* tex_coords = (Coord2f*)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof(Coord2f), 3 * surface->no_triangles());

		// triangles loop
		for (int i = 0, k = 0; i < surface->no_triangles(); ++i)
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

		diff = current_material->get_texture(current_material->kDiffuseMapSlot)->get_bilinear_texel(tex_coord__.u, tex_coord__.v);
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

Color4f Raytracer::shader(RTCRayHit ray_hit, float depth, float ior)
{
		
	if (depth >= 8 && ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		return {0, 0, 0, 1};
	}


	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		current_material = surfaces_[ray_hit.hit.geomID]->get_material();

		// what did we hit 
		if (surfaces_[ray_hit.hit.geomID]->get_material()->get_name() == "green_plastic_transparent" ||
			surfaces_[ray_hit.hit.geomID]->get_material()->get_name() == "wire_214229166")
		{
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

			if (!(acosf(cos_01) > 0 && (1 - pow((n1 / n2), 2) * (1 - pow(cos_01, 2))) > 0)) 
			{
				return shader(ray_hit, depth + 1, n2);
			}

			float cos_02 = sqrt(1 - pow((n1 / n2), 2) * (1 - pow(cos_01, 2)));
			Vector3 refracted_dir{};
			refracted_dir.x = (n1 / n2) * d.x + ((n1 / n2) * cos_01 - cos_02) * n.x;
			refracted_dir.y = (n1 / n2) * d.y + ((n1 / n2) * cos_01 - cos_02) * n.y;
			refracted_dir.z = (n1 / n2) * d.z + ((n1 / n2) * cos_01 - cos_02) * n.z;
			refracted_dir.Normalize();

			Vector3 reflected_dir{};
			reflected_dir.x = (2 * v.DotProduct(n)) * n.x - v.x;
			reflected_dir.y = (2 * v.DotProduct(n)) * n.y - v.y;
			reflected_dir.z = (2 * v.DotProduct(n)) * n.z - v.z;
			reflected_dir.Normalize();

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

			Color4f r = shader(refl_ray_hit, depth + 1);
			Color4f t = shader(refr_ray_hit, depth + 1, n2);

			
			//float dd = normal.DotProduct(d);
			//if (dd >= -0.2 && dd <= 0.2)
			//{
			//	return {0, 0, 1, 1};
			//}

			Vector3 att_coef = surfaces_[ray_hit.hit.geomID]->get_material()->emission;
			att_coef = { 0.4f, 0.001f, 0.4f };
			//att_coef = { 0.005f, 0.001f, 0.005f };

			//att_coef = { 0.0f, 0.0f, 0.0f };

			Vector3 length{};
			if (ior == 1.0f)
			{
				length.x = ray_hit.ray.org_x - p.x;
				length.y = ray_hit.ray.org_y - p.y;
				length.z = ray_hit.ray.org_z - p.z;
			}
			else
			{
				length.x = 0.0f;
				length.y = 0.0f;
				length.z = 0.0f;
			}
			
			float e = 2.71828182846f;
			Vector3 attenuation{};
			//attenuation.x = pow(e, att_coef.x * length.L2Norm());
			//attenuation.y = pow(e, att_coef.y * length.L2Norm());
			//attenuation.z = pow(e, att_coef.z * length.L2Norm());

			attenuation.x = exp(-1.0f * att_coef.x * length.L2Norm());
			attenuation.y = exp(-1.0f * att_coef.y * length.L2Norm());
			attenuation.z = exp(-1.0f * att_coef.z * length.L2Norm());
			//Color4f color = add_color(multiply_color(t, T), multiply_color(r, R));
			//Color4f color = multiply_color(t, T);

			//R = 1;
			//T = 0;
			Color4f color{ 0, 0, 0, 1 };
			color.r = t.r * T * attenuation.x + r.r * R;// * attenuation.x;
			color.g = t.g * T * attenuation.y + r.g * R;// * attenuation.y;
			color.b = t.b * T * attenuation.z + r.b * R;// * attenuation.z;
			return color;
			//color.r = color.r * attenuation.x;
			//color.g = color.g * attenuation.y;
			//color.b = color.b * attenuation.z;

			//return color;

			// ambient
			//color.r += 0.03f;
			//color.g += 0.03f;
			//color.b += 0.03f;
			//return multiply_color(color, attenuation);

			
		}
		else //if (surfaces_[ray_hit.hit.geomID]->get_material()->get_name() == "black_plastic")
		{
			float n1 = 1.0f;
			float n2 = 1.0f;
			Vector3 d = { ray_hit.ray.dir_x,ray_hit.ray.dir_y,ray_hit.ray.dir_z };
			Vector3 v = -d;
			RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
			Vector3 normal{};
			rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
				RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
			normal.Normalize();

			//if ((normal.x * v.x + normal.y * v.y + normal.z * v.z) > 0.0f) {
			//	normal.x *= -1;
			//	normal.y *= -1;
			//	normal.z *= -1;
			//}
			//if (normal.DotProduct(v) < 0)
			//	normal = -normal;
			if (normal.DotProduct(camera_.view_direction) < 0)
				normal *= -1;

			Vector3 n = normal;
			float cos_01 = n.DotProduct(v);
			if (cos_01 < 0)
				cos_01 = -normal.DotProduct(v);

			float cos_02 = sqrt(1 - pow((n1 / n2), 2) * (1 - pow(cos_01, 2)));
			Vector3 refracted_dir{};
			refracted_dir.x = (n1 / n2) * d.x + ((n1 / n2) * cos_01 - cos_02) * n.x;
			refracted_dir.y = (n1 / n2) * d.y + ((n1 / n2) * cos_01 - cos_02) * n.y;
			refracted_dir.z = (n1 / n2) * d.z + ((n1 / n2) * cos_01 - cos_02) * n.z;
			refracted_dir.Normalize();

			Vector3 reflected_dir{};
			reflected_dir.x = (2 * v.DotProduct(n)) * n.x - v.x;
			reflected_dir.y = (2 * v.DotProduct(n)) * n.y - v.y;
			reflected_dir.z = (2 * v.DotProduct(n)) * n.z - v.z;
			reflected_dir.Normalize();

			float Rs = pow((n2 * cos_02 - n1 * cos_02) / (n2 * cos_02 + n1 * cos_01), 2);
			float Rp = pow((n2 * cos_01 - n1 * cos_02) / (n2 * cos_01 + n1 * cos_02), 2);
			float R = (Rs + Rp) / 2;	// refl
			float T = 1 - R;			// refr

			Vector3 p{};
			p.x = ray_hit.ray.org_x + (ray_hit.ray.dir_x * ray_hit.ray.tfar);
			p.y = ray_hit.ray.org_y + (ray_hit.ray.dir_y * ray_hit.ray.tfar);
			p.z = ray_hit.ray.org_z + (ray_hit.ray.dir_z * ray_hit.ray.tfar);

			Vector3 light_dir{};
			light_dir.x = light_position.x - p.x;
			light_dir.y = light_position.y - p.y;
			light_dir.z = light_position.z - p.z;
			light_dir.Normalize();

			RTCRay refl_ray = generate_ray(p, reflected_dir);

			RTCRayHit refl_ray_hit = generate_ray_hit(refl_ray);
			Color4f reflected_color = shader(refl_ray_hit, depth + 1);

			Color4f blinn_phong = calc_blinn_phong(ray_hit);
			Vector3 diffuse = current_material->diffuse;
			diffuse.x = gamma_util::c_linear(diffuse.x);
			diffuse.y = gamma_util::c_linear(diffuse.z);
			diffuse.z = gamma_util::c_linear(diffuse.z);
			if (current_material->get_texture(current_material->kDiffuseMapSlot) != nullptr)
			{
				float u = ray_hit.hit.u;
				float v = ray_hit.hit.v;
				Coord2f tex_coord__;
				rtcInterpolate0(geometry, ray_hit.hit.primID, u, v,
					RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord__.u, 2);

				tex_coord__.v = abs(1 - tex_coord__.v);
				Color3f diff{ 0, 0, 0 };
				diff = current_material->get_texture(current_material->kDiffuseMapSlot)->get_bilinear_texel(tex_coord__.u, tex_coord__.v);
				diffuse.x = diff.r;
				diffuse.y = diff.g;
				diffuse.z = diff.b;
			}

			// shadow check
			float shadow = 1.0f;
			RTCRay shadow_ray = generate_ray(p, light_dir);
			RTCRayHit shadow_ray_hit = generate_ray_hit(shadow_ray);
			if (shadow_ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
			{
				shadow = 0.3f;
			}

			// phong	
			Vector3 specular = current_material->specular;

			Color4f final_color{ 0, 0, 0, 1 };
			float factor = max(normal.DotProduct(light_dir), 0.0f);
			float cos__ = min(max(v.DotProduct(reflected_dir), 0.0f), 1.0f);
			float spec_ = pow(cos__, current_material->shininess);
			//float spec_ = pow(max(normal.DotProduct(reflected_dir), 0.0), current_material->shininess);
			Vector3 ambient = current_material->ambient;

			
			if (shadow == 1.0f)
			{
				final_color.r += diffuse.x * factor + specular.x * spec_;
				final_color.g += diffuse.y * factor + specular.y * spec_;
				final_color.b += diffuse.z * factor + specular.z * spec_;
			}
			else
			{
				final_color.r += diffuse.x * factor;
				final_color.g += diffuse.y * factor;
				final_color.b += diffuse.z * factor;
			}

			//float factor = normal.DotProduct(light_dir);


			//final_color.r = diffuse.x * factor + specular.x * spec_;
			//final_color.g = diffuse.y * factor + specular.y * spec_;
			//final_color.b = diffuse.z * factor + specular.z * spec_;
			float reflectivity = current_material->reflectivity;
			//reflectivity = 0.97;
			
			final_color.r += reflected_color.r * (1 - reflectivity);
			final_color.g += reflected_color.g * (1 - reflectivity);
			final_color.b += reflected_color.b * (1 - reflectivity);

			final_color.r = final_color.r * shadow;
			final_color.g = final_color.g * shadow;
			final_color.b = final_color.b * shadow;

			add_color(final_color, { 0.03, 0.03, 0.03, 1 });

			return final_color;
			
		}
		
	}
	// background
	Vector3 dir{ -ray_hit.ray.dir_x, -ray_hit.ray.dir_y, ray_hit.ray.dir_z };
	const float theta = acos(dir.z);
	const float phi = atan2f(dir.y, dir.x) + M_PI;
	const float u = 1.0f - phi * 0.5f * M_1_PI;
	const float v = theta * M_1_PI;
	Color3f bg = background->get_bilinear_texel(u, v);
	//bg.r = gamma_util::c_linear(bg.r);
	//bg.g = gamma_util::c_linear(bg.g);
	//bg.b = gamma_util::c_linear(bg.b);
	return Color4f{ bg.r, bg.g, bg.b, 1.0f };

	return Color4f{ 1, 1, 1, 1 };

}

// hemisphere sampling
inline Vector3 sample_hemisphere(Vector3 normal)
{
	float xi1 = Random::Float();
	float xi2 = Random::Float();

	float x = 2 * cos(2 * M_PI * xi1) * sqrt(xi2 * (1 - xi2));
	float y = 2 * sin(2 * M_PI * xi1) * sqrt(xi2 * (1 - xi2));
	float z = 1 - (2 * xi2);

	Vector3 result{x,y,z};

	if (result.DotProduct(normal) < 0)
		result *= -1;
	
	return result;
}

Color4f Raytracer::shader_BRDF(RTCRayHit ray_hit, float depth, float ior)
{
	if (depth > 10) 
		return Color4f{ 0.0f, 0.0f, 0.0f, 1.0f };

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		current_material = surfaces_[ray_hit.hit.geomID]->get_material();

		//if (current_material->get_name() == "white_lambert" || 
		//	current_material->get_name() == "red_lambert" || 
		//	current_material->get_name() == "green_lambert")
		{
			RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
			Vector3 omega_o = { -ray_hit.ray.dir_x, -ray_hit.ray.dir_y, -ray_hit.ray.dir_z };
			Vector3 position{};
			position.x = ray_hit.ray.org_x * ray_hit.ray.tfar;
			position.y = ray_hit.ray.org_y * ray_hit.ray.tfar;
			position.z = ray_hit.ray.org_z * ray_hit.ray.tfar;

			Vector3 normal{};
			rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
				RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
			normal.Normalize();
			if ((normal.x * ray_hit.ray.dir_x + normal.y * ray_hit.ray.dir_y + normal.z * ray_hit.ray.dir_z) > 0.0f) {
				normal.x *= -1;
				normal.y *= -1;
				normal.z *= -1;
			}
			if (normal.DotProduct(omega_o) < 0)
				normal = -normal;
			
			// BRDF
			Vector3 emission = current_material->emission;
			Vector3 albedo = current_material->diffuse;
			Color4f L_e = { emission.x, emission.y, emission.z, 1.0f }; //get_emission(p, omega_o);
			
			if (L_e.r != 0.0f && L_e.g != 0.0f && L_e.b != 0.0f)
				return L_e;
			Vector3 omega_i = sample_hemisphere(normal);
			//omega_i.Normalize();
			RTCRay ray = generate_ray(position, omega_i);
			RTCRayHit reflected_ray_hit = generate_ray_hit(ray);

			float pdf = 2 * M_PI;//1 / 2 * M_PI;
			Color4f L_i = shader_BRDF(reflected_ray_hit, depth + 1);
			Color4f f_r{ albedo.x / M_PI, albedo.y / M_PI, albedo.z / M_PI, 1.0f };
			
			//Color4f L_r = L_i * f_r * (omega_i * normal) / pdf;
			//Color4f L_r = multiply_color(multiply_color(multiply_color(L_i, f_r), normal.DotProduct(omega_i)), pdf);
			Color4f L_r{};
			L_r.r = (L_i.r * f_r.r * omega_i.DotProduct(normal) * pdf);
			L_r.g = (L_i.g * f_r.g * omega_i.DotProduct(normal) * pdf);
			L_r.b = (L_i.b * f_r.b * omega_i.DotProduct(normal) * pdf);
			L_r.a = 1.0f;
			return L_r;
		}
		//return Color4f{ 1.0f, 0.0f, 1.0f };
	}
	else // return background
	{
		return Color4f{ 0, 0, 0, 1 };
		//return Color4f{ 0.5f, 0.5f, 0.5f, 1.0f };
		// background
		Vector3 dir{ -ray_hit.ray.dir_x, -ray_hit.ray.dir_y, ray_hit.ray.dir_z };
		const float theta = acos(dir.z);
		const float phi = atan2f(dir.y, dir.x) + M_PI;
		const float u = 1.0f - phi * 0.5f * M_1_PI;
		const float v = theta * M_1_PI;
		Color3f bg = background->get_texel(u, v);

		Color4f result{ bg.r, bg.g, bg.b, 1.0f };
		return result;
	}

	return Color4f{1.0f, 1.0f, 1.0f, 1.0f};
}

// 640, 480
Color4f Raytracer::get_pixel(const int x, const int y, const float t)
{

	Color4f final_color{ 0, 0, 0, 1 };
	
	// UNIFORM SAMPLING
	//float samples_in_direction = 2;
	//for (int y_dir = 0; y_dir < samples_in_direction; y_dir++)
	//	for (int x_dir = 0; x_dir < samples_in_direction; x_dir++)
	//	{
	//		float x_step = x_dir / samples_in_direction;
	//		float y_step = y_dir / samples_in_direction;
	//
	//		RTCRay ray;
	//		ray = camera_.GenerateRay(x+x_step, y+y_dir);
	//		RTCRayHit ray_hit = generate_ray_hit(ray);
	//		Color4f tmp = shader(ray_hit, 0);
	//		final_color.r += tmp.r;
	//		final_color.g += tmp.g;
	//		final_color.b += tmp.b;
	//	}
#if DOF
	float fr;
	float aperture = 1;

	int depth_of_field = 1;
	int sample_amount = 10;

	for (int i = 0; i < sample_amount; i++)
	{
		float x_bonus = Random::Float();
		float y_bonus = Random::Float();

		RTCRay ray;
		ray = camera_.GenerateRay(x + x_bonus, y + y_bonus);
		RTCRayHit org = generate_ray_hit(ray);
		RTCRayHit ray_hit;

		Color4f final_color_{0, 0, 0, 1};
		
		for (int j = 0; j < depth_of_field; j++)
		{
			Vector3 shift{ Random::Float() * 2 - 1, Random::Float() * 2 - 1, 0 };
			shift.Normalize();

			// shift the camera
			Vector3 new_camera = camera_.view_from_;
			new_camera.x = 0;
			new_camera.y = shift.x * aperture * 0.2f;
			new_camera.z = shift.y * aperture * 0.2f;

			//new_camera = camera_.M_c_w_ * new_camera;

			// hit point
			Vector3 p{ 0, 0, 0 };
			// check -> mb generovat novy dir podle nove kamery?
			//p.x = org.hit.Ng_x;
			//p.x = org.hit.Ng_y;
			//p.x = org.hit.Ng_z;

			p.x = org.ray.dir_x * org.ray.tfar;
			p.y = org.ray.dir_y * org.ray.tfar;
			p.z = org.ray.dir_z * org.ray.tfar;
			fr = org.ray.tfar;

			//float new_focal = fr * tan(camera_.fov_y_ / 2);
			//camera_.f_y_ = new_focal;

			Vector3 direction{ 0,0,0 };
			direction.x = p.x - new_camera.x;
			direction.y = p.y - new_camera.y;
			direction.z = p.z - new_camera.z;

			//direction.x = p.x - (org.ray.org_x + new_camera.x);
			//direction.y = p.y - (org.ray.org_y + new_camera.y);
			//direction.z = p.z - (org.ray.org_z + new_camera.z);


			//direction.x = p.x - org.ray.org_x;
			//direction.y = p.y - org.ray.org_y;
			//direction.z = p.z - org.ray.org_z;

			direction.Normalize();

			ray = camera_.GenerateRay(x + x_bonus, y + y_bonus);
			ray.org_x = org.ray.org_x + new_camera.x;
			ray.org_y = org.ray.org_y + new_camera.y;
			ray.org_z = org.ray.org_z + new_camera.z;

			ray.dir_x = direction.x;
			ray.dir_y = direction.y;
			ray.dir_z = direction.z;


			//Camera camera = Camera(640, 480, camera_.fov_y_, new_camera, Vector3(0, 0, 35));
			//ray = camera.GenerateRay(x + x_bonus, y + y_bonus);
			//
			//ray.org_x = new_camera.x;
			//ray.org_y = new_camera.y;
			//ray.org_z = new_camera.z;
			//
			//ray.dir_x = direction.x;
			//ray.dir_y = direction.y;
			//ray.dir_z = direction.z;


			ray_hit = generate_ray_hit(ray);


			Color4f tmp = shader(ray_hit, 0);
			final_color_.r += tmp.r;
			final_color_.g += tmp.g;
			final_color_.b += tmp.b;
		}

		final_color_.r = final_color_.r / depth_of_field;
		final_color_.g = final_color_.g / depth_of_field;
		final_color_.b = final_color_.b / depth_of_field;
		final_color_.a = 1.0f;

		final_color.r += final_color_.r;
		final_color.g += final_color_.g;
		final_color.b += final_color_.b;

	}

	final_color.r = final_color.r / depth_of_field;
	final_color.g = final_color.g / depth_of_field;
	final_color.b = final_color.b / depth_of_field;
	final_color.r = gamma_util::c_srgb(final_color.r);
	final_color.g = gamma_util::c_srgb(final_color.g);
	final_color.b = gamma_util::c_srgb(final_color.b);
	final_color.a = 1.0f;

	return final_color;


#else
	
	int sample_amount = 10;
	for (int i = 0; i < sample_amount; i++)
	{
		float x_bonus = Random::Float();
		float y_bonus = Random::Float();

		RTCRay ray;
		ray = camera_.GenerateRay(x + x_bonus, y + y_bonus);
		RTCRayHit ray_hit = generate_ray_hit(ray);
		
		Color4f tmp = shader(ray_hit, 0);
		final_color.r += tmp.r;
		final_color.g += tmp.g;
		final_color.b += tmp.b;
	}

	final_color.r = final_color.r / sample_amount;
	final_color.g = final_color.g / sample_amount;
	final_color.b = final_color.b / sample_amount;
	final_color.a = 1.0f;
	final_color.r = gamma_util::c_srgb(final_color.r);
	final_color.g = gamma_util::c_srgb(final_color.g);
	final_color.b = gamma_util::c_srgb(final_color.b);

	return final_color;
#endif

	Color4f sampled{0, 0, 0, 1};
	float number_of_samples = 20;
	for (int i = 0; i < number_of_samples; i++)
	{
		RTCRay ray;
		ray = camera_.GenerateRay(x, y);
		RTCRayHit ray_hit = generate_ray_hit(ray);
		Color4f tmp = shader_BRDF(ray_hit, 0);
		sampled.r += tmp.r;
		sampled.g += tmp.g;
		sampled.b += tmp.b;
	}
	sampled.r = sampled.r / number_of_samples;
	sampled.g = sampled.g / number_of_samples;
	sampled.b = sampled.b / number_of_samples;

	return sampled;
}

int Raytracer::Ui()
{
	static float f = 0.0f;
	static int counter = 0;

	// Use a Begin/End pair to created a named window
	ImGui::Begin("Ray Tracer Params");

	ImGui::Text("Surfaces = %d", surfaces_.size());
	ImGui::Text("Materials = %d", materials_.size());
	ImGui::Separator();
	ImGui::Checkbox("Vsync", &vsync_);

	//ImGui::Checkbox( "Demo Window", &show_demo_window ); // Edit bools storing our window open/close state
	//ImGui::Checkbox( "Another Window", &show_another_window );

	ImGui::SliderFloat("float", &f, 0.0f, 1.0f); // Edit 1 float using a slider from 0.0f to 1.0f    
	//ImGui::ColorEdit3( "clear color", ( float* )&clear_color ); // Edit 3 floats representing a color

	// Buttons return true when clicked (most widgets return true when edited/activated)
	if (ImGui::Button("Button"))
		counter++;
	ImGui::SameLine();
	ImGui::Text("counter = %d", counter);

	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
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