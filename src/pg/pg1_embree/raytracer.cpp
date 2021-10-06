#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"

#define SPECULAR_STRENGTH 0.5

Raytracer::Raytracer( const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config ) : SimpleGuiDX11( width, height )
{
	InitDeviceAndScene( config );

	camera_ = Camera( width, height, fov_y, view_from, view_at );
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
			Triangle & triangle = surface->get_triangle( i );

			// vertices loop
			for ( int j = 0; j < 3; ++j, ++k )
			{
				const Vertex & vertex = triangle.vertex( j );

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

		rtcCommitGeometry( mesh );
		unsigned int geom_id = rtcAttachGeometry( scene_, mesh );
		rtcReleaseGeometry( mesh );
	} // end of surfaces loop

	rtcCommitScene( scene_ );
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


	return Vector3{diff.r, diff.g, diff.b};
}

Vector3 Raytracer::calc_normal(RTCRayHit ray_hit) 
{
	RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);

	Vector3 normal_{};
	rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
		RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal_.x, 3);
	normal_.Normalize();
	return normal_;
}

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
	Vertex A = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(0);
	Vertex B = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(1);
	Vertex C = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(2);


	// TODO(ondra): hit points are _inside_ the faces..
	Vector3 light_dir = light_position - frag_pos;
	light_dir.Normalize();
	float value = calc_normal(ray_hit).DotProduct(light_dir);
	if (value >= 0)
		frag_pos += calc_normal(ray_hit);
	else
		frag_pos += calc_normal(ray_hit).Abs();

	
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
	// TODO(ondra): tohle je duvod proc chci abs, ne max..
	//diffuse_strength = max(diffuse_strength, 0.0f);
	float diffuse_strength = light_dir.DotProduct(normal);
	diffuse_strength = abs(diffuse_strength);

	diff.x = diff.x * diffuse_strength;
	diff.y = diff.y * diffuse_strength;
	diff.z = diff.z * diffuse_strength;

	// reflecton
	Vector3 reflect_dir = reflect(camera_.view_direction, normal);
	Vector3 halfway_dir = camera_.view_direction + light_dir;
	halfway_dir.Normalize();
	float spec_ = pow(max(normal.DotProduct(halfway_dir), 0.0), current_material->shininess);
	//float spec_ = pow(max(normal_.DotProduct(reflect_dir), 0.0), current_material->shininess);
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
	ray.tnear = FLT_MIN; // start of ray segment
	
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

// 640, 480
Color4f Raytracer::get_pixel( const int x, const int y, const float t )
{
	// TODO generate primary ray and perform ray cast on the scene
	RTCRay ray = camera_.GenerateRay(x, y);
	
	RTCRayHit ray_hit = generate_ray_hit(ray);

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		Color4f blinn_phong = calc_blinn_phong(ray_hit);
		Vector3 frag_pos = get_fragment_position(ray_hit);
		bool test = generate_shadow_ray(frag_pos, light_position);
		
		Color4f final_color{0,0,0,1};
		
		final_color.r = blinn_phong.r;
		final_color.g = blinn_phong.g;
		final_color.b = blinn_phong.b;

		float shadow_factor = 0.35f;
		if (test)
			return Color4f
			{
				final_color.r * shadow_factor, 
				final_color.g * shadow_factor, 
				final_color.b * shadow_factor, 1
			};


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