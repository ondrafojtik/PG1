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


// 640, 480
Color4f Raytracer::get_pixel( const int x, const int y, const float t )
{
	// TODO generate primary ray and perform ray cast on the scene
	RTCRay ray = camera_.GenerateRay(x, y);
	
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

	
	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		current_material = surfaces_[ray_hit.hit.geomID]->get_material();

		Color3f diff{ 0, 0, 0 };
		Color3f spec{ 0, 0, 0 };
		Color3f normal{ 0, 0, 0 };
		Color3f oppacity{ 0, 0, 0 };


		if (current_material->get_texture(current_material->kDiffuseMapSlot) != nullptr)
		{
			Vertex A = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(0);
			Vertex B = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(1);
			Vertex C = surfaces_[ray_hit.hit.geomID]->get_triangle(0).vertex(2);

			float u = ray_hit.hit.u;
			float v = ray_hit.hit.v;
			float w = 1 - u - v;
			float u_ = (A.texture_coords[0].u * w) + (B.texture_coords[0].u * u) + (C.texture_coords[0].u * v);
			float v_ = (A.texture_coords[0].v * w) + (B.texture_coords[0].v * u) + (C.texture_coords[0].v * v);
		
			diff = current_material->get_texture(current_material->kDiffuseMapSlot)->get_texel(u_, v_);
			return Color4f{ diff.r, diff.g, diff.b, 1.0f };
			return Color4f{ 0.0f, 0.0f, 1.0f, 1.0f };
		}
		else
			diff = { current_material->diffuse.x, current_material->diffuse.y, current_material->diffuse.z };

		if (current_material->get_texture(current_material->kSpecularMapSlot) != nullptr)
			spec = current_material->get_texture(current_material->kSpecularMapSlot)->get_texel(ray_hit.hit.u, ray_hit.hit.v);
		else
			spec = { current_material->specular.x, current_material->specular.y, current_material->specular.z };


		Vector3 normal_{};
		normal_.x = ray_hit.hit.Ng_x;
		normal_.y = ray_hit.hit.Ng_y;
		normal_.z = ray_hit.hit.Ng_z;
		normal_.Normalize();
		normal = { normal_.x, normal_.y, normal_.z };

		Vector3 light_position{ -200, 150, -120 };
		Vector3 frag_pos
		{
			ray_hit.ray.org_x + (ray_hit.ray.dir_x * ray_hit.ray.tfar),
			ray_hit.ray.org_y + (ray_hit.ray.dir_y * ray_hit.ray.tfar),
			ray_hit.ray.org_z + (ray_hit.ray.dir_z * ray_hit.ray.tfar)
		};
		Vector3 light_direction
		{
			frag_pos.x - light_position.x,
			frag_pos.y - light_position.y,
			frag_pos.z - light_position.z
		};
		light_direction.Normalize();

		/*
		DIFFERENT WAY OF CALCULATING THE HIT POINT

		Vector3 point_hit
		{
			(1.0f - u - v) * A.position.x + u * B.position.x + v * C.position.x,
			(1.0f - u - v) * A.position.y + u * B.position.y + v * C.position.y,
			(1.0f - u - v) * A.position.z + u * B.position.z + v * C.position.z
		}
		
		*/


		// diffuse
		float diffuse_strength = light_direction.DotProduct(normal_);
		diffuse_strength = max(diffuse_strength, 0.0f);
		diff.r = diff.r * diffuse_strength;
		diff.g = diff.g * diffuse_strength;
		diff.b = diff.b * diffuse_strength;


		// reflecton
		Vector3 reflect_dir = camera_.view_direction - 2 * (camera_.view_direction.DotProduct(normal_) * normal_);
		Vector3 halfway_dir = camera_.view_direction + light_direction;
		halfway_dir.Normalize();
		float spec_ = pow(max(normal_.DotProduct(halfway_dir), 0.0), current_material->shininess);
		//float spec_ = pow(max(normal_.DotProduct(reflect_dir), 0.0), current_material->shininess);
		spec.r = SPECULAR_STRENGTH * spec_;
		spec.g = SPECULAR_STRENGTH * spec_;
		spec.b = SPECULAR_STRENGTH * spec_;

		//relection = direction - 2 (dot(direction, normal))normal;

		Color3f final_color{ 0, 0, 0 };
		final_color.r = diff.r + spec.r;
		final_color.g = diff.g + spec.g;
		final_color.b = diff.b + spec.b;


		return Color4f{final_color.r, final_color.g, final_color.b, 1.0f};
	}

	float x_ = (float)x;
	float y_ = (float)y;
	return Color4f{ x_/ 640.0f, 0.0f, y_ / 480.0f, 1.0f };
	
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
