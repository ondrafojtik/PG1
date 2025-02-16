#include "stdafx.h"
#include "tutorials.h"
#include "raytracer.h"
#include "structs.h"
#include "texture.h"
#include "mymath.h"

#include <OpenImageDenoise/oidn.h>

/* error reporting function */
void error_handler( void * user_ptr, const RTCError code, const char * str )
{
	if ( code != RTC_ERROR_NONE )
	{
		std::string descr = str ? ": " + std::string( str ) : "";

		switch ( code )
		{
		case RTC_ERROR_UNKNOWN: throw std::runtime_error( "RTC_ERROR_UNKNOWN" + descr );
		case RTC_ERROR_INVALID_ARGUMENT: throw std::runtime_error( "RTC_ERROR_INVALID_ARGUMENT" + descr ); break;
		case RTC_ERROR_INVALID_OPERATION: throw std::runtime_error( "RTC_ERROR_INVALID_OPERATION" + descr ); break;
		case RTC_ERROR_OUT_OF_MEMORY: throw std::runtime_error( "RTC_ERROR_OUT_OF_MEMORY" + descr ); break;
		case RTC_ERROR_UNSUPPORTED_CPU: throw std::runtime_error( "RTC_ERROR_UNSUPPORTED_CPU" + descr ); break;
		case RTC_ERROR_CANCELLED: throw std::runtime_error( "RTC_ERROR_CANCELLED" + descr ); break;
		default: throw std::runtime_error( "invalid error code" + descr ); break;
		}
	}
}

/* adds a single triangle to the scene */
unsigned int add_triangle( const RTCDevice device, RTCScene scene )
{
	// geometries are objects that represent an array of primitives of the same type, so lets create a triangle
	RTCGeometry mesh = rtcNewGeometry( device, RTC_GEOMETRY_TYPE_TRIANGLE );

	// and depending on the geometry type, different buffers must be bound (typically, vertex and index buffer is required)

	// set vertices in the newly created buffer
	Vertex3f * vertices = ( Vertex3f * )rtcSetNewGeometryBuffer(
		mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof( Vertex3f ), 3 );
	vertices[0].x = 0; vertices[0].y = 0; vertices[0].z = 0;
	vertices[1].x = 2; vertices[1].y = 0; vertices[1].z = 0;
	vertices[2].x = 0; vertices[2].y = 3; vertices[2].z = 0;

	// set triangle indices
	Triangle3ui * triangles = ( Triangle3ui * )rtcSetNewGeometryBuffer( mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof( Triangle3ui ), 1 );
	triangles[0].v0 = 0; triangles[0].v1 = 1; triangles[0].v2 = 2;

	// see also rtcSetSharedGeometryBuffer, rtcSetGeometryBuffer		

	/*
	The parametrization of a triangle uses the first vertex p0 as base point, the vector (p1 - p0) as u-direction and the vector (p2 - p0) as v-direction.
	Thus vertex attributes t0, t1, t2 can be linearly interpolated over the triangle the following way:

	t_uv = (1-u-v)*t0 + u*t1 + v*t2	= t0 + u*(t1-t0) + v*(t2-t0)
	*/

	// sets the number of slots (vertexAttributeCount parameter) for vertex attribute buffers (RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE)
	rtcSetGeometryVertexAttributeCount( mesh, 2 );

	// set vertex normals
	Normal3f * normals = ( Normal3f * )rtcSetNewGeometryBuffer( mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3, sizeof( Normal3f ), 3 );
	normals[0].x = 0; normals[0].y = 0; normals[0].z = 1;
	normals[1].x = 0; normals[1].y = 0; normals[1].z = 1;
	normals[2].x = 0; normals[2].y = 0; normals[2].z = 1;

	// set texture coordinates
	Coord2f * tex_coords = ( Coord2f * )rtcSetNewGeometryBuffer( mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2, sizeof( Coord2f ), 3 );
	tex_coords[0].u = 0; tex_coords[0].v = 1;
	tex_coords[1].u = 1; tex_coords[1].v = 1;
	tex_coords[2].u = 0; tex_coords[2].v = 0;

	// changes to the geometry must be always committed
	rtcCommitGeometry( mesh );

	// geometries can be attached to a single scene
	unsigned int geom_id = rtcAttachGeometry( scene, mesh );
	// release geometry handle
	rtcReleaseGeometry( mesh );

	return geom_id;
}

/* generate a single ray and get the closest intersection with the scene */
int generate_and_trace_ray( RTCScene & scene )
{
	// setup a primary ray
	RTCRay ray;
	ray.org_x = 0.1f; // ray origin
	ray.org_y = 0.2f;
	ray.org_z = 2.0f;
	ray.tnear = FLT_MIN; // start of ray segment

	ray.dir_x = 0.0f; // ray direction
	ray.dir_y = 0.0f;
	ray.dir_z = -1.0f;
	ray.time = 0.0f; // time of this ray for motion blur

	ray.tfar = FLT_MAX; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	// setup a hit
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
	rtcInitIntersectContext( &context );
	rtcIntersect1( scene, &context, &ray_hit );

	if ( ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID )
	{
		// we hit something
		RTCGeometry geometry = rtcGetGeometry( scene, ray_hit.hit.geomID );
		Normal3f normal;
		// get interpolated normal
		rtcInterpolate0( geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3 );
		// and texture coordinates
		Coord2f tex_coord;
		rtcInterpolate0( geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord.u, 2 );

		printf( "normal = (%0.3f, %0.3f, %0.3f)\n", normal.x, normal.y, normal.z );
		printf( "tex_coord = (%0.3f, %0.3f)\n", tex_coord.u, tex_coord.v );
	}

	return EXIT_SUCCESS;
}

/* simple tutorial how to find the closest intersection with a single triangle using embree */
int tutorial_1( const char * config )
{
	RTCDevice device = rtcNewDevice( config );
	error_handler( nullptr, rtcGetDeviceError( device ), "Unable to create a new device.\n" );
	rtcSetDeviceErrorFunction( device, error_handler, nullptr );

	ssize_t triangle_supported = rtcGetDeviceProperty( device, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED );

	// create a new scene bound to the specified device
	RTCScene scene = rtcNewScene( device );

	// add a single triangle geometry to the scene
	unsigned int triangle_geom_id = add_triangle( device, scene );

	// commit changes to scene
	rtcCommitScene( scene );

	generate_and_trace_ray( scene );

	// release scene and detach with all geometries
	rtcReleaseScene( scene );

	rtcReleaseDevice( device );

	return EXIT_SUCCESS;
}

std::vector<float> to_float_buffer(Texture* albedo, int normal)
{
	std::vector<float> color_buffer;
	for (int i = 0; i < albedo->height(); i++)
	{
		for (int j = 0; i < albedo->width(); j++)
			for (int k = 0; k < 3; k++)
			{
				float value = albedo->data_[(i * albedo->width() + j) * 3 + k] / 255.0f;
				if (normal == 1)
					value = (value * 2) - 1;

				color_buffer.push_back(value);
			}
	}
	return color_buffer;
}

/* texture loading and texel access */
int tutorial_2()
{
	// create texture
	//Texture texture( "../../../data/test4.png" );
	//Color3f texel = texture.get_texel( ( 1.0f / texture.width() ) * 2.5f, 0.0f );
	//printf( "(r = %0.3f, g = %0.3f, b = %0.3f)\n", texel.r, texel.g, texel.b );

	// denoise
	Texture albedo("../../../data/denoise/albedo_100spp.png");
	Texture normal("../../../data/denoise/normal_100spp.png");
	Texture color("../../../data/denoise/color_100spp.png");

	auto albedo_ = to_float_buffer(&albedo, 0);
	auto normal_ = to_float_buffer(&normal, 1);
	auto color_ = to_float_buffer(&color, 0);

	//float* outputPtr = new float(albedo.width() * albedo.height()*3);
	std::vector<float> outputPtr;

#if 1
	// Create an Intel Open Image Denoise device
	OIDNDevice device = oidnNewDevice(OIDN_DEVICE_TYPE_DEFAULT);
	oidnCommitDevice(device);

	// Create a filter for denoising a beauty (color) image using optional auxiliary images too
	OIDNFilter filter = oidnNewFilter(device, "RT"); // generic ray tracing filter
	oidnSetSharedFilterImage(filter, "color", color_.data(),
		OIDN_FORMAT_FLOAT3, color.width(), color.height(), 0, 0, 0); // beauty
	oidnSetSharedFilterImage(filter, "albedo", albedo_.data(),
		OIDN_FORMAT_FLOAT3, albedo.width(), albedo.height(), 0, 0, 0); // auxiliary
	oidnSetSharedFilterImage(filter, "normal", normal_.data(),
		OIDN_FORMAT_FLOAT3, normal.width(), normal.height(), 0, 0, 0); // auxiliary
	oidnSetSharedFilterImage(filter, "output", outputPtr.data(),
		OIDN_FORMAT_FLOAT3, albedo.width(), albedo.height(), 0, 0, 0); // denoised beauty
	oidnSetFilter1b(filter, "hdr", false); // beauty image is HDR
	oidnCommitFilter(filter);

	// Filter the image
	oidnExecuteFilter(filter);

	// Check for errors
	const char* errorMessage;
	if (oidnGetDeviceError(device, &errorMessage) != OIDN_ERROR_NONE)
		printf("Error: %s\n", errorMessage);

	// Cleanup
	oidnReleaseFilter(filter);
	oidnReleaseDevice(device);

	FILE* file = fopen("../../../data/denoise/test.png", "wt");
	fprintf(file, "PF\n");
	fprintf(file, "%d %d\n", albedo.width(), albedo.height());
	fprintf(file, "-1\n");
	file = fopen("../../../data/denoise/test.png", "a+b");
	fwrite(outputPtr.data(), sizeof(float) * 3, albedo.width() * albedo.height(), file);
//	fclose(file);

#endif

	return EXIT_SUCCESS;
}

/* raytracer mainloop */
int tutorial_3( const std::string file_name, const char * config )
{
	//SimpleGuiDX11 gui( 640, 480 );
	//gui.MainLoop();

	//original!
#if LEGO_SPACESHIP
	Raytracer raytracer( 640, 480, deg2rad( 45.0 ),
		Vector3(110, -130, 120), Vector3( 0, 0, 35 ), config );
	// 110, -130, 120
	// 130, -100, 90

	
	// 100, -100, 75
	// 130, -100, 90	-- fajn
	// 100, -150, 120	-- best ?
	// 190, -120, 140
#endif

#if REFRACTION_SPHERE
	Raytracer raytracer(640, 480, deg2rad(45.0),
		Vector3(0, -3.0f, 0), Vector3(0, 0, 0), config);
#endif

#if CORNELL_BOX
	Raytracer raytracer(640, 480, deg2rad(45.0),
		Vector3(0, -650, 250), Vector3(0, 0, 250), config);
#endif

	raytracer.LoadScene( file_name );
	raytracer.MainLoop();

	return EXIT_SUCCESS;
}
