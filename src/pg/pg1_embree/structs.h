#pragma once

struct Vertex3f { float x, y, z; }; // a single vertex position structure matching certain format

using Normal3f = Vertex3f; // a single vertex normal structure matching certain format

struct Coord2f { float u, v; }; // texture coord structure

struct Triangle3ui { unsigned int v0, v1, v2; }; // indicies of a single triangle, the struct must match certain format, e.g. RTC_FORMAT_UINT3

struct RTC_ALIGN( 16 ) Color4f
{
	struct { float r, g, b, a; }; // a = 1 means that the pixel is opaque
};

struct Color3f { 
	float r, g, b; 
	
	Color3f operator*(const float val)
	{
		return Color3f{ this->r * val, this->g * val, this->b * val };
	}

	Color3f operator/(const float val)
	{
		return Color3f{ this->r / val, this->g / val, this->b / val };
	}

	Color3f operator+(const Color3f& val)
	{
		return Color3f{ this->r + val.r, this->g + val.g, this->b + val.b };
	}

};
