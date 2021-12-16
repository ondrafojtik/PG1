#include "stdafx.h"
#include "texture.h"

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

Texture::Texture( const char * file_name )
{
	// image format
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
	// pointer to the image, once loaded
	FIBITMAP * dib =  nullptr;
	// pointer to the image data
	BYTE * bits = nullptr;

	// check the file signature and deduce its format
	fif = FreeImage_GetFileType( file_name, 0 );
	// if still unknown, try to guess the file format from the file extension
	if ( fif == FIF_UNKNOWN )
	{
		fif = FreeImage_GetFIFFromFilename( file_name );
	}
	// if known
	if ( fif != FIF_UNKNOWN )
	{
		// check that the plugin has reading capabilities and load the file
		if ( FreeImage_FIFSupportsReading( fif ) )
		{
			dib = FreeImage_Load( fif, file_name );
		}
		// if the image loaded
		if ( dib )
		{
			// retrieve the image data
			bits = FreeImage_GetBits( dib );
			//FreeImage_ConvertToRawBits()
			// get the image width and height
			width_ = int( FreeImage_GetWidth( dib ) );
			height_ = int( FreeImage_GetHeight( dib ) );

			// if each of these is ok
			if ( ( bits != 0 ) && ( width_ != 0 ) && ( height_ != 0 ) )
			{				
				// texture loaded
				scan_width_ = FreeImage_GetPitch( dib ); // in bytes
				pixel_size_ = FreeImage_GetBPP( dib ) / 8; // in bytes

				data_ = new BYTE[scan_width_ * height_]; // BGR(A) format
				FreeImage_ConvertToRawBits( data_, dib, scan_width_, pixel_size_ * 8, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE );
			}

			FreeImage_Unload( dib );
			bits = nullptr;
		}
	}	
}

Texture::~Texture()
{	
	if ( data_ )
	{
		// free FreeImage's copy of the data
		delete[] data_;
		data_ = nullptr;
		
		width_ = 0;
		height_ = 0;
	}
}

Color3f Texture::get_texel( const float u, const float v ) const
{
	//assert( ( u >= 0.0f && u <= 1.0f ) && ( v >= 0.0f && v <= 1.0f ) );	
	
	const int x = max( 0, min( width_ - 1, int( u * width_ ) ) );
	const int y = max( 0, min( height_ - 1, int( v * height_ ) ) );

	const int offset = y * scan_width_ + x * pixel_size_;
	const float b = data_[offset] / 255.0f;
	const float g = data_[offset + 1] / 255.0f;
	const float r = data_[offset + 2] / 255.0f;
	
	return Color3f{ r, g, b };
}

Color3f Texture::get_bilinear_texel(const float u, const float v) const
{
	int x = max(0, min(width_ - 1, floor(u * width_)));
	int y = max(0, min(height_ - 1, floor(v * height_)));

	int x0 = x;
	int y0 = y;

	int x1 = min(width_ - 1, x0 + 1);
	int y1 = min(height_ - 1, y0 + 1);

	int p1_o = y0 * scan_width_ + x0 * pixel_size_;
	int p2_o = y0 * scan_width_ + x1 * pixel_size_;
	int p3_o = y1 * scan_width_ + x1 * pixel_size_;
	int p4_o = y1 * scan_width_ + x0 * pixel_size_;

	float x__ = x - x0;
	float y__ = y - y0;


	Color3f result = 
		Color3f(
		(Color3f{ (float)(data_[p1_o + 2]), (float)(data_[p1_o + 1]), (float)(data_[p1_o]) } * (1 - x__)	* (1 - y__)) +
		(Color3f{ (float)(data_[p2_o + 2]), (float)(data_[p2_o + 1]), (float)(data_[p2_o]) } * x__			* (1 - y__)) +
		(Color3f{ (float)(data_[p3_o + 2]), (float)(data_[p3_o + 1]), (float)(data_[p3_o]) } * x__			* y__) +
		(Color3f{ (float)(data_[p4_o + 2]), (float)(data_[p4_o + 1]), (float)(data_[p4_o]) } * (1 - x__)	* y__));

	result = result / 255;
	result.r = c_linear(result.r);
	result.g = c_linear(result.g);
	result.b = c_linear(result.b);

	return result;
}

int Texture::width() const
{
	return width_;
}

int Texture::height() const
{
	return height_;
}
