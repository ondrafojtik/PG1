#pragma once

#define LEGO_SPACESHIP 1
#define REFRACTION_SPHERE 0
#define DOF 0



#include "targetver.h"

#define _CRT_SECURE_NO_WARNINGS

#define _USE_MATH_DEFINES // for C++
#include <cmath>

// std libs
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <chrono>
#include <mutex>
#include <thread>
#include <atomic>
#include <tchar.h>
#include <vector>
#include <map>
#include <random>
#include <functional>
#include <iostream>

#include "Random.h"
#include <xmmintrin.h>
#include <pmmintrin.h>

// Intel Embree 3.6.1 (high performance ray tracing kernels)
#include <embree3/rtcore.h>

// Dear ImGui (bloat-free graphical user interface library for C++)
#include <imgui.h>
#include "imgui_impl_win32.h"
#include "imgui_impl_dx11.h"
#include <d3d11.h>
#define DIRECTINPUT_VERSION 0x0800
#include <dinput.h>



#define CORNELL_BOX 0

