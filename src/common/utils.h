#pragma once

#include <chrono>
#include <iostream>
#include <vector>
#include <sstream>

#include <refresh/compression/lib/zstd_wrapper.h>

using namespace std;
using namespace std::chrono;

#define COMPACT_ENCODING

void show_timings(vector<pair<string, time_point<high_resolution_clock>>>& timings);

constexpr uint32_t no_bytes(uint32_t x)
{
	if (x < 256)
		return 1;
	if (x < 256 * 256)
		return 2;
	if (x < 256 * 256 * 256)
		return 3;
	return 4;
}

template<typename T>
void save_int_lsb(FILE *f, T x, int n_bytes)
{
	switch (n_bytes)
	{
	case 8: putc(x & 0xff, f);		x >>= 8;
	case 7: putc(x & 0xff, f);		x >>= 8;
	case 6: putc(x & 0xff, f);		x >>= 8;
	case 5: putc(x & 0xff, f);		x >>= 8;
	case 4: putc(x & 0xff, f);		x >>= 8;
	case 3: putc(x & 0xff, f);		x >>= 8;
	case 2: putc(x & 0xff, f);		x >>= 8;
	case 1: putc(x & 0xff, f);
	}
}

template<typename T>
void save_int_lsb(refresh::zstd_file &f, T x, int n_bytes)
{
	switch (n_bytes)
	{
	case 8: f.put(x & 0xff);		x >>= 8;
	case 7: f.put(x & 0xff);		x >>= 8;
	case 6: f.put(x & 0xff);		x >>= 8;
	case 5: f.put(x & 0xff);		x >>= 8;
	case 4: f.put(x & 0xff);		x >>= 8;
	case 3: f.put(x & 0xff);		x >>= 8;
	case 2: f.put(x & 0xff);		x >>= 8;
	case 1: f.put(x & 0xff);
	}
}

/*template<typename T>
void append_int_msb(vector<uint8_t>& v, T x, int n_bytes)
{
	switch (n_bytes)
	{
	case 8: v.emplace_back((x >> 56) & 0xff);
	case 7: v.emplace_back((x >> 48) & 0xff);
	case 6: v.emplace_back((x >> 40) & 0xff);
	case 5: v.emplace_back((x >> 32) & 0xff);
	case 4: v.emplace_back((x >> 24) & 0xff);
	case 3: v.emplace_back((x >> 16) & 0xff);
	case 2: v.emplace_back((x >> 8) & 0xff);
	case 1: v.emplace_back(x & 0xff);
	}
}*/

void encode_shared_prefix(vector<uint8_t>& v_data, vector<uint8_t>& v_prev, vector<uint8_t>& v_curr);


// EOF