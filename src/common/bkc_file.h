#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include <cinttypes>
#include <thread>
#include <mutex>

using namespace std;

#include <types/satc_data.h>
#include "defs.h"

struct bkc_record_t {
	uint64_t sample_id;
	uint64_t barcode;
	uint64_t leader;
	uint64_t follower;
	uint64_t count;
	
	bkc_record_t() :
		sample_id(0),
		barcode(0),
		leader(0),
		follower(0),
		count(0)
	{}

	bkc_record_t(uint64_t _sample_id, uint64_t _barcode, uint64_t _leader, uint64_t _follower, uint64_t _count) :
		sample_id(_sample_id),
		barcode(_barcode),
		leader(_leader),
		follower(_follower),
		count(_count) 
	{}

	bkc_record_t(const bkc_record_t&) = default;
	bkc_record_t(bkc_record_t&&) = default;

	bkc_record_t& operator=(const bkc_record_t&) = default;
	bkc_record_t& operator=(bkc_record_t&&) = default;
};

class CBKCFile
{
	enum class open_mode_t {none, reading, writing};
	open_mode_t open_mode{ open_mode_t::none };

	const int BUFFER_SIZE = 16 << 20;

	mutex mtx;

	refresh::zstd_file f{ 9 };

	uint8_t sample_id_size_in_bytes;
	uint8_t barcode_size_in_bytes;
	uint8_t leader_size_in_bytes;
	uint8_t follower_size_in_bytes;
	uint8_t counter_size_in_bytes;
	uint8_t barcode_len_in_symbols;
	uint8_t leader_len_in_symbols;
	uint8_t gap_len_in_symbols;
	uint8_t follower_len_in_symbols;
	Header::ordering_t ordering{ Header::default_ordering };

	uint32_t zstd_level;
	output_format_t output_format{ output_format_t::unknown };

	vector<uint8_t> rec_prev, rec_curr;
	uint32_t rec_len;

	void save_header();

//	void add_record(uint64_t sample_id, uint64_t barcode, uint64_t leader, uint64_t follower, uint64_t count);

public:
	CBKCFile() = default;
	~CBKCFile()
	{
		Close();
	};

	void SetParams(uint8_t _sample_id_size_in_bytes, uint8_t _barcode_size_in_bytes, uint8_t _leader_size_in_bytes, uint8_t _follower_size_in_bytes, uint8_t _counter_size_in_bytes,
		uint8_t _barcode_len_in_symbols, uint8_t _leader_len_in_symbols, uint8_t _gap_len_in_symbols, uint8_t _follower_len_in_symbols, uint32_t _zstd_level);

	bool Create(const string& file_name, const output_format_t _output_format);
	bool Open(const string& file_name);
	bool Close();

//	void AddRecord(uint64_t sample_id, uint64_t barcode, uint64_t leader, uint64_t follower, uint64_t count);
//	void AddRecords(vector<bkc_record_t> &records);
	void AddPacked(vector<uint8_t>& packed);
	bool GetRecord(uint64_t &sample_id, uint64_t &barcode, uint64_t &leader, uint64_t &follower, uint64_t &count);

	void GetLens(uint8_t& _barcode_len_in_symbols, uint8_t& _leader_len_in_symbols, uint8_t& _follower_len_in_symbols, uint8_t& _counter_size_in_bytes);
};

// EOF