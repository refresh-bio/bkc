#pragma once

#include <cstdio>
#include <string>
#include <zlib-ng/zlib.h>
#include <vector>
#include <array>

#include <refresh/memory_chunk/lib/memory_chunk.h>

using namespace std;
using namespace refresh;

// *********************************************************************************************
class CFastXReader
{
protected:
	const size_t BUFFER_SIZE = 64 << 20;
	const size_t ZLIB_BUFFER_SIZE = 32 << 20;

	FILE* in = nullptr;
	gzFile gz_in = nullptr;
	bool is_gzipped = false;
	bool is_fastq = true;
	int rec_lines;

	string file_name;

	array<int, 4> eol_positions;
	vector<char> internal_buffer;

	void close();

	bool is_gzipped_name(const string& fn)
	{
		return fn.size() > 3 && fn.substr(fn.size() - 3, 3) == ".gz";
	}

	void find_last_eols(memory_chunk<char>& mc);

public:
	CFastXReader(bool is_fastq) : 
		is_fastq(is_fastq) ,
		rec_lines(is_fastq ? 4 : 2)
	{}
	~CFastXReader() {
		close();
	}

	bool Open(const string& _file_name);
	void Close();

	bool ReadBlock(memory_chunk<char>& mc);
	bool Eof();
};

// *********************************************************************************************
struct read_desc_t
{
	memory_chunk<char>::iterator header;
	memory_chunk<char>::iterator bases;
	memory_chunk<char>::iterator plus;
	memory_chunk<char>::iterator quality;

	read_desc_t() :
		header(nullptr),
		bases(nullptr),
		plus(nullptr),
		quality(nullptr)
	{}

	read_desc_t(memory_chunk<char>::iterator _header, memory_chunk<char>::iterator _bases, memory_chunk<char>::iterator _plus, memory_chunk<char>::iterator _quality) :
		header(_header),
		bases(_bases),
		plus(_plus),
		quality(_quality)
	{}

	read_desc_t(const read_desc_t&) = default;
};

// *********************************************************************************************
class CReadReader
{
	memory_chunk<char> block;
	memory_chunk<char>::iterator p_start;

	bool is_eob;
	bool is_fastq;
	int rec_lines;
	char first_symbol;

	bool find_eol(memory_chunk<char>::iterator& iter);
	void skip_eols(memory_chunk<char>::iterator& iter);

public:
	CReadReader(bool is_fastq) : 
		is_eob(false),
		is_fastq(is_fastq),
		rec_lines(is_fastq ? 4 : 2),
		first_symbol(is_fastq ? '@' : '>')
	{};
	void Assign(memory_chunk<char>& _block);
	bool GetRead(read_desc_t &read_desc);
	bool Eob();
	bool ShrinkBlock();
};


// EOF
