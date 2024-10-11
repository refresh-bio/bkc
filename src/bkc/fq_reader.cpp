#include <iostream>
#include <cstring>
#include "fq_reader.h"

// *********************************************************************************************
void CFastXReader::close()
{
	if (in)
	{
		fclose(in);
		in = nullptr;
	}
	else if (gz_in)
	{
		gzclose(gz_in);
		gz_in = nullptr;
	}

	file_name.clear();

	is_gzipped = false;
}

// *********************************************************************************************
bool CFastXReader::Open(const string &_file_name)
{
	if (in || gz_in)
		close();

	if (is_gzipped_name(_file_name))
	{
		gz_in = gzopen(_file_name.c_str(), "rb");

		if (!gz_in)
			return false;

		gzbuffer(gz_in, ZLIB_BUFFER_SIZE);
		is_gzipped = true;
	}
	else
	{
		in = fopen(_file_name.c_str(), "rb");

		if (!in)
			return false;

		setvbuf(in, nullptr, _IOFBF, BUFFER_SIZE);
		is_gzipped = false;
	}

	file_name = _file_name;

	return true;
}

// *********************************************************************************************
void CFastXReader::Close()
{
	close();
}

// *********************************************************************************************
bool CFastXReader::ReadBlock(memory_chunk<char>& mc)
{
	if (!in && !gz_in)
		return false;

	mc.resize(internal_buffer.size());
	memcpy(mc.data(), internal_buffer.data(), internal_buffer.size());
	internal_buffer.clear();

	size_t to_read = mc.capacity() - mc.size();
	size_t readed;
	size_t filled = mc.size();

	mc.resize(mc.capacity());

	if (is_gzipped)
		readed = gzread(gz_in, mc.data() + filled, to_read);
	else
		readed = fread(mc.data() + filled, 1, to_read, in);

	mc.resize(filled + readed);

	find_last_eols(mc);

	if (eol_positions[rec_lines-1] < 0 && mc.size() != mc.capacity())
		return false;		// No file end but impossible to find any complete read!

	if (eol_positions[rec_lines - 1] + 1 != (int) mc.size())
		internal_buffer.assign(mc.data() + eol_positions[rec_lines - 1] + 1, mc.data() + mc.size());

	mc.resize(eol_positions[rec_lines - 1] + 1);

	return true;
}

// *********************************************************************************************
bool CFastXReader::Eof()
{
	if (in)
		return internal_buffer.empty() && feof(in);

	if (gz_in)
		return internal_buffer.empty() && gzeof(gz_in);

	return false;
}

// *********************************************************************************************
void CFastXReader::find_last_eols(memory_chunk<char>& mc)
{
	fill_n(eol_positions.begin(), rec_lines, -1);

	int eol_cnt = 0;

	char* p = mc.data();

	for (int i = 0; i < (int)mc.size(); ++i, ++p)
		if (*p == '\n')
			eol_positions[eol_cnt++ % rec_lines] = i;
}

// *********************************************************************************************
//
// *********************************************************************************************

// *********************************************************************************************
void CReadReader::Assign(memory_chunk<char>& _block) 
{ 
	block = _block.get_copy(); 
	is_eob = false;
	p_start = block.data();
}

// *********************************************************************************************
bool CReadReader::GetRead(read_desc_t& read_desc)
{
	if (is_eob)
		return false;

	if (p_start == block.end())
	{
		is_eob = true;
		return false;
	}

	if (*p_start != first_symbol)
		return false;

	auto p = p_start;
//	auto p_end = block.end();

	read_desc.header = p;

	if (!find_eol(p))
		return false;

	skip_eols(p);

	read_desc.bases = p;

	if (!find_eol(p))
		return false;

	skip_eols(p);

	if (is_fastq)
	{
		read_desc.plus = p;
		if (!find_eol(p))
			return false;

		skip_eols(p);

		read_desc.quality = p;

		if (!find_eol(p))
			return false;

		skip_eols(p);
	}

	p_start = p;

	return true;
}

// *********************************************************************************************
bool CReadReader::Eob()
{
	return is_eob;
}

// *********************************************************************************************
bool CReadReader::ShrinkBlock()
{
	copy(p_start, block.end(), block.begin());
	block.resize(block.end() - p_start);

	return true;
}

// *********************************************************************************************
bool CReadReader::find_eol(memory_chunk<char>::iterator& iter)
{
	auto p_end = block.end();

	if (iter == p_end)
		return false;

	for (; iter != p_end; ++iter)
		if (*iter == '\n' || *iter == '\r')
			break;

	return true;
}

// *********************************************************************************************
// Skip EOLs and change them into '\0' for easier parsing
void CReadReader::skip_eols(memory_chunk<char>::iterator& iter)
{
	while (iter != block.end() && *iter != 0)
		if (*iter == '\n' || *iter == '\r')
			*iter++ = 0;
		else
			break;
}

// EOF
