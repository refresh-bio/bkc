#include "bkc_file.h"
#include "utils.h"

// *********************************************************************************************
void CBKCFile::SetParams(uint8_t _sample_id_size_in_bytes, uint8_t _barcode_size_in_bytes, uint8_t _leader_size_in_bytes, uint8_t _follower_size_in_bytes, uint8_t _counter_size_in_bytes,
	uint8_t _barcode_len_in_symbols, uint8_t _leader_len_in_symbols, uint8_t _gap_len_in_symbols, uint8_t _follower_len_in_symbols, uint32_t _zstd_level)
{
	lock_guard<mutex> lck(mtx);

	sample_id_size_in_bytes = _sample_id_size_in_bytes;
	barcode_size_in_bytes = _barcode_size_in_bytes;
	leader_size_in_bytes = _leader_size_in_bytes;
	follower_size_in_bytes = _follower_size_in_bytes;
	counter_size_in_bytes = _counter_size_in_bytes;
	barcode_len_in_symbols = _barcode_len_in_symbols;
	leader_len_in_symbols = _leader_len_in_symbols;
	gap_len_in_symbols = _gap_len_in_symbols;
	follower_len_in_symbols = _follower_len_in_symbols;

	zstd_level = _zstd_level;
}

// *********************************************************************************************
bool CBKCFile::Create(const string& file_name, const output_format_t _output_format)
{
	lock_guard<mutex> lck(mtx);

	if (open_mode != open_mode_t::none)
		return false;

	output_format = _output_format;

	if (f.is_opened())
		f.close();

	f.set_compression_level(zstd_level);
	f.set_io_buffer_size(BUFFER_SIZE);

	switch (output_format)
	{
	case output_format_t::bkc:
		f.open_writing(file_name.c_str(),
			[&](FILE* f) {
				uint8_t to_save[] = { 'B', 'K', 'C', 1, 1, 0,
					(uint8_t)ordering,
					sample_id_size_in_bytes,
					barcode_size_in_bytes,
					leader_size_in_bytes,
					follower_size_in_bytes,
					counter_size_in_bytes,
					barcode_len_in_symbols,
					leader_len_in_symbols,
					follower_len_in_symbols,
					gap_len_in_symbols };
		return fwrite(to_save, 1, sizeof(to_save), f) == sizeof(to_save);
			});
		break;
	case output_format_t::splash:
		f.open_writing(file_name.c_str());
		break;
	default:
		open_mode = open_mode_t::none;
		return false;
	}

	if (!f.is_opened())
		return false;

	if(output_format == output_format_t::splash)
		save_header();

	open_mode = open_mode_t::writing;

	return true;
}

// *********************************************************************************************
bool CBKCFile::Open(const string& file_name)
{
	lock_guard<mutex> lck(mtx);

	if (open_mode != open_mode_t::none)
		return false;

	if (f.is_opened())
		f.close();

/*	switch (output_format)
	{
	case output_format_t::bkc:*/
		f.open_reading(file_name.c_str(),
			[&](FILE* f) {
				if (getc(f) != 'B')	return false;
				if (getc(f) != 'K')	return false;
				if (getc(f) != 'C')	return false;
				if (getc(f) != 1)	return false;
				if (getc(f) != 1)	return false;
				if (getc(f) != 0)	return false;

				ordering = (Header::ordering_t)getc(f);
				sample_id_size_in_bytes = getc(f);
				barcode_size_in_bytes = getc(f);
				leader_size_in_bytes = getc(f);
				follower_size_in_bytes = getc(f);
				counter_size_in_bytes = getc(f);
				barcode_len_in_symbols = getc(f);
				leader_len_in_symbols = getc(f);
				follower_len_in_symbols = getc(f);
				gap_len_in_symbols = getc(f);

				if (feof(f))
					return false;

				return true;
			});
/*		break;
	case output_format_t::splash:
		f.open_writing(file_name.c_str());
		break;
	default:
		open_mode = open_mode_t::none;
		return false;
	}*/

	if (!f.is_opened())
		return false;

//	if (output_format == output_format_t::splash)
//		save_header();

	open_mode = open_mode_t::reading;

	rec_len = sample_id_size_in_bytes + barcode_size_in_bytes + leader_size_in_bytes + follower_size_in_bytes + counter_size_in_bytes;

	return true;
}

// *********************************************************************************************
bool CBKCFile::Close()
{
	f.close();
	open_mode = open_mode_t::none;

	return true;
}

// *********************************************************************************************
void CBKCFile::save_header()
{

	if (!f.is_opened_for_writing())
		return;

	save_int_lsb(f, sample_id_size_in_bytes, 1);
	save_int_lsb(f, barcode_size_in_bytes, 1);
	save_int_lsb(f, leader_size_in_bytes, 1);
	save_int_lsb(f, follower_size_in_bytes, 1);
	save_int_lsb(f, counter_size_in_bytes, 1);
	save_int_lsb(f, barcode_len_in_symbols, 1);
	save_int_lsb(f, leader_len_in_symbols, 1);
	save_int_lsb(f, follower_len_in_symbols, 1);
	save_int_lsb(f, gap_len_in_symbols, 1);
}

// *********************************************************************************************
/*void CBKCFile::add_record(uint64_t sample_id, uint64_t barcode, uint64_t leader, uint64_t follower, uint64_t count)
{
#ifdef COMPACT_ENCODING
	rec_curr.clear();

	append_int_msb(rec_curr, sample_id, sample_id_size_in_bytes);
	append_int_msb(rec_curr, barcode, barcode_size_in_bytes);
	append_int_msb(rec_curr, leader, leader_size_in_bytes);
	append_int_msb(rec_curr, follower, follower_size_in_bytes);
	append_int_msb(rec_curr, count, counter_size_in_bytes);

	auto rec_compact = encode_shared_prefix(rec_prev, rec_curr);

	fwrite(rec_compact.data(), 1, rec_compact.size(), f);
	rec_prev = move(rec_curr);
#else
	save_int_lsb(f, sample_id, sample_id_size_in_bytes);
	save_int_lsb(f, barcode, barcode_size_in_bytes);
	save_int_lsb(f, leader, leader_size_in_bytes);
	save_int_lsb(f, follower, follower_size_in_bytes);
	save_int_lsb(f, count, counter_size_in_bytes);
#endif
}

// *********************************************************************************************
void CBKCFile::AddRecord(uint64_t sample_id, uint64_t barcode, uint64_t leader, uint64_t follower, uint64_t count)
{
	lock_guard<mutex> lck(mtx);

	add_record(sample_id, barcode, leader, follower, count);
}

// *********************************************************************************************
void CBKCFile::AddRecords(vector<bkc_record_t>& records)
{
	lock_guard<mutex> lck(mtx);

	for (auto& r : records)
		add_record(r.sample_id, r.barcode, r.leader, r.follower, r.count);
}
*/
// *********************************************************************************************
void CBKCFile::AddPacked(vector<uint8_t>& packed)
{
	lock_guard<mutex> lck(mtx);

	f.write((char*) packed.data(), packed.size());
}

// *********************************************************************************************
bool CBKCFile::GetRecord(uint64_t& sample_id, uint64_t& barcode, uint64_t& leader, uint64_t& follower, uint64_t& count)
{
	char no_same_symbols;

	if (!f.get(no_same_symbols))
		return false;

	rec_curr.assign(rec_prev.begin(), rec_prev.begin() + no_same_symbols);
	rec_curr.resize(rec_len);
	if (f.read((char*)rec_curr.data() + no_same_symbols, rec_len - no_same_symbols) != rec_len - no_same_symbols)
		return false;

	auto p = rec_curr.begin();

	switch (ordering)
	{
	case Header::default_ordering:
		load_int_msb(p, sample_id, sample_id_size_in_bytes);
		load_int_msb(p, barcode, barcode_size_in_bytes);
		load_int_msb(p, leader, leader_size_in_bytes);
		load_int_msb(p, follower, follower_size_in_bytes);
		load_int_msb(p, count, counter_size_in_bytes);
		break;
	default:
		return false;
	}

	swap(rec_curr, rec_prev);

	return true;
}

// *********************************************************************************************
void CBKCFile::GetLens(uint8_t& _barcode_len_in_symbols, uint8_t& _leader_len_in_symbols, uint8_t& _follower_len_in_symbols, uint8_t& _counter_size_in_bytes)
{
	_barcode_len_in_symbols = barcode_len_in_symbols;
	_leader_len_in_symbols = leader_len_in_symbols;
	_follower_len_in_symbols = follower_len_in_symbols;
	_counter_size_in_bytes = counter_size_in_bytes;
}

// EOF
