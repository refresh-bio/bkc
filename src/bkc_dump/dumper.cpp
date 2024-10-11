#include "dumper.h"

#include <string>

#include <refresh/conversions/lib/conversions.h>
#include "../common/bkc_file.h"

// *********************************************************************************************
bool CDumper::SetParams(const CParams& _params)
{
	params = _params;

	return true;
}

// *********************************************************************************************
bool CDumper::Dump()
{
	FILE* f_out = fopen(params.output_file_name.c_str(), "wb");

	if (!f_out)
	{
		cerr << "Cannot open: " << params.output_file_name << endl;
		return false;
	}

	setvbuf(f_out, nullptr, _IOFBF, 16 << 20);

	char line[1024];

	for (uint32_t i = 0; i < params.no_splits.get(); ++i)
	{
		CBKCFile bkc_file;
		string fn = params.input_file_name;

		uint8_t barcode_len_in_symbols;
		uint8_t leader_len_in_symbols;
		uint8_t follower_len_in_symbols;
		uint8_t counter_size_in_bytes;

		if (params.no_splits.get() > 1)
			fn += "." + to_string(i);

		if (!bkc_file.Open(fn))
		{
			cerr << "Cannot open: " << fn << endl;
			fclose(f_out);
			return false;
		}

		bkc_file.GetLens(barcode_len_in_symbols, leader_len_in_symbols, follower_len_in_symbols, counter_size_in_bytes);

		uint64_t sample_id;
		cbc_t cbc;
		leader_t leader;
		follower_t follower;
		uint64_t counter;

		while (bkc_file.GetRecord(sample_id, cbc, leader, follower, counter))
		{
			char* p = line;

			p += refresh::int_to_pchar(sample_id, p, '\t');
			p += refresh::kmer_to_pchar(cbc, p, barcode_len_in_symbols, false, '\t');
			p += refresh::kmer_to_pchar(leader, p, leader_len_in_symbols, false, '\t');
			if (follower_len_in_symbols)
				p += refresh::kmer_to_pchar(follower, p, follower_len_in_symbols, false, '\t');
			p += refresh::int_to_pchar(counter, p, '\n');

			fwrite(line, 1, p - line, f_out);
		}
	}

	fclose(f_out);

	return true;
}
