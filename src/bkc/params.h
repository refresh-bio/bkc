#pragma once

#include <cstdint>
#include <iostream>
#include <vector>

#include "../common/defs.h"
#include <types/common_types.h>

using namespace std;

// *********************************************************************************************
inline output_format_t output_format_from_string(const std::string& str) {
	if (str == "bkc" || str == "BKC")
		return output_format_t::bkc;
	else if (str == "splash" || str == "SPLASH")
		return output_format_t::splash;
	else
	{
		return output_format_t::unknown;
	}
}

// *********************************************************************************************
inline counting_mode_t counting_mode_from_string(const std::string& str) {
	if (str == "single")
		return counting_mode_t::single;
	else if (str == "pair")
		return counting_mode_t::pair;
	else if (str == "filter")
		return counting_mode_t::filter;
	else
	{
		return counting_mode_t::unknown;
	}
}

// *********************************************************************************************
inline std::string to_string(output_format_t output_format) {
	switch (output_format) {
	case output_format_t::bkc:
		return "bkc";
	case output_format_t::splash:
	default:
		return "unknown";
	}
}

// *********************************************************************************************
inline std::string to_string(export_filtered_input_t export_filtered_input) {
	switch (export_filtered_input) {
	case export_filtered_input_t::none:
		return "none";
	case export_filtered_input_t::first:
		return "first";
	case export_filtered_input_t::second:
		return "second";
	case export_filtered_input_t::both:
		return "both";
	default:
		return "unknown";
	}
}

// *********************************************************************************************
inline string technology_str(technology_t technology)
{
	if (technology == technology_t::ten_x)
		return "10x";
	else if (technology == technology_t::visium)
		return "visium";
	else
		return "unknown";
}

// *********************************************************************************************
struct CParams
{
	counting_mode_t counting_mode = counting_mode_t::single;
	param_t<uint32_t> cbc_len{ 10, 16, 16 };
	param_t<uint32_t> umi_len{ 8, 16, 12 };
	param_t<uint32_t> leader_len{ 1, 31, 27 };
	param_t<uint32_t> follower_len{ 0, 31, 0 };
	param_t<uint32_t> gap_len{ 0, ~0u, 0 };
	param_t<uint32_t> soft_cbc_umi_len_limit{ 0, 1'000'000'000, 0 };
	param_t<uint32_t> no_splits{ 1, 256, 1 };
	param_t<uint32_t> no_threads{ 0, 256, 8 };
	param_t<uint32_t> max_count{ 1, ~0u, 65535 };
	param_t<uint32_t> zstd_level{ 0, 19, 6 };
	bool canonical_mode{ false };
	uint32_t sample_id{ 0 };
	vector<string> cbc_file_names;
	vector<string> read_file_names;
	string out_file_name{ "./results.bkc" };
	param_t<uint32_t> poly_ACGT_len{ 0, 31, 0 };
	string artifacts;
	bool apply_filter_illumina_adapters{ false };
	param_t<uint32_t> verbosity_level{ 0, 2, 0 };
	param_t<uint32_t> rare_leader_thr{ 0, 255, 5 };
	bool apply_cbc_correction{ false };
	param_t<uint32_t> cbc_filtering_thr{ 0, ~0u, 0 };			// auto
	technology_t technology{ technology_t::ten_x };
	vector<string> predefined_cbc;
	bool export_cbc_logs{ false };
	string predefined_cbc_fn;
	string cbc_log_file_name;
	export_filtered_input_t export_filtered_input { export_filtered_input_t::none };
	string filtered_input_path{};
	bool allow_strange_cbc_umi_reads{ false };
	input_format_t input_format{ input_format_t::fastq };
	output_format_t output_format {output_format_t::bkc};
};
