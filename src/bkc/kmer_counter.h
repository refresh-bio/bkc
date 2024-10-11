#pragma once

#include <atomic>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <thread>
#include <memory>
#include <cinttypes>
#include <chrono>
#include <memory>

#include <refresh/parallel_queues/lib/parallel-queues.h>
#include <refresh/memory_chunk/lib/memory_chunk.h>
#include <refresh/allocators/lib/memory_monotonic.h>

#include "memory_pool.h"
#include "../common/utils.h"
#include "../common/bkc_file.h"
#include "params.h"

#include <filters/poly_ACGT_filter.h>
#include <filters/artifacts_filter.h>
#include <types/base_coding.h>
#include <types/common_types.h>

using namespace std;
using namespace std::chrono;
using namespace refresh;

struct leader_follower_t {
	leader_t leader;
	follower_t follower;

	leader_follower_t() = default;
	leader_follower_t(leader_t _leader, follower_t _follower) :
		leader(_leader), follower(_follower) {}
	leader_follower_t(const leader_follower_t&) = default;
	leader_follower_t(leader_follower_t&&) = default;
	leader_follower_t& operator=(const leader_follower_t&) = default;
	leader_follower_t& operator=(leader_follower_t&&) = default;

	bool operator<(const leader_follower_t& rhs) const {
		if (this->leader != rhs.leader)
			return this->leader < rhs.leader;
		return this->follower < rhs.follower;
	}

	bool operator==(const leader_follower_t& rhs) const 
	{
		return this->leader == rhs.leader && this->follower == rhs.follower;
	}
};

struct leader_follower_count_t : public leader_follower_t
{
	uint64_t count;

	leader_follower_count_t() = default;
	leader_follower_count_t(leader_t _leader, follower_t _follower, uint64_t _count) :
		leader_follower_t(_leader, _follower),
		count(_count)
	{}
	leader_follower_count_t(const leader_follower_count_t&) = default;
	leader_follower_count_t(const leader_follower_t& x) : leader_follower_t(x), count(1) {};
	leader_follower_count_t(leader_follower_count_t&&) = default;
	leader_follower_count_t& operator=(const leader_follower_count_t&) = default;
	leader_follower_count_t& operator=(leader_follower_count_t&&) = default;

	bool operator<(const leader_follower_count_t& rhs) {
		if (this->leader != rhs.leader)
			return this->leader < rhs.leader;
		return this->follower < rhs.follower;
	}

	bool operator==(const leader_follower_t& rhs)
	{
		return this->leader == rhs.leader && this->follower == rhs.follower;
	}

	bool equal_lf(const leader_follower_t& rhs)
	{
		return this->leader == rhs.leader && this->follower == rhs.follower;
	}
};

struct kmer_count_t
{
	kmer_t kmer;
	uint64_t count;

	kmer_count_t() = default;
	kmer_count_t(kmer_t _kmer, uint64_t _count) :
		kmer(_kmer),
		count(_count)
	{}
	kmer_count_t(const kmer_count_t&) = default;
	kmer_count_t(const kmer_t& x) : kmer(x), count(1) {};
	kmer_count_t(kmer_count_t&&) = default;
	kmer_count_t& operator=(const kmer_count_t&) = default;
	kmer_count_t& operator=(kmer_count_t&&) = default;

	bool operator<(const kmer_count_t& rhs) {
		return this->kmer < rhs.kmer;
	}

	bool operator==(const kmer_t& rhs)
	{
		return this->kmer == rhs;
	}

	bool equal_lf(const kmer_t& rhs)
	{
		return this->kmer == rhs;
	}
};

class CBarcodedCounter
{
	const size_t no_blocks_in_queue = 3;
	const size_t no_chunks_per_file = 3;
	const size_t chunk_size = 64 << 20;
	const int gz_filtered_file_buffer_size = 16 << 20;
//	const bool filtered_input_in_FASTA = false;
	bool filtered_input_in_FASTA = false;
	bool canonical_mode = false;

	vector<pair<string, time_point<high_resolution_clock>>> times;

	vector<string> file_names;
	int no_threads = 1;
	int no_reading_threads = 0;

	string out_file_name = "./results.bkc";

	counting_mode_t counting_mode = counting_mode_t::unknown;
	uint32_t cbc_len = 16;
	uint32_t umi_len = 12;
	uint32_t leader_len = 27;
	uint32_t follower_len = 27;
	uint32_t soft_cbc_umi_len_limit = 0;
	bool allow_strange_cbc_umi_reads = false;
	uint32_t gap_len = 0;
	uint64_t min_leader_count = 1;
	uint32_t no_splits = 1;
	uint32_t rare_leader_thr = 0;
	uint32_t max_count = 65535;
	uint32_t zstd_level = 6;
	uint32_t verbosity_level = 0;
	bool apply_cbc_correction = false;
	uint32_t cbc_filtering_thr = 0;
	bool export_cbc_logs = false;
	string cbc_log_file_name;
	export_filtered_input_t export_filtered_input = export_filtered_input_t::none;
	string filtered_input_path;
	input_format_t input_format;
	output_format_t output_format;

	vector<string> cbc_file_names;
	vector<string> read_file_names;

	uint8_t sample_id_size_in_bytes;
	uint8_t barcode_size_in_bytes;
	uint8_t leader_size_in_bytes;
	uint8_t follower_size_in_bytes;
	uint8_t counter_size_in_bytes;

	uint64_t sample_id;

	atomic<uint64_t> no_sample_reads;

	vector<thread> reading_threads;
	vector<thread> counting_threads;
	vector<thread> reads_loading_threads;
	vector<thread> reads_exporting_threads;

	vector<unique_ptr<CMemoryPool<char>>> memory_pools;
	vector<unique_ptr<parallel_queue<pair<int, memory_chunk<char>>>>> block_queues;

	unique_ptr<parallel_queue<pair<int, string>>> fn_queue;

	using umi_t = uint64_t;
	using readfid_t = uint64_t;
	using umi_readfid_t = pair<umi_t, readfid_t>;

	vector<unordered_map<cbc_t, vector<umi_readfid_t>, refresh::MurMur64Hash>> cbc_dict;
	vector<pair<uint64_t, cbc_t>> cbc_vec, cbc_for_corr_vec;
	unordered_map<cbc_t, vector<vector<umi_readfid_t>>, refresh::MurMur64Hash> global_cbc_umi_dict;
	unordered_map<cbc_t, vector<readfid_t>, refresh::MurMur64Hash> global_cbc_dict;
	unordered_map<cbc_t, cbc_t, refresh::MurMur64Hash> correction_map;
	unordered_map<cbc_t, uint64_t, refresh::MurMur64Hash> cbc_stats;

	unordered_set<cbc_t> predefined_cbc;

	vector<vector<bool>> valid_reads;
	vector<unique_ptr<memory_monotonic_safe>> mma;
	vector<vector<uint8_t*>> sample_reads;
	vector<uint64_t> file_no_reads;
	vector<uint64_t> file_no_reads_after_cleanup;

	const int max_records_in_buffer = 128 << 10;
//	const int max_records_in_buffer = 2048 << 10;
	vector<shared_ptr<CBKCFile>> bkc_files;

	vector<unordered_map<uint64_t, uint64_t, refresh::MurMur64Hash>> leader_counts;
	unordered_set<uint64_t> valid_leaders;

	atomic<uint64_t> total_no_kmer_leaders_counts;
	atomic<uint64_t> sum_kmer_leaders_counts;

	atomic<uint64_t> total_no_kmer_pair_counts;
	atomic<uint64_t> sum_kmer_pair_counts;

	atomic<uint64_t> total_no_kmer_counts;
	atomic<uint64_t> sum_kmer_counts;

	atomic<uint64_t> a_total_read_len;
	atomic<uint64_t> a_total_no_reads;

	PolyACGTFilter poly_ACGT_filter;
	ArtifactsFilter artifacts_filter;
	BaseCoding4 base_coding4;
	BaseCoding3 base_coding3;

	uint64_t encode_read_id(uint64_t file_id, uint64_t read_no)
	{
		return (file_id << 32) + read_no;
	}

	pair<uint64_t, uint64_t> decode_read_id(uint64_t x)
	{
		return make_pair(x >> 32, x & 0xffffffffull);
	}

	template <typename T> void clear_vec(T&v)
	{
		v.clear();
		v.shrink_to_fit();
		T{}.swap(v);
	}

	std::string get_dedup_file_name(const std::string& input_path, const uint32_t id);

	void join_threads(vector<thread>& threads);
	void start_reading_threads();
	void start_counting_threads();
	void start_reads_loading_threads();
	void start_reads_exporting_threads();

	void init_queues_and_pools();
	void reinit_queues();

	bool init_bkc_files();

	void merge_cbc_dict();
	void sort_cbc_dict();
	void gather_cbc_stats();

	double calc_dist(vector<pair<uint64_t, uint64_t>>::iterator p, vector<pair<uint64_t, uint64_t>>::iterator q);
	int find_split(vector<pair<uint64_t, uint64_t>>& arr);
	void find_trusted_thr();
	void find_predefined_cbc();

	void find_CBC_corrections();
	void remove_non_trusted_CBC();
	void remove_duplicated_UMI();

	void create_valid_reads_lists();

	void list_cbc_dict(const string& suffix);

	string kmer_to_string(uint64_t kmer, int len);

	void enumerate_kmer_leaders_from_read(uint8_t* bases, vector<leader_t>& kmer_leaders);
	void enumerate_kmer_leaders_for_cbc(cbc_t cbc, vector<leader_t>& kmer_leaders);
	void count_leaders();
	void determine_valid_leaders();

	void enumerate_kmer_pairs_from_read(uint8_t* bases, vector<leader_follower_t>& kmer_pairs);
	void enumerate_kmers_from_read(uint8_t* bases, vector<kmer_t>& kmers);

	void enumerate_kmer_pairs_for_cbc(cbc_t cbc, vector<leader_follower_t>& kmer_pairs);
	void enumerate_kmers_for_cbc(cbc_t cbc, vector<kmer_t>& kmers);

	void sort_and_gather_kmer_pairs_for_cbc(vector<leader_follower_t>& kmer_pairs, vector<leader_follower_count_t>& kmer_pair_counts);
	void sort_and_gather_kmers_for_cbc(vector<kmer_t>& kmers, vector<kmer_count_t>& kmer_counts);

	void filter_rare_leader_sample_cbc(vector<leader_follower_count_t>& kmer_pair_counts);
	void filter_rare_kmer_sample_cbc(vector<kmer_count_t>& kmer_counts);

	void store_kmer_pairs(cbc_t cbc, vector<leader_follower_count_t>& kmer_pair_counts, vector<vector<bkc_record_t>> &record_buffers);
	void store_kmers(cbc_t cbc, vector<kmer_count_t>& kmer_pair_counts, vector<vector<bkc_record_t>> &record_buffers);

	void count_kmer_pairs();
	void count_kmers();

	void pack_records(vector<bkc_record_t>& records, vector<uint8_t>& packed_buffer);

	void set_CBC_file_names();
	void set_read_file_names();

public:
	CBarcodedCounter() = default;

	void SetParams(const CParams& params);

	void ShowTimings() {
		if(verbosity_level > 0)
			show_timings(times);
	}

	bool ProcessCBC();
	bool ProcessExportFilteredCBCReads();
	bool ProcessReads();
};

// EOF
