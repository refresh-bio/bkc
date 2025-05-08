#include "kmer_counter.h"
#include "fq_reader.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <atomic>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <random>
#include <filesystem>

#include <refresh/hash_tables/lib/murmur_hash.h>
#include <refresh/sort/lib/pdqsort_par.h>
#include <filters/illumina_adapters_static.h>
#include <types/kmer.h>

#include <mimalloc.h>

//#define AGGRESIVE_MEMORY_SAVING
#define USE_READ_COMPRESSION

// *********************************************************************************************
void CBarcodedCounter::join_threads(vector<thread>& threads)
{
	for (auto& t : threads)
		t.join();
}

// *********************************************************************************************
void CBarcodedCounter::SetParams(const CParams& params)
{
	no_threads = params.no_threads.get();
	cbc_len = params.cbc_len.get();
	umi_len = params.umi_len.get();
	soft_cbc_umi_len_limit = params.soft_cbc_umi_len_limit.get();
	allow_strange_cbc_umi_reads = params.allow_strange_cbc_umi_reads;
	counting_mode = params.counting_mode;
	leader_len = params.leader_len.get();
	follower_len = params.follower_len.get();
	gap_len = params.gap_len.get();
	min_leader_count = 1;					// currently fixed
	rare_leader_thr = params.rare_leader_thr.get();
	no_splits = params.no_splits.get();
	max_count = params.max_count.get();
	zstd_level = params.zstd_level.get();

	if (counting_mode == counting_mode_t::single)
		follower_len = 0;

	sample_id_size_in_bytes = no_bytes(params.sample_id);
	barcode_size_in_bytes = (cbc_len + 3) / 4;
	leader_size_in_bytes = (leader_len + 3) / 4;
	follower_size_in_bytes = (follower_len + 3) / 4;
	counter_size_in_bytes = no_bytes(max_count);

	sample_id = params.sample_id;
	canonical_mode = params.canonical_mode;
	out_file_name = params.out_file_name;

	poly_ACGT_filter = PolyACGTFilter(params.poly_ACGT_len.get());

	artifacts_filter = ArtifactsFilter(params.artifacts);
	if (params.apply_filter_illumina_adapters)
		artifacts_filter.Add(12, IlluminaAdaptersStatic::Get12Mers());

	verbosity_level = params.verbosity_level.get();

	apply_cbc_correction = params.apply_cbc_correction;
	cbc_filtering_thr = params.cbc_filtering_thr.get();

	input_format = params.input_format;
	filtered_input_in_FASTA = input_format == input_format_t::fasta;
	output_format = params.output_format;

	predefined_cbc.clear();

	for (const auto& s : params.predefined_cbc)
		predefined_cbc.insert(base_coding4.encode_bases_2b(s));

	if (params.export_cbc_logs)
	{
		export_cbc_logs = true;
		cbc_log_file_name = params.cbc_log_file_name;
	}

	if (params.export_filtered_input != export_filtered_input_t::none)
	{
		export_filtered_input = params.export_filtered_input;
		filtered_input_path = params.filtered_input_path;

		if (filtered_input_path.empty())
			filtered_input_path = ".";

		if (filtered_input_path.back() != '\\' && filtered_input_path.back() != '/')
			filtered_input_path.push_back('/');
	}

	cbc_file_names = params.cbc_file_names;
	read_file_names = params.read_file_names;
}

// *********************************************************************************************
void CBarcodedCounter::set_CBC_file_names()
{
	file_names = cbc_file_names;

	fn_queue = make_unique<parallel_queue<pair<int, string>>>(cbc_file_names.size());

	for(int i = 0; i < (int) cbc_file_names.size(); ++i)
		fn_queue->push(make_pair(i, cbc_file_names[i]));

	fn_queue->mark_completed();
}

// *********************************************************************************************
void CBarcodedCounter::set_read_file_names()
{
	file_names = read_file_names;

	fn_queue = make_unique<parallel_queue<pair<int, string>>>(read_file_names.size());

	for (int i = 0; i < (int) read_file_names.size(); ++i)
		fn_queue->push(make_pair(i, read_file_names[i]));

	fn_queue->mark_completed();
}

// *********************************************************************************************
void CBarcodedCounter::start_reading_threads()
{
	reading_threads.clear();
	reading_threads.reserve(no_reading_threads);

	for (int i = 0; i < no_reading_threads; ++i)
	{
		reading_threads.push_back(thread([&, i] {
			int thread_id = i;
			pair<int, string> id_fn;
			CFastXReader fqx(input_format == input_format_t::fastq);
			memory_chunk<char> mc;

			while (fn_queue->pop(id_fn))
			{
				if(verbosity_level >= 2)
					std::cerr << "Reading thread " + to_string(thread_id) + " opens: " + id_fn.second + "\n";

				if (fqx.Open(id_fn.second))
				{
					if (verbosity_level >= 2)
						std::cerr << "File " + id_fn.second + " opened\n";
				}
				else
				{
					std::cerr << "Error: File " + id_fn.second + " cannot be opened\n";
					exit(1);
				}

				while (!fqx.Eof())
				{
					memory_pools[thread_id]->Pop(mc);

					if (!fqx.ReadBlock(mc))
					{
						memory_pools[thread_id]->Push(mc);
						break;
					}
					else
					{
//						cerr << "Reading thread " + to_string(thread_id) + " loaded block of size: " + to_string(mc.size()) + "\n";
					}

					block_queues[thread_id]->push(make_pair(id_fn.first, move(mc)));
				}
			}

			if (verbosity_level >= 2)
				std::cerr << "Reading thread " + to_string(thread_id) + " completed\n";

			block_queues[thread_id]->mark_completed();
		}));
	}
}

// *********************************************************************************************
void CBarcodedCounter::start_counting_threads()
{
	counting_threads.clear();
	counting_threads.reserve(no_reading_threads);

	file_no_reads.clear();
	file_no_reads.resize(file_names.size(), 0);

	no_sample_reads = 0;

	for (int i = 0; i < no_reading_threads; ++i)
	{
		counting_threads.push_back(thread([&, i] {
			int thread_id = i;

			pair<int, memory_chunk<char>> id_mc;
			CReadReader read_reader(input_format == input_format_t::fastq);
			read_desc_t read_desc;

			auto& my_block_queue = block_queues[thread_id];
			auto& my_memory_pool = memory_pools[thread_id];
//			auto& my_cbc_dict = cbc_dict[thread_id];

			unordered_map<cbc_t, vector<umi_readfid_t>, refresh::MurMur64Hash> my_cbc_dict;
			my_cbc_dict.max_load_factor(0.8);

			BaseCoding4 bc4;

			int total_no_reads = 0;

			int file_id = -1;
			int file_read_id = 0;

			while (my_block_queue->pop(id_mc))
			{
				if (id_mc.first != file_id)
				{
					file_id = id_mc.first;
					file_read_id = 0;
				}

//				cerr << "Counting thread " + to_string(thread_id) + " got block of size : " + to_string(id_mc.second.size()) + "\n";

				read_reader.Assign(id_mc.second);

				int no_reads = 0;

				while (read_reader.GetRead(read_desc))
				{
					auto read_len = strlen(read_desc.bases);
						
					if (read_len < cbc_len + umi_len || read_len > cbc_len + umi_len + soft_cbc_umi_len_limit)
					{
						std::cerr << "Strange read: " + string(read_desc.header) + " " + string(read_desc.bases) + "\n";
						if (!allow_strange_cbc_umi_reads)
							exit(1);
						++no_reads;
						++file_read_id;
						continue;
					}

					cbc_t cbc = bc4.encode_bases_2b(read_desc.bases, read_desc.bases + cbc_len);
					umi_t umi = bc4.encode_bases_2b(read_desc.bases + cbc_len, read_desc.bases + cbc_len + umi_len);

					if (cbc != ~0ull && umi != ~0ull)
						my_cbc_dict[cbc].emplace_back(umi, encode_read_id(id_mc.first, file_read_id++));
					else
						file_read_id++;

					++no_reads;
				}

				total_no_reads += no_reads;

				file_no_reads[file_id] += no_reads;

//				cerr << "Counting thread " + to_string(thread_id) + " found " + to_string(no_reads) + " reads in block\n";

				my_memory_pool->Push(id_mc.second);
			}

			cbc_dict[thread_id] = move(my_cbc_dict);

			if (verbosity_level >= 2)
				std::cerr << "Counting thread " + to_string(thread_id) + " found " + to_string(total_no_reads) + " reads in total and completed\n";

			no_sample_reads += total_no_reads;

			}));
	}
}

// *********************************************************************************************
std::string CBarcodedCounter::get_dedup_file_name(const std::string& input_path, const uint32_t id)
{
	std::filesystem::path path(input_path);

	//remove .gz if present
	if (path.extension() == ".gz")
		path = path.stem();
	//to remove ".fastq" or ".fasta"
	path = path.stem();

	std::filesystem::path out(filtered_input_path);

	out /= path.filename();
	out += ".dedup"s + (filtered_input_in_FASTA ? ".fasta" : ".fastq") + ".gz";

	return out.string();
}

// *********************************************************************************************
void CBarcodedCounter::start_reads_exporting_threads()
{
	reads_exporting_threads.clear();
	reads_exporting_threads.reserve(no_threads);

	no_sample_reads = 0;

	a_total_read_len = 0;
	a_total_no_reads = 0;

	for (int i = 0; i < no_reading_threads; ++i)
	{
		reads_exporting_threads.push_back(thread([&, i] {
			int thread_id = i;

			pair<int, memory_chunk<char>> id_mc;
			CReadReader read_reader(input_format == input_format_t::fastq);
			read_desc_t read_desc;

			auto& my_block_queue = block_queues[thread_id];
			auto& my_memory_pool = memory_pools[thread_id];

			int total_no_reads = 0;

			int file_id = -1;
			int file_read_id = 0;			// read id after relabelling
			int file_read_id_raw = 0;		// read id before relabelling

			uint64_t my_total_read_len = 0;
			uint64_t my_total_no_reads = 0;

			gzFile filtered_file = nullptr;

			while (my_block_queue->pop(id_mc))
			{
				if (id_mc.first != file_id)
				{
					file_id = id_mc.first;
					file_read_id = 0;
					file_read_id_raw = 0;

					if (filtered_file)
						gzclose(filtered_file);

					string out_fn = get_dedup_file_name(file_names[file_id], 1);

					filtered_file = gzopen(out_fn.c_str(), "w3");
					if (filtered_file)
						gzbuffer(filtered_file, gz_filtered_file_buffer_size);
					else
						cerr << "Cannot create filtered file: " << file_names[file_id] << endl;
				}
				//				cerr << "Counting thread " + to_string(thread_id) + " got block of size : " + to_string(id_mc.second.size()) + "\n";

				read_reader.Assign(id_mc.second);

				int no_reads = 0;

				while (read_reader.GetRead(read_desc))
				{
					if (!valid_reads[file_id][file_read_id_raw])
					{
						++file_read_id_raw;
						++no_reads;
						continue;
					}

					int read_len = strlen(read_desc.bases);

					if (filtered_file)
					{
						gzputs(filtered_file, read_desc.header);
						gzputc(filtered_file, '\n');
						gzputs(filtered_file, read_desc.bases);
						gzputc(filtered_file, '\n');

						if (!filtered_input_in_FASTA)
						{
							gzputs(filtered_file, read_desc.plus);
							gzputc(filtered_file, '\n');
							gzputs(filtered_file, read_desc.quality);
							gzputc(filtered_file, '\n');
						}
					}

					my_total_no_reads++;
					my_total_read_len += read_len;

					++file_read_id;
					++file_read_id_raw;
					++no_reads;
				}

				total_no_reads += no_reads;

//				cerr << "Reads loading thread " + to_string(thread_id) + " found " + to_string(no_reads) + " reads in block\n";

				my_memory_pool->Push(id_mc.second);
			}

			if (filtered_file)	
				gzclose(filtered_file);

			if (verbosity_level >= 2)
				std::cerr << "Reads exporting thread " + to_string(thread_id) + " found " + to_string(total_no_reads) + " reads in total and completed\n";

//			no_sample_reads += total_no_reads;

//			a_total_no_reads += my_total_no_reads;
//			a_total_read_len += my_total_read_len;

			}));
	}
}

// *********************************************************************************************
void CBarcodedCounter::start_reads_loading_threads()
{
	reads_loading_threads.clear();
	reads_loading_threads.reserve(no_threads);

	no_sample_reads = 0;

	a_total_read_len = 0;
	a_total_no_reads = 0;

	mma.reserve(file_names.size());
	for (int i = 0; i < (int) file_names.size(); ++i)
		mma.push_back(move(make_unique<memory_monotonic_safe>(16 << 20, 1)));

	sample_reads.resize(file_names.size());
	for (int i = 0; i < (int) sample_reads.size(); ++i)
		sample_reads[i].resize(file_no_reads_after_cleanup[i], nullptr);

	for (int i = 0; i < no_reading_threads; ++i)
	{
		reading_threads.push_back(thread([&, i] {
			int thread_id = i;

			pair<int, memory_chunk<char>> id_mc;
			CReadReader read_reader(input_format == input_format_t::fastq);
			read_desc_t read_desc;

			auto& my_block_queue = block_queues[thread_id];
			auto& my_memory_pool = memory_pools[thread_id];
			auto& my_mma = mma[thread_id];

			int total_no_reads = 0;

			int file_id = -1;
			int file_read_id = 0;			// read id after relabelling
			int file_read_id_raw = 0;		// read id before relabelling

			uint64_t my_total_read_len = 0;
			uint64_t my_total_no_reads = 0;

			gzFile filtered_file = nullptr;

			while (my_block_queue->pop(id_mc))
			{
				if (id_mc.first != file_id)
				{
					file_id = id_mc.first;
					file_read_id = 0;
					file_read_id_raw = 0;

					if (((uint32_t) export_filtered_input) & (uint32_t) export_filtered_input_t::second)
					{
						if (filtered_file)
							gzclose(filtered_file);

						string out_fn = get_dedup_file_name(file_names[file_id], 2);

						filtered_file = gzopen(out_fn.c_str(), "w3");
						if (filtered_file)
							gzbuffer(filtered_file, gz_filtered_file_buffer_size);
						else
							cerr << "Cannot create filtered file: " << file_names[file_id] << endl;
					}
				}
				//				cerr << "Counting thread " + to_string(thread_id) + " got block of size : " + to_string(id_mc.second.size()) + "\n";

				read_reader.Assign(id_mc.second);

				int no_reads = 0;

				while (read_reader.GetRead(read_desc))
				{
					if (!valid_reads[file_id][file_read_id_raw])
					{
						++file_read_id_raw;
						++no_reads;
						continue;
					}

					int read_len = strlen(read_desc.bases);

					if (filtered_file)
					{
						gzputs(filtered_file, read_desc.header);
						gzputc(filtered_file, '\n');
						gzputs(filtered_file, read_desc.bases);
						gzputc(filtered_file, '\n');

						if (!filtered_input_in_FASTA)
						{
							gzputs(filtered_file, read_desc.plus);
							gzputc(filtered_file, '\n');
							gzputs(filtered_file, read_desc.quality);
							gzputc(filtered_file, '\n');
						}
					}

					my_total_no_reads++;
					my_total_read_len += read_len;

#ifdef USE_READ_COMPRESSION
					size_t pred_len = (read_len + 1 + 2) / 3;
					uint8_t *p = (uint8_t*)(my_mma->allocate(pred_len));
					size_t enc_len = base_coding3.encode_bases(read_desc.bases, read_len, p);

					if (pred_len != enc_len)
						std::cerr << to_string(read_len) + "   -   "  + to_string(pred_len) + " : " + to_string(enc_len) + "\n";

#else
					uint8_t* p = (uint8_t*)(my_mma->allocate(read_len + 1));
					memcpy(p, read_desc.bases, read_len + 1);			// !!! Add compression of reads (at least 2 bases -> 1 byte) here
#endif

					sample_reads[file_id][file_read_id] = p;

					++file_read_id;
					++file_read_id_raw;
					++no_reads;
				}

				total_no_reads += no_reads;

//				cerr << "Reads loading thread " + to_string(thread_id) + " found " + to_string(no_reads) + " reads in block\n";

				my_memory_pool->Push(id_mc.second);
			}

			if (filtered_file)	
				gzclose(filtered_file);

			if (verbosity_level >= 2)
				std::cerr << "Reads loading thread " + to_string(thread_id) + " found " + to_string(total_no_reads) + " reads in total and completed\n";

			no_sample_reads += total_no_reads;

			a_total_no_reads += my_total_no_reads;
			a_total_read_len += my_total_read_len;

			}));
	}
}

// *********************************************************************************************
void CBarcodedCounter::init_queues_and_pools()
{
	block_queues.clear();

	for(int i = 0; i < no_reading_threads; ++i)
		block_queues.emplace_back(make_unique<parallel_queue<pair<int, memory_chunk<char>>>>(no_blocks_in_queue));

	memory_pools.clear();

	for (int i = 0; i < no_reading_threads; ++i)
		memory_pools.emplace_back(make_unique<CMemoryPool<char>>(no_chunks_per_file, chunk_size));

	cbc_dict.resize(no_reading_threads);
}

// *********************************************************************************************
void CBarcodedCounter::reinit_queues()
{
	block_queues.clear();

	for(int i = 0; i < no_reading_threads; ++i)
		block_queues.emplace_back(make_unique<parallel_queue<pair<int, memory_chunk<char>>>>(no_blocks_in_queue));
}

// *********************************************************************************************
bool CBarcodedCounter::init_bkc_files()
{
	bkc_files.clear();
	bkc_files.reserve(no_splits);

	bool r = true;

	for (uint32_t i = 0; i < no_splits; ++i)
	{
		bkc_files.push_back(make_shared<CBKCFile>());
		bkc_files.back()->SetParams(sample_id_size_in_bytes, barcode_size_in_bytes, leader_size_in_bytes, follower_size_in_bytes, counter_size_in_bytes,
			cbc_len, leader_len, gap_len, follower_len, zstd_level);

		if(no_splits == 1)
			r &= bkc_files.back()->Create(out_file_name, output_format);
		else
			r &= bkc_files.back()->Create(out_file_name + "." + to_string(i), output_format);
	}

	if (!r)
		bkc_files.clear();
		
	return r;
}

// *********************************************************************************************
void CBarcodedCounter::merge_cbc_dict()
{
	vector<thread> threads;

	vector<cbc_t> cbcs;
	atomic_int id{ 0 };

	for (int i = 1; i < (int) cbc_dict.size(); ++i)
		for (auto& x : cbc_dict[i])
			if (!cbc_dict[0].count(x.first))
				cbc_dict[0][x.first] = {};

	cbcs.reserve(cbc_dict[0].size());

	for (auto& x : cbc_dict[0])
		cbcs.emplace_back(x.first);

	threads.reserve(no_threads);

	for (int t = 0; t < no_threads; ++t)
		threads.emplace_back([&] {
		int curr_id;

		vector<vector<umi_readfid_t>*> ptrs(cbc_dict.size(), nullptr);

		while (true)
		{
			curr_id = id.fetch_add(1);
			if (curr_id >= (int)cbcs.size())
				break;

			cbc_t curr_cbc = cbcs[curr_id];
			size_t no_items = 0;

			for (int i = 0; i < (int) cbc_dict.size(); ++i)
			{
				auto p = cbc_dict[i].find(curr_cbc);

				if (p == cbc_dict[i].end())
					ptrs[i] = nullptr;
				else
				{
					ptrs[i] = &(p->second);
					no_items += p->second.size();
				}
			}

			ptrs[0]->reserve(no_items);

			for (int i = 1; i < (int) cbc_dict.size(); ++i)
			{
				if (!ptrs[i])
					continue;

				ptrs[0]->insert(ptrs[0]->end(), ptrs[i]->begin(), ptrs[i]->end());
				clear_vec(*ptrs[i]);
			}
		}
			});

	join_threads(threads);

	for (int i = 1; i < (int) cbc_dict.size(); ++i)
		cbc_dict[i].clear();
}

// *********************************************************************************************
void CBarcodedCounter::sort_cbc_dict()
{
	atomic_int id{ 0 };

	vector<thread> threads;

	vector<vector<umi_readfid_t>*> umis;

	umis.reserve(cbc_dict[0].size());
	for (auto& x : cbc_dict[0])
		umis.emplace_back(&(x.second));

	threads.reserve(no_threads);

	for (int i = 0; i < no_threads; ++i)
		threads.emplace_back([&] {
		int curr_id = -1;
		while (true)
		{
			curr_id = id.fetch_add(1);
			if (curr_id >= (int)umis.size())
				break;

			stable_sort(umis[curr_id]->begin(), umis[curr_id]->end());
		}
			});

	join_threads(threads);
}

// *********************************************************************************************
void CBarcodedCounter::gather_cbc_stats()
{
	cbc_stats.max_load_factor(0.8);
	cbc_stats.reserve(cbc_dict[0].size());

	for (size_t i = 0; i < cbc_dict.size(); ++i)
		for (const auto& x : cbc_dict[i])
			cbc_stats[x.first] += x.second.size();
}

// *********************************************************************************************
void CBarcodedCounter::find_CBC_corrections()
{
	unordered_map<cbc_t, vector<cbc_t>> candidate_corrections;

	for (auto x : cbc_vec)
	{
		for (uint32_t i = 0; i < cbc_len; ++i)
		{
			auto xc = x.second;
			uint32_t shift = 2 * i;
			cbc_t mask = 3ull << shift;

			for (uint32_t j = 0; j < 4; ++j)
			{
				xc = (xc & ~mask) + (((cbc_t)j) << shift);
				candidate_corrections[xc].emplace_back(x.second);
			}
		}
	}

	for (auto x : cbc_for_corr_vec)
	{
		auto p = candidate_corrections.find(x.second);

		if (p != candidate_corrections.end() && p->second.size() == 1)
			correction_map[p->first] = p->second.front();
	}

	clear_vec(cbc_for_corr_vec);
}

// *********************************************************************************************
void CBarcodedCounter::remove_non_trusted_CBC()
{
	unordered_set<cbc_t> trusted_CBC;

	trusted_CBC.max_load_factor(0.8);

	for (auto x : cbc_vec)
		trusted_CBC.insert(x.second);

	global_cbc_umi_dict.max_load_factor(0.8);

	for (size_t i = 0; i < cbc_dict.size(); ++i)
	{
		for (auto& x : cbc_dict[i])
		{
			cbc_t cbc = x.first;

			if (apply_cbc_correction)
			{
				auto p = correction_map.find(cbc);
				if (p != correction_map.end())
					cbc = p->second;
			}

			if (trusted_CBC.count(cbc))
			{
				auto& r = global_cbc_umi_dict[cbc];
				if (r.empty())
					r.reserve(cbc_dict.size() - i);
				r.emplace_back(move(x.second));
			}

			clear_vec(x.second);
		}

		cbc_dict[i].clear();
	}

	clear_vec(cbc_dict);
}

// *********************************************************************************************
void CBarcodedCounter::remove_duplicated_UMI()
{
	atomic_int id{ 0 };
	atomic_uint64_t total_no_after_removal { 0 };
	atomic<uint64_t> a_total_no_reads_before_UMI_cleaning{ 0 };

	vector<thread> threads;

	vector<cbc_t> v_cbc;

	global_cbc_dict.reserve(global_cbc_umi_dict.size());
	v_cbc.reserve(global_cbc_umi_dict.size());

	cbc_vec.clear();

	global_cbc_dict.max_load_factor(0.8);

	for (auto& x : global_cbc_umi_dict)
	{
		v_cbc.emplace_back(x.first);
		global_cbc_dict.emplace(x.first, vector<uint64_t>());
		cbc_vec.emplace_back(0, x.first);
	}

	threads.reserve(no_threads);

	for (int i_thread = 0; i_thread < no_threads; ++i_thread)
		threads.emplace_back([&id, &v_cbc, this, &a_total_no_reads_before_UMI_cleaning, &total_no_after_removal] {
		int curr_id = -1;
		
		vector<vector<umi_readfid_t>::iterator> p_src;

		while (true)
		{
			curr_id = id.fetch_add(1);
			if (curr_id >= (int)v_cbc.size())
				break;

			mt19937_64 mt(v_cbc[curr_id]);

			p_src.clear();

			size_t cur_cbc_size = 0;
			const umi_readfid_t umi_guard{ ~0ull, ~0ull };

			auto& gd_src = global_cbc_umi_dict[v_cbc[curr_id]];
			auto& gd_dest = global_cbc_dict[v_cbc[curr_id]];

			for (auto& x : gd_src)
			{
				cur_cbc_size += x.size();
				std::stable_sort(x.begin(), x.end());

				// Add guard records
				x.emplace_back(umi_guard);

				p_src.emplace_back(x.begin());
			}

			a_total_no_reads_before_UMI_cleaning += cur_cbc_size;

			size_t no_streams = gd_src.size();

			auto min_umi = *p_src[0];
			size_t min_umi_id = 0;

			for (size_t i = 1; i < no_streams; ++i)
			{
				if (*p_src[i] < min_umi)
				{
					min_umi = *p_src[i];
					min_umi_id = i;
				}
			}

			auto prev_umi = min_umi;
			gd_dest.emplace_back(p_src[min_umi_id]->second);
			++p_src[min_umi_id];
			size_t no_same_umi = 1;

			while (true)
			{
				// Find min UMI
				min_umi = *p_src[0];
				size_t min_umi_id = 0;

				for (size_t i = 1; i < no_streams; ++i)
				{
					if (*p_src[i] < min_umi)
					{
						min_umi = *p_src[i];
						min_umi_id = i;
					}
				}

				if (min_umi.first != prev_umi.first)
				{
					if (no_same_umi > 1)
					{
						size_t id_to_preserve = mt() % no_same_umi;
						gd_dest[gd_dest.size() - no_same_umi] = gd_dest[gd_dest.size() - id_to_preserve - 1];
						gd_dest.resize(gd_dest.size() - no_same_umi + 1);
					}
					no_same_umi = 0;
				}

				if (min_umi == umi_guard)
					break;

				prev_umi = min_umi;
				gd_dest.emplace_back(p_src[min_umi_id]->second);
				++p_src[min_umi_id];
				++no_same_umi;
			}

			total_no_after_removal += gd_dest.size();

			cbc_vec[curr_id].first = gd_dest.size();

			gd_dest.shrink_to_fit();

			clear_vec(gd_src);
		}
		});

	join_threads(threads);

	global_cbc_umi_dict.clear();

	stable_sort(cbc_vec.begin(), cbc_vec.end(), greater<pair<uint64_t, cbc_t>>());

	if (cbc_filtering_thr)
	{
		while (!cbc_vec.empty())
		{
			if (cbc_vec.back().first < cbc_filtering_thr)
			{
				total_no_after_removal -= cbc_vec.back().first;
				cbc_vec.pop_back();
			}
			else
				break;
		}
	}

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. of reads before UMI cleaning: " + to_string(a_total_no_reads_before_UMI_cleaning) + "\n";
		std::cerr << "Total no. of reads after UMI cleaning: " + to_string(total_no_after_removal) + "\n";
	}
}

// *********************************************************************************************
void CBarcodedCounter::create_valid_reads_lists()
{
	uint64_t no_valid_reads = 0;

	valid_reads.clear();
	valid_reads.resize(file_names.size());

	vector<vector<uint32_t>> relabelling(file_names.size());
	file_no_reads_after_cleanup.resize(file_names.size(), 0);

	for (int i = 0; i < (int)file_names.size(); ++i)
	{
		relabelling[i].resize(file_no_reads[i], 0);
		valid_reads[i].resize(file_no_reads[i], false);
	}

	uint64_t file_id;
	uint64_t read_id;

	// Count and mark valid reads
	for (auto& x : global_cbc_dict)
	{
		for (auto y : x.second)
		{
			tie(file_id, read_id) = decode_read_id(y);
			valid_reads[file_id][read_id] = true;
			file_no_reads_after_cleanup[file_id]++;
		}

		no_valid_reads += x.second.size();
	}

	// Determine new read ids mapping
	for (size_t file_id = 0; file_id < file_names.size(); ++file_id)
	{
		uint32_t cnt = 0;

		for (size_t i = 0; i < file_no_reads[file_id]; ++i)
			if (valid_reads[file_id][i])
				relabelling[file_id][i] = cnt++;
	}

	// Change read ids
	for (auto& x : global_cbc_dict)
	{
		for (auto& y : x.second)
		{
			tie(file_id, read_id) = decode_read_id(y);
			y = encode_read_id(file_id, relabelling[file_id][read_id]);
		}
	}

	if (verbosity_level >= 2)
		cout << "No. valid reads: " + to_string(no_valid_reads) + " of " + to_string(no_sample_reads) + " reads in sample\n";
}

// *********************************************************************************************
void CBarcodedCounter::list_cbc_dict(const string &suffix)
{
	if (!export_cbc_logs)
		return;

	ofstream ofs(cbc_log_file_name + suffix);

	uint64_t cum_sum = 0;

	ofs << "No. CBCs: " << cbc_vec.size() << endl;
//	ofs << "No. non-unique CBCs: " << count_if(cbc_vec.begin(), cbc_vec.end(), [](auto x) {return x.first > 1; }) << endl << endl;

	for (auto& x : cbc_vec)
	{
		ofs << base_coding4.decode_bases_2b(x.second, cbc_len) + " " + to_string(x.first) + "   cum: " + to_string(cum_sum) + "\n";
		cum_sum += x.first;
	}

	ofs << "Total no. of reads: " << cum_sum << endl;
}

// *********************************************************************************************
// Warning: p should be before q
double CBarcodedCounter::calc_dist(vector<pair<uint64_t, uint64_t>>::iterator p, vector<pair<uint64_t, uint64_t>>::iterator q)
{
	double dx = (double)(q - p);
	double dy = (double)(p->first) - (double)(q->first);

	return sqrt(dx * dx + dy * dy);
}

// *********************************************************************************************
int CBarcodedCounter::find_split(vector<pair<uint64_t, uint64_t>>& arr)
{
	size_t size = arr.size();
	double best_area = 0;
	int best_split = 0;

	if (size < 3)
		return 0;

	double dist_ac = calc_dist(arr.begin(), arr.begin() + size - 1);

	for (int i = 1; i < (int) size-1; ++i)
	{
		double dist_ab = calc_dist(arr.begin(), arr.begin() + i);
		double dist_bc = calc_dist(arr.begin() + i, arr.begin() + size - 1);

		double sum_div2 = (dist_ab + dist_ac + dist_bc) / 2;

		double area = sqrt(sum_div2 * (sum_div2 - dist_ab) * (sum_div2 - dist_ac) * (sum_div2 - dist_bc));
		if (area > best_area)
		{
			best_area = area;
			best_split = i;
		}
	}

	return best_split;
}

// *********************************************************************************************
void CBarcodedCounter::find_trusted_thr()
{
	cbc_vec.reserve(cbc_stats.size());

	for (auto& x : cbc_stats)
		cbc_vec.emplace_back(x.second, x.first);

	cbc_stats.clear();

	stable_sort(cbc_vec.begin(), cbc_vec.end(), greater<pair<uint64_t, cbc_t>>());

	vector<pair<uint64_t, cbc_t>> cbc_sum(cbc_vec);

	uint64_t tot = 0;

	// Calculate cummulative sums
	for (int i = 0; i < (int) cbc_sum.size(); ++i)
	{
		cbc_sum[i].first += tot;
		tot = cbc_sum[i].first;
	}

	int best_split = 0;

	if (cbc_filtering_thr == 0)
	{
		const int no_iters = 100;
		for (int i = 0; i < no_iters; ++i)
		{
			int curr_split = find_split(cbc_sum);
			if (curr_split == best_split)
				break;

			best_split = curr_split;

			if (best_split * 3 < (int)cbc_sum.size())
				cbc_sum.resize(best_split * 3);
		}
	}
	else
	{
		best_split = -1;
		for (size_t i = 0; i < cbc_vec.size(); ++i)
		{
			if (cbc_vec[i].first < cbc_filtering_thr)
			{
				best_split = i;
				break;
			}
		}

		if (best_split < 0)
			best_split = cbc_vec.size();
	}

	if (best_split <= 0)
		best_split = 1;

	if (apply_cbc_correction)
		cbc_for_corr_vec.assign(cbc_vec.begin() + best_split, cbc_vec.end());

	if (verbosity_level >= 2)
	{
		std::cerr << "No. of unique CBCs: " + to_string(cbc_vec.size()) + "\n";
		std::cerr << "No. of CBCs after filtering: " + to_string(best_split) + "\n";
	}

	cbc_vec.resize(best_split);
}

// *********************************************************************************************
void CBarcodedCounter::find_predefined_cbc()
{
	cbc_vec.reserve(predefined_cbc.size());

	for (auto& x : cbc_stats)
		if (predefined_cbc.count(x.first))
			cbc_vec.emplace_back(x.second, x.first);

	cbc_stats.clear();

	stable_sort(cbc_vec.begin(), cbc_vec.end(), greater<pair<uint64_t, cbc_t>>());

	if (verbosity_level >= 2)
	{
		std::cerr << "No. of unique CBCs: " + to_string(cbc_vec.size()) + "\n";
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_leaders_from_read(uint8_t* bases, vector<leader_t>& kmer_leaders)
{
	CKmer leader(leader_len, kmer_mode_t::direct);
	CKmer follower(follower_len, kmer_mode_t::direct);

	int read_len = strlen((char*)bases);
	int follower_start_pos = leader_len + gap_len;

	if (leader_len + gap_len + follower_len > (uint32_t) read_len)
		return;

	for(uint32_t i = 0; i < leader_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			leader.insert(symbol);
		else
			leader.Reset();
	}

	for(uint32_t i = follower_start_pos; i < follower_start_pos + follower_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			follower.insert(symbol);
		else
			follower.Reset();
	}

	// leader and follower contain almost complete k-mers (without last symbols)

	for (int i = follower_start_pos + follower_len - 1; i < read_len; ++i)
	{
		uint64_t t_symbol = dna_code(bases[i]);
		uint64_t a_symbol = dna_code(bases[i - follower_len - gap_len]);

		if (t_symbol < 4)
			follower.insert(t_symbol);
		else
			follower.Reset();

		if (a_symbol < 4)
			leader.insert(a_symbol);
		else
			leader.Reset();

		if (leader.is_full() && follower.is_full())
			kmer_leaders.emplace_back(leader.data_dir());
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_pairs_from_read(uint8_t* bases, vector<leader_follower_t>& kmer_pairs)
{
	CKmer leader(leader_len, kmer_mode_t::direct);
	CKmer follower(follower_len, kmer_mode_t::direct);

	int read_len = (int) strlen((char*)bases);
	int follower_start_pos = leader_len + gap_len;

	if (leader_len + gap_len + follower_len > (uint32_t) read_len)
		return;

	for(uint32_t i = 0; i < leader_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			leader.insert(symbol);
		else
			leader.Reset();
	}

	for(uint32_t i = follower_start_pos; i < follower_start_pos + follower_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			follower.insert(symbol);
		else
			follower.Reset();
	}

	// leader and follower contain almost complete k-mers (without last symbols)

	for (int i = follower_start_pos + follower_len - 1; i < read_len; ++i)
	{
		uint64_t t_symbol = dna_code(bases[i]);
		uint64_t a_symbol = dna_code(bases[i - follower_len - gap_len]);

		if (t_symbol < 4)
			follower.insert(t_symbol);
		else
			follower.Reset();

		if (a_symbol < 4)
			leader.insert(a_symbol);
		else
			leader.Reset();

		if (leader.is_full() && follower.is_full() && (min_leader_count <= 1 || valid_leaders.count(leader.data())))
			kmer_pairs.emplace_back(leader.data_aligned_dir(), follower.data_aligned_dir());
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmers_from_read(uint8_t* bases, vector<kmer_t>& kmers)
{
//	CKmer kmer(leader_len, kmer_mode_t::direct);			// !!! TODO - add support for canonical
	CKmer kmer(leader_len, canonical_mode ? kmer_mode_t::canonical : kmer_mode_t::direct);

	int read_len = (int) strlen((char*)bases);
	
	if (leader_len > (uint32_t) read_len)
		return;

	for(uint32_t i = 0; i < leader_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			kmer.insert(symbol);
		else
			kmer.Reset();
	}

	// kmer contains almost complete k-mer (without last symbol)

	for (int i = leader_len - 1; i < read_len; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);

		if (symbol < 4)
			kmer.insert(symbol);
		else
			kmer.Reset();

		if (kmer.is_full() && (min_leader_count <= 1 || valid_leaders.count(kmer.data())))
//			kmers.emplace_back(kmer.data_aligned_dir());
			kmers.emplace_back(kmer.data_aligned());
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_leaders_for_cbc(cbc_t cbc, vector<leader_t>& kmer_leaders)
{
	kmer_leaders.clear();

	uint64_t file_id;
	uint64_t read_id;

#ifdef USE_READ_COMPRESSION
	vector<uint8_t> decompressed_read;
#endif

	for (auto x : global_cbc_dict[cbc])
	{
		tie(file_id, read_id) = decode_read_id(x);

#ifdef USE_READ_COMPRESSION
		base_coding3.decode_bases(sample_reads[file_id][read_id], decompressed_read);
		enumerate_kmer_leaders_from_read(decompressed_read.data(), kmer_leaders);
#else
		enumerate_kmer_leaders_from_read(sample_reads[file_id][read_id], kmer_leaders);
#endif
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_pairs_for_cbc(cbc_t cbc, vector<leader_follower_t>& kmer_pairs)
{
	kmer_pairs.clear();

	uint64_t file_id;
	uint64_t read_id;

#ifdef USE_READ_COMPRESSION
	vector<uint8_t> decompressed_read;
#endif

	for (auto x : global_cbc_dict[cbc])
	{
		tie(file_id, read_id) = decode_read_id(x);

#ifdef USE_READ_COMPRESSION
		base_coding3.decode_bases(sample_reads[file_id][read_id], decompressed_read);

		enumerate_kmer_pairs_from_read(decompressed_read.data(), kmer_pairs);
#else
		enumerate_kmer_pairs_from_read(sample_reads[file_id][read_id], kmer_pairs);
#endif
	}

#ifdef AGGRESIVE_MEMORY_SAVING
	kmer_pairs.shrink_to_fit();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmers_for_cbc(cbc_t cbc, vector<kmer_t>& kmers)
{
	kmers.clear();

	uint64_t file_id;
	uint64_t read_id;

#ifdef USE_READ_COMPRESSION
	vector<uint8_t> decompressed_read;
#endif

	for (auto x : global_cbc_dict[cbc])
	{
		tie(file_id, read_id) = decode_read_id(x);

#ifdef USE_READ_COMPRESSION
		base_coding3.decode_bases(sample_reads[file_id][read_id], decompressed_read);

		enumerate_kmers_from_read(decompressed_read.data(), kmers);
#else
		enumerate_kmers_from_read(sample_reads[file_id][read_id], kmers);
#endif
	}

#ifdef AGGRESIVE_MEMORY_SAVING
	kmer_pairs.shrink_to_fit();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::sort_and_gather_kmer_pairs_for_cbc(vector<leader_follower_t>& kmer_pairs, vector<leader_follower_count_t>& kmer_pair_counts)
{
//	std::sort(kmer_pairs.begin(), kmer_pairs.end());
	refresh::sort::pdqsort(kmer_pairs.begin(), kmer_pairs.end());

	kmer_pair_counts.clear();

	if (kmer_pairs.empty())
		return;

	kmer_pair_counts.emplace_back(kmer_pairs.front());

	for (int i = 1; i < (int) kmer_pairs.size(); ++i)
		if (kmer_pair_counts.back().equal_lf(kmer_pairs[i]))
			kmer_pair_counts.back().count++;
		else
		{
			if (poly_ACGT_filter.IsPolyACGT(kmer_pair_counts.back().leader, leader_len) ||
					artifacts_filter.ContainsArtifact(kmer_pair_counts.back().leader, leader_len))
				kmer_pair_counts.pop_back();
			kmer_pair_counts.emplace_back(kmer_pairs[i]);
		}
	if (poly_ACGT_filter.IsPolyACGT(kmer_pair_counts.back().leader, leader_len) ||
			artifacts_filter.ContainsArtifact(kmer_pair_counts.back().leader, leader_len))
		kmer_pair_counts.pop_back();

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmer_pairs);
	kmer_pair_counts.shrink_to_fit();
#else
	kmer_pairs.clear();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::sort_and_gather_kmers_for_cbc(vector<kmer_t>& kmers, vector<kmer_count_t>& kmer_counts)
{
//	std::sort(kmer_pairs.begin(), kmer_pairs.end());
	refresh::sort::pdqsort(kmers.begin(), kmers.end());

	kmer_counts.clear();

	if (kmers.empty())
		return;

	kmer_counts.emplace_back(kmers.front());

	for (int i = 1; i < (int) kmers.size(); ++i)
		if (kmer_counts.back().equal_lf(kmers[i]))
			kmer_counts.back().count++;
		else
		{
			if (poly_ACGT_filter.IsPolyACGT(kmer_counts.back().kmer, leader_len) ||
					artifacts_filter.ContainsArtifact(kmer_counts.back().kmer, leader_len))
				kmer_counts.pop_back();
			kmer_counts.emplace_back(kmers[i]);
		}
	if (poly_ACGT_filter.IsPolyACGT(kmer_counts.back().kmer, leader_len) ||
			artifacts_filter.ContainsArtifact(kmer_counts.back().kmer, leader_len))
		kmer_counts.pop_back();

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmers);
	kmer_counts.shrink_to_fit();
#else
	kmers.clear();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::filter_rare_leader_sample_cbc(vector<leader_follower_count_t>& kmer_pair_counts)
{
	if (rare_leader_thr < 1)
		return;

	auto p_leader_begin = kmer_pair_counts.begin();
	auto p_end = kmer_pair_counts.end();
	auto p_dest = p_leader_begin;

	while (p_leader_begin != p_end)
	{
		auto leader = p_leader_begin->leader;
		auto p = p_leader_begin;
		uint64_t count = 0;

		for (; p != p_end && p->leader == leader; ++p)
			count += p->count;

		if (count > rare_leader_thr)
		{
			if (p_leader_begin == p_dest)
				p_dest = p;
			else
				p_dest = copy(p_leader_begin, p, p_dest);
		}

		p_leader_begin = p;
	}

	kmer_pair_counts.erase(p_dest, p_end);

#ifdef AGGRESIVE_MEMORY_SAVING
	kmer_pair_counts.shrink_to_fit();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::filter_rare_kmer_sample_cbc(vector<kmer_count_t>& kmer_counts)
{
	if (rare_leader_thr < 1)
		return;

	auto p_leader_begin = kmer_counts.begin();
	auto p_end = kmer_counts.end();
	auto p_dest = p_leader_begin;

	while (p_leader_begin != p_end)
	{
		auto leader = p_leader_begin->kmer;
		auto p = p_leader_begin;
		uint64_t count = 0;

		for (; p != p_end && p->kmer == leader; ++p)
			count += p->count;

		if (count > rare_leader_thr)
		{
			if (p_leader_begin == p_dest)
				p_dest = p;
			else
				p_dest = copy(p_leader_begin, p, p_dest);
		}

		p_leader_begin = p;
	}

	kmer_counts.erase(p_dest, p_end);

#ifdef AGGRESIVE_MEMORY_SAVING
	kmer_counts.shrink_to_fit();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::count_leaders()
{
	atomic_int id{ 0 };
	atomic_uint64_t total_no_after_removal{ 0 };

	vector<thread> threads;
	vector<cbc_t> cbcs;

	total_no_kmer_leaders_counts = 0;
	sum_kmer_leaders_counts = 0;

	cbcs.reserve(global_cbc_umi_dict.size());
	for (auto& x : global_cbc_umi_dict)
		cbcs.emplace_back(x.first);

	leader_counts.resize(no_threads);

	threads.reserve(no_threads);

	for (int i = 0; i < no_threads; ++i)
		threads.emplace_back([&, i] {
		int thread_id = i;
		int curr_id = -1;

		vector<leader_t> kmer_leaders;

		auto& my_leader_counts = leader_counts[thread_id];

		while (true)
		{
			std::cerr << to_string(curr_id) + "\n";

			curr_id = id.fetch_add(1);
			if (curr_id >= (int)cbcs.size())
				break;

			enumerate_kmer_leaders_for_cbc(cbcs[curr_id], kmer_leaders);
			for (auto x : kmer_leaders)
				my_leader_counts[x] += 1;
		}
			});

	join_threads(threads);

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. k-mer leaders counts: " << total_no_kmer_leaders_counts << endl;
		std::cerr << "Sum of k-mer leaders counts: " << sum_kmer_leaders_counts << endl;
	}
}

// *********************************************************************************************
void CBarcodedCounter::determine_valid_leaders()
{
	// !!! Consider parallelization

	uint64_t total_no_considered_leaders_threads = 0;

	for (auto& x : leader_counts)
		total_no_considered_leaders_threads += x.size();

	for (int i = 1; i < (int)leader_counts.size(); ++i)
		for (auto x : leader_counts[i])
			leader_counts[0][x.first] += x.second;

	for (auto x : leader_counts[0])
		if (x.second >= min_leader_count)
			valid_leaders.emplace(x.first);

	clear_vec(leader_counts);

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. of considered leaders in threads: " << total_no_considered_leaders_threads << endl;
		std::cerr << "No. valid leaders: " << valid_leaders.size() << endl;
	}
}

// *********************************************************************************************
void CBarcodedCounter::count_kmer_pairs()
{
	atomic_int id{ 0 };
	atomic_uint64_t total_no_after_removal{ 0 };

	vector<thread> threads;
	vector<cbc_t> cbcs;
	
	total_no_kmer_pair_counts = 0;
	sum_kmer_pair_counts = 0;

	cbcs.reserve(global_cbc_dict.size());
	for (auto& x : global_cbc_dict)
		cbcs.emplace_back(x.first);

	threads.reserve(no_threads);

	for (int i = 0; i < no_threads; ++i)
		threads.emplace_back([&] {
		int curr_id = -1;

		vector<leader_follower_t> kmer_pairs;
		vector<leader_follower_count_t> kmer_pair_counts;

		vector<vector<bkc_record_t>> record_buffers;

		vector<uint8_t> packed_buffer;

		record_buffers.resize(no_splits);

//		zstd_in_memory zim{ (int) zstd_level };
//		vector<uint8_t> zstd_working_space;

		while (true)
		{
			curr_id = id.fetch_add(1);
			if (curr_id >= (int)cbcs.size())
				break;

			enumerate_kmer_pairs_for_cbc(cbcs[curr_id], kmer_pairs);
			sort_and_gather_kmer_pairs_for_cbc(kmer_pairs, kmer_pair_counts);
			filter_rare_leader_sample_cbc(kmer_pair_counts);
			store_kmer_pairs(cbcs[curr_id], kmer_pair_counts, record_buffers);

			for (uint32_t i = 0; i < no_splits; ++i)
				if ((int) record_buffers[i].size() >= max_records_in_buffer)
				{
					pack_records(record_buffers[i], packed_buffer);
					bkc_files[i]->AddPacked(packed_buffer);

/*					zstd_working_space.resize(packed_buffer.size() + zim.get_overhead(packed_buffer.size()));
					auto packed_size = zim.compress(packed_buffer.data(), packed_buffer.size(), zstd_working_space.data(), zstd_working_space.size(), (int) zstd_level);
					zstd_working_space.resize(packed_size);
					bkc_files[i]->AddPacked(zstd_working_space);*/

					record_buffers[i].clear();
				}
		}

		for (uint32_t i = 0; i < no_splits; ++i)
		{
			pack_records(record_buffers[i], packed_buffer);
			bkc_files[i]->AddPacked(packed_buffer);
		}
		});

	join_threads(threads);

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. k-mer pair counts: " << total_no_kmer_pair_counts << endl;
		std::cerr << "Sum of k-mer pair counts: " << sum_kmer_pair_counts << endl;
	}
}

// *********************************************************************************************
void CBarcodedCounter::count_kmers()
{
	atomic_int id{ 0 };
	atomic_uint64_t total_no_after_removal{ 0 };

	vector<thread> threads;
	vector<cbc_t> cbcs;
	
	total_no_kmer_counts = 0;
	sum_kmer_counts = 0;

	cbcs.reserve(global_cbc_dict.size());
	for (auto& x : global_cbc_dict)
		cbcs.emplace_back(x.first);

	threads.reserve(no_threads);

	for (int i = 0; i < no_threads; ++i)
		threads.emplace_back([&] {
		int curr_id = -1;

		vector<kmer_t> kmers;
		vector<kmer_count_t> kmer_counts;

		vector<vector<bkc_record_t>> record_buffers;

		vector<uint8_t> packed_buffer;

		record_buffers.resize(no_splits);

//		zstd_in_memory zim{ (int) zstd_level };
//		vector<uint8_t> zstd_working_space;

		while (true)
		{
			curr_id = id.fetch_add(1);
			if (curr_id >= (int)cbcs.size())
				break;

			enumerate_kmers_for_cbc(cbcs[curr_id], kmers);
			sort_and_gather_kmers_for_cbc(kmers, kmer_counts);
			filter_rare_kmer_sample_cbc(kmer_counts);
			store_kmers(cbcs[curr_id], kmer_counts, record_buffers);

			for (uint32_t i = 0; i < no_splits; ++i)
				if ((int) record_buffers[i].size() >= max_records_in_buffer)
				{
					pack_records(record_buffers[i], packed_buffer);
					bkc_files[i]->AddPacked(packed_buffer);

/*					zstd_working_space.resize(packed_buffer.size() + zim.get_overhead(packed_buffer.size()));
					auto packed_size = zim.compress(packed_buffer.data(), packed_buffer.size(), zstd_working_space.data(), zstd_working_space.size(), (int) zstd_level);
					zstd_working_space.resize(packed_size);
					bkc_files[i]->AddPacked(zstd_working_space);*/

					record_buffers[i].clear();
				}
		}

		for (uint32_t i = 0; i < no_splits; ++i)
		{
			pack_records(record_buffers[i], packed_buffer);
			bkc_files[i]->AddPacked(packed_buffer);
		}
		});

	join_threads(threads);

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. k-mers: " << total_no_kmer_counts << endl;
		std::cerr << "Sum of k-mer counts: " << sum_kmer_counts << endl;
	}
}

// *********************************************************************************************
void CBarcodedCounter::pack_records(vector<bkc_record_t>& records, vector<uint8_t>& packed_buffer)
{
	vector<uint8_t> rec_prev, rec_curr;

#ifdef COMPACT_ENCODING
	packed_buffer.clear();

	for (auto& x : records)
	{
		rec_curr.clear();

		append_int_msb(rec_curr, x.sample_id, sample_id_size_in_bytes);
		append_int_msb(rec_curr, x.barcode, barcode_size_in_bytes);
		append_int_msb(rec_curr, x.leader, leader_size_in_bytes);
		append_int_msb(rec_curr, x.follower, follower_size_in_bytes);
		append_int_msb(rec_curr, x.count, counter_size_in_bytes);

		encode_shared_prefix(packed_buffer, rec_prev, rec_curr);

		swap(rec_prev, rec_curr);
	}
#else
	save_int_lsb(sample_id, sample_id_size_in_bytes);
	save_int_lsb(barcode, barcode_size_in_bytes);
	save_int_lsb(leader, leader_size_in_bytes);
	save_int_lsb(follower, follower_size_in_bytes);
	save_int_lsb(count, counter_size_in_bytes);
#endif
}

// *********************************************************************************************
void CBarcodedCounter::store_kmer_pairs(cbc_t cbc, vector<leader_follower_count_t>& kmer_pair_counts, vector<vector<bkc_record_t>>& record_buffers)
{
	total_no_kmer_pair_counts += kmer_pair_counts.size();

	string cbc_str = base_coding4.decode_bases_2b(cbc, cbc_len);

	refresh::MurMur64Hash mh;

	uint64_t sum = 0;
	for (const auto& x : kmer_pair_counts)
	{
		uint64_t h = mh(x.leader) % no_splits;

		record_buffers[h].emplace_back(sample_id, cbc, x.leader, x.follower, x.count);
		sum += x.count;
	}
	
	sum_kmer_pair_counts += sum;

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmer_pair_counts);
#else
	kmer_pair_counts.clear();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::store_kmers(cbc_t cbc, vector<kmer_count_t>& kmer_counts, vector<vector<bkc_record_t>>& record_buffers)
{
	total_no_kmer_counts += kmer_counts.size();

	string cbc_str = base_coding4.decode_bases_2b(cbc, cbc_len);

	refresh::MurMur64Hash mh;

	uint64_t sum = 0;
	for (const auto& x : kmer_counts)
	{
		uint64_t h = mh(x.kmer) % no_splits;

		record_buffers[h].emplace_back(sample_id, cbc, x.kmer, 0, x.count);
		sum += x.count;
	}
	
	sum_kmer_pair_counts += sum;

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmer_counts);
#else
	kmer_counts.clear();
#endif
}

// *********************************************************************************************
string CBarcodedCounter::kmer_to_string(uint64_t kmer, int len)
{
	string str;

	for (int i = 0; i < len; ++i)
	{
		auto c = kmer & 3;
		str.push_back("ACGT"[c]);
		kmer >>= 2;
	}

	reverse(str.begin(), str.end());

	return str;
}

// *********************************************************************************************
bool CBarcodedCounter::ProcessCBC()
{
	set_CBC_file_names();

	times.emplace_back("", high_resolution_clock::now());

	if (!no_threads || file_names.empty())
		return false;

	no_reading_threads = max(min(no_threads / 2, (int)file_names.size()), 1);

	init_queues_and_pools();

	start_reading_threads();
	start_counting_threads();

	join_threads(reading_threads);
	join_threads(counting_threads);

	times.emplace_back("Reading and counting", high_resolution_clock::now());

	if (verbosity_level >= 1)
		std::cerr << "Gathering CBC statistics\n";
	gather_cbc_stats();
	mi_collect(true);
	times.emplace_back("Gathering CBC statistics", high_resolution_clock::now());

	if (!predefined_cbc.empty())
	{
		if (verbosity_level >= 1)
			std::cerr << "Looking for predefined CBC in data\n";

		find_predefined_cbc();

		mi_collect(true);
		times.emplace_back("Looking for predefined CBC in data", high_resolution_clock::now());
	}
	else
	{
		if (verbosity_level >= 1)
			std::cerr << "Looking for trusted threshold\n";
		find_trusted_thr();		// !!! This can be parallelized
		mi_collect(true);
		times.emplace_back("Looking for trusted threshold", high_resolution_clock::now());

		if (apply_cbc_correction)
		{
			if (verbosity_level >= 1)
				std::cerr << "CBCs correction\n";
			find_CBC_corrections();
			mi_collect(true);
			times.emplace_back("CBCs correction", high_resolution_clock::now());
		}
	}

	if (verbosity_level >= 1)
		std::cerr << "Removing non-trusted CBCs\n";
	remove_non_trusted_CBC();
	mi_collect(true);
	list_cbc_dict(".before_umi_deduplication");
	times.emplace_back("Removing non-trusted CBCs", high_resolution_clock::now());

	if (verbosity_level >= 1)
		std::cerr << "Removing duplicated UMIs\n";
	remove_duplicated_UMI();
	list_cbc_dict(".after_umi_deduplication");
	mi_collect(true);
	times.emplace_back("Removing duplicated UMIs", high_resolution_clock::now());

	if (verbosity_level >= 1)
		std::cerr << "Creating valid reads list\n";
	create_valid_reads_lists();
	mi_collect(true);
	times.emplace_back("Creating valid reads list", high_resolution_clock::now());

	return true;
}

// *********************************************************************************************
bool CBarcodedCounter::ProcessExportFilteredCBCReads()
{
	set_CBC_file_names();

	if (!no_threads || file_names.empty())
		return false;

	no_reading_threads = max(min(no_threads / 2, (int)file_names.size()), 1);

	if (verbosity_level >= 1)
		std::cerr << "Reads loading\n";

	reinit_queues();

	init_bkc_files();

	start_reading_threads();
	start_reads_exporting_threads();

	join_threads(reading_threads);
	join_threads(reads_exporting_threads);
	mi_collect(true);

	times.emplace_back("CBC reads filtering", high_resolution_clock::now());

	return true;
}

// *********************************************************************************************
bool CBarcodedCounter::ProcessExportFilteredReads()
{
	set_read_file_names();

	if (!no_threads || file_names.empty())
		return false;

	no_reading_threads = max(min(no_threads / 2, (int)file_names.size()), 1);

	if (verbosity_level >= 1)
		std::cerr << "Reads loading\n";

	reinit_queues();

	init_bkc_files();

	start_reading_threads();
	start_reads_exporting_threads();

	join_threads(reading_threads);
	join_threads(reads_exporting_threads);
	mi_collect(true);

	times.emplace_back("CBC reads filtering", high_resolution_clock::now());

	return true;
}

// *********************************************************************************************
bool CBarcodedCounter::ProcessReads()
{
	set_read_file_names();

	if (!no_threads || file_names.empty())
		return false;

	no_reading_threads = max(min(no_threads / 2, (int)file_names.size()), 1);

	if (verbosity_level >= 1)
		std::cerr << "Reads loading\n";

	reinit_queues();

	init_bkc_files();

	start_reading_threads();
	start_reads_loading_threads();

	join_threads(reading_threads);
	join_threads(reads_loading_threads);
	mi_collect(true);

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. of loaded reads: " << a_total_no_reads << endl;
		std::cerr << "Total len of loaded reads: " << a_total_read_len << endl;
	}

	times.emplace_back("Reads loading", high_resolution_clock::now());

#if 0			// Currently not used
	if (min_leader_count > 1)
	{
		if (verbosity_level >= 1)
			std::cerr << "Enumerating and counting leader k-mers\n";
		count_leaders();
		times.emplace_back("Enumerating and counting leader k-mers", high_resolution_clock::now());

		if (verbosity_level >= 1)
			std::cerr << "Determine valid leaders\n";
		determine_valid_leaders();
		times.emplace_back("Determine valid leaders", high_resolution_clock::now());
	}
#endif

	if (counting_mode == counting_mode_t::single)
	{
		if (verbosity_level >= 1)
			std::cerr << "Enumerating and counting k-mers\n";
		count_kmers();
		mi_collect(true);
		times.emplace_back("Enumerating and counting k-mers", high_resolution_clock::now());
	}
	else
	{
		if (verbosity_level >= 1)
			std::cerr << "Enumerating and counting leader-follower pairs\n";
		count_kmer_pairs();
		mi_collect(true);
		times.emplace_back("Enumerating and counting leader-follower pairs", high_resolution_clock::now());
	}

	mma.clear();

	bkc_files.clear();

	return true;
}

// EOF
