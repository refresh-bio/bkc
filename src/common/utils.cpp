#include "utils.h"

void show_timings(vector<pair<string, time_point<high_resolution_clock>>>& timings)
{
	cerr << "***** Timings\n";

	for (size_t i = 1; i < timings.size(); ++i)
	{
		duration<double> diff = timings[i].second - timings[i-1].second;

		cerr << timings[i].first << ": " << diff.count() << "s\n";
	}

	duration<double> diff = timings.back().second - timings.front().second;

	cerr << "Total: " << diff.count() << "s\n";
}

void encode_shared_prefix(vector<uint8_t> &v_data, vector<uint8_t>& v_prev, vector<uint8_t>& v_curr)
{
	size_t i;
	size_t len = min(v_prev.size(), v_curr.size());

	for (i = 0; i < len; ++i)
		if (v_curr[i] != v_prev[i])
			break;

	v_data.emplace_back((uint8_t)i);
	v_data.insert(v_data.end(), v_curr.begin() + i, v_curr.end());
}
