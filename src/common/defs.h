#pragma once

#include <cinttypes>
#include <sstream>
#include <string>

using namespace std;

// *********************************************************************************************
enum class technology_t { unknown, ten_x, visium };
enum class counting_mode_t { unknown, single, pair };
enum class output_format_t { unknown, bkc, splash };
enum class export_filtered_input_t { none = 0, first = 1, second = 2, both = 3 };

const string BKC_VERSION = "1.0.0";
const string BKC_DATE = "2024-10-11";

using kmer_t = uint64_t;
using leader_t = uint64_t;
using follower_t = uint64_t;
using cbc_t = uint64_t;

// *********************************************************************************************
template<typename T> class param_t
{
	T min_value;
	T max_value;
	T def_value;
	T value;

public:
	param_t(T min_value, T max_value, T def_value) :
		min_value(min_value), max_value(max_value), def_value(def_value), value(def_value)
	{}

	bool set(T new_value)
	{
		if (new_value < min_value || new_value > max_value)
			return false;
		value = new_value;
		return true;
	}

	T get()		const
	{
		return value;
	}

	string str()	const
	{
		stringstream s;
		s << "(default: " << def_value << ", min: " << min_value << ", max: " << max_value << ")";

		return s.str();
	}
};

