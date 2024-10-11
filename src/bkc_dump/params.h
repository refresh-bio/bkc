#pragma once

#include "../common/defs.h"

// *********************************************************************************************
struct CParams
{
	param_t<uint32_t> no_splits{ 1, 256, 1 };
	param_t<uint32_t> no_threads{ 0, 256, 8 };
	string input_file_name{  };
	string output_file_name{ "./results.txt" };
};

