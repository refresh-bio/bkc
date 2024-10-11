#include <iostream>
#include "params.h"
#include "dumper.h"

void usage();
bool parse_args(int argc, char** argv);


CParams params;


// *********************************************************************************************
void usage()
{
    cerr
        << "BKC-dump: dumping k-mers or k-mer pairs from BKC files (v.1.0.0)\n";
    cerr
        << "Usage:\n"
        << "    bxc-dump [options]\n"
        << "    --input_name <file_name> - BKC file name\n"
//        << "    --n_threads <int> - no. threads " << params.no_threads.str() << endl
        << "    --output_name <file_name> - output file name (default: " << params.output_file_name << ")\n"
        << "    --n_splits <int> - no. splits " << params.no_splits.str() << endl;
}

// *********************************************************************************************
bool parse_args(int argc, char** argv)
{
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == "--n_splits"s && i + 1 < argc)
		{
			if (!params.no_splits.set(atoi(argv[++i])))
			{
				cerr << "Incorrect value for n_splits: " << argv[i] << endl;
				return false;
			}
		}
/*		else if (argv[i] == "--n_threads"s && i + 1 < argc)
		{
			if (!params.no_threads.set(atoi(argv[++i])))
			{
				cerr << "Incorrect value for n_threads: " << argv[i] << endl;
				return false;
			}
		}*/
		else if (argv[i] == "--input_name"s && i + 1 < argc)
		{
			params.input_file_name = argv[++i];
		}
		else if (argv[i] == "--output_name"s && i + 1 < argc)
		{
			params.output_file_name = argv[++i];
		}
		else
		{
			cerr << "Unknown parameter: " << argv[i] << endl;
			return false;
		}
	}

	if (params.input_file_name.empty())
		return false;

	return true;
}

// *********************************************************************************************
int main(int argc, char **argv)
{
	if (!parse_args(argc, argv))
	{
		usage();
		return 1;
	}

	CDumper dumper;

	dumper.SetParams(params);
	dumper.Dump();

    return 0;
}

