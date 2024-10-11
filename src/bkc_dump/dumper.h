#pragma once
#include "params.h"

class CDumper
{
	CParams params;

//	bool open_file()

public:
	CDumper() = default;

	bool SetParams(const CParams& _params);
	bool Dump();
};