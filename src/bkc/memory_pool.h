#pragma once

#include <mutex>
#include <condition_variable>
#include <vector>
#include <stack>
#include <refresh/memory_chunk/lib/memory_chunk.h>

using namespace std;
using namespace refresh;

template<typename T> class CMemoryPool
{
	size_t max_no_chunks;
	size_t chunk_size;

	stack<memory_chunk<T>, vector<memory_chunk<T>>> chunks;

	mutex mtx;
	condition_variable cv;

	void deallocate()
	{
		while (!chunks.empty())
		{
			delete[] chunks.top().data();
			chunks.pop();
		}
	}

	void allocate()
	{
		for (size_t i = 0; i < max_no_chunks; ++i)
			chunks.push(memory_chunk<T>(new T[chunk_size], chunk_size));
	}

public:
	CMemoryPool() :
		max_no_chunks(0),
		chunk_size(0)
	{}

	CMemoryPool(size_t _max_no_chunks, size_t _chunk_size) :
		max_no_chunks(_max_no_chunks),
		chunk_size(_chunk_size)
	{
		allocate();
	}

	CMemoryPool(const CMemoryPool<T>&) = delete;
	
	CMemoryPool(CMemoryPool<T>&& x)
	{
		lock_guard<mutex> lck(x.mtx);

		chunks = move(x.chunks);
		max_no_chunks = x.max_no_chunks;
		chunk_size = x.chunk_size;

		x.max_no_chunks = 0;
		x.chunk_size = 0;
	}

	~CMemoryPool()
	{
		deallocate();
	}

	void Resize(size_t _max_no_chunks, size_t _chunk_size)
	{
		lock_guard<mutex> lck(mtx);

		deallocate();

		max_no_chunks = _max_no_chunks;
		chunk_size = _chunk_size;

		allocate();
	}

	void Clear()
	{
		lock_guard<mutex> lck(mtx);

		deallocate();
	}

	size_t Capacity()
	{
		lock_guard<mutex> lck(mtx);

		return max_no_chunks;
	}

	size_t AvailableChunks()
	{
		lock_guard<mutex> lck(mtx);

		return chunks.size();
	}

	void Pop(memory_chunk<T>& mc)
	{
		unique_lock<mutex> lck(mtx);

		cv.wait(lck, [&] {return !chunks.empty(); });

		mc = move(chunks.top());
		chunks.pop();
	}

	void Push(memory_chunk<T>& mc)
	{
		lock_guard<mutex> lck(mtx);

		chunks.push(move(mc));

		cv.notify_one();
	}
};

// EOF