/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_OBJECT_POOL_H_
#define INFOMAP_OBJECT_POOL_H_

#include <cassert>
#include <cstddef>
#include <memory>
#include <new>
#include <utility>
#include <vector>

namespace infomap {

/**
 * Address-stable object pool.
 *
 * Storage is a vector of raw byte chunks of geometrically growing capacity.
 * Chunks are never reallocated, so any pointer returned by alloc() stays valid
 * until that object is free()d (or the pool is destroyed). Objects are
 * placement-constructed on alloc and their destructor is run on free; freed
 * slots are recycled via an explicit free list.
 *
 * Teardown contract: by the time the pool is destroyed, every object must have
 * been free()d already (the owner reclaims its whole object graph first). The
 * destructor therefore only releases raw memory and asserts the pool is empty;
 * it does NOT run object destructors. A non-empty pool at teardown is a leak
 * (caught by the assert in debug and by LeakSanitizer otherwise).
 *
 * Thread safety: NOT thread-safe. Infomap relies on this: each InfomapBase owns
 * its own pools, which gives every parallel partition/trial worker a private,
 * lock-free allocator. Do NOT share one pool across instances that may run
 * concurrently (e.g. parallel sub-Infomap recursion) — concurrent alloc/free
 * races on m_chunks/m_free and corrupts the heap. Making it shareable would
 * require locking every alloc/free, serializing the hot path and negating the
 * point of the per-instance design.
 */
template <typename T>
class ObjectPool {
  // Chunk storage comes from operator new[], which only guarantees alignment
  // for ordinary types. The slot arithmetic below relies on that guarantee.
  static_assert(alignof(T) <= alignof(std::max_align_t),
                "ObjectPool does not support over-aligned types");

  // A small initial chunk keeps the many tiny per-module sub-Infomap instances
  // (each with its own node and edge pool) close to their real footprint —
  // recursive partitioning spawns one pool pair per refined module, most of
  // them only a handful of nodes. A larger cap keeps the big main-instance pool
  // fast by amortizing chunk allocation. Tuned against test/bench/native_benchmark:
  // matches or beats pre-pool peak RSS while improving large-graph runtime.
  static constexpr std::size_t kInitialChunk = 4;
  static constexpr std::size_t kMaxChunk = 1024;

  struct Chunk {
    std::unique_ptr<unsigned char[]> storage;
    std::size_t capacity = 0;
    std::size_t used = 0;

    explicit Chunk(std::size_t cap)
        : storage(new unsigned char[cap * sizeof(T)]), capacity(cap) {}

    T* slot(std::size_t i) noexcept
    {
      return reinterpret_cast<T*>(storage.get() + i * sizeof(T));
    }
  };

  std::vector<Chunk> m_chunks;
  std::vector<T*> m_free;
  std::size_t m_liveCount = 0;

  T* nextRawSlot()
  {
    if (m_chunks.empty() || m_chunks.back().used == m_chunks.back().capacity) {
      std::size_t cap = m_chunks.empty()
          ? kInitialChunk
          : (m_chunks.back().capacity < kMaxChunk ? m_chunks.back().capacity * 2 : kMaxChunk);
      m_chunks.emplace_back(cap);
    }
    Chunk& c = m_chunks.back();
    return c.slot(c.used++);
  }

public:
  ObjectPool() = default;
  ObjectPool(const ObjectPool&) = delete;
  ObjectPool& operator=(const ObjectPool&) = delete;
  ObjectPool(ObjectPool&&) = delete;
  ObjectPool& operator=(ObjectPool&&) = delete;

  ~ObjectPool()
  {
    // All objects must already be freed by the owner's teardown.
    assert(m_liveCount == 0 && "ObjectPool destroyed with live objects (leak)");
  }

  template <typename... Args>
  T* alloc(Args&&... args)
  {
    T* slot;
    if (!m_free.empty()) {
      slot = m_free.back();
      m_free.pop_back();
    } else {
      slot = nextRawSlot();
    }
    T* obj = ::new (static_cast<void*>(slot)) T(std::forward<Args>(args)...);
    ++m_liveCount;
    return obj;
  }

  void free(T* obj) noexcept
  {
    if (obj == nullptr)
      return;
    assert(m_liveCount > 0 && "ObjectPool::free without a matching alloc (double free?)");
    obj->~T();
    m_free.push_back(obj);
    --m_liveCount;
  }

  // Ensure room for n more objects without geometric ramp-up, for call sites
  // that know the bulk count up front (leaf nodes = network size, sub-network
  // clones = parent child-degree). One exactly-sized chunk replaces the
  // 4,8,16,... ramp + tail slack, which is the bulk of pool over-reservation —
  // most impactful for the swarm of small sub-Infomap pools.
  void reserve(std::size_t n)
  {
    std::size_t available = m_free.size();
    if (!m_chunks.empty())
      available += m_chunks.back().capacity - m_chunks.back().used;
    if (available >= n)
      return;
    // Appending a chunk strands the current back chunk's unused tail
    // (allocations only come from the back chunk), so size the new chunk
    // against the free list alone to actually guarantee room for n.
    m_chunks.emplace_back(n - m_free.size());
  }

  std::size_t liveCount() const noexcept { return m_liveCount; }
  std::size_t chunkCount() const noexcept { return m_chunks.size(); }
};

} // namespace infomap

#endif // INFOMAP_OBJECT_POOL_H_
