#ifndef CCEH_H_
#define CCEH_H_

#ifndef UTIL_PERSIST_H_
#include "util/persist.h"
#endif

#ifndef UTIL_RWLOCK_H_
#include "util/rwlock.h"
#endif

#ifndef UTIL_TIMER_H_
#include "util/timer.h"
#endif

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <shared_mutex>
#include <mutex>
#include <atomic>
#include <vector>

#define INPLACE
#define USE_AEP

#ifdef QUARTZ
#include "../../quartz/src/lib/pmalloc.h"
#endif

#ifdef USE_AEP
#ifndef AEP_H_
#include "aep.h"
#endif
#endif

#ifdef USE_LIBVMEMALLOC
#include <libvmmalloc.h>
#endif
/*
 * When recovering, we need to clflush the data in the back_ into the corresponding segment
 * if which doesn't hold the data. So the number of clflush instructions is linear to the data
 * that we lost, which won't be bigger than the size of ALL CPU cache lines. If the size of ALL CPU
 * caches is 16MB, the number of pairs that we lost won't be bigger than 400M. And in our test, we
 * just need 0.68s to recover that.
 */
#define TEST_RECOVERY

// DramCache use these macros
#define ATOMIC_LOAD(ptr)            __atomic_load_n (ptr, __ATOMIC_RELAXED)
#define ATOMIC_STORE(ptr, val)      __atomic_store_n (ptr, val, __ATOMIC_RELAXED)
#define TO_DRAM_CACHE_ENTRIES(v)    ((unsigned long*)((v>>64)&0xfffffffffffffff0))
#define TO_DRAM_CACHE_DEPTH(v)      ((unsigned long)v&0x0fffffffffffffff)
#define TO_DRAM_CACHE_ENTRIES_AND_DEPTH(e, d)   ((((__int128)e)<<64)|d)

#define TO_DRAM_CACHE_ENTRIES_VERSION(v)        ((unsigned long)((v>>64)&0xf))
#define TO_DRAM_CACHE_DEPTH_VERSION(v)          ((unsigned long)v>>60)
#define CHECK_VERSION(v)                        ((TO_DRAM_CACHE_ENTRIES_VERSION(v))==(TO_DRAM_CACHE_DEPTH_VERSION(v)))

#define IS_EMPTY(key)        (key==INVALID_KEY)

// radix tree use these macros
const unsigned long kSegmentFlag = 0x0;
const unsigned long kRadixTreeNodeFlag = 0x1;

#define IS_SEGMENT(x)           ((x&(0x1))==kSegmentFlag)
#define IS_RADIX_TREE_NODE(x)   ((x&kRadixTreeNodeFlag)==kRadixTreeNodeFlag)
#define TO_SEGMENT(x)           ((Segment*)x)
#define TO_RADIX_TREE_NODE(x)   ((RadixTreeNode*)(x&(~kRadixTreeNodeFlag)))

typedef size_t Key_t;
typedef const char* Value_t;
const unsigned long INVALID_KEY = 0xffffffffffffffff;
const unsigned long KEY_SENTINEL = 0xfffffffffffffffe;
//INVALID_PATTERN means this segment won't be used anymore.
const unsigned long INVALID_PATTERN = 0xffffffffffffffff;
const Value_t NONE = 0x0;

struct Pair {
public:
    Key_t key_;
    Value_t value_;
    Pair(void): key_{INVALID_KEY}{};
};
//statistic information
extern unsigned long clflush_seg_cnt_;
extern unsigned long clflush_dir_cnt_;
extern unsigned long recovery_leaf_node_cnt_;
extern unsigned long recovery_clflush_cnt_;     // not thread safe

// ndeh means Non-Doubling Extendible Hashing
const unsigned int kRadixTreeNodeSingleSlotSize = sizeof(unsigned long);
const unsigned int kDramCacheEntrySize = kRadixTreeNodeSingleSlotSize;
const unsigned int kRadixTreeNodeSlotsShift = 3;
const unsigned int kRadixTreeNodeSlotMark = (~((unsigned long)(~0x0)<<kRadixTreeNodeSlotsShift));
const unsigned int kRadixTreeNodeSlotsNum = (1<<kRadixTreeNodeSlotsShift);
const unsigned int kRadixTreeNodeSlotsSize = kRadixTreeNodeSlotsNum*kRadixTreeNodeSingleSlotSize;
//const int kCacheLineSize = 64; // is already defined in util/persist.h
const unsigned int kPairSize = sizeof(Pair);
const unsigned int kNumPairPerCacheLine = kCacheLineSize/kPairSize;

const unsigned int kSegmentBits = 8;
const unsigned int kSegmentSize = (1 << kSegmentBits) * kCacheLineSize;
const unsigned int kSegmentSlotsNum = (kSegmentSize/kPairSize);

const unsigned int kBucketSize = kCacheLineSize;
const unsigned int kNumPairPerBucket = kNumPairPerCacheLine;
const unsigned int kBucketMask = (1 << kSegmentBits)-1;

const unsigned int kHashKeyBits = 64;
const unsigned int kNumBucketToLinearProbe = 4;

const unsigned int kSegmentSplittingFlag = -1;

/*
 * kNumBucketToLinearProbe  kSegmentSlotsNum    max-loadfactor  insert-ops
 * 4                        32                  74%
 * 6                        48                  82%
 */

class RadixTreeNode;

enum CCEHStatus {
    Ok = 0,
    ErrorNotFindSpace,
    SegmentSplitting,
    InsertRetry,
    NotFindKey,
    WrongPattern,
};

class Segment {
public:
    int local_depth_;
    std::shared_mutex shared_mutex_;
    std::shared_mutex *shared_mutex_p; 

#ifdef INPLACE
    //COW do not need is_splitting_, and use pattern_ to check invalid.
    //just a flag, do not need persist.
    volatile bool is_splitting_;
#endif
    //can be calcaulted from depth, do not need persist
    volatile unsigned long pattern_;
    Pair pairs_[kSegmentSlotsNum];
    Segment(int local_depth): local_depth_(local_depth) {
         shared_mutex_p = new std::shared_mutex();
        // the all keys in pairs are initialized with INVALID(-1).
#ifdef INPLACE
        is_splitting_ = false;
#endif
    }
    CCEHStatus Insert(Key_t key, Value_t value, unsigned long segment_pairs_idx, unsigned long key_hashed);
    CCEHStatus Insert4Split(Key_t key, Value_t value, unsigned long segment_pairs_idx);
    CCEHStatus Insert4Recovery(Key_t key, Value_t value, unsigned long segment_pairs_idx);
    Segment** Split(RadixTreeNode *parent, unsigned long key_hashed);
    void SetParent(RadixTreeNode *parent);
    int GetSize();
    CCEHStatus Delete(Key_t key, unsigned long segment_pairs_idx, unsigned long key_hashed);
    void Clear();
    static Segment* NewSegment(int local_depth) {
        Segment *seg;
#ifdef QUARTZ
        seg = (Segment *)pmalloc(sizeof(Segment));
        seg->minus_depth_ = 0;
        seg->parent_ = parent;
        for(int i = 0; i < kSegmentSlotsNum + kSegmentStashSlotsNum; i++)
            seg->pairs_[i].sed_crossed_kv_ = INVALID_CROSSED_KV;

#ifdef INPLACE
        seg->is_splitting_ = false;
#endif

#else
#ifdef USE_AEP
    seg = (Segment *)vmem_malloc(sizeof(Segment));
	if(seg == NULL) {
		std::cout << "vmem malloc" << std::endl;
		exit(1);
	}
	new (seg) Segment(local_depth);
#else
#ifdef USE_LIBVMEMALLOC
     seg = (Segment *)memalign(kCacheLineSize,sizeof(Segment));
	if(seg == NULL) {
		std::cout << "memalign" << std::endl;
		exit(1);
	}
	new (seg) Segment(local_depth);

    #else
        seg = new Segment(local_depth);
#endif
#endif
#endif
        return seg;
    }
private:

};

class DramCache {
public:
    // we store entries_ and depth_ in __int128 to be stored and loaded atomically,
    // and the first 8 bytes is entries_ and the last 8 bytes is depth_.
    // the cache capacity is pow(2, depth_).
    // The depth_ starts from 1, and it equals to (max(depth_)*3-min(minus_depth_)), where the minus_depth_
    // belongs to the segments whose parents' depth_ equals to maximum depth_ of all RadixTreeNode.
    alignas(16) volatile __int128 entries_and_depth_;
    DramCache() {
        unsigned long* entries;
        unsigned long depth = 1;

#ifdef USE_AEP
       entries = (unsigned long*)vmem_malloc( pow(2, depth)*kDramCacheEntrySize);
       if (entries == NULL){
           std::cout << "vmem malloc" << std::endl;
		   exit(1);
       }
      
#else
#ifdef USE_LIBVMEMALLOC
       entries = (unsigned long*)memalign(kCacheLineSize,  pow(2, depth)*kDramCacheEntrySize);
       if(entries == NULL) {
		std::cout << "vmem malloc" << std::endl;
		exit(1);
	   }
     
#else
        if(posix_memalign((void**)&entries, 64, pow(2, depth)*kDramCacheEntrySize) != 0)
            exit(-1);
#endif
#endif
        ATOMIC_STORE(&entries_and_depth_, TO_DRAM_CACHE_ENTRIES_AND_DEPTH(entries, depth));
    }
    ~DramCache(void) {
        ;
    }
    void SetNewSegments(Segment* s[2], int start_entry, int interval);
    unsigned long Double(void);
    unsigned long GetOccupiedSpace();
private:

};

class CCEH {
public:
    int dir_doubling_cnt_;                  // thread-safe
    Timer timer;
    double breakdown_seg_split_;            // not thread-safe
    double breakdown_maintenance_dir_;      // not thread-safe
    CCEH() {
        dram_cache_ = new DramCache();
        Segment *seg = nullptr;
        for(int i=0; i<2; i++) {
            seg = Segment::NewSegment(1);
            seg->pattern_ = i;
            TO_DRAM_CACHE_ENTRIES(dram_cache_->entries_and_depth_)[i] = (unsigned long)seg;
        }
	mutex_p = new std::mutex();
        dir_doubling_cnt_ = 0;
        breakdown_seg_split_ = 0;
        breakdown_maintenance_dir_ = 0;
    }
    ~CCEH() {
        if(dram_cache_ != nullptr)
            delete dram_cache_;
    }
    void Insert(Key_t key, Value_t value);
    bool Delete(Key_t key);
    Value_t Get(Key_t key);
    void Recovery(int nr_threads);
    unsigned long StatisticGetDirectoryOccupiedSpace();
    double StatisticGetLoadfactor();
private:
    DramCache *dram_cache_;
    std::mutex mutex_;
    std::mutex *mutex_p; 
};


#endif  // EXTENDIBLE_PTR_H_
