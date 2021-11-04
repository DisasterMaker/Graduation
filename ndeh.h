#ifndef NDEH_H_
#define NDEH_H_

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
#define TEST_RECOVERY

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
#define RECOVER_BACK

// DramCache use these macros
#define ATOMIC_LOAD(ptr)            __atomic_load_n (ptr, __ATOMIC_RELAXED)
#define ATOMIC_STORE(ptr, val)      __atomic_store_n (ptr, val, __ATOMIC_RELAXED)
#define TO_DRAM_CACHE_ENTRIES(v)    ((unsigned long*)((v>>64)&0xfffffffffffffff0))
#define TO_DRAM_CACHE_DEPTH(v)      ((unsigned long)v&0x0fffffffffffffff)
#define TO_DRAM_CACHE_ENTRIES_AND_DEPTH(e, d)   ((((__int128)e)<<64)|d)

#define TO_DRAM_CACHE_ENTRIES_VERSION(v)        ((unsigned long)((v>>64)&0xf))
#define TO_DRAM_CACHE_DEPTH_VERSION(v)          ((unsigned long)v>>60)
#define CHECK_VERSION(v)                        ((TO_DRAM_CACHE_ENTRIES_VERSION(v))==(TO_DRAM_CACHE_DEPTH_VERSION(v)))

// Segment use these macros
#define FIR_FOUR_BYTES(x)               (((unsigned long)x>>32)&0xffffffff)
#define LAST_FOUR_BYTES(x)              ((unsigned long)x&0xffffffff)
#define TO_FIR_CROSSED_KV(key, value)   ((FIR_FOUR_BYTES(key)<<32)|FIR_FOUR_BYTES(value))
#define TO_SED_CROSSED_KV(key, value)   ((LAST_FOUR_BYTES(key)<<32)|LAST_FOUR_BYTES(value))          
#define TO_KEY(fir_crossed_kv, sed_crossed_kv)      ((FIR_FOUR_BYTES(fir_crossed_kv)<<32)| \
                                                        FIR_FOUR_BYTES(sed_crossed_kv))
#define TO_VALUE(fir_crossed_kv, sed_crossed_kv)    ((LAST_FOUR_BYTES(fir_crossed_kv)<<32)| \
                                                        LAST_FOUR_BYTES(sed_crossed_kv))
#define IS_EMPTY(sed_crossed_kv)        (sed_crossed_kv==INVALID_CROSSED_KV)

// radix tree use these macros
const unsigned long kSegmentFlag = 0x0;
const unsigned long kRadixTreeNodeFlag = 0x1;

#define IS_SEGMENT(x)           ((x&(0x1))==kSegmentFlag)
#define IS_RADIX_TREE_NODE(x)   ((x&kRadixTreeNodeFlag)==kRadixTreeNodeFlag)
#define TO_SEGMENT(x)           ((Segment*)x)
#define TO_RADIX_TREE_NODE(x)   ((RadixTreeNode*)(x&(~kRadixTreeNodeFlag)))

typedef size_t Key_t;
typedef const char* Value_t;
const unsigned long INVALID_CROSSED_KV = 0xffffffffffffffff;
const unsigned long CROSSED_KV_SENTINEL = 0xfffffffffffffffe;
//INVALID_PATTERN means this segment won't be used anymore.
const unsigned long INVALID_PATTERN = 0xffffffffffffffff;
const Value_t NONE = 0x0;

struct Pair {
public:
    unsigned long fir_crossed_kv_;
    unsigned long sed_crossed_kv_;
    Pair(void): fir_crossed_kv_{INVALID_CROSSED_KV}, \
    sed_crossed_kv_{INVALID_CROSSED_KV}{};
};
//statistic information
extern unsigned long clflush_seg_cnt_;
extern unsigned long clflush_dir_cnt_;
extern unsigned long recovery_leaf_node_cnt_;
extern unsigned long recovery_clflush_cnt_;     // not thread safe

// ndeh means Non-Doubling Extendible Hashing
namespace ndeh {
const unsigned int kRadixTreeNodeSingleSlotSize = sizeof(unsigned long);
const unsigned int kDramCacheEntrySize = kRadixTreeNodeSingleSlotSize;
const unsigned int kRadixTreeNodeSlotsShift = 3;
const unsigned int kRadixTreeNodeSlotMark = (~((unsigned long)(~0x0)<<kRadixTreeNodeSlotsShift));
const unsigned int kRadixTreeNodeSlotsNum = (1<<kRadixTreeNodeSlotsShift);
const unsigned int kRadixTreeNodeSlotsSize = kRadixTreeNodeSlotsNum*kRadixTreeNodeSingleSlotSize;
//const int kCacheLineSize = 64; // is already defined in util/persist.h
const unsigned int kPairSize = sizeof(Pair);
const unsigned int kNumPairPerCacheLine = kCacheLineSize/kPairSize;

const unsigned int kSegmentBits =8;
const unsigned int kSegmentSize = (1 << kSegmentBits) * kCacheLineSize;
const unsigned int kSegmentSlotsNum = (kSegmentSize/kPairSize);

const unsigned int kBucketSize = kCacheLineSize;
const unsigned int kNumPairPerBucket = kNumPairPerCacheLine;
const unsigned int kBucketMask = (1 << kSegmentBits)-1;

const unsigned int kHashKeyBits = 64;
const unsigned int kNumBucketToLinearProbe = 4;

const unsigned int kSegmentSplittingFlag = -1;
const unsigned int kSegmentStashSlotsNum = 16;

/*
 * kNumBucketToLinearProbe  kSegmentSlotsNum    max-loadfactor  insert-ops
 * 4                        32                  74%
 * 6                        48                  82%
 */

class RadixTreeNode;
enum NdehStatus {
    Ok = 0,
    ErrorNotFindSpace,
    SegmentSplitting,
    InsertRetry,
    NotFindKey,
    WrongPattern,
};

class Segment {
public:
    //need persist
    int minus_depth_;
    //can be get from recover, do not need persist
    RadixTreeNode *parent_;
    RwLock rw_lock_;
    std::shared_mutex shared_mutex_;
    std::shared_mutex *shared_mutex_p; 

#ifdef INPLACE
    //COW do not need is_splitting_, and use pattern_ to check invalid.
    //just a flag, do not need persist.
    volatile bool is_splitting_;
#endif
    //can be calcaulted from depth, do not need persist
    volatile unsigned long pattern_;
#ifdef INPLACE
    // used when recoverying.
    unsigned long saved_pattern_;
#endif
    Pair pairs_[kSegmentSlotsNum + kSegmentStashSlotsNum];
    Segment(RadixTreeNode *parent): minus_depth_(0), \
        parent_(parent) {
        shared_mutex_p = new std::shared_mutex();
        // the all keys in pairs are initialized with INVALID(-1).
#ifdef INPLACE
        is_splitting_ = false;
#endif
    }
    NdehStatus Insert(Key_t key, Value_t value, unsigned long segment_pairs_idx, unsigned long key_hashed);
    NdehStatus Insert4Split(Key_t key, Value_t value, unsigned long segment_pairs_idx);
    NdehStatus Insert4Recovery(Key_t key, Value_t value, unsigned long segment_pairs_idx);
    Segment** Split(RadixTreeNode *parent, unsigned long key_hashed);
    void SetParent(RadixTreeNode *parent);
    int GetSize();
    NdehStatus Delete(Key_t key, unsigned long segment_pairs_idx, unsigned long key_hashed);
    void Clear();
    static Segment* NewSegment(RadixTreeNode *parent) {
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
	new (seg) Segment(parent);
#else
#ifdef USE_LIBVMEMALLOC
     seg = (Segment *)memalign(kCacheLineSize,sizeof(Segment));
	if(seg == NULL) {
		std::cout << "memalign" << std::endl;
		exit(1);
	}
	new (seg) Segment(parent);
    #else
	seg = new Segment(parent);
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
    DramCache(int depth = 1) {
        unsigned long *entries;
        if(posix_memalign((void**)&entries, 64, pow(2, depth)*kDramCacheEntrySize) != 0) {
            std::cout << "allocate memory failed, size: " << pow(2, depth)*kDramCacheEntrySize << std::endl;
            exit(-1);
        }
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

// we use radix tree as directory in PM
class RadixTreeNode {
public:
    // the elements in slots are all pointers, if the last bit is 1, means
    // this is a leaf node, and the pointer is pointing to a segment
    unsigned long slots_[kRadixTreeNodeSlotsNum];
    int depth_;
    Segment *back_;
    int recovery_used_;
    const static unsigned long INVALID_SLOT = 0x0;
    RadixTreeNode(int depth) {
        depth_ = depth;
        back_ = nullptr;
    }
    ~RadixTreeNode() {
        if(slots_ != nullptr)
            free(slots_);
    }
    void SetNewSegments(Segment* s[2], int start_slot, bool does_need_mfence);
    static RadixTreeNode *NewRadixTreeNode(int depth) {
        RadixTreeNode *radix_tree_node;
#ifdef QUARTZ
        radix_tree_node = (RadixTreeNode *)pmalloc(sizeof(RadixTreeNode));
        radix_tree_node->depth_ = depth;
        radix_tree_node->back_ = nullptr;
#else
#ifdef USE_AEP
	radix_tree_node = (RadixTreeNode *)vmem_malloc(sizeof(RadixTreeNode));
	if(radix_tree_node == NULL) {
		std::cout << "vmem malloc" << std::endl;
		exit(1);
	}
	new (radix_tree_node) RadixTreeNode(depth);
#else
#ifdef USE_LIBVMEMALLOC
    radix_tree_node = (RadixTreeNode *)memalign(kCacheLineSize,sizeof(RadixTreeNode));
	if(radix_tree_node == NULL) {
		std::cout << "vmem malloc" << std::endl;
		exit(1);
	}
	new (radix_tree_node) RadixTreeNode(depth);
    #else
        radix_tree_node = new RadixTreeNode(depth);
#endif

#endif

#endif

        return radix_tree_node;
    }
    bool IsLeafNode();
    void Recovery();
private:
    
};

// not optimized too early
class RadixTree {
public:
    RadixTree(DramCache *dram_cache) {
        root_ = RadixTreeNode::NewRadixTreeNode(1);
        dram_cache_ = dram_cache;
        Segment *seg = nullptr;
        for(int i=0; i<2; i++) {
            seg = Segment::NewSegment(root_);
            seg->minus_depth_ = kRadixTreeNodeSlotsShift-1;
            seg->pattern_ = i;
            for(int j=0; j<kRadixTreeNodeSlotsNum/2; j++)
                root_->slots_[i*kRadixTreeNodeSlotsNum/2 + j] = (unsigned long)seg;
            TO_DRAM_CACHE_ENTRIES(dram_cache_->entries_and_depth_)[i] = (unsigned long)seg;
        }
    }
    void Recovery(int nr_threads);
    static void RevoveryThread(std::vector<RadixTreeNode*> *nodes_needed_recovery, int start, int end);
    double GetLoadfactor();
    double rebuild_dir_time_;
    double recovery_rt_dir_time_;
private:
    RadixTreeNode *root_;  
    DramCache *dram_cache_;
};

class Ndeh {
public:
    int dir_doubling_cnt_;                  // thread-safe
    Timer timer;
    double breakdown_seg_split_;            // not thread-safe
    double breakdown_maintenance_dir_;      // not thread-safe
    Ndeh() {
        dram_cache_ = new DramCache();
        radix_tree_ = new RadixTree(dram_cache_);
        dir_doubling_cnt_ = 0;
        breakdown_seg_split_ = 0;
        breakdown_maintenance_dir_ = 0;
    }
    ~Ndeh() {
        if(dram_cache_ != nullptr)
            delete dram_cache_;
        if(radix_tree_ != nullptr)
            delete radix_tree_;
    }
    void Insert(Key_t key, Value_t value);
    bool Delete(Key_t key);
    Value_t Get(Key_t key);
    void Recovery(int nr_threads);
    unsigned long StatisticGetDirectoryOccupiedSpace();
    double StatisticGetLoadfactor();
public:
    DramCache *dram_cache_;
    RadixTree *radix_tree_;      // used as directory
    std::mutex mutex_;
};

} // namespace ndeh

#endif
