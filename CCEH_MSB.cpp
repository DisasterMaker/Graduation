#ifndef CCEH_H_
#include "CCEH.h"
#endif

#ifndef UTIL_HASH_H_
#include "util/hash.h"
#endif

#include <iostream>
#include <queue>
#include <cstdio>
#include <vector>
#include <thread>
#include <sys/time.h>

unsigned long clflush_seg_cnt_ = 0;
unsigned long clflush_dir_cnt_ = 0;
unsigned long recovery_leaf_node_cnt_ = 0;
unsigned long recovery_clflush_cnt_ = 0;

// we do not need lw_lock_ to lock, because which won't cause unconsistency
void DramCache::SetNewSegments(Segment** s, int start_entry, int interval) {
    __int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&entries_and_depth_);
    unsigned long *entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    int tmp_start_entry;
    for(int i=1; i>=0; i--) {
        tmp_start_entry = start_entry + pow(2, interval) * i;
        for(int j=0; j < pow(2, interval); j++) {
            entries[tmp_start_entry+j] = (unsigned long)s[i];
        }
        clflush((char*)&entries[tmp_start_entry], pow(2, interval)*sizeof(unsigned long));
    }
}

unsigned long DramCache::Double(void) {
    __int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&entries_and_depth_);
    unsigned long *old_entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    unsigned long old_depth = TO_DRAM_CACHE_DEPTH(dram_cache_entries_and_depth);
    unsigned long version = TO_DRAM_CACHE_ENTRIES_VERSION(dram_cache_entries_and_depth);
    unsigned long * new_entries;
    new_entries = (unsigned long*)vmem_malloc(pow(2, old_depth+1)*kDramCacheEntrySize);
       if (new_entries == NULL){
           std::cout << "vmem malloc" << std::endl;
		   exit(1);
       }
    /* if(__glibc_unlikely(posix_memalign((void**)&new_entries, 64, pow(2, old_depth+1)*kDramCacheEntrySize)!=0)) {
        std::cout << "error when use posix_memalign get memory" << std::endl;
        exit(-1);
    }*/
    for(int i=0; i<pow(2, old_depth); i++) {
        for(int j=0; j<2; j++) {
            new_entries[i*2+j] = old_entries[i];
        }
    }
    clflush((char*)new_entries, pow(2, old_depth+1)*kDramCacheEntrySize);
    ATOMIC_STORE(&entries_and_depth_, TO_DRAM_CACHE_ENTRIES_AND_DEPTH(new_entries|((version+1)%16), \
        (old_depth+1)|(((version+1)%16)<<60)));
    clflush((char*)&entries_and_depth_, sizeof(__int128));
    return old_depth+1;
}

inline unsigned long DramCache::GetOccupiedSpace() {
    return pow(2, TO_DRAM_CACHE_DEPTH(ATOMIC_LOAD(&entries_and_depth_)))*kDramCacheEntrySize;
}

#ifdef INPLACE
CCEHStatus Segment::Insert(Key_t key, Value_t value, \
    unsigned long segment_pairs_idx, unsigned long key_hashed) {
    if(shared_mutex_p->try_lock_shared() == false)
        return SegmentSplitting;
    //check if the pattern is right, it may be wrong when a thread get a wrong seg with old index.
    //this check must be protected by rwlock.
    if((key_hashed>>(kHashKeyBits-(local_depth_))) != pattern_) {
        shared_mutex_p->unlock_shared();
        return WrongPattern;
    }
    int slot = 0;
    unsigned long invalid_key = INVALID_KEY;
    unsigned long expected_key;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        expected_key = pairs_[slot].key_;
        if((h(&expected_key, sizeof(Key_t))>> \
            (kHashKeyBits-(local_depth_)))!=pattern_ \
            && expected_key != KEY_SENTINEL \
            && (CAS(&pairs_[slot].key_, &expected_key, KEY_SENTINEL))) {
            pairs_[slot].value_ = value;
            mfence();
            pairs_[slot].key_ = key;
            shared_mutex_p->unlock_shared();
            return Ok;
        }
        if((IS_EMPTY(pairs_[slot].key_)==true) \
            && (CAS(&pairs_[slot].key_, &invalid_key, KEY_SENTINEL))){
            pairs_[slot].value_ = value;
            mfence();
            pairs_[slot].key_ = key;
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    shared_mutex_p->unlock_shared();
    return ErrorNotFindSpace;
}
#else
CCEHStatus Segment::Insert(Key_t key, Value_t value, \
    unsigned long segment_pairs_idx, unsigned long key_hashed) {
    if(shared_mutex_p->try_lock_shared() == false)
        return SegmentSplitting;
    //check if the pattern is right, it may be wrong when a thread get a wrong seg with old index.
    //this check must be protected by rwlock.
    if((key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_))) != pattern_) {
        shared_mutex_p->unlock_shared();
        return WrongPattern;
    }
    int slot = 0;
    unsigned long invalid_crossed_kv = INVALID_CROSSED_KV;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if((IS_EMPTY(pairs_[slot].sed_crossed_kv_)==true) \
            && (CAS(&pairs_[slot].sed_crossed_kv_, &invalid_crossed_kv, CROSSED_KV_SENTINEL))){
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    //find space in stash slots
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i + kSegmentSlotsNum;
        if((IS_EMPTY(pairs_[slot].sed_crossed_kv_)==true) \
            && (CAS(&pairs_[slot].sed_crossed_kv_, &invalid_crossed_kv, CROSSED_KV_SENTINEL))){
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    shared_mutex_p->unlock_shared();
    return ErrorNotFindSpace;
}
#endif

CCEHStatus Segment::Insert4Split(Key_t key, Value_t value, \
    unsigned long segment_pairs_idx) {
    int slot = 0;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(IS_EMPTY(pairs_[slot].key_)==true){
            pairs_[slot].key_ = key;
            pairs_[slot].value_ = value;
            return Ok;
        }
    }
    return ErrorNotFindSpace;
}

void Segment::Clear() {
    for(int i=0; i<kSegmentSlotsNum; i++) {
        if(IS_EMPTY(pairs_[i].key_) == false) {
            pairs_[i].key_ = INVALID_KEY;
        }
    }
}

#ifdef INPLACE
Segment** Segment::Split(RadixTreeNode *parent, unsigned long key_hashed) {
    shared_mutex_p->lock();
    //must be protected by write lock.
    if(is_splitting_ == true) {
        shared_mutex_p->unlock();
        return nullptr;
    }
    //pattern checking is also a must, because the parameter may be wrong.
    if((key_hashed>>(kHashKeyBits-(local_depth_))) \
        != pattern_) {
        shared_mutex_p->unlock();
        return nullptr;
    }
    unsigned long old_pattern = pattern_;
    int old_local_depth = local_depth_;
    Segment** split = new Segment*[2];
    split[0] = this;
    split[0]->is_splitting_ = true;
    split[1] = Segment::NewSegment(0);
    split[1]->is_splitting_ = true;

    //update s[0], s[1] patterns
    split[0]->pattern_ = (key_hashed>>(kHashKeyBits-old_local_depth)) << 1;
    split[1]->pattern_ = ((key_hashed>>(kHashKeyBits-old_local_depth)) << 1) + 1;

    // update local_depth_
    split[0]->local_depth_ = old_local_depth+1;
    split[1]->local_depth_ = old_local_depth+1;

    unsigned long tmp_key_hashed = 0;
    CCEHStatus reval;
    for(int i = 0; i<kSegmentSlotsNum; i++) {
        tmp_key_hashed = h(&pairs_[i].key_, sizeof(Key_t));
        if(tmp_key_hashed&((unsigned long)1<<(kHashKeyBits-old_local_depth-1)))
            reval = split[1]->Insert4Split(pairs_[i].key_, pairs_[i].value_, (tmp_key_hashed&kBucketMask)*kNumPairPerBucket);
        if(__glibc_unlikely(reval == ErrorNotFindSpace)) {
            split[1]->Clear();
            goto direct_insert;
        }
    }

split_finished:
    clflush((char*)split[1], sizeof(Segment));
    clflush_seg_cnt_++;
    clflush((char*)&local_depth_, sizeof(int));
    shared_mutex_p->unlock();
    return split;

direct_insert:
    for(int i = 0; i<kSegmentSlotsNum; i++) {
        tmp_key_hashed = h(&pairs_[i].key_, sizeof(Key_t));
        if(tmp_key_hashed&((unsigned long)1<<(kHashKeyBits-old_local_depth-1))) {
            split[1]->pairs_[i].key_ = pairs_[i].key_;
            split[1]->pairs_[i].value_ = pairs_[i].value_;
        }
    }
    goto split_finished;
}
#else
Segment** Segment::Split(RadixTreeNode *parent, unsigned long key_hashed) {
    shared_mutex_p->lock();
    //pattern checking is also a must, because the parameter may be wrong.
    if((key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_))) \
        != pattern_) {
        shared_mutex_p->unlock();
        return nullptr;
    }
    Segment** split = new Segment*[2];
    split[0] = Segment::NewSegment(parent);
    split[1] = Segment::NewSegment(parent);
    //update s[0], s[1] patterns
    this->pattern_ = INVALID_PATTERN;   //means we won't use this segment
    split[0]->pattern_ = (key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_+1)+1)) << 1;
    split[1]->pattern_ = ((key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_+1)+1)) << 1) + 1;
    //update minus_depth_
    split[0]->minus_depth_ = (minus_depth_-1+kRadixTreeNodeSlotsShift)%kRadixTreeNodeSlotsShift;
    split[1]->minus_depth_ = (minus_depth_-1+kRadixTreeNodeSlotsShift)%kRadixTreeNodeSlotsShift;

    unsigned long tmp_key_hashed = 0;
    Key_t tmp_key;
    unsigned long last_four_bytes_key;
    Value_t tmp_value;
    CCEHStatus reval;
    for(int i = 0; i<kSegmentSlotsNum + kSegmentStashSlotsNum; i++) {
        // we do not need to move invalid key-value to new segments
        if(IS_EMPTY(pairs_[i].sed_crossed_kv_) == true)
            continue;
        tmp_key = TO_KEY(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
        last_four_bytes_key = LAST_FOUR_BYTES(tmp_key);
        tmp_value = (Value_t)TO_VALUE(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
        tmp_key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        if(tmp_key_hashed&((unsigned long)1<<(kHashKeyBits-(kRadixTreeNodeSlotsShift*parent_->depth_-minus_depth_)-1)))
            reval = split[1]->Insert4Split(tmp_key, tmp_value, (tmp_key_hashed&kBucketMask)*kNumPairPerBucket);
        else
            reval = split[0]->Insert4Split(tmp_key, tmp_value, (tmp_key_hashed&kBucketMask)*kNumPairPerBucket);
        if(__glibc_unlikely(reval == ErrorNotFindSpace)) {
            split[1]->Clear();
            split[0]->Clear();
            goto direct_insert;
        }
    }

split_finished:
    clflush((char*)split[0], sizeof(Segment));
    clflush((char*)split[1], sizeof(Segment));
    shared_mutex_p->unlock();
    return split;

direct_insert:
    for(int i = 0; i<kSegmentSlotsNum + kSegmentStashSlotsNum; i++) {
        // we do not need to move invalid key-value to new segments
        if(IS_EMPTY(pairs_[i].sed_crossed_kv_) == true)
            continue;
        tmp_key = TO_KEY(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
        last_four_bytes_key = LAST_FOUR_BYTES(tmp_key);
        key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        if(key_hashed&((unsigned long)1<<(kHashKeyBits-(kRadixTreeNodeSlotsShift*parent_->depth_-minus_depth_)-1))) {
            split[1]->pairs_[i].fir_crossed_kv_ = pairs_[i].fir_crossed_kv_;
            split[1]->pairs_[i].sed_crossed_kv_ = pairs_[i].sed_crossed_kv_;
        }
        else {
            split[0]->pairs_[i].fir_crossed_kv_ = pairs_[i].fir_crossed_kv_;
            split[0]->pairs_[i].sed_crossed_kv_ = pairs_[i].sed_crossed_kv_;
        }
    }
    goto split_finished;
}
#endif

static double mysecond() {
	    struct timeval tp;
	        int i;
		    i = gettimeofday(&tp, NULL);
		        return ( (double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

void CCEH::Insert(Key_t key, Value_t value) {
retry:
    unsigned long key_hashed = h(&key, sizeof(Key_t));
    unsigned long bucket_idx = (key_hashed&kBucketMask);
    unsigned long segment_pairs_idx = bucket_idx*kNumPairPerBucket;

    //__int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&dram_cache_->entries_and_depth_);
    __int128 dram_cache_entries_and_depth = dram_cache_->entries_and_depth_;
    if(!CHECK_VERSION(dram_cache_entries_and_depth))
        goto retry;
    unsigned long *dram_cache_entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    unsigned long dram_cache_depth = TO_DRAM_CACHE_DEPTH(dram_cache_entries_and_depth);
    unsigned long segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-dram_cache_depth);
    Segment *seg = (Segment*)dram_cache_entries[segment_idx_in_dram_cache];
    CCEHStatus status = seg->Insert(key, value, segment_pairs_idx, key_hashed);
    // split segment
    if(status == ErrorNotFindSpace) {
        // print current load factor
        //std::cout << StatisticGetLoadfactor() << std::endl;
        timer.Start();
        Segment **s = seg->Split(nullptr, key_hashed);
        timer.Stop();
        breakdown_seg_split_ += timer.GetSeconds();
        if(s == nullptr)
            goto retry;
	
        //timer.Start();
	double tmp_a = mysecond();
        mutex_p->lock();
        //reget the entries and depth, becasuse it may has been changed.
        dram_cache_depth = TO_DRAM_CACHE_DEPTH(ATOMIC_LOAD(&dram_cache_->entries_and_depth_));
        segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-dram_cache_depth);

        unsigned long seg_depth = s[1]->local_depth_-1;
        if(seg_depth < dram_cache_depth) {
            //updating dram cache first won't cause unconsistency.
            dram_cache_->SetNewSegments(s, segment_idx_in_dram_cache&((~0x0)<<(dram_cache_depth-seg_depth)), \
            dram_cache_depth-(seg_depth+1));
        }else {
            // double the dram cache
            unsigned long new_dram_cache_depth = dram_cache_->Double();
            dir_doubling_cnt_++;
            unsigned long new_segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-new_dram_cache_depth);
            dram_cache_->SetNewSegments(s, \
                new_segment_idx_in_dram_cache&((~0x0)<<(new_dram_cache_depth-seg_depth)), \
                new_dram_cache_depth-(seg_depth+1));
        }
        //timer.Stop();
        breakdown_maintenance_dir_ += mysecond() - tmp_a;;
        mutex_p->unlock();
#ifdef INPLACE
        // COW won't use seg anymore, so do nothing.
        s[0]->is_splitting_ = false;
        s[1]->is_splitting_ = false;
#endif
        goto retry;
    }else if(status == SegmentSplitting || status == WrongPattern) {
        goto retry;
    }
}

unsigned long CCEH::StatisticGetDirectoryOccupiedSpace() {
    return dram_cache_->GetOccupiedSpace();
}

double CCEH::StatisticGetLoadfactor() {
    //return radix_tree_->GetLoadfactor();
    return 0;
}

Value_t CCEH::Get(Key_t key) {
    unsigned long key_hashed = h(&key, sizeof(Key_t));
    unsigned long bucket_idx = (key_hashed&kBucketMask);
    unsigned long segment_pairs_idx = bucket_idx*kNumPairPerBucket;
retry:
    //__int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&dram_cache_->entries_and_depth_);
    __int128 dram_cache_entries_and_depth = dram_cache_->entries_and_depth_;
    if(!CHECK_VERSION(dram_cache_entries_and_depth))
        goto retry;
    unsigned long *dram_cache_entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    unsigned long dram_cache_depth = TO_DRAM_CACHE_DEPTH(dram_cache_entries_and_depth);
    unsigned long segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-dram_cache_depth);
    Segment *seg = (Segment*)dram_cache_entries[segment_idx_in_dram_cache];
#ifdef INPLACE
    if(seg->shared_mutex_p->try_lock_shared() == false)
        goto retry;
    if((key_hashed>>(kHashKeyBits-seg->local_depth_)) \
        != seg->pattern_) {
        seg->shared_mutex_p->unlock_shared();
        goto retry;
    }
#endif
    int slot = 0;
    // the direction of insert and get is opposite, because of update
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(seg->pairs_[slot].key_ == key) {
#ifdef INPLACE
            seg->shared_mutex_p->unlock_shared();
#endif
            return seg->pairs_[slot].value_;
        }
    }
#ifdef INPLACE
    seg->shared_mutex_p->unlock_shared();
#endif
    return NONE;
}
/*
 * When deleting, we just set the secondary crossed-kv to be INVALID_CROSSED_KV.
 * TODO: we also need to delete the key-value in back_ if not nullptr.
 */
CCEHStatus Segment::Delete(Key_t key, unsigned long segment_pairs_idx, unsigned long key_hashed) {
    if(shared_mutex_p->try_lock_shared() == false)
        return SegmentSplitting;
    if((key_hashed>>(kHashKeyBits-local_depth_)) != pattern_) {
        shared_mutex_p->unlock_shared();
        return WrongPattern;
    }
    int slot = 0;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(pairs_[slot].key_ == key){
            pairs_[slot].key_ = INVALID_KEY;
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    shared_mutex_p->unlock_shared();
    return NotFindKey;
}
/*
 * TODO: shrink the radix tree in NVM and the flat directory in DRAM
 */
bool CCEH::Delete(Key_t key) {
    unsigned long key_hashed = h(&key, sizeof(Key_t));
    unsigned long bucket_idx = (key_hashed&kBucketMask);
    unsigned long segment_pairs_idx = bucket_idx*kNumPairPerBucket;
retry:
    //__int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&dram_cache_->entries_and_depth_);
    __int128 dram_cache_entries_and_depth = dram_cache_->entries_and_depth_;
    if(!CHECK_VERSION(dram_cache_entries_and_depth))
        goto retry;
    unsigned long *dram_cache_entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    unsigned long dram_cache_depth = TO_DRAM_CACHE_DEPTH(dram_cache_entries_and_depth);
    unsigned long segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-dram_cache_depth);
    Segment *seg = (Segment*)dram_cache_entries[segment_idx_in_dram_cache];
    CCEHStatus status = seg->Delete(key, segment_pairs_idx, key_hashed);
    if(status == SegmentSplitting)
        goto retry;
    else
        return status == Ok ? true: false;
}

void CCEH::Recovery(int nr_threads) {
    //radix_tree_->Recovery(nr_threads);
}

