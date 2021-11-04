#ifndef NDEH_H_
#include "ndeh.h"
#endif

#ifndef UTIL_HASH_H_
#include "util/hash.h"
#endif

#include <iostream>
#include <queue>
#include <cstdio>
#include <vector>
#include <thread>

unsigned long clflush_seg_cnt_ = 0;
unsigned long clflush_dir_cnt_ = 0;
unsigned long recovery_leaf_node_cnt_ = 0;
unsigned long recovery_clflush_cnt_ = 0;

namespace ndeh {

// we do not need lw_lock_ to lock, because which won't cause unconsistency
void DramCache::SetNewSegments(Segment** s, int start_entry, int interval) {
    __int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&entries_and_depth_);
    unsigned long *entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    int tmp_start_entry;
    for(int i=0; i<2; i++) {
        tmp_start_entry = start_entry + pow(2, interval) * i;
        for(int j=0; j < pow(2, interval); j++) {
            entries[tmp_start_entry+j] = (unsigned long)s[i];
        }
    }
}

unsigned long DramCache::Double(void) {
    __int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&entries_and_depth_);
    unsigned long *old_entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    unsigned long old_depth = TO_DRAM_CACHE_DEPTH(dram_cache_entries_and_depth);
    unsigned long version = TO_DRAM_CACHE_ENTRIES_VERSION(dram_cache_entries_and_depth);
    unsigned long * new_entries;
    if(__glibc_unlikely(posix_memalign((void**)&new_entries, 64, pow(2, old_depth+1)*kDramCacheEntrySize)!=0)) {
        std::cout << "error when use posix_memalign get memory" << std::endl;
        exit(-1);
    }
    for(int i=0; i<pow(2, old_depth); i++) {
        for(int j=0; j<2; j++) {
            new_entries[i*2+j] = old_entries[i];
        }
    }
    ATOMIC_STORE(&entries_and_depth_, TO_DRAM_CACHE_ENTRIES_AND_DEPTH(new_entries|((version+1)%16), \
        (old_depth+1)|(((version+1)%16)<<60)));
    return old_depth+1;
}

inline unsigned long DramCache::GetOccupiedSpace() {
    return pow(2, TO_DRAM_CACHE_DEPTH(ATOMIC_LOAD(&entries_and_depth_)))*kDramCacheEntrySize;
}

//get minus_depth_ from s[1], because s[0] maybe changed with INPLACE
void RadixTreeNode::SetNewSegments(Segment** s, int start_slot, bool does_need_mfence) {
    // update from right to left
    int tmp_start_slot;
    for(int i=1; i>=0; i--) {
        tmp_start_slot = start_slot + pow(2, s[1]->minus_depth_) * i;
        for(int j=0; j < pow(2, s[1]->minus_depth_); j++) {
            slots_[tmp_start_slot+j] = (unsigned long)s[i];
        }
        if(does_need_mfence == true)
            mfence();
    }
}
/*
 * TODO: to make sure what a leaf node is, for now, if the slots in a node
 * point some segments, we treate it as a leaf node.
 */
bool RadixTreeNode::IsLeafNode() {
    for(int i=0; i<pow(2, kRadixTreeNodeSlotsShift); i++) {
        if(IS_SEGMENT(slots_[i]))
            return true;
    }
    return false;
}
/*
 * TODO: consider the situation that we can not find space when inserting.
 */
NdehStatus Segment::Insert4Recovery(Key_t key, Value_t value, unsigned long segment_pairs_idx) {
    int slot;
#ifdef TEST_RECOVERY
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(TO_KEY(pairs_[slot].fir_crossed_kv_, pairs_[slot].sed_crossed_kv_) == key) {
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            clflush((char*)&pairs_[slot].fir_crossed_kv_, kPairSize);
            recovery_clflush_cnt_++;
            return Ok;
        }
    }
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i+kSegmentSlotsNum;
        if(TO_KEY(pairs_[slot].fir_crossed_kv_, pairs_[slot].sed_crossed_kv_) == key) {
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            clflush((char*)&pairs_[slot].fir_crossed_kv_, kPairSize);
            recovery_clflush_cnt_++;
            return Ok;
        }
    }

#else
    // try find key first
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(TO_KEY(pairs_[slot].fir_crossed_kv_, pairs_[slot].sed_crossed_kv_) == key)
            return Ok;
    }
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i+kSegmentSlotsNum;
        if(TO_KEY(pairs_[slot].fir_crossed_kv_, pairs_[slot].sed_crossed_kv_) == key)
            return Ok;
    }
    //if not found, insert key and value
    unsigned long tmp_key_hashed;
    unsigned long last_four_bytes_key;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        last_four_bytes_key = FIR_FOUR_BYTES(pairs_[i].sed_crossed_kv_);
        tmp_key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        if(IS_EMPTY(pairs_[slot].sed_crossed_kv_) \
            || ((tmp_key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_))) != pattern_)) {
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            clflush((char*)&pairs_[slot].fir_crossed_kv_, kPairSize);
            recovery_clflush_cnt_++;
            return Ok;
        }
    }
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i+kSegmentSlotsNum;
        last_four_bytes_key = FIR_FOUR_BYTES(pairs_[i].sed_crossed_kv_);
        tmp_key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        if(IS_EMPTY(pairs_[slot].sed_crossed_kv_) \
            || ((tmp_key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_))) != pattern_)) {
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            recovery_clflush_cnt_++;
            clflush((char*)&pairs_[slot].fir_crossed_kv_, kPairSize);
            return Ok;
        }
    }
#endif
    return Ok;
}

void RadixTreeNode::Recovery() {
    //segment recovery
    int buddy;
    for(int i=0; i<pow(2, kRadixTreeNodeSlotsShift); i++) {
        if(IS_SEGMENT(slots_[i])) {
            buddy = i + pow(2, TO_SEGMENT(slots_[i])->minus_depth_);
            for(int j=buddy-1; i<j; j--) {
                if(TO_SEGMENT(slots_[j])->minus_depth_ != TO_SEGMENT(slots_[i])->minus_depth_)
                    slots_[j] = slots_[i];
            }
            i = buddy - 1;
        }
    }

// just INPLACE mode has back_
#ifdef RECOVERY_BACK
#ifdef INPLACE
    if(back_ == nullptr)
        return;
    //data recovery
    Key_t tmp_key;
    Value_t tmp_value;
    unsigned long last_four_bytes_key;
    unsigned long tmp_key_hashed;
    int seg_pos;
    for(int i=0; i<kSegmentSlotsNum+kSegmentStashSlotsNum; i++) {
        if(IS_EMPTY(back_->pairs_[i].sed_crossed_kv_) == true)
            continue;
        tmp_key = TO_KEY(back_->pairs_[i].fir_crossed_kv_, back_->pairs_[i].sed_crossed_kv_);
        last_four_bytes_key = LAST_FOUR_BYTES(tmp_key);
        tmp_key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        //check whether the key belongs to this segment.
        if((tmp_key_hashed>>(kHashKeyBits-(depth_-1)*kRadixTreeNodeSlotsShift))!= back_->saved_pattern_)
           continue;
        tmp_value = (Value_t)TO_VALUE(back_->pairs_[i].fir_crossed_kv_, back_->pairs_[i].sed_crossed_kv_);
        seg_pos = (tmp_key_hashed>>(kHashKeyBits - depth_*kRadixTreeNodeSlotsShift))%8;
        if(IS_SEGMENT(slots_[seg_pos])) {
            TO_SEGMENT(slots_[seg_pos])->Insert4Recovery(tmp_key, tmp_value, \
                (tmp_key_hashed&kBucketMask)*kNumPairPerBucket);
        }
    }

#ifndef TEST_RECOVERY
    back_ = nullptr;
#endif

#endif
#endif

}

double RadixTree::GetLoadfactor() {
    unsigned long capacity = 0;
    unsigned long key_value_cnt = 0;
    // BFS
    std::queue<unsigned long> queue;
    queue.push((unsigned long)root_|kRadixTreeNodeFlag);
    unsigned long front_value;
    RadixTreeNode *radix_tree_node;
    unsigned long tmp;
    while(queue.empty() == false) {
        front_value = queue.front();
        queue.pop();
        if(IS_RADIX_TREE_NODE(front_value)) {
            radix_tree_node = TO_RADIX_TREE_NODE(front_value);
            for(int i=0; i<pow(2, kRadixTreeNodeSlotsShift); i++) {
                queue.push(radix_tree_node->slots_[i]);
                if(IS_SEGMENT(radix_tree_node->slots_[i]))
                    i += pow(2, TO_SEGMENT(radix_tree_node->slots_[i])->minus_depth_)-1;
            }
        }else {
            capacity += kSegmentSlotsNum + kSegmentStashSlotsNum;
            key_value_cnt += TO_SEGMENT(front_value)->GetSize();
        }
    }
    return key_value_cnt*1.0/capacity;
}

void RadixTree::RevoveryThread(std::vector<RadixTreeNode*> *nodes_needed_recovery, int start, int end) {
    for(int i=start; i<end; i++)
        ((*nodes_needed_recovery)[i])->Recovery();
}

void RadixTree::Recovery(int nr_threads) {
    // BFS
    // the radix tree depth is the max depth of leaf nodes
    int radix_tree_depth = 0;
    std::vector<RadixTreeNode*> nodes_needed_recovery;
    std::queue<RadixTreeNode*> queue;
    root_->depth_ = 1;
    root_->recovery_used_ = 0;
    queue.push(root_);
    RadixTreeNode* front_value;
    unsigned long tmp;
    Timer timer;
    timer.Start();
    while(queue.empty() == false) {
        front_value = queue.front();
        queue.pop();
        if(front_value->IsLeafNode() == true) {
            recovery_leaf_node_cnt_++;
            nodes_needed_recovery.push_back(front_value);
            if(front_value->depth_ > radix_tree_depth)
                radix_tree_depth = front_value->depth_;
            continue;
        }
        for(int i=0; i<pow(2, kRadixTreeNodeSlotsShift); i++) {
            if(IS_RADIX_TREE_NODE(front_value->slots_[i])) {
                TO_RADIX_TREE_NODE(front_value->slots_[i])->depth_=front_value->depth_+1;
                TO_RADIX_TREE_NODE(front_value->slots_[i])->recovery_used_ = \
                    front_value->recovery_used_*pow(2, kRadixTreeNodeSlotsShift)+i;
                queue.push(TO_RADIX_TREE_NODE(front_value->slots_[i]));
            }
        }
    }
    //we use multiple threads to recover
    std::thread *threads[nr_threads];
    for(int i=0; i<nr_threads-1; i++) {
        threads[i] = new std::thread(RevoveryThread, &nodes_needed_recovery, \
            i*(nodes_needed_recovery.size()/nr_threads), \
            (i+1)*(nodes_needed_recovery.size()/nr_threads));
    }
    threads[nr_threads-1] = new std::thread(RevoveryThread, &nodes_needed_recovery, \
        (nr_threads-1)*(nodes_needed_recovery.size()/nr_threads), \
        nodes_needed_recovery.size());
    for(int i=0; i<nr_threads; i++)
        threads[i]->join();

    timer.Stop();
    recovery_rt_dir_time_ =  timer.GetSeconds();
    static bool is_rebuild_directory_in_dram = false;
    if(is_rebuild_directory_in_dram == false) {
        is_rebuild_directory_in_dram = true;
        Timer timer;
        timer.Start();
        // alloc memory for directory in DRAM
        dram_cache_ = new DramCache(radix_tree_depth*kRadixTreeNodeSlotsShift);
        //rebuild directory in DRAM
        for(int i=0; i<nodes_needed_recovery.size(); i++) {
            int node_start = nodes_needed_recovery[i]->recovery_used_ * \
                pow(pow(2, kRadixTreeNodeSlotsShift), radix_tree_depth-nodes_needed_recovery[i]->depth_+1);
            for(int j=0; j<pow(2, kRadixTreeNodeSlotsShift); j++) {
                if(IS_SEGMENT(nodes_needed_recovery[i]->slots_[j])) {
                    int start = node_start + j * pow(pow(2, kRadixTreeNodeSlotsShift), \
                        radix_tree_depth-nodes_needed_recovery[i]->depth_);
                    for(int k=0; k < pow(pow(2, kRadixTreeNodeSlotsShift), \
                        radix_tree_depth-nodes_needed_recovery[i]->recovery_used_); k++) {
                        TO_DRAM_CACHE_ENTRIES(dram_cache_->entries_and_depth_)[k+start] = \
                        (unsigned long)TO_SEGMENT(nodes_needed_recovery[i]->slots_[j]);
                    }
                }
            }
        }
        timer.Stop();
        rebuild_dir_time_ = timer.GetSeconds();
    }
}

#ifdef INPLACE
NdehStatus Segment::Insert(Key_t key, Value_t value, \
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
    unsigned long expected_sed_crossed_kv;
    unsigned long last_four_bytes_key;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        expected_sed_crossed_kv = pairs_[slot].sed_crossed_kv_;
        last_four_bytes_key = FIR_FOUR_BYTES(expected_sed_crossed_kv);
        if((h(&last_four_bytes_key, sizeof(Key_t))>> \
            (kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_)))!=pattern_ \
            && expected_sed_crossed_kv != CROSSED_KV_SENTINEL \
            && (CAS(&pairs_[slot].sed_crossed_kv_, &expected_sed_crossed_kv, CROSSED_KV_SENTINEL))) {
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            shared_mutex_p->unlock_shared();
            return Ok;
        }
        if((IS_EMPTY(pairs_[slot].sed_crossed_kv_)==true) \
            && (CAS(&pairs_[slot].sed_crossed_kv_, &invalid_crossed_kv, CROSSED_KV_SENTINEL))){
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    //find space in stash slots
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i + kSegmentSlotsNum;
        expected_sed_crossed_kv = pairs_[slot].sed_crossed_kv_;
        last_four_bytes_key = FIR_FOUR_BYTES(expected_sed_crossed_kv);
        if((h(&last_four_bytes_key, sizeof(Key_t))>> \
            (kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_)))!=pattern_ \
            && expected_sed_crossed_kv != CROSSED_KV_SENTINEL \
            && (CAS(&pairs_[slot].sed_crossed_kv_, &expected_sed_crossed_kv, CROSSED_KV_SENTINEL))) {
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            shared_mutex_p->unlock_shared();
            return Ok;
        }
        if((IS_EMPTY(pairs_[slot].sed_crossed_kv_)==true) \
            && (CAS(&pairs_[slot].sed_crossed_kv_, &invalid_crossed_kv, CROSSED_KV_SENTINEL))){
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    shared_mutex_p->unlock_shared();
    return ErrorNotFindSpace;
}
#else
NdehStatus Segment::Insert(Key_t key, Value_t value, \
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

NdehStatus Segment::Insert4Split(Key_t key, Value_t value, \
    unsigned long segment_pairs_idx) {
    int slot = 0;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(IS_EMPTY(pairs_[slot].sed_crossed_kv_)==true){
            pairs_[slot].fir_crossed_kv_ = TO_FIR_CROSSED_KV(key, value);
            pairs_[slot].sed_crossed_kv_ = TO_SED_CROSSED_KV(key, value);
            return Ok;
        }
    }
    return ErrorNotFindSpace;
}

inline void Segment::SetParent(RadixTreeNode *parent) {
    parent_ = parent;
}

void Segment::Clear() {
    for(int i=0; i<kSegmentSlotsNum+kSegmentStashSlotsNum; i++) {
        if(IS_EMPTY(pairs_[i].sed_crossed_kv_) == false) {
            pairs_[i].fir_crossed_kv_ = INVALID_CROSSED_KV;
            pairs_[i].sed_crossed_kv_ = INVALID_CROSSED_KV;
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
    if((key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_))) \
        != pattern_) {
        shared_mutex_p->unlock();
        return nullptr;
    }
    unsigned long old_pattern = pattern_;
    int old_minus_depth = minus_depth_;
    int return_seg_cnt = old_minus_depth == 0 ? 3 : 2;
    Segment** split = new Segment*[return_seg_cnt];
    if(return_seg_cnt == 3) {
        split[0] = Segment::NewSegment(parent);
        split[0]->is_splitting_ = true;
        split[1] = Segment::NewSegment(parent);
        split[1]->is_splitting_ = true;
        split[2] = this;
        split[2]->is_splitting_ = true;
        split[2]->saved_pattern_ = split[2]->pattern_;
        split[2]->pattern_ = INVALID_PATTERN;       //means we won't use this segment anymore
    }else {
        split[0] = this;
        split[0]->is_splitting_ = true;
        split[1] = Segment::NewSegment(parent);
        split[1]->is_splitting_ = true;
    }
    //update s[0], s[1] patterns
    split[0]->pattern_ = (key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_+1)+1)) << 1;
    split[1]->pattern_ = ((key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_+1)+1)) << 1) + 1;

    //update minus_depth_
    split[0]->minus_depth_ = (old_minus_depth-1+kRadixTreeNodeSlotsShift)%kRadixTreeNodeSlotsShift;
    split[1]->minus_depth_ = (old_minus_depth-1+kRadixTreeNodeSlotsShift)%kRadixTreeNodeSlotsShift;

    unsigned long tmp_key_hashed = 0;
    Key_t tmp_key;
    unsigned long last_four_bytes_key;
    Value_t tmp_value;
    NdehStatus reval;
    for(int i = 0; i<kSegmentSlotsNum+kSegmentStashSlotsNum; i++) {
        // we do not need to move invalid key-value to new segments
        if(IS_EMPTY(pairs_[i].sed_crossed_kv_) == true)
            continue;
        tmp_key = TO_KEY(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
        last_four_bytes_key = LAST_FOUR_BYTES(tmp_key);
        tmp_key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        if((tmp_key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-old_minus_depth))) \
            != old_pattern)
           continue;
        tmp_value = (Value_t)TO_VALUE(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
        if(tmp_key_hashed&((unsigned long)1<<(kHashKeyBits-(kRadixTreeNodeSlotsShift*parent_->depth_-old_minus_depth)-1)))
            reval = split[1]->Insert4Split(tmp_key, tmp_value, (tmp_key_hashed&kBucketMask)*kNumPairPerBucket);
        else if(return_seg_cnt == 3)
            reval = split[0]->Insert4Split(tmp_key, tmp_value, (tmp_key_hashed&kBucketMask)*kNumPairPerBucket);
        if(__glibc_unlikely(reval == ErrorNotFindSpace)) {
            split[1]->Clear();
            if(return_seg_cnt == 3)
                split[0]->Clear();
            goto direct_insert;
        }
    }

split_finished:
    if(return_seg_cnt == 3) {
        clflush_seg_cnt_++;
        clflush((char*)split[2], sizeof(Segment));
    }else if(parent_->back_ == nullptr) {
            clflush_seg_cnt_++;
            clflush((char*)split[1], sizeof(Segment));
            clflush((char*)&split[0]->minus_depth_, sizeof(int));
    }
    //shared_mutex_p->unlock();
    return split;

direct_insert:
    for(int i = 0; i<kSegmentSlotsNum+kSegmentStashSlotsNum; i++) {
        // we do not need to move invalid key-value to new segments
        if(IS_EMPTY(pairs_[i].sed_crossed_kv_) == true)
            continue;
        tmp_key = TO_KEY(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
        last_four_bytes_key = LAST_FOUR_BYTES(tmp_key);
        tmp_key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
        if((tmp_key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-old_minus_depth))) \
            != old_pattern)
           continue;
        if(tmp_key_hashed&((unsigned long)1<<(kHashKeyBits-(kRadixTreeNodeSlotsShift*parent_->depth_-old_minus_depth)-1))) {
            split[1]->pairs_[i].fir_crossed_kv_ = pairs_[i].fir_crossed_kv_;
            split[1]->pairs_[i].sed_crossed_kv_ = pairs_[i].sed_crossed_kv_;
        }else if(return_seg_cnt == 3) {
            split[0]->pairs_[i].fir_crossed_kv_ = pairs_[i].fir_crossed_kv_;
            split[0]->pairs_[i].sed_crossed_kv_ = pairs_[i].sed_crossed_kv_;
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
    NdehStatus reval;
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

int Segment::GetSize() {
    int cnt = 0;
    Key_t tmp_key;
    unsigned long last_four_bytes_key;
    unsigned long tmp_pattern;
    unsigned long key_hashed;
    for(int i=0; i<kSegmentSlotsNum+kSegmentStashSlotsNum; i++) {
        if(IS_EMPTY(pairs_[i].sed_crossed_kv_)==false) {
            tmp_key = TO_KEY(pairs_[i].fir_crossed_kv_, pairs_[i].sed_crossed_kv_);
            last_four_bytes_key = LAST_FOUR_BYTES(tmp_key);
            key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
            tmp_pattern = key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_));
            if(tmp_pattern==pattern_)
                cnt++;
        }
    }
    return cnt;
}

void Ndeh::Insert(Key_t key, Value_t value) {
retry:
    unsigned long last_four_bytes_key = LAST_FOUR_BYTES(key);
    unsigned long key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
    unsigned long bucket_idx = (key_hashed&kBucketMask);
    unsigned long segment_pairs_idx = bucket_idx*kNumPairPerBucket;

    __int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&dram_cache_->entries_and_depth_);
    /*__int128 dram_cache_entries_and_depth = dram_cache_->entries_and_depth_;
    if(!CHECK_VERSION(dram_cache_entries_and_depth))
        goto retry;*/
    unsigned long *dram_cache_entries = TO_DRAM_CACHE_ENTRIES(dram_cache_entries_and_depth);
    unsigned long dram_cache_depth = TO_DRAM_CACHE_DEPTH(dram_cache_entries_and_depth);
    unsigned long segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-dram_cache_depth);
    Segment *seg = (Segment*)dram_cache_entries[segment_idx_in_dram_cache];
    NdehStatus status = seg->Insert(key, value, segment_pairs_idx, key_hashed);
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

        timer.Start();
        mutex_.lock();
        //reget the entries and depth, becasuse it may has been changed.
        dram_cache_depth = TO_DRAM_CACHE_DEPTH(ATOMIC_LOAD(&dram_cache_->entries_and_depth_));
        segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-dram_cache_depth);
        
        int old_seg_minus_depth = (s[1]->minus_depth_+1)%kRadixTreeNodeSlotsShift;
        unsigned long seg_depth = seg->parent_->depth_*kRadixTreeNodeSlotsShift-old_seg_minus_depth;
        unsigned long segment_idx_in_radix_tree = key_hashed>>(kHashKeyBits-seg_depth);
        if(seg_depth < dram_cache_depth) {
            // need to allocate a new radix tree node
            if(old_seg_minus_depth == 0) {
                RadixTreeNode *new_node = new RadixTreeNode(seg->parent_->depth_+1);
                unsigned long parent_slot_idx = segment_idx_in_radix_tree&kRadixTreeNodeSlotMark;
                new_node->back_ = s[2];
                for(int i=0; i<2; i++)
                    s[i]->SetParent(new_node);
                // update redix tree node
                new_node->SetNewSegments(s, 0, false);
                clflush((char*)new_node, sizeof(RadixTreeNode));
                // the clflush operation of back_ shouldn't be added to clflush_dir_cnt_
                clflush_dir_cnt_++;
                // update the parent's slots
                seg->parent_->slots_[parent_slot_idx] = (unsigned long)new_node|kRadixTreeNodeFlag;
                clflush((char*)(seg->parent_->slots_), sizeof(unsigned long));
                clflush_dir_cnt_++;
            }else {
                for(int i=0; i<2; i++)
                    s[i]->SetParent(seg->parent_);
                //update radix tree node
                seg->parent_->SetNewSegments(s, (segment_idx_in_radix_tree<<old_seg_minus_depth)&kRadixTreeNodeSlotMark, true);
                clflush((char*)seg->parent_->slots_, sizeof(kCacheLineSize));
                clflush_dir_cnt_++;
            }
            //updating dram cache first won't cause unconsistency.
            dram_cache_->SetNewSegments(s, segment_idx_in_dram_cache&((~0x0)<<(dram_cache_depth-seg_depth)), \
            dram_cache_depth-(seg_depth+1));
        }else {
            // need to allocate a new radix tree node
            if(old_seg_minus_depth == 0) {
                RadixTreeNode *new_node = new RadixTreeNode(seg->parent_->depth_+1);
                unsigned long parent_slot_idx = segment_idx_in_radix_tree&kRadixTreeNodeSlotMark;
                new_node->back_ = s[2];
                for(int i=0; i<2; i++)
                    s[i]->SetParent(new_node);
                // update redix tree node
                new_node->SetNewSegments(s, 0, false);
                clflush((char*)new_node, sizeof(RadixTreeNode));
                // the clflush operation of back_ shouldn't be added to clflush_dir_cnt_
                clflush_dir_cnt_++;
                // update the parent's slots
                seg->parent_->slots_[parent_slot_idx] = (unsigned long)new_node|kRadixTreeNodeFlag;
                clflush((char*)(seg->parent_->slots_), sizeof(unsigned long));
                clflush_dir_cnt_++;
            }else {
                for(int i=0; i<2; i++)
                    s[i]->SetParent(seg->parent_);
                //update radix tree node
                seg->parent_->SetNewSegments(s, (segment_idx_in_radix_tree<<old_seg_minus_depth)&kRadixTreeNodeSlotMark, true);
                clflush((char*)seg->parent_->slots_, sizeof(kCacheLineSize));
                clflush_dir_cnt_++;
            }
            // double the dram cache
            unsigned long new_dram_cache_depth = dram_cache_->Double();
            dir_doubling_cnt_++;
            unsigned long new_segment_idx_in_dram_cache = key_hashed>>(kHashKeyBits-new_dram_cache_depth);
            dram_cache_->SetNewSegments(s, \
                new_segment_idx_in_dram_cache&((~0x0)<<(new_dram_cache_depth-seg_depth)), \
                new_dram_cache_depth-(seg_depth+1));
        }
        mutex_.unlock();
#ifdef INPLACE
        // COW won't use seg anymore, so do nothing.
        s[0]->is_splitting_ = false;
        s[1]->is_splitting_ = false;
        seg->shared_mutex_p->unlock();
#endif
        timer.Stop();
        breakdown_maintenance_dir_ += timer.GetSeconds();
        goto retry;
    }else if(status == SegmentSplitting || status == WrongPattern) {
        goto retry;
    }
}

unsigned long Ndeh::StatisticGetDirectoryOccupiedSpace() {
    return dram_cache_->GetOccupiedSpace();
}

double Ndeh::StatisticGetLoadfactor() {
    return radix_tree_->GetLoadfactor();
}

// here is a bug, when delete the segment that is seached.
Value_t Ndeh::Get(Key_t key) {
    unsigned long last_four_bytes_key = LAST_FOUR_BYTES(key);
    unsigned long key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
    unsigned long bucket_idx = (key_hashed&kBucketMask);
    unsigned long segment_pairs_idx = bucket_idx*kNumPairPerBucket;
retry:
   // __int128 dram_cache_entries_and_depth = ATOMIC_LOAD(&dram_cache_->entries_and_depth_);
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
    if((key_hashed>>(kHashKeyBits-(seg->parent_->depth_*kRadixTreeNodeSlotsShift-seg->minus_depth_))) \
        != seg->pattern_) {
        seg->shared_mutex_p->unlock_shared();
        goto retry;
    }
#endif
    int slot = 0;
    // the direction of insert and get is opposite, because of update
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(TO_KEY(seg->pairs_[slot].fir_crossed_kv_, seg->pairs_[slot].sed_crossed_kv_) == key) {
#ifdef INPLACE
            seg->shared_mutex_p->unlock_shared();
#endif
            return (Value_t)TO_VALUE(seg->pairs_[slot].fir_crossed_kv_, seg->pairs_[slot].sed_crossed_kv_);
        }
    }
    // stash 
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i+kSegmentSlotsNum;
        if(TO_KEY(seg->pairs_[slot].fir_crossed_kv_, seg->pairs_[slot].sed_crossed_kv_) == key) {
#ifdef INPLACE
            seg->shared_mutex_p->unlock_shared();
#endif
            return (Value_t)TO_VALUE(seg->pairs_[slot].fir_crossed_kv_, seg->pairs_[slot].sed_crossed_kv_);
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
NdehStatus Segment::Delete(Key_t key, unsigned long segment_pairs_idx, unsigned long key_hashed) {
    if(shared_mutex_p->try_lock_shared() == false)
        return SegmentSplitting;
    if((key_hashed>>(kHashKeyBits-(parent_->depth_*kRadixTreeNodeSlotsShift-minus_depth_))) != pattern_) {
        shared_mutex_p->unlock_shared();
        return WrongPattern;
    }
    int slot = 0;
    for(int i=0; i<kNumBucketToLinearProbe*kNumPairPerBucket; i++) {
        slot = (i+segment_pairs_idx)%kSegmentSlotsNum;
        if(TO_KEY(pairs_[slot].fir_crossed_kv_, pairs_[slot].sed_crossed_kv_) == key){
            pairs_[slot].sed_crossed_kv_ = INVALID_CROSSED_KV;
            shared_mutex_p->unlock_shared();
            return Ok;
        }
    }
    for(int i=0; i<kSegmentStashSlotsNum; i++) {
        slot = i+kSegmentSlotsNum;
        if(TO_KEY(pairs_[slot].fir_crossed_kv_, pairs_[slot].sed_crossed_kv_) == key){
            pairs_[slot].sed_crossed_kv_ = INVALID_CROSSED_KV;
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
bool Ndeh::Delete(Key_t key) {
    unsigned long last_four_bytes_key = LAST_FOUR_BYTES(key);
    unsigned long key_hashed = h(&last_four_bytes_key, sizeof(Key_t));
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
    NdehStatus status = seg->Delete(key, segment_pairs_idx, key_hashed);
    if(status == SegmentSplitting)
        goto retry;
    else
        return status == Ok ? true: false; 
}

void Ndeh::Recovery(int nr_threads) {
    radix_tree_->Recovery(nr_threads);
}

} // namespace ndeh
