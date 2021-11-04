#ifndef NDEH_H_
#include "ndeh.h"
#endif

#ifndef PRINT_CS_H_
#include "printcs.h"
#endif

#include <iostream>
#include <thread>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <fstream>

#include<time.h>
#define random(x) (rand()%x)

#ifndef CONFIG_SCAS_COUNTER
#define CONFIG_SCAS_COUNTER
#endif

#define NR_OPERATIONS   160000000
#define INPUT_FILE      "data"       // using input_gen.cpp to generate data file

#define NR_THREADS      16        // must be equal or less nr_cpus
#define NR_RECOVERY_THREADS     1   //the default value

double *times;
unsigned long *directory_size;

uint64_t kWriteLatencyInNS = 0;
uint64_t clflushCount = 0;

double mysecond() {
    struct timeval tp;
    int i;
    i = gettimeofday(&tp, NULL);
    return ( (double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

void insert(ndeh::Ndeh *ndeh, Key_t *keys, Value_t *values, int cur_nr_threads) {
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        ndeh->Insert(keys[i], values[i]);
    }
}

void get(ndeh::Ndeh *ndeh, Key_t *keys, Value_t *values, int cur_nr_threads) {
    Value_t tmp_value;
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        tmp_value = ndeh->Get(keys[i]);
        if(__glibc_unlikely(tmp_value != values[i])){
            std::cout << "get failed, key: " << (unsigned long)keys[i] << \
            " value: " << (unsigned long)values[i] << \
            " wrong value: " << (unsigned long)tmp_value << std::endl;
        }
    }
}

void _delete(ndeh::Ndeh *ndeh, Key_t *keys, Value_t *values, int cur_nr_threads) {
    bool is_delete_success;
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        is_delete_success = ndeh->Delete(keys[i]);
        if(__glibc_unlikely(is_delete_success == false)){
            std::cout << "delele failed, key: " << (unsigned long)keys[i] << std::endl;
        }
    }
}

void print_result(int thread_num, ndeh::Ndeh *ndeh) {
    double total_time = times[1] - times[0];
    std::cout << "thread num: " << thread_num;
#ifdef INPLACE
    std::cout<< ", mode: INPLACE";
#else
    std::cout << ", mode: COW";
#endif
    std::cout <<", total time: " << total_time << \
    ", OPS: " << NR_OPERATIONS*1.0/(total_time) << \
    ", clflush segment count: " << clflush_seg_cnt_ << \
    ", clflush directory count: " << clflush_dir_cnt_ << \
    ", load factor: " << ndeh->StatisticGetLoadfactor() << std::endl;
}

#ifdef TEST_RECOVERY
void print_recovery_result(int nr_threads) {
    double total_time = times[1] - times[0];
    std::cout << "thread num: " << nr_threads;
#ifdef INPLACE
    std::cout<< ", mode: INPLACE";
#else
    std::cout << "mode: COW";
#endif
    std::cout <<", total time: " << total_time << \
    ", leaf node count: " << recovery_leaf_node_cnt_ << \
    ", clflush count: " << recovery_clflush_cnt_ << std::endl;
}
#endif

void clear_cache() {
    const size_t size = 1024*1024*512;
    int* dummy = new int[size];
    for (int i=100; i<size; i++) {
        dummy[i] = i;
    }
    for (int i=100;i<size-(1024*1024);i++) {
        dummy[i] = dummy[i-rand()%100] + dummy[i+rand()%100];
    }
    delete[] dummy;
}

int main(int argc, char *argv[]) {
    // register_printcs_with_signal(SIGSEGV);
    // get keys and values from file
    init_vmp();
    if (argv[1]==NULL) {
    std::cout << "please input thread num" << std::endl;
    return 1;
    }
    Key_t *keys = new Key_t[NR_OPERATIONS];
    Value_t *values = new Value_t[NR_OPERATIONS];
    std::ifstream ifs;
    ifs.open(INPUT_FILE);
    unsigned long tmp;
    for(long long int i=0; i<NR_OPERATIONS; ++i) {
        ifs >> keys[i];
        ifs >> tmp;
        values[i] = (Value_t)tmp;
    }

    int nr_cpus = get_nprocs_conf();
    std::thread *threads[nr_cpus];
    int thread_idx;
    thread_idx = atoi(argv[1]);

    times = new double [2];
    directory_size = new unsigned long[nr_cpus];
    ndeh::Ndeh *ndeh;
    //test insert
    std::cout << "test ndeh insert: " << std::endl;
        ndeh = new ndeh::Ndeh();
        times[0] = mysecond();
        for(long long int i=0; i<thread_idx; i++)
            threads[i] = new std::thread(insert, ndeh, keys+i*NR_OPERATIONS/(thread_idx),
                 values + i*NR_OPERATIONS/(thread_idx), thread_idx);
        for(int i=0; i<thread_idx; i++)
            threads[i]->join();
        times[1] = mysecond();
        print_result(thread_idx, ndeh);
        delete ndeh; 
    
    // test get
        std::cout << "test ndeh get: " << std::endl;
        ndeh = new ndeh::Ndeh();
        // insert data first
        for(long long int i=0; i<thread_idx; i++)
            threads[i] = new std::thread(insert, ndeh, keys+i*NR_OPERATIONS/(thread_idx),
                 values + i*NR_OPERATIONS/(thread_idx), thread_idx);
        for(int i=0; i<thread_idx; i++)
            threads[i]->join();

        times[0] = mysecond();
        for(long long int i=0; i<thread_idx; i++)
            threads[i] = new std::thread(get, ndeh, keys+i*NR_OPERATIONS/(thread_idx),
                values + i*NR_OPERATIONS/(thread_idx), thread_idx);
        for(int i=0; i<thread_idx; i++)
            threads[i]->join();
        times[1] = mysecond();
        print_result(thread_idx, ndeh);
        delete ndeh;
    
/*    // test _delete
      std::cout << "test _delete: " << std::endl;
    for(int thread_idx=0; thread_idx<NR_THREADS; thread_idx++) {
        ndeh = new ndeh::Ndeh();
        // insert data first
        for(long long int i=0; i<=thread_idx; i++)
            threads[i] = new std::thread(insert, ndeh, keys+i*NR_OPERATIONS/(thread_idx+1),
                 values + i*NR_OPERATIONS/(thread_idx+1), thread_idx+1);
        for(int i=0; i<=thread_idx; i++)
            threads[i]->join();

        times[0] = mysecond();
        for(long long int i=0; i<=thread_idx; i++)
            threads[i] = new std::thread(_delete, ndeh, keys+i*NR_OPERATIONS/(thread_idx+1),
                values + i*NR_OPERATIONS/(thread_idx+1), thread_idx+1);
        for(int i=0; i<=thread_idx; i++)
            threads[i]->join();
        times[1] = mysecond();
        print_result(thread_idx, ndeh);
        delete ndeh;
    }
#ifdef TEST_RECOVERY
    // test recovery
    std::cout << "test recovery: " << std::endl;
    // insert data first
    ndeh = new ndeh::Ndeh();
    for(int i=0; i<NR_OPERATIONS; i++) {
        ndeh->Insert(keys[i], values[i]);
    }
    //note: the data maybe cached
    for(int i=1; i<=NR_RECOVERY_THREADS; i++) {
        recovery_leaf_node_cnt_ = 0;
        recovery_clflush_cnt_ = 0;
        times[0] = mysecond();
        ndeh->Recovery(i);
        times[1] = mysecond();
        print_recovery_result(i);
    }
    std::cout << "rebuild dir in DRAM: " <<  ndeh->radix_tree_->rebuild_dir_time_ << " " << \
    "recovery rt dir:  " << ndeh->radix_tree_->recovery_rt_dir_time_ << std::endl;
    Value_t tmp_value;
    for(int i=0; i<NR_OPERATIONS; i++) {
        tmp_value = ndeh->Get(keys[i]);
        if(__glibc_unlikely(tmp_value != values[i])){
            std::cout << "get failed, key: " << (unsigned long)keys[i] << \
            " value: " << (unsigned long)values[i] << \
            " wrong value: " << (unsigned long)tmp_value << std::endl;
        }
    }
#endif
*/
    return 0;
}
