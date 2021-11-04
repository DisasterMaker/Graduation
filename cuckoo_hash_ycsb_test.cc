#include <libcuckoo/cuckoohash_map.hh>

#include <iostream>
#include <thread>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <fstream>
#include <time.h>
#include "../../util/pair.h"

#ifndef CONFIG_SCAS_COUNTER
#define CONFIG_SCAS_COUNTER
#endif

#define NR_LOAD         80000000       // the number of operations in load_ycsb
#define NR_OPERATIONS   80000000    // the data you generated must be greater or equal to this
#define LOAD_YCSB      "/home/eric/test/others/ndeh-multi-thread/ndeh/src/rw-50-50-load-0.8"      //  ycsb load file
#define RUN_YCSB       "/home/eric/test/others/ndeh-multi-thread/ndeh/src/rw-50-50-run-0.8"       // ycsb workload file

#define NR_THREADS      16        // must be equal or less nr_cpus


double *times;
uint64_t kWriteLatencyInNS = 0;
uint64_t clflushCount = 0;

double mysecond() {
    struct timeval tp;
    int i;
    i = gettimeofday(&tp, NULL);
    return ( (double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

void insert(cuckoohash_map<Key_t, Value_t> *cuckoo_hash, Key_t *keys, Value_t *values, int cur_nr_threads) {
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        cuckoo_hash->insert(keys[i], values[i]);
    }
}

void get(cuckoohash_map<Key_t, Value_t> *cuckoo_hash, Key_t *keys, Value_t *values, int cur_nr_threads) {
    Value_t tmp_value;
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        tmp_value = cuckoo_hash->find(keys[i]);
        if(__glibc_unlikely(tmp_value != values[i])){
             std::cout << "get failed, key: " << (unsigned long)keys[i] << \
             " value: " << (unsigned long)values[i] << \
             " wrong value: " << (unsigned long)tmp_value << std::endl;
        }
    }
}

void cuckoperation(cuckoohash_map<Key_t, Value_t> *cuckoo_hash, Key_t *operationtype, Key_t *keys, Value_t *values, int cur_nr_threads){
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        Value_t tmp;
        if(operationtype[i] == 1)
        cuckoo_hash->insert(keys[i], values[i]);
        else
        {
            tmp = cuckoo_hash->find(keys[i]);
            if(__glibc_unlikely(tmp != values[i])){
            std::cout << "get failed, key: " << (unsigned long)keys[i] << \
            " value: " << (unsigned long)values[i] << \
            " wrong value: " << (unsigned long)tmp << std::endl;
            }
        }
    }
}

void print_result(int thread_num, cuckoohash_map<Key_t, Value_t> *cuckoo_hash) {
    double total_time = times[1] - times[0];
    std::cout << "thread num: " << thread_num+1;
    std::cout << ", total time: " << total_time << \
    ", OPS: " << NR_OPERATIONS*1.0/(total_time) <<  std::endl;
}

int main(int argc, char* argv[]) {

    // register_printcs_with_signal(SIGSEGV);
    // get keys and values from file
    //init_vmp();
    if (argv[1]==NULL) {
    std::cout << "please input thread num" << std::endl;
    return 1;
    }
    Key_t *loadkeys = new Key_t[NR_OPERATIONS];
    Key_t *operationtypes = new Key_t[NR_OPERATIONS];
    Value_t *loadvalues = new Value_t[NR_OPERATIONS];
    std::string tmp_string;
    std::ifstream ifs;
    ifs.open(LOAD_YCSB);
    unsigned long tmp ;
    for(int i=0; i<NR_OPERATIONS; ++i) {
        ifs >> tmp_string;
        ifs >> loadkeys[i];
        loadvalues[i]=(Value_t)loadkeys[i];
    }
    ifs.close();
    
    //get operation types and keys from ycsb_workload
    Key_t *runkeys = new Key_t[NR_OPERATIONS];
    Value_t *runvalues = new Value_t[NR_OPERATIONS];
    std::string operation;
    std::string insert_op = "insert";
    std::string read_op = "read";
    ifs.open(RUN_YCSB);
     for(int i=0; i<NR_OPERATIONS; ++i) {
        ifs >> operation;
        if(operation.compare(insert_op) == 0)
           operationtypes[i]=1;
        else if (operation.compare(read_op) == 0)
        {
             operationtypes[i] = 2;
        }
        else
        {
              perror("undefined operation type!");
        }
        ifs >> runkeys[i];
        runvalues[i]=(Value_t)runkeys[i];
    }
    ifs.close();


    int nr_cpus = get_nprocs_conf();
    std::thread *threads[nr_cpus];
    int thread_idx;
    thread_idx = atoi(argv[1]);
    times = new double [2];
    
    cuckoohash_map<Key_t, Value_t> *cuckoo_hash;
     //test insert
    std::cout << "test insert: " << std::endl;
    cuckoo_hash = new cuckoohash_map<Key_t, Value_t> ();
    //insert load_ycsb
    for(int i=0; i<NR_OPERATIONS; i++) {
            cuckoo_hash->insert(loadkeys[i], loadvalues[i]);
    }
    times[0] = mysecond();
    for(long long int i=0; i<thread_idx; i++)
        threads[i] = new std::thread(cuckoperation, cuckoo_hash, operationtypes+i*NR_OPERATIONS/thread_idx ,runkeys+i*NR_OPERATIONS/thread_idx,
                runvalues + i*NR_OPERATIONS/thread_idx, thread_idx);
    for(int i=0; i<thread_idx; i++)
        threads[i]->join();
    times[1] = mysecond();
    print_result(thread_idx, cuckoo_hash);
    delete cuckoo_hash;

    return 0; 
}
