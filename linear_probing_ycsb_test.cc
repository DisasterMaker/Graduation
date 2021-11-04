#ifndef LINEAR_HASH_H_
#include "linear_probing.h"
#endif

#include <iostream>
#include <thread>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <fstream>
#include <time.h>
#include <string>

#define NR_LOAD         100000000        // the number of operations in load_ycsb
#define NR_OPERATIONS   16000000    // the data you generated must be greater or equal to this
#define LOAD_YCSB      "mix-load-1"      //  ycsb load file
#define RUN_YCSB       "mix-run-0.16"        // ycsb workload file


#define NR_THREADS      4       // must be equal or less nr_cpus

double *times;
uint64_t kWriteLatencyInNS = 0;
uint64_t clflushCount = 0;


double mysecond() {
    struct timeval tp;
    int i;
    i = gettimeofday(&tp, NULL);
    return ( (double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

void insert(LinearProbingHash *linear_probing, Key_t *keys, Value_t *values, int cur_nr_threads) {
   /// double start = mysecond();
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        linear_probing->Insert(keys[i], values[i]);
    }
    //times[thread_id] = mysecond() - start;
}

void get(LinearProbingHash *linear_probing, Key_t *keys, Value_t *values, int cur_nr_threads) {
    Value_t tmp_value;
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        tmp_value = linear_probing->Get(keys[i]);
        if(__glibc_unlikely(tmp_value != values[i])){
             std::cout << "get failed, key: " << (unsigned long)keys[i] << \
             " value: " << (unsigned long)values[i] << \
             " wrong value: " << (unsigned long)tmp_value << std::endl;
        }
    }
}

void lpoperation(LinearProbingHash *linear_probing,Key_t *operationtype, Key_t *keys, Value_t *values, int cur_nr_threads){
    for(int i=0; i<NR_OPERATIONS/cur_nr_threads; i++) {
        Value_t tmp;
        if(operationtype[i] == 1)
        linear_probing->Insert(keys[i], values[i]);
        else
        {
            tmp = linear_probing->Get(keys[i]);
            if(__glibc_unlikely(tmp != values[i])){
             std::cout << "get failed, key: " << (unsigned long)keys[i] << \
             " value: " << (unsigned long)values[i] << \
             " wrong value: " << (unsigned long)tmp << std::endl;
            }
        }
    }
}

void print_result(int thread_num, LinearProbingHash *linear_probing) {
    double total_time = times[1] - times[0];
    std::cout <<RUN_YCSB <<": thread num: " << thread_num;
    std::cout << ", total time: " << total_time << \
    ", OPS: " << NR_OPERATIONS*1.0/(total_time) << \
    ", load factor: " << linear_probing->Utilization() << std::endl;

}


int main(int argc, char* argv[])
{
    //get operation type and keys from load_ycsb
    init_vmp();
    if (argv[1]==NULL) {
    std::cout << "please input thread num" << std::endl;
    return 1;
    }
    Key_t *loadkeys = new Key_t[NR_LOAD];
    Key_t *operationtypes = new Key_t[NR_OPERATIONS];
    Value_t *loadvalues = new Value_t[NR_LOAD];
    std::string tmp_string;
    std::ifstream ifs;
    ifs.open(LOAD_YCSB);
    unsigned long tmp ;
    for(int i=0; i<NR_LOAD; ++i) {
        ifs >> tmp_string;
        ifs >> loadkeys[i];
        loadvalues[i]=(Value_t)loadkeys[i];
    }
    ifs.close();
   // delete[] operationtypes;

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
    LinearProbingHash *linear_probing;

    times = new double [2];
    int nr_cpus = get_nprocs_conf();
    std::thread *threads[nr_cpus];
    int thread_idx;
    thread_idx = atoi(argv[1]);
    //insert load_ycsb
    std::cout << "begin load phase: " << std::endl;
       //int thread_idx=1;
        linear_probing = new LinearProbingHash(16);
        // for(int i=0; i<=thread_idx; i++) {
        //     threads[i] = new std::thread(insert, linear_probing,  loadkeys+i*NR_OPERATIONS/(thread_idx+1),
        //                   loadvalues + i*NR_OPERATIONS/(thread_idx+1),thread_idx+1);
        // }
        // for(int i=0; i<=thread_idx; i++) {
        //     threads[i]->join();
        // }
        // single thread is ok, we don't need multiple thread when loading
        for(int i=0; i<NR_LOAD; i++) {
            linear_probing->Insert(loadkeys[i], loadvalues[i]);
        }
        //run ycsb
        std::cout<<"run ycsb test"<<std::endl;
        times[0] = mysecond();
        for(long long int i=0; i<thread_idx; i++)
            threads[i] = new std::thread(lpoperation, linear_probing, operationtypes+i*NR_OPERATIONS/(thread_idx) ,runkeys+i*NR_OPERATIONS/(thread_idx),
                runvalues + i*NR_OPERATIONS/(thread_idx), thread_idx);
        for(int i=0; i<thread_idx; i++)
            threads[i]->join();
        times[1] = mysecond();
        print_result(thread_idx, linear_probing);
        //delete linear_probing;

    return 0;
}
