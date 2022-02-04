#pragma once

#include "gettime.h"
#include "utilities.h"

using namespace std;

struct DummyProfiler{
    void incFindMinTime(double t){
    }
    void incHeapifyTime(double t){
    }
    void incFindMinNum(){
    }

    void incHeapifyNum(){
    }
    void incHeapifySize(size_t i){
    }
    void setReadTime(double t){
    }

    void setInitTime(double t){
    }
    void incVTime(double t){
    }

    void incInsertTime(double t){
    }
    void incUpdTime(double t){
    }

    void report(){
    }
};

struct Profiler{

    double readFileTime=0;
    double initTime=0;
    double getVTime=0;
    double insertTime=0;
    double updateTime=0;
    size_t findminNum=0;
    size_t heapifyNum=0;
    size_t heapifySize=0;
    double findMinTime=0;
    double heapifyTime=0;

    Profiler(){

    }

    void setReadTime(double t){
        readFileTime = t;
    }

    void setInitTime(double t){
        initTime = t;
    }

    void incVTime(double t){
        getVTime+=t;
    }

    void incInsertTime(double t){
        insertTime+=t;
    }
    void incUpdTime(double t){
        updateTime+=t;
    }

    void incFindMinNum(){
        pbbs::write_add(&findminNum, 1);
    }
    void incHeapifyNum(){
        pbbs::write_add(&heapifyNum, 1);
    }
    void incHeapifySize(size_t i){
        pbbs::write_add(&heapifySize, i);
    }
    void incFindMinTime(double t){
        pbbs::write_add(&findMinTime, t);
    }
    void incHeapifyTime(double t){
        pbbs::write_add(&heapifyTime, t);
    }

    void report(){
        // cout << "read: " << readFileTime << endl;
        cout << "init: " << initTime << endl;
        cout << "get_v_list: " << getVTime << endl;
        cout << "insert: " << insertTime << endl;
        cout << "findminNum: " << findminNum << endl;
        cout << "heapifyNum: " << heapifyNum << endl;
        cout << "heapifySize: " << heapifySize << endl;

    }
};