#ifndef HEAP_H
#define HEAP_H

#include "util.h"
#include "globals.h"

#include <iostream>
#include <vector>
#include <queue>

#include <utility>
#include <functional>
#include <unordered_map>
#include <map>

#include <random>

using namespace std;

typedef pair<double,int> eipair; // energy index pair


// print pair
template<typename A, typename B>
void print_pair(const std::pair<A,B> &p)
{
    cout<<"("<<p.first <<","<<p.second<<") ";
}

// return the queue as vector
template<typename A, typename B>
void print_pairs(vector<pair<A,B>> v) {
    for(const auto& p:v){
        print_pair(p);
    }
}

class HeapMap{
    public:

        /* 
        min heap of (du,idx) pair, only insert du<0
         when element is delete/modified, do not modify pq
         check valid element at popping
         use idxmap to keep track
         */
        priority_queue<eipair, vector<eipair>, greater<eipair>> pq;

        GridData grid;
        
        /*
        a map that keeps track of latest valid idx-du pair, where idx is interface points
        du maybe positive
        */
        unordered_map<int,double> idxmap;
        

        /*
        reserve memory. Seems not helpful 
        https://stackoverflow.com/questions/29235978/how-to-preallocatereserve-a-priority-queuevector
        */
        void reserve(int n){
            std::vector<eipair> container;
            container.reserve(n);
            std::priority_queue<eipair, std::vector<eipair>, greater<eipair>> pq2 (greater<eipair>(), move(container));
            pq = pq2;
            idxmap.reserve(n);
        }
        /*
        insert element (idx, du), if du is updated, we insert a new copy to the queue, and update map[idx]=du 
        */ 

        
        void insert(int idx, double du){
            if(du<0){
                pq.emplace(du,idx);
            }
            idxmap[idx] = du;
            
        }

        /*
        remove element, some index may get out of interface
        */
        void remove(int idx){
            idxmap.erase(idx);
        }

        /*
        check if empty
        */
        bool empty(){
            return pq.empty();
        }



        /*
        check if ind is valid: in the map, also the value in queue is the same as that in the map
        */
        bool valid(eipair x){
            // auto [u,idx] = x;
            double u = x.first;
            double idx = x.second;
            auto search = idxmap.find(idx);
            if (search != idxmap.end()){
                // if idxmap contain
                if(  abs(idxmap[idx] - u) < 1e-12 ){
                    return true;
                }
            }
            return false;
        }


        bool inmap(int idx){
            return idxmap.find(idx)!=idxmap.end();
        }
        

        /*
        return the (idx, du) with lowest energy in the the idxmap
        look at next element, if not valid, remove.
        */

        void trim(){
            while( !pq.empty()){
                eipair temp = pq.top();
                if (!valid(temp)){
                    pq.pop();
                }else{
                    return;
                }
            }
            return;
        }

        eipair pop(){
            while( !pq.empty()){
                eipair temp = pq.top();
                pq.pop();
                if (valid(temp)){
                    remove(temp.second);
                    trim();
                    return temp;        
                }
            }
            std::cout<<"pop empty queue"<<endl;
            return {0.0,-1};
        }

        /*
        return element (idx, du) with negative du as (du,idx), should be the same as heap
        */

        vector<eipair> GetNegMap(){
            vector<eipair> y;
            for(auto const& x:idxmap){
                if (x.second<0){
                    y.push_back(make_pair(x.second,x.first));    
                }
            }
            sort(y.begin(),y.end());
            return y;
        }

        
        void PrintMap(){
            vector<eipair> y = GetNegMap();
            print_pairs(y);
        }


        /*
        return the queue as vector of (du,idx) pairs
        */
        vector<eipair> GetQueue(){
            HeapMap tmp = *this;
            vector<eipair> y;
            while(!tmp.empty()) {
                y.push_back(tmp.pop());
            }

            if (y[0].second==-1){
                y.erase(y.begin());
            }

            return y;
        }

        void PrintQueue(){
            vector<eipair> y = GetQueue();
            print_pairs(y);
        }
};



// return the queue as vector
template<typename T>
vector<eipair> get_queue(T q) {
    vector<eipair> y;
    while(!q.empty()) {
        y.push_back(q.pop());
    }
    return y;
}


// convert elements in a map to a vector
// sort by key
template<template <typename...> class Map, typename K, typename V>
vector<pair<K,V>> get_map(const Map<K,V> &q){
    vector<pair<K,V>> y;
    for(auto const& x:q){
        y.push_back(make_pair(x.first, x.second));
    }
    sort(y.begin(),y.end());
    return y;
}



// flip a pair
template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

// flip a map
//https://stackoverflow.com/questions/33315837/making-a-template-work-with-both-stdmap-and-stdunordered-map
template<template <typename...> class Map, typename K, typename V>
Map<V, K> flip_map(const Map<K, V>& map)
{
    Map<V, K> result;
    for (const auto& p : map) result.emplace(p.second, p.first);
    return result;
}


// print the queue
template<typename T>
void print_queue(T q) { // NB: pass by value so the print uses a copy
    while(!q.empty()) {
        std::cout << q.top() << ' ';
        q.pop();
    }
    std::cout << '\n';
}

// print the queue
template<typename T>
void print_map(T &q) { // NB: pass by value so the print uses a copy
    for(auto const& x : q){
        std::cout <<"(" <<x.first<<":"<<x.second<<") ";
    }
    std::cout << '\n';
}



#endif