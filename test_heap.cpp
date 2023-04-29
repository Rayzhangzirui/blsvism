#include "heap2.h"


#include <random>
std::random_device dev;
std::mt19937 rng(dev());

#include <gtest/gtest.h>


using namespace std;

ostream& operator<<(ostream& stream, const eipair& x){
    return stream<<"("<<x.first<<","<<x.second<<")";
}


// the queue match the map
void match(HeapMap &hm){
    vector<eipair> from_map = hm.GetNegMap();
    vector<eipair> from_heap = hm.GetQueue();

    EXPECT_EQ(from_map.size(),from_heap.size());
    for(int i = 0; i < from_map.size(); i ++){
        // eipair ei = hm.pop();
        EXPECT_EQ(from_map[i].first, from_heap[i].first);
        EXPECT_EQ(from_map[i].second,from_heap[i].second);
    }

    cout<<"map:"<<endl;
    hm.PrintMap();
    cout<<endl;
    cout<<"queue:"<<endl;
    hm.PrintQueue();
    cout<<endl;
}


// simple
TEST(heap, simple){

    
    HeapMap hm;

    hm.insert(1,-1.0);
    hm.insert(2,-2.0);
    hm.insert(3,-3.0);
    match(hm);

    hm.insert(3,-4.0);
    match(hm);

    hm.insert(1,1.0);
    hm.insert(2,1.0);
    hm.insert(3,1.0);
    match(hm);
    
}


// insert element into heap, the popped element should be in order
TEST(heap, sort){

    int N = 50; // number of element
    int m = 100; // number of sample

    uniform_real_distribution<double> rndreal(-1.0, 1.0);
    uniform_int_distribution<int> rndint(0,N-1);

    HeapMap hm;

    for(int i = 0; i < m; i++){
        double u = rndreal(rng);
        hm.insert( rndint(rng) , u );
    }

    match(hm);
    
}


/*
repeatedly insert the same key, should pop the one in the map
*/ 
TEST(heap, insert){


    int N = 1; // range of index
    int m = 10; // number of random insertion

    uniform_real_distribution<double> rndreal(-1,1);
    uniform_int_distribution<int> rndint(0,N-1);

    HeapMap hm;

    for (int i = 0; i < m; i++){
        // randomly insert/update element
        int idx = rndint(rng);
        double u = rndreal(rng);

        hm.insert(idx,u);
    }

    match(hm);

}

/*
repeatedly insert key-value pair, randomly remove key
*/ 
TEST(heap, random){

    int N = 50; // range of index
    int m = 100; // number of samples

    uniform_real_distribution<double> rndreal(-1.0, 1.0);
    uniform_int_distribution<int> rndint(0,N-1);
    bernoulli_distribution rndb(0.5);

    map<int,double> idx_u_map;

    HeapMap hm;

    for (int i = 0; i < m; i++){
        // randomly insert/update element
        int idx = rndint(rng);
        double u = rndreal(rng);
        
        idx_u_map[idx] = u;        
        hm.insert(idx,u);

        if(rndb(rng) && idx_u_map.size() > 0){
            // randomly remove element
            auto it = idx_u_map.begin();
            std::advance(it, rand() % idx_u_map.size());
            idx_u_map.erase(it->first);

            hm.remove(it->first);
        }
    }

    match(hm);

}

/*
N = range of index. m = number of samples
*/
void TestSpeed(int N, int m){
    uniform_real_distribution<double> rndreal(-1.0, 1.0);
    uniform_int_distribution<int> rndint(0,N-1);
    bernoulli_distribution rndb(0.5); // frequency of remove

    HeapMap hm;

    for (int i = 0; i < m; i++){
        // randomly insert/update element
        int idx = rndint(rng);
        double u = rndreal(rng);
        
        hm.insert(idx,u);

        if(rndb(rng)){
            // randomly remove element
            hm.remove(idx);
        }
    }

}

/*
repeatedly insert key-value pair, randomly remove key
*/ 
TEST(heap, speed){
    int N = 200*200;
    int m = N*10;
    TestSpeed(N,m);
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
