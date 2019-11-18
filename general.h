#ifndef PROJECT_GENERAL_H
#define PROJECT_GENERAL_H

#include <iostream>
#include <fstream>
#include <ctime>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<vector<size_t>> umat;
typedef vector<size_t> uvec;
typedef unsigned long long uword;

struct comb
{
    int val;
    int counter;
    comb() : val(0), counter(1) {};
    comb(int val) : val(val), counter(1) {};
    comb(int val, int counter) : val(val), counter(counter) {};

    comb operator^(const comb & x)
    {
        comb res(val ^ x.val, x.counter + counter);
        return res;
    }
};

enum flagtype {
    CAST_TO_ONES,
    RETURN_INDICES,
    REVERSE
};

const size_t n = 6;
const size_t k = 20;
const uword linear_function_count = uword(pow(2, n + 1));
const uword number_of_sets = uword(pow(2, n));
const uword number_of_function_values = pow(2, k - 1);
const uword max_number = uword(pow(2, number_of_sets));

template <typename T>
ostream & operator<<(ostream & os, const vector<T> & v);

template <typename T>
ofstream & operator<<(ofstream & os, const vector<T> & v);

#endif
