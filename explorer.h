#ifndef PROJECT_EXPLORER_H
#define PROJECT_EXPLORER_H

#include "general.h"

struct Explorer
{
    umat L_VALUES, F_VALUES, BINS_WITH_1, coins;
    uvec all_alphas, number_of_coins, sets_with_1, F, lcomb, indices5, indices6;
    uword F_VALUES_SIZE, q_start, q, q_max; //=463828
    ofstream outfile;
    size_t ITER_SIZE;

    Explorer(uword q_start, uword q, uword q_max, string filename, size_t r);

    umat initL();
    umat initF();

    uword parityCheck(uword x);
    uvec toBinary(uword num, size_t count, flagtype flag);
    uvec toBinary(uword num, size_t count, flagtype flag, size_t & size);

    void printResults(uvec FUNCTION, uvec ALL_ALPHAS);
    void printCoins();

    bool linCheck(uword & linfunc, size_t & new_set);
    bool modify0(size_t & new_set);
    bool modify1(size_t & new_set);
    bool modify2(size_t & new_set);
    bool modify3(size_t & new_set);
    bool modify4(size_t & new_set);
    bool modify5(size_t & new_set);
    bool modify6(size_t & new_set);


    bool explore();

    bool checkF(uvec FUNCTION, uvec ALL_ALPHAS);
    double calcPk(uword q_start);
    uword howMany(uvec ALL_ALPHAS);
};

#endif
