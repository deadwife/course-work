#include "explorer.h"

vector<uword> initQ(size_t CHILDNUM)
{
    vector<uword> qs(CHILDNUM);
    qs[0] = 750000343139;
    qs[1] = 69999002147;
    qs[2] = 800000344721;
    qs[3] = 900000344471;
    qs[4] = 999990000193;
    return qs;
}

vector<string> initFiles(size_t CHILDNUM)
{
    vector<string> filenames(CHILDNUM);
    filenames[0] = "../results/child0.txt";
    filenames[1] = "../results/child1.txt";
    filenames[2] = "../results/child2.txt";
    filenames[3] = "../results/child3.txt";
    filenames[4] = "../results/child4.txt";
    return filenames;
}

int main()
{
    uword q_start;
    size_t pid = getpid();
    string filename;
    size_t CHILDNUM = 5;
    vector<string> filenames = initFiles(CHILDNUM);
    vector<uword> qs = initQ(CHILDNUM);
	//father process setup
    filename = "../results/father.txt";
    q_start = 169999001801;

    //children setup
    for(size_t i = 0; i < CHILDNUM; ++i) {
        if(!fork()) {
            q_start = qs[i];
            filename = filenames[i];
            break;
        }
    }
    uword s0 = (2 * number_of_function_values - 1) / q_start + 1;
    uword s1 = (max_number - max_number / (2 * number_of_function_values)) / q_start;
    uword q_max = q_start * s1;

	uword q = q_start * s0;
    size_t r = 17;

    Explorer explorer(q_start, q, q_max, filename, r);
	//main searching executes here
    explorer.explore();

    if(pid == getpid()) {
        for(size_t i = 0; i < CHILDNUM; ++i) {
            cout << "Waiting for child number " << i << "...\n ";
            wait(NULL);
        }
    }
    else exit(0);
    return 0;
}
