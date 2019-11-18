#include "explorer.h"

template <typename T>
ostream & operator<<(ostream & os, const vector<T> & v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        os << v[i];
        if (i != v.size() - 1)
            os << ",";
    }
    return os;
}

template <typename T>
ofstream & operator<<(ofstream & os, const vector<T> & v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        os << v[i];
        if (i != v.size() - 1)
            os << ",";
    }
    return os;
}

Explorer::Explorer(uword q_start, uword q, uword q_max, string filename, size_t r) :
                    BINS_WITH_1(number_of_sets, uvec(n + 1)),
                    coins(linear_function_count, uvec(n + 1, 0)),
                    all_alphas(k),
                    number_of_coins(linear_function_count, 0),
                    sets_with_1(number_of_sets),
                    F(k),
                    lcomb(linear_function_count - 1),
                    indices5(11),
                    indices6(26),
                    q_start(q_start),
                    q(q),
                    q_max(q_max),
                    ITER_SIZE(r)
{
    L_VALUES = initL();
    F_VALUES = initF();
    F_VALUES_SIZE = F_VALUES.size();
    outfile.open(filename, ios::trunc);
    for(uword i = 0; i < number_of_sets; ++i)
    {
        size_t cur_set = i | (1 << n);
        sets_with_1[i] = cur_set;
        BINS_WITH_1[i] = toBinary(cur_set, n + 1, CAST_TO_ONES);
    }
    indices5 = {4, 8, 10, 12, 16, 18, 20, 22, 24, 26, 28};
    indices6 = {4, 8, 10, 12, 16, 18, 20, 22, 24, 26, 28, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58};
}

uword Explorer::parityCheck(uword x) {
    return __builtin_parity(x);
}

umat Explorer::initL() {
    umat L(linear_function_count, uvec(number_of_sets));
    for(size_t i = 0; i < linear_function_count; ++i) {
        for(size_t j = 0; j < number_of_sets; ++j) {
            size_t a = (i & j) | (i & (1 << n));
            L[i][j] = parityCheck(a);
        }
    }
    return L;
}

umat Explorer::initF() {
    umat F(number_of_function_values, uvec(k));
    uvec f(k);
    size_t number_of_ones, number_of_zeros, counter = 0;
    for(size_t i = linear_function_count - 1; i < number_of_function_values; ++i) {
        f = toBinary(i, k, CAST_TO_ONES, number_of_ones);
        number_of_zeros = k - number_of_ones;
        if(number_of_ones >= (n + 1) && number_of_zeros >= (n + 1)) {
            F[counter] = f;
            ++counter;
        }
    }
    auto first = F.begin();
    return umat(first, first + counter);
}

uvec Explorer::toBinary(uword num, size_t count, flagtype flag)
{
    uvec binary_representation(count);
    uword one = 1;
    switch(flag)
    {
        case CAST_TO_ONES:
        {
            for(uword i = 0; i < count; ++i)
            {
                binary_representation[i] = uword(bool(num & (one << i)));
            }
            break;
        }
        case RETURN_INDICES:
        {
            uword result;
            uword counter = 0;
            for(uword i = 0; i < count; ++i)
            {
                result = num & (one << i);
                if(result)
                {
                    binary_representation[counter] = i;
                    ++counter;
                }
            }
            auto first = binary_representation.begin();
            return uvec(first, first + counter);
        }
        case REVERSE:
        {
            for(uword i = 0; i < count; ++i)
            {
                binary_representation[count - 1 - i] = uword(bool(num & (one << i)));
            }
            break;
        }
        default:
        {
            cout << "Error occured in function toBinary..." << endl << flush;
        }
    }
    return binary_representation;
}

uvec Explorer::toBinary(uword num, size_t count, flagtype flag, size_t & size)
{
    uvec binary_representation(count);
    uword one = 1, result, counter = 0;
    switch(flag) {
        case CAST_TO_ONES: {
            for(uword i = 0; i < count; ++i) {
                result = uword(bool(num & (one << i)));
                binary_representation[i] = result;
                if(result) ++counter;
            }
            size = counter;
            break;
        }
        case RETURN_INDICES: {
            for(uword i = 0; i < count; ++i) {
                result = num & (one << i);
                if(result) {
                    binary_representation[counter] = i;
                    ++counter;
                }
            }
            size = counter;
            auto first = binary_representation.begin();
            return uvec(first, first + counter);
        }
        default: {
            cout << "Error occured in function toBinary..." << endl << flush;
        }
    }
    return binary_representation;
}


void Explorer::printResults(uvec FUNCTION, uvec ALL_ALPHAS)
{
    outfile << "Found universal function of " << n << " variables!\nF: [" << FUNCTION << "]\n";
    outfile << "D(f): [" << ALL_ALPHAS << "]\n";
    outfile << "q_start = " << q_start << ", current_q = " << q << endl << endl;
    printCoins();
    cout << "Found!" << endl << flush;
}

void Explorer::printCoins()
{
    outfile << "Coins:\n " << endl;
    for(uword i = 0; i < linear_function_count; ++i)
    {
        uvec bin = toBinary(i, n + 1, CAST_TO_ONES);
        bool flag = false;
        for(int z = n - 1; z >= 0; --z)
        {
            if(bin[z] == 1)
            {
                if(flag) outfile << " + ";
                flag = true;
                outfile << "x" << n - z;
            }
        }
        if(i == number_of_sets) outfile << "1";
        else if(bin[n] == 1) outfile << " + 1";
        if(i == 0) outfile << "0";
        outfile << ":\n";
        for(uword j = 0; j < n + 1; ++j)
        {
            outfile << "F(" << coins[i][j] << ")=F(" << toBinary(coins[i][j], n, REVERSE) << ")=" << L_VALUES[i][coins[i][j]]<< "\t";
        }
        outfile << endl << endl;
    }
}


bool Explorer::explore()
{
    bool flag = true;
    size_t q_size;
    time_t start_time = clock();
    while(q <= q_max) {
        all_alphas = toBinary(q, number_of_sets, RETURN_INDICES, q_size);
        if(q_size != k) {
            q += q_start;
            continue;
        }

        srand(time(0));
        uword random_f = rand() % ITER_SIZE;
        for(size_t i = random_f; i < F_VALUES_SIZE; i += ITER_SIZE) {
            F = F_VALUES[i];
            fill(number_of_coins.begin(), number_of_coins.end(), 0);
            flag = true;
            size_t set;
            size_t new_set_with_1;
            for(uword linfunc = 0; linfunc < linear_function_count; ++linfunc) {
                for(size_t f_index = 0; f_index < k; ++f_index) {
                    set = all_alphas[f_index];
                    new_set_with_1 = sets_with_1[set];
                    if(L_VALUES[linfunc][set] == F[f_index] && linCheck(linfunc, new_set_with_1)) {
                        ++number_of_coins[linfunc];
                    }
                }
                if(number_of_coins[linfunc] != (n + 1)) {
                    flag = false;
                    break;
                }
            }
            if(flag) {
                if(checkF(F, all_alphas)) return true;
                else return false;
            }
        }
        q += q_start;
    }
    outfile << "Nothing found with q = " << q_start << endl << flush;
    outfile << "Time taken: " << (clock() - start_time) / 1000000.0f << endl << flush;
    return false;
}

bool Explorer::checkF(uvec FUNCTION, uvec ALL_ALPHAS) {
    fill(number_of_coins.begin(), number_of_coins.end(), 0);
    bool flag = true;
    size_t set;
    size_t new_set;
    for(uword linfunc = 0; linfunc < linear_function_count; ++linfunc) {
        for(size_t f_index = 0; f_index < k; ++f_index) {
            set = ALL_ALPHAS[f_index];
            new_set = sets_with_1[set];
            if(L_VALUES[linfunc][set] == FUNCTION[f_index] && linCheck(linfunc, new_set)) {
                coins[linfunc][number_of_coins[linfunc]] = set;
                ++number_of_coins[linfunc];
            }
        }
        if(number_of_coins[linfunc] != (n + 1)) {
            flag = false;
            break;
        }
    }
    if(flag) {
        printResults(FUNCTION, ALL_ALPHAS);
        return true;
    }
    cout << "Invalid function!" << endl << flush;
    return false;
}

uword Explorer::howMany(uvec ALL_ALPHAS)
{
    uword fcounter = 0;
    for(uword i = 0; i < F_VALUES_SIZE; ++i)
    {
        uvec FUNCTION = F_VALUES[i];
        fill(number_of_coins.begin(), number_of_coins.end(), 0);
        bool flag = true;
        size_t set;
        size_t new_set;
        for(uword linfunc = 0; linfunc < linear_function_count; ++linfunc)
        {
            for(size_t f_index = 0; f_index < k; ++f_index)
            {
                set = ALL_ALPHAS[f_index];
                new_set = sets_with_1[set];
                if(L_VALUES[linfunc][set] == FUNCTION[f_index] && linCheck(linfunc, new_set))
                {
                    coins[linfunc][number_of_coins[linfunc]] = set;
                    ++number_of_coins[linfunc];
                }
            }
            if(number_of_coins[linfunc] != (n + 1))
            {
                flag = false;
                break;
            }
        }
        if(flag) ++fcounter;
    }
    return fcounter;
}

double Explorer::calcPk(uword q_start_test)
{
    q_start = q_start_test;
    uword s0 = (2 * number_of_function_values - 1) / q_start + 1;
    uword q = q_start * s0;
    uword s1 = (max_number - max_number / (2 * number_of_function_values)) / q_start;
    uword q_max = q_start * s1;

    bool flag = true;
    size_t q_size;
    uword sk = 0;

    while(q <= q_max)
    {
        q += q_start;
        all_alphas = toBinary(q, number_of_sets, CAST_TO_ONES, q_size);
        if (q_size != k) continue;
        ++sk;
    }
    outfile << sk / double(s1 - s0) << flush;
    return sk / double(s1 - s0);

}

bool Explorer::linCheck(uword & linfunc, size_t & new_set)
{
    switch(number_of_coins[linfunc])
    {
        case 0: {
            return modify0(new_set);
        }
        case 1: {
	        return modify1(new_set);
        }
        case 2: {
            return modify2(new_set);
        }
        case 3: {
            return modify3(new_set);
        }
        case 4: {
            return modify4(new_set);
        }
        case 5: {
	        return modify5(new_set);
        }
        case 6: {
            return modify6(new_set);
        }
        case n + 1: {
            return false;
        }
        default: {
            cout << "DEFAULT! ERROR..." << endl << flush;
            return false;
        }
    }
}

bool Explorer::modify0(size_t &new_set)
{
	lcomb[0] = new_set;
	return true;
}

bool Explorer::modify1(size_t &new_set)
{
	lcomb[1] = new_set ^ lcomb[0];
	lcomb[2] = new_set;
	return true;
}

bool Explorer::modify2(size_t &new_set)
{
	lcomb[3] = new_set ^ lcomb[0];
	lcomb[4] = new_set ^ lcomb[1];
	lcomb[5] = new_set ^ lcomb[2];
	lcomb[6] = new_set;
	return true;
}

bool Explorer::modify3(size_t &new_set)
{
	int sum = new_set ^ lcomb[4];
	if(sum)
	{
		lcomb[7] = new_set ^ lcomb[0];
		lcomb[8] = new_set ^ lcomb[1];
		lcomb[9] = new_set ^ lcomb[2];
		lcomb[10] = new_set ^ lcomb[3];
		lcomb[11] = sum;
		lcomb[12] = new_set ^ lcomb[5];
		lcomb[13] = new_set ^ lcomb[6];
		lcomb[14] = new_set;
		return true;
	}
	return false;
}

bool Explorer::modify4(size_t &new_set)
{
	int sum1 = new_set ^ lcomb[4];
	int sum2 = new_set ^ lcomb[8];
	if(sum1 && sum2)
	{
		int sum3 = new_set ^ lcomb[10];
		int sum4 = new_set ^ lcomb[12];
		if(sum3 && sum4)
		{
			lcomb[15] = new_set ^ lcomb[0];
			lcomb[16] = new_set ^ lcomb[1];
			lcomb[17] = new_set ^ lcomb[2];
			lcomb[18] = new_set ^ lcomb[3];
			lcomb[19] = sum1;
			lcomb[20] = new_set ^ lcomb[5];
			lcomb[21] = new_set ^ lcomb[6];
			lcomb[22] = new_set ^ lcomb[7];
			lcomb[23] = sum2;
			lcomb[24] = new_set ^ lcomb[9];
			lcomb[25] = sum3;
			lcomb[26] = new_set ^ lcomb[11];
			lcomb[27] = sum4;
			lcomb[28] = new_set ^ lcomb[13];
			lcomb[29] = new_set ^ lcomb[14];
			lcomb[30] = new_set;
			return true;
		}
	}
	return false;
}

bool Explorer::modify5(size_t &new_set)
{
	uvec sum(11);
	for(size_t i = 0; i < 11; ++i)
	{
		sum[i] = new_set ^ lcomb[indices5[i]];
		if(sum[i] == 0) return false;
	}
	lcomb[31] = new_set ^ lcomb[0];
	lcomb[32] = new_set ^ lcomb[1];
	lcomb[33] = new_set ^ lcomb[2];
	lcomb[34] = new_set ^ lcomb[3];
	lcomb[35] = sum[0];
	lcomb[36] = new_set ^ lcomb[5];
	lcomb[37] = new_set ^ lcomb[6];
	lcomb[38] = new_set ^ lcomb[7];
	lcomb[39] = sum[1];
	lcomb[40] = new_set ^ lcomb[9];
	lcomb[41] = sum[2];
	lcomb[42] = new_set ^ lcomb[11];
	lcomb[43] = sum[3];
	lcomb[44] = new_set ^ lcomb[13];
	lcomb[45] = new_set ^ lcomb[14];
	lcomb[46] = new_set ^ lcomb[15];
	lcomb[47] = sum[4];
	lcomb[48] = new_set ^ lcomb[17];
	lcomb[49] = sum[5];
	lcomb[50] = new_set ^ lcomb[19];
	lcomb[51] = sum[6];
	lcomb[52] = new_set ^ lcomb[21];
	lcomb[53] = sum[7];
	lcomb[54] = new_set ^ lcomb[23];
	lcomb[55] = sum[8];
	lcomb[56] = new_set ^ lcomb[25];
	lcomb[57] = sum[9];
	lcomb[58] = new_set ^ lcomb[27];
	lcomb[59] = sum[10];
	lcomb[60] = new_set ^ lcomb[29];
	lcomb[61] = new_set ^ lcomb[30];
	lcomb[62] = new_set;
	return true;
}

bool Explorer::modify6(size_t &new_set)
{
	uvec sum(26);
	for(size_t i = 0; i < 26; ++i)
	{
		sum[i] = new_set ^ lcomb[indices6[i]];
		if(sum[i] == 0) return false;
	}
	lcomb[63] = new_set ^ lcomb[0];
	lcomb[64] = new_set ^ lcomb[1];
	lcomb[65] = new_set ^ lcomb[2];
	lcomb[66] = new_set ^ lcomb[3];
	lcomb[67] = sum[0];
	lcomb[68] = new_set ^ lcomb[5];
	lcomb[69] = new_set ^ lcomb[6];
	lcomb[70] = new_set ^ lcomb[7];
	lcomb[71] = sum[1];
	lcomb[72] = new_set ^ lcomb[9];
	lcomb[73] = sum[2];
	lcomb[74] = new_set ^ lcomb[11];
	lcomb[75] = sum[3];
	lcomb[76] = new_set ^ lcomb[13];
	lcomb[77] = new_set ^ lcomb[14];
	lcomb[78] = new_set ^ lcomb[15];
	lcomb[79] = sum[4];
	lcomb[80] = new_set ^ lcomb[17];
	lcomb[81] = sum[5];
	lcomb[82] = new_set ^ lcomb[19];
	lcomb[83] = sum[6];
	lcomb[84] = new_set ^ lcomb[21];
	lcomb[85] = sum[7];
	lcomb[86] = new_set ^ lcomb[23];
	lcomb[87] = sum[8];
	lcomb[88] = new_set ^ lcomb[25];
	lcomb[89] = sum[9];
	lcomb[90] = new_set ^ lcomb[27];
	lcomb[91] = sum[10];
	lcomb[92] = new_set ^ lcomb[29];
	lcomb[93] = new_set ^ lcomb[30];
	lcomb[94] = new_set ^ lcomb[31];
	lcomb[95] = sum[11];
	lcomb[96] = new_set ^ lcomb[33];
	lcomb[97] = sum[12];
	lcomb[98] = new_set ^ lcomb[35];
	lcomb[99] = sum[13];
	lcomb[100] = new_set ^ lcomb[37];
	lcomb[101] = sum[14];
	lcomb[102] = new_set ^ lcomb[39];
	lcomb[103] = sum[15];
	lcomb[104] = new_set ^ lcomb[41];
	lcomb[105] = sum[16];
	lcomb[106] = new_set ^ lcomb[43];
	lcomb[107] = sum[17];
	lcomb[108] = new_set ^ lcomb[45];
	lcomb[109] = sum[18];
	lcomb[110] = new_set ^ lcomb[47];
	lcomb[111] = sum[19];
	lcomb[112] = new_set ^ lcomb[49];
	lcomb[113] = sum[20];
	lcomb[114] = new_set ^ lcomb[51];
	lcomb[115] = sum[21];
	lcomb[116] = new_set ^ lcomb[53];
	lcomb[117] = sum[22];
	lcomb[118] = new_set ^ lcomb[55];
	lcomb[119] = sum[23];
	lcomb[120] = new_set ^ lcomb[57];
	lcomb[121] = sum[24];
	lcomb[122] = new_set ^ lcomb[59];
	lcomb[123] = sum[25];
	lcomb[124] = new_set ^ lcomb[61];
	lcomb[125] = new_set ^ lcomb[62];
	lcomb[126] = new_set;
	return true;
}
