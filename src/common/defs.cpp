#include "defs.h"

const u1_t* get_dna_encode_table()
{
    static bool table_is_filled = false;
    static u1_t table[256];
    if (table_is_filled) return table;
    std::fill(table, table + 256, 16);
	
#define __encode_one_nucl(L, U, v) \
		do { \
        int l = (L); \
        int u = (U); \
        table[l] = table[u] = (v); \
    } while (0)

    __encode_one_nucl('-', '-', 15);
    __encode_one_nucl('a', 'A', 0);
    __encode_one_nucl('c', 'C', 1);
    __encode_one_nucl('m', 'M', 6);
    __encode_one_nucl('g', 'G', 2);
    __encode_one_nucl('r', 'R', 4);
    __encode_one_nucl('s', 'S', 9);
    __encode_one_nucl('v', 'V', 13);
    __encode_one_nucl('t', 'T', 3);
    __encode_one_nucl('w', 'W', 8);
    __encode_one_nucl('y', 'Y', 5);
    __encode_one_nucl('h', 'H', 12);
    __encode_one_nucl('k', 'K', 7);
    __encode_one_nucl('d', 'D', 11);
    __encode_one_nucl('b', 'B', 10);
    __encode_one_nucl('n', 'N', 14);
	
	table_is_filled = true;
	return table;
}

const char*
get_dna_decode_table()
{
    static const char dt[] = { 
                              'A', 'C', 'G', 'T','-' 
                             };
    return dt;
}

const u1_t* 
get_dna_complement_table()
{
	static const u1_t ct[] = { 3, 2, 1, 0 };
	return ct;
}
