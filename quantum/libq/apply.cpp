//
// Created by jewoo on 2022-11-06.
//

#include <cstdio>
#include <cstring>
#include <cmath>

#include "libq.h"

namespace libq {

    static inline unsigned int hash64(state_t key, int width) {
        unsigned int k32 = (key & 0xFFFFFFFF) ^ (key >> 32);
        k32 *= 0x9e370001UL;
        k32 = k32 >> (32 - width);
        return k32;
    }

    state_t get_state(state_t a, qureg *reg) {
        unsigned int i = hash64(a, reg->hash_w);
        while (reg->hash[i]) {
            if (reg->state[reg->hash[i] - 1] == a) {
                return reg->hash[i] - 1;
            }
            i++;
            if (i == static_cast<unsigned int>((1 << reg->hash_w))) {
                break;
            }
        }
        return -1;
    }

    void libq_add_hash(state_t a, int pos, qureg *reg) {
        int mark = 0;

        int i = hash64(a, reg->hash_w);
        while (reg->hash[i]) {
            i++;
            if (i == (1 << reg->hash_w)) {
                if (!mark) {
                    i = 0;
                    mark = 1;
                } else {
                    // TODO(rhundt): Handle full hashtable.
                }
            }
        }
        reg->hash[i] = pos + 1;
        if (reg->hash_cashing && reg->hits < HASH_CASH_SIZE) {
            reg->hash_hits[reg->hits] = i;
            reg->hits += 1;
        }
    }

    void libq_reconstruct_hash(qureg *reg) {

    }
}