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
        reg->hash_computes += 1;
        if (reg->hash_cashing && reg->hits < HASH_CASH_SIZE) {
            for (int i = 0; i < reg->hits; ++i) {
                reg->hash[reg->hash_hits[i]] = 0;
                reg->hash_hits[i] = 0;
            }
            reg->hits = 0;
        } else {
            memset(reg->hash, 0, (1 << reg->hash_w) * sizeof(int));
            memset(reg->hash_hits, 0, reg->hits * sizeof(int));
            reg->hits = 0;
            // TODO(rhundt): Apparently the compiler doesn't convert this loop,
            // Investigate.
//            for (int i = 0; i < (i << reg->hash_w); ++i) {
//                reg->hash[i] = 0;
//            }
        }
        for (int i = 0; i < reg->size; ++i) {
            libq_add_hash(reg->state[i], i, reg);
        }
    }

    void libq_gate1(int target, cmplx m[4], qureg *reg) {
        int add_size{0};

        libq_reconstruct_hash(reg);

        // calculae the number of basis states to be added.
        for (int i = 0; i < reg->size; ++i) {
            /* determine whether XORed basis state already exists. */
            if (get_state(reg->state[i] ^ (static_cast<state_t>(1) << target), reg) == static_cast<state_t>(-1)) {
                add_size++;
            }

            /* allocate memory for the new basis states */
            if (add_size) {
                reg->state = static_cast<state_t *>(
                        realloc(reg->state, (reg->size + add_size) * sizeof(state_t)));
                reg->amplitude = static_cast<cmplx *>(
                        realloc(reg->amplitude, (reg->size + add_size) * sizeof(cmplx)));

                memset(&reg->state[reg->size], 0, add_size * sizeof(int));
                memset(&reg->amplitude[reg->size], 0, add_size * sizeof(cmplx));

                if (reg->size + add_size > reg->max_size) { reg->max_size = reg->size + add_size; }
            }
        }
        char *done = static_cast<char *>(calloc(reg->size + add_size,
                                                sizeof(char)));
        int next_state = reg->size;
        float limit = (1.0 / (static_cast<state_t>(1) << reg->width)) * 1e-6;

        /* perform the actual matix multiplication */
        for (int i = 0; i < reg->size; ++i) {
            if (!done[i]) {
                /* determine if the target of the basis state is set */
                int is_set = reg->state[i] & (static_cast<state_t>(1) << target);
                int xor_index =
                        get_state(reg->state[i] ^ (static_cast<state_t>(1) << target), reg);
                cmplx tnot = xor_index >= 0 ? reg->amplitude[xor_index] : 0;
                cmplx t = reg->amplitude[i];

                if (is_set) {
                    reg->amplitude[i] = m[2] * tnot + m[3] * t;
                } else {
                    reg->amplitude[i] = m[0] * t + m[1] * tnot;
                }

                if (xor_index >= 0) {
                    if (is_set) {
                        reg->amplitude[xor_index] = m[0] * tnot + m[1] * t;
                    } else {
                        reg->amplitude[xor_index] = m[2] * t + m[3] * tnot;
                    }
                } else { /* new basis state will be created */
                    if (std::abs(m[1]) == 0.0 && is_set) { break; }
                    if (std::abs(m[2]) == 0.0 && !is_set) { break; }

                    reg->state[next_state] =
                            reg->state[i] ^ (static_cast<state_t>(1) << target);
                    reg->amplitude[next_state] = is_set ? m[1] * t : m[2] * t;
                    next_state += 1;
                }
                if (xor_index >= 0) {
                    done[xor_index] = 1;
                }
            }
        }

        reg->size += add_size;
        free(done);

        /* remove basis states with extremply small amplitude. */
        if (reg->hash_w) {
            int dec_size{0};
            for (int i = 0, j = 0; i < reg->size; ++i) {
                if (probability(reg->amplitude[i]) < limit) {
                    j++;
                    dec_size++;
                } else if (j) {
                    reg->state[i - j] = reg->state[i];
                    reg->amplitude[i - j] = reg->amplitude[i];
                }
            }

            if (dec_size) {
                reg->size -= dec_size;

                // These seem redundant.
//                reg->amplitude = static_cast<cmplx *>(
//                        realloc(reg->amplitude, reg->size * sizeof(cmplx)));
//                reg->state = static_cast<state_t *>(
//                        realloc(reg->state, reg->size * sizeof(state_t)));
            }
        }

        if (reg->size > (1 << (reg->hash_w - 1))) {
            fprintf(stderr,
                    "Warning: inefficient hash table (size %i vs hash %i)\n",
                    reg->size, 1 << reg->hash_w);
        }
    }
}