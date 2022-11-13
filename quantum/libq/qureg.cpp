//
// Created by jewoo on 2022-11-13.
//
#include <iostream>

#include "libq.h"

namespace libq {

    float probability(cmplx ampl) {
        return ampl.real() * ampl.real() + ampl.imag() * ampl.imag();
    }

    qureg *new_qureg(state_t init_val, int width) {
        qureg *reg = new qureg;

        reg->width = width;
        reg->size = 1;
        reg->max_size = 0;
        reg->hash_computes = 0;
        reg->hash_w = width + 2;


        /* allocate memory for 1 basis state */
        reg->state = static_cast<state_t * >(calloc(1, sizeof(state_t)));
        reg->amplitude = static_cast<cmplx *>(calloc(1, sizeof(cmplx)));

        /* Allocate the hash table */
        reg->hash = static_cast<int *>(calloc(1 << reg->hash_w, sizeof(int)));

        reg->hash_cashing = true;
        reg->hash_hits = nullptr;
        reg->hits = 0;

        if (reg->hash_cashing) {
            reg->hash_hits = static_cast<int *>(calloc(HASH_CASH_SIZE, sizeof(int)));
        }

        /* Initialize the quantum register */
        reg->state[0] = init_val;
        reg->amplitude[0] = 1;

        return reg;
    }

    void delete_qureg(qureg *reg) {
        if (reg->hash) {
            free(reg->hash);
            reg->hash = nullptr;
        }
        if (reg->hash_hits) {
            free(reg->hash_hits);
            reg->hash_hits = nullptr;
            reg->hits = 0;
        }
        free(reg->amplitude);
        reg->amplitude = nullptr;

        if (reg->state) {
            free(reg->state);
            reg->state = nullptr;
        }

        delete reg;
    }

    void print_qureg(qureg *reg) {
        printf("States with non-zero probability:\n");
        for (int i = 0; i < reg->size; ++i) {
            printf("  % f %+fi|%llu> (%e) (|", reg->amplitude[i].real(),
                   reg->amplitude[i].imag(), reg->state[i],
                   probability(reg->amplitude[i]));
            for (int j = reg->width - 1; j >= 0; --j) {
                if (j % 4 == 3) {
                    printf(" ");
                }
                printf("%i", (((static_cast<state_t>(1) << j) & reg->state[i]) > 0));
            }
            printf(">)\n");
        }
    }

    void print_qureg_stats(qureg *reg) {
        printf("# of qubits            : %d\n", reg->width);
        printf("# of hash computes     : %d\n", reg->hash_computes);
        printf("Maximum # if states: %d, theroretical: %d, %.3f%%\n", reg->max_size,
               2 << reg->width,
               100.0 * reg->max_size / (2 << reg->width));
    }
}