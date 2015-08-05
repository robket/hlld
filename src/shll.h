#include <stdint.h>
#include "bitmap.h"

#ifndef SHLL_H
#define SHLL_H

// Ensure precision in a sane bound
#define HLL_MIN_PRECISION 4      // 16 registers
#define HLL_MAX_PRECISION 18     // 262,144 registers

typedef struct fpm {
	uint32_t timestamp;
	uint8_t r_value;
	struct fpm *next;
} fpm_t; // FuturePossibleMaximum

typedef struct fpm_list {
	fpm_t *first; // First is oldest and largest
} fpm_list_t;

typedef struct {
    unsigned char precision;
    uint32_t max_t;
    fpm_list_t *registers;
    shlld_bitmap *bm;
} shll_t;

/**
 * Initializes a new SHLL
 * @arg precision The digits of precision to use
 * @arg h The HLL to initialize
 * @return 0 on success
 */
int shll_init(unsigned char precision, shll_t *h);

/**
 * Initializes a new HLL from a bitmap
 * @arg precision The digits of precision to use
 * @arg bm The bitmap to use
 * @arg h The HLL to initialize
 * @return 0 on success
 */
int shll_init_from_bitmap(unsigned char precision, shlld_bitmap *bm, shll_t *h);

/**
 * Destroys an hll. Closes the bitmap, but does not free it.
 * @return 0 on success
 */
int shll_destroy(shll_t *h);

/**
 * Adds a new key to the HLL
 * @arg h The hll to add to
 * @arg key The key to add
 */
void shll_add(shll_t *h, char *key, uint32_t timestamp);

/**
 * Adds a new hash to the HLL
 * @arg h The hll to add to
 * @arg hash The hash to add
 */
void shll_add_hash(shll_t *h, uint64_t hash, uint32_t timestamp);

/**
 * Estimates the cardinality of the HLL
 * @arg h The hll to query
 * @return An estimate of the cardinality
 */
double shll_size(shll_t *h, uint32_t start_time);

/**
 * Computes the minimum digits of precision
 * needed to hit a target error.
 * @arg error The target error rate
 * @return The number of digits needed, or
 * negative on error.
 */
int shll_precision_for_error(double err);

/**
 * Computes the upper bound on variance given
 * a precision
 * @arg prec The precision to use
 * @return The expected variance in the count,
 * or zero on error.
 */
double shll_error_for_precision(int prec);

/**
 * Computes the bytes required for a HLL of the
 * given precision.
 * @arg prec The precision to use
 * @return The bytes required or 0 on error.
 */
uint64_t shll_bytes_for_precision(int prec);

#endif
