/*
 * Based on the Google paper
 * "HyperLogLog in Practice: Algorithmic Engineering of a
State of The Art Cardinality Estimation Algorithm"
 *
 * We implement a HyperLogLog using 6 bits for register,
 * and a 64bit hash function. For our needs, we always use
 * a dense representation and avoid the sparse/dense conversions.
 *
 */
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "shll.h"
#include "hll_constants.h"

#define REG_WIDTH 6     // Bits per register
#define INT_WIDTH 32    // Bits in an int

#define NUM_REG(precision) ((1 << precision))
#define INT_CEIL(num, denom) (((num) + (denom) - 1) / (denom))

// Link the external murmur hash in
extern void MurmurHash3_x64_128(const void * key, const int len, const uint32_t seed, void *out);


/**
 * Initializes a new HLL
 * @arg precision The digits of precision to use
 * @arg h The HLL to initialize
 * @return 0 on success
 */
int shll_init(unsigned char precision, shll_t *h) {
    // Ensure the precision is somewhat sane
    if (precision < HLL_MIN_PRECISION || precision > HLL_MAX_PRECISION)
        return -1;

    // Store precision
    h->precision = precision;

    // Determine how many registers are needed
    int reg = NUM_REG(precision);

    // Allocate and zero out the registers
    h->bm = NULL;
    h->registers = calloc(reg, sizeof(fpm_list_t));
    if (!h->registers) return -1;
    return 0;
}

static uint32_t to_uint32(uint8_t* address) {
    uint32_t result = 0;
    result += *address << 24;
    result += *(address + 1) << 16;
    result += *(address + 2) << 8;
    result += *(address + 3);
    return result;
}

static void to_uint8(uint32_t value, uint8_t* address) {
    *address = (uint8_t) value >> 24;
    *(address + 1) = (uint8_t) value >> 16;
    *(address + 2) = (uint8_t) value >> 8;
    *(address + 3) = (uint8_t) value;
}

static fpm_t * insert_new_fpm(fpm_list_t *fpm_list, fpm_t *last_fpm, uint32_t timestamp, uint8_t r_value) {
    // Insert new value after last_fpm
    fpm_t *new_fpm = calloc(1, sizeof(fpm_t));
    if (last_fpm != NULL) {
        last_fpm->next = new_fpm;
    } else {
        fpm_list->first = new_fpm;
    }
    new_fpm->r_value = r_value;
    new_fpm->timestamp = timestamp;
    new_fpm->next = NULL;
    fpm_list->size++;
    return new_fpm;
}

/**
 * Initializes a new HLL from a bitmap
 * @arg precision The digits of precision to use
 * @arg bm The bitmap to use
 * @arg h The HLL to initialize
 * @return 0 on success
 */
int shll_init_from_bitmap(unsigned char precision, shlld_bitmap *bm, shll_t *h) {
    // Ensure the precision is somewhat sane
    if (precision < HLL_MIN_PRECISION || precision > HLL_MAX_PRECISION)
        return -1;

    // Store precision
    h->precision = precision;

    // Determine how many registers are needed
    int reg = NUM_REG(precision);

    // Allocate and zero out the registers
    h->registers = calloc(reg, sizeof(fpm_list_t));
    if (!h->registers) return -1;
    int mmap_index = 0;
    for (int i = 0; i < reg; i++) {
        uint8_t size = *(bm->mmap +mmap_index++);
        fpm_t *last_fpm = NULL;
        for (int j = 0; j < size; j++) {
            uint8_t r_value = *(bm->mmap + mmap_index++);
            uint32_t timestamp = to_uint32((bm->mmap + mmap_index));
            mmap_index += 4;
            last_fpm = insert_new_fpm((h->registers + i), last_fpm, timestamp, r_value);
        }
    }
    return 0;
}

void shll_to_bitmap(shll_t *h, shlld_bitmap *bitmap) {
    int reg = NUM_REG(h->precision);
    // Get total count of stored r-values
    int total_r_values = 0;
    for (int i = 0; i < reg; i++) {
        total_r_values += h->registers[i].size;
    }
    // Need reg * uint8 + total_r_values * (uint8 + uint32)
    // = (reg * + total_r_values * (1 + 4)) * uint8

    int size = (reg + total_r_values * 5) * sizeof(uint8_t);
    uint8_t* mmap = malloc(size);
    int mmap_index = 0;
    for (int i = 0; i < reg; i++) {
        uint8_t number_of_r_values = h->registers[i].size;
        *(mmap + mmap_index++) = number_of_r_values;
        fpm_t *node = h->registers[i].first;
        for (int j = 0; j < number_of_r_values; j++) {
            *(mmap + mmap_index++) = node->r_value;
            to_uint8(node->timestamp, (mmap + mmap_index));
            mmap_index += 4;
        }
    }
    bitmap->mmap = mmap;
    bitmap->size = size;
}


/**
 * Destroys an hll. Closes the bitmap, but does not free it.
 * @return 0 on success
 */
int shll_destroy(shll_t *h) {
    if (h->registers) {
        int reg = NUM_REG(h->precision);
        for (int i = 0; i < reg; i++) {
            fpm_list_t fpm_list = h->registers[i];
            fpm_t *node = fpm_list.first;
            fpm_list.first = NULL;
            while (node != NULL) {
                fpm_t *next_node = node->next;
                free(node);
                node = next_node;
                fpm_list.size--;
            }
        }
        free(h->registers);
        h->registers = NULL;
    }
    // TODO figure out what happens here when we get this from bitmap
    //if (h->bm) {
    //    bitmap_close(h->bm);
    //    h->bm = NULL;
    //}
    return 0;
}

static int get_fpm_register(shll_t *h, int idx, uint32_t start_time) {
    fpm_list_t fpm_list = (h->registers)[idx];
    fpm_t *node = fpm_list.first;

    while (node != NULL && node->timestamp < start_time) {
        node = node->next;
    }
    if (node != NULL) {
        return node->r_value;
    } else {
        return 0;
    }
}

static void set_fpm_register(shll_t *h, int idx, uint32_t timestamp, uint8_t r_value) {
    // Get pointer to list
    fpm_list_t fpm_list = h->registers[idx];
    fpm_t *node = fpm_list.first;

    // Delete nodes older than timestamp - max_t from beginning
    while (node != NULL && node->timestamp < timestamp - h->max_t) {
        fpm_t *next_node = node->next;
        free(node);
        node = next_node;
        fpm_list.size--;
    }
    fpm_list.first = node;
    fpm_t *previous_node = node;

    // Find first node with value smaller than/equal to new value
    while (node != NULL && node->r_value > r_value) {
        previous_node = node;
        node = node->next;
    }

    insert_new_fpm(&fpm_list, previous_node, timestamp, r_value);
    // Delete old element and subsequent elements
    while (node != NULL) {
        fpm_t *next_node = node->next;
        free(node);
        node = next_node;
        fpm_list.size--;
    }
}

/**
 * Adds a new key to the HLL
 * @arg h The hll to add to
 * @arg key The key to add
 */
void shll_add(shll_t *h, char *key, uint32_t timestamp) {
    // Compute the hash value of the key
    uint64_t out[2];
    MurmurHash3_x64_128(key, strlen(key), 0, &out);

    // Add the hashed value
    shll_add_hash(h, out[1], timestamp);
}

/**
 * Adds a new hash to the HLL
 * @arg h The hll to add to
 * @arg hash The hash to add
 */
void shll_add_hash(shll_t *h, uint64_t hash, uint32_t timestamp) {
    // Determine the index using the first p bits
    int idx = hash >> (64 - h->precision);

    // Shift out the index bits
    hash = hash << h->precision | (1 << (h->precision -1));

    // Determine the count of leading zeros
    int leading = __builtin_clzll(hash) + 1;

    // Update the register if the new value is larger
    set_fpm_register(h, idx, timestamp, leading);
}

/*
 * Returns the bias correctors from the
 * hyperloglog paper
 */
static double alpha(unsigned char precision) {
    switch (precision) {
        case 4:
            return 0.673;
        case 5:
            return 0.697;
        case 6:
            return 0.709;
        default:
            return 0.7213 / (1 + 1.079 / NUM_REG(precision));
    }
}

/*
 * Computes the raw cardinality estimate
 */
static double raw_estimate(shll_t *h, int *num_zero, uint32_t start_time) {
    unsigned char precision = h->precision;
    int num_reg = NUM_REG(precision);
    double multi = alpha(precision) * num_reg * num_reg;

    int reg_val;
    double inv_sum = 0;
    for (int i=0; i < num_reg; i++) {
        reg_val = get_fpm_register(h, i, start_time);
        inv_sum += pow(2.0, -1 * reg_val);
        if (!reg_val) *num_zero += 1;
    }
    return multi * (1.0 / inv_sum);
}

/*
 * Estimates cardinality using a linear counting.
 * Used when some registers still have a zero value.
 */
static double linear_count(shll_t *h, int num_zero) {
    int registers = NUM_REG(h->precision);
    return registers *
        log((double)registers / (double)num_zero);
}

/**
 * Binary searches for the nearest matching index
 * @return The matching index, or closest match
 */
static int binary_search(double val, int num, const double *array) {
    int low=0, mid, high=num-1;
    while (low < high) {
        mid = (low + high) / 2;
        if (val > array[mid]) {
            low = mid + 1;
        } else if (val == array[mid]) {
            return mid;
        } else {
            high = mid - 1;
        }
    }
    return low;
}

/**
 * Interpolates the bias estimate using the
 * empirical data collected by Google, from the
 * paper mentioned above.
 */
static double bias_estimate(shll_t *h, double raw_est) {
    // Determine the samples available
    int samples;
    int precision = h->precision;
    switch (precision) {
        case 4:
            samples = 80;
            break;
        case 5:
            samples = 160;
            break;
        default:
            samples = 200;
            break;
    }

    // Get the proper arrays based on precision
    double *estimates = *(rawEstimateData+(precision-4));
    double *biases = *(biasData+(precision-4));

    // Get the matching biases
    int idx = binary_search(raw_est, samples, estimates);
    if (idx == 0)
        return biases[0];
    else if (idx == samples)
        return biases[samples-1];
    else
        return (biases[idx] + biases[idx-1]) / 2;
}

/**
 * Estimates the cardinality of the HLL
 * @arg h The hll to query
 * @return An estimate of the cardinality
 */
double shll_size(shll_t *h, uint32_t start_time) {
    int num_zero = 0;
    double raw_est = raw_estimate(h, &num_zero, start_time);

    // Check if we need to apply bias correction
    int num_reg = NUM_REG(h->precision);
    if (raw_est <= 5 * num_reg) {
        raw_est -= bias_estimate(h, raw_est);
    }

    // Check if linear counting should be used
    double alt_est;
    if (num_zero) {
        alt_est = linear_count(h, num_zero);
    } else {
        alt_est = raw_est;
    }

    // Determine which estimate to use
    if (alt_est <= switchThreshold[h->precision-4]) {
        return alt_est;
    } else {
        return raw_est;
    }
}


/**
 * Computes the minimum number of registers
 * needed to hit a target error.
 * @arg error The target error rate
 * @return The number of registers needed, or
 * negative on error.
 */
int shll_precision_for_error(double err) {
    // Check that the error bound is sane
    if (err >= 1 || err <= 0)
        return -1;

    /*
     * Error of HLL is 1.04 / sqrt(m)
     * m is given by 2^p, so solve for p,
     * and use the ceiling.
     */
    double p = log2(pow(1.04 / err, 2));
    return ceil(p);
}

/**
 * Computes the upper bound on variance given
 * a precision
 * @arg prec The precision to use
 * @return The expected variance in the count,
 * or zero on error.
 */
double shll_error_for_precision(int prec) {
    // Check that the error bound is sane
    if (prec < HLL_MIN_PRECISION || prec > HLL_MAX_PRECISION)
        return 0;

    /*
     * Error of HLL is 1.04 / sqrt(m)
     * m is given by 2^p
     */
    int registers = pow(2, prec);
    return 1.04 / sqrt(registers);
}
