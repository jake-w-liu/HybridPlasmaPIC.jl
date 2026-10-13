/* Standalone MPFR shared-cache diagnostic. See docs/MPFR_UPSTREAM_REPORT.md.
 * A passing fresh process does not certify thread safety or correct rounding. */
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpfr.h>

/* Ascending requests grow the cache from 128 to 8192 bits across 64 calls. */
#define WORKERS 4
#define STEPS 16
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t condition = PTHREAD_COND_INITIALIZER;
static int ready = 0;
static int go = 0;

struct result {
    int id;
    int failures;
};

static void require_pthread(int error, const char *operation)
{
    if (error != 0) {
        fprintf(stderr, "%s: %s\n", operation, strerror(error));
        /* A failed synchronization operation makes the whole probe unusable. */
        exit(EXIT_FAILURE);
    }
}

static void *worker(void *input)
{
    struct result *out = input;
    require_pthread(pthread_mutex_lock(&mutex), "worker mutex lock");
    ++ready;
    require_pthread(pthread_cond_broadcast(&condition), "worker readiness broadcast");
    while (!go)
        require_pthread(pthread_cond_wait(&condition, &mutex), "worker start wait");
    require_pthread(pthread_mutex_unlock(&mutex), "worker mutex unlock");

    /* libm's double-precision reference does not use MPFR's constant cache. */
    const double expected = acos(-1.0);
    for (int step = 0; step < STEPS; ++step) {
        mpfr_prec_t precision = 128 * (WORKERS * step + out->id + 1);
        mpfr_t value;
        mpfr_init2(value, precision);
        mpfr_const_pi(value, MPFR_RNDN);
        double observed = mpfr_get_d(value, MPFR_RNDN);
        if (!isfinite(observed) || fabs(observed - expected) > 1e-12) {
            fprintf(stderr, "BAD_PI precision=%ld observed=%.17g expected=%.17g\n",
                    (long)precision, observed, expected);
            ++out->failures;
        }
        mpfr_clear(value);
    }
    mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
    return NULL;
}

int main(void)
{
    printf("mpfr=%s header=%s gmp=%s tls=%d shared=%d patches=%s\n",
           mpfr_get_version(), MPFR_VERSION_STRING, gmp_version,
           mpfr_buildopt_tls_p(), mpfr_buildopt_sharedcache_p(), mpfr_get_patches());
    fflush(stdout);
    if (!mpfr_buildopt_tls_p()) {
        fprintf(stderr, "This probe requires MPFR thread-local state.\n");
        return EXIT_FAILURE;
    }
    pthread_t jobs[WORKERS];
    struct result results[WORKERS];
    for (int i = 0; i < WORKERS; ++i) {
        results[i] = (struct result){i, 0};
        require_pthread(pthread_create(&jobs[i], NULL, worker, &results[i]), "thread creation");
    }
    require_pthread(pthread_mutex_lock(&mutex), "coordinator mutex lock");
    while (ready < WORKERS)
        require_pthread(pthread_cond_wait(&condition, &mutex), "coordinator readiness wait");
    go = 1;
    require_pthread(pthread_cond_broadcast(&condition), "coordinator start broadcast");
    require_pthread(pthread_mutex_unlock(&mutex), "coordinator mutex unlock");
    int failures = 0;
    for (int i = 0; i < WORKERS; ++i) {
        require_pthread(pthread_join(jobs[i], NULL), "thread join");
        failures += results[i].failures;
    }
    mpfr_free_cache();
    require_pthread(pthread_cond_destroy(&condition), "condition cleanup");
    require_pthread(pthread_mutex_destroy(&mutex), "mutex cleanup");
    printf("constants=%d failures=%d\n", WORKERS * STEPS, failures);
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
