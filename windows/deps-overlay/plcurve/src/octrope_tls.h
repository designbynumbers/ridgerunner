/*
 * Thread-local storage helpers for octrope reentrancy.
 * Enables concurrent octrope calls (e.g. parallel line-search) and
 * OpenMP workers inside octrope_struts after the master broadcasts
 * by_oct / num_edges into each worker's TLS.
 */
#ifndef OCTROPE_TLS_H
#define OCTROPE_TLS_H

#if defined(_MSC_VER)
#  define OCTROPE_TLS __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
#  define OCTROPE_TLS __thread
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  define OCTROPE_TLS _Thread_local
#else
#  define OCTROPE_TLS
#endif

#endif /* OCTROPE_TLS_H */
