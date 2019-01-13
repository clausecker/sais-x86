/*
 * sais.c for sais-lite
 * Copyright (c) 2008--2010 Yuta Mori All Rights Reserved.
 * Copyright (c) 2018--2019 Robert Clausecker All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "sais.h"

#ifndef UCHAR_SIZE
# define UCHAR_SIZE 256
#endif
#ifndef MINBUCKETSIZE
# define MINBUCKETSIZE 256
#endif

#define sais_index_type int
#define SAIS_LMSSORT2_LIMIT 0x3fffffff

static void
getBucketStarts(const sais_index_type *C, sais_index_type *B, sais_index_type k)
{
	sais_index_type i, sum = 0;

	for (i = 0; i < k; ++i) {
		sum += C[i];
		B[i] = sum - C[i];
	}
}

static void
getBucketEnds(const sais_index_type *C, sais_index_type *B, sais_index_type k)
{
	sais_index_type i, sum = 0;

	for (i = 0; i < k; ++i) {
		sum += C[i];
		B[i] = sum;
	}
}

static sais_index_type
LMSpostproc2(sais_index_type *restrict SA, sais_index_type n, sais_index_type m)
{
	sais_index_type i, j, d, name;

	/* compact all the sorted LMS substrings into the first m items of SA */
	assert(0 < n);
	for (i = 0, name = 0; (j = SA[i]) < 0; ++i) {
		j = ~j;

		if (n <= j)
			name += 1;

		SA[i] = j;
		assert(i + 1 < n);
	}

	if (i < m)
		for (d = i++; ; ++i) {
			assert(i < n);
			if ((j = SA[i]) < 0) {
				j = ~j;
				if (n <= j)
					name += 1;

				SA[d++] = j;
				SA[i] = 0;

				if (d == m)
					break;
			}
		}

	if(name < m) {
		/* store the lexicographic names */
		for (i = m - 1, d = name + 1; 0 <= i; --i) {
      			if(n <= (j = SA[i])) {
				j -= n;
				--d;
			}

			SA[m + (j >> 1)] = d;
		}
	} else {
		/* unset flags */
		for (i = 0; i < m; ++i) {
			if (n <= (j = SA[i])) {
				j -= n;
				SA[i] = j;
			}
		}
	}

	return (name);
}

/* instantiate sais-chr.c for char and sais_index_type */
#define TSIZE(ident) ident ## _idx
#define CHR sais_index_type
#include "sais-chr.orig.c"
#undef CHR
#undef TSIZE

#define TSIZE(ident) ident ## _chr
#define CHR unsigned char
#include "sais-chr.orig.c"
#undef CHR
#undef TSIZE

int
sais(const unsigned char *restrict T, int *restrict SA, int n) {
	if (T == NULL || SA == NULL || n < 0)
		return (-1);

	if (n <= 1) {
		if (n == 1)
			SA[0] = 0;

		return (0);
	}

	return (sais_main_chr(T, SA, 0, n, UCHAR_SIZE));
}
