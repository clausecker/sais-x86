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

#ifndef TSIZE
# error sais-chr.c must be included from sais.c, do not compile alone!
#endif

extern sais_index_type
TSIZE(sais_main)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type fs, sais_index_type n, sais_index_type k);

/* find the start or end of each bucket */
static void
TSIZE(getCounts)(const CHR *restrict T, sais_index_type *restrict C, sais_index_type n, sais_index_type k)
{
	sais_index_type i;

	for (i = 0; i < k; ++i)
		C[i] = 0;

	for (i = 0; i < n; ++i)
		++C[T[i]];
}

/* sort all type LMS suffixes */
static void
TSIZE(LMSsort1)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type *C, sais_index_type *B,
    sais_index_type n, sais_index_type k)
{
	sais_index_type *b, i, j;
	sais_index_type c0, c1;

	/* compute SAl */
	if (C == B)
		TSIZE(getCounts)(T, C, n, k);

	getBucketStarts(C, B, k);
	j = n - 1;
	b = SA + B[c1 = T[j]];
	--j;
	*b++ = T[j] < c1 ? ~j : j;

	for (i = 0; i < n; ++i) {
		if (0 < (j = SA[i])) {
			assert(T[j] >= T[j + 1]);
      			if ((c0 = T[j]) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}

			assert(i < (b - SA));
			--j;
      			*b++ = T[j] < c1 ? ~j : j;
			SA[i] = 0;
		} else if (j < 0)
			SA[i] = ~j;
	}

	/* compute SAs */
	if(C == B)
		TSIZE(getCounts)(T, C, n, k);

	getBucketEnds(C, B, k);
	for (i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
		if(0 < (j = SA[i])) {
			assert(T[j] <= T[j + 1]);
			if ((c0 = T[j]) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}

			assert((b - SA) <= i);
			--j;
			*--b = (T[j] > c1) ? ~(j + 1) : j;
			SA[i] = 0;
		}
	}
}

static sais_index_type
TSIZE(LMSpostproc1)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type n, sais_index_type m)
{
	sais_index_type i, j, p, q, plen, qlen, name;
	sais_index_type c0, c1;

	/* compact all the sorted substrings into the first m items of SA
	   2*m must be not larger than n (proveable) */
	assert(0 < n);
	for (i = 0; (p = SA[i]) < 0; ++i) {
		SA[i] = ~p;
		assert(i + 1 < n);
	}

	if(i < m)
		for (j = i++; ; ++i) {
			assert(i < n);
			if ((p = SA[i]) < 0) {
				SA[j++] = ~p;
				SA[i] = 0;
				if (j == m)
					break;
			}
		}

	/* store the length of all substrings */
	i = n - 1;
	j = n - 1;
	c0 = T[n - 1];


	for (;;) {
		do c1 = c0;
		while (0 <= --i && (c0 = T[i]) >= c1);

		if (0 > i)
			break;

		do c1 = c0;
		while(0 <= --i && (c0 = T[i]) <= c1);

		if (0 > i)
			break;

		SA[m + (i + 1 >> 1)] = j - i;
		j = i + 1;
	}

	/* find the lexicographic names of all substrings */
	for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
 		p = SA[i];
		plen = SA[m + (p >> 1)];

		if (plen == qlen && q + plen < n) {
			for (j = 0; j < plen && T[p + j] == T[q + j]; ++j)
				;

			if (j == plen)
				goto p_equals_q;
		}

		++name;
		q = p;
		qlen = plen;

	p_equals_q:
		SA[m + (p >> 1)] = name;
	}

	return (name);
}

static void
TSIZE(LMSsort2)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type *restrict C, sais_index_type *restrict B, sais_index_type *restrict D,
    sais_index_type n, sais_index_type k)
{
	sais_index_type *b, i, j, t, d;
	sais_index_type c0, c1;

	assert(C != B);

	/* compute SAl */
	getBucketStarts(C, B, k);
	j = n - 1;
	b = SA + B[c1 = T[j]];
	--j;
	t = T[j] < c1;
	j += n;
	*b++ = t & 1 ? ~j : j;
	for (i = 0, d = 0; i < n; ++i) {
		if (0 < (j = SA[i])) {
			if (n <= j) {
				d += 1;
				j -= n;
			}

			assert(T[j] >= T[j + 1]);
			if ((c0 = T[j]) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}

			assert(i < b - SA);
			--j;
			t = c0 << 1 | T[j] < c1;

			if (D[t] != d) {
				j += n;
				D[t] = d;
			}

			*b++ = t & 1 ? ~j : j;
			SA[i] = 0;
		} else if (j < 0)
			SA[i] = ~j;
	}

	for (i = n - 1; 0 <= i; --i)
		if (0 < SA[i] && SA[i] < n) {
			SA[i] += n;
			for (j = i - 1; SA[j] < n; --j)
				;

			SA[j] -= n;
			i = j;
		}

	/* compute SAs */
	getBucketEnds(C, B, k);

	for (i = n - 1, d += 1, b = SA + B[c1 = 0]; 0 <= i; --i)
		if (0 < (j = SA[i])) {
			if (n <= j) {
				d += 1;
				j -= n;
			}

			assert(T[j] <= T[j + 1]);
			if ((c0 = T[j]) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}

			assert((b - SA) <= i);
			--j;
			t = c0 << 1 | T[j] > c1;
			if (D[t] != d) {
				j += n;
				D[t] = d;
			}

			*--b = t & 1 ? ~(j + 1) : j;
			SA[i] = 0;
		}
}

/* compute SA */
static void
TSIZE(induceSA)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type *C, sais_index_type *B,
    sais_index_type n, sais_index_type k)
{
	sais_index_type *b, i, j;
	sais_index_type c0, c1;

	/* compute SAl */
	if (C == B)
		TSIZE(getCounts)(T, C, n, k);

	getBucketStarts(C, B, k);
	j = n - 1;
	b = SA + B[c1 = T[j]];
	*b++ = 0 < j && T[j - 1] < c1 ? ~j : j;

	for (i = 0; i < n; ++i) {
		j = SA[i];
		SA[i] = ~j;

		if (0 < j) {
			--j;
			assert(T[j] >= T[j + 1]);
			if ((c0 = T[j]) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}

			assert(i < b - SA);
			*b++ = 0 < j && T[j - 1] < c1 ? ~j : j;
		}
	}

	/* compute SAs */
	if (C == B)
		TSIZE(getCounts)(T, C, n, k);

	getBucketEnds(C, B, k);
	for (i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i)
		if (0 < (j = SA[i])) {
			--j;
			assert(T[j] <= T[j + 1]);
			if ((c0 = T[j]) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}

			assert((b - SA) <= i);
			*--b = j == 0 || T[j - 1] > c1 ? ~j : j;
		} else
			SA[i] = ~j;
}

/* stage 1: reduce the problem by at least 1/2, sort all the LMS-substrings */
static sais_index_type
TSIZE(sais_stage1)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type *C, sais_index_type *B, sais_index_type n, sais_index_type k,
    sais_index_type *restrict name, int flags)
{
	sais_index_type *D, *b;
	sais_index_type i, j, m, c0, c1, t;

	TSIZE(getCounts)(T, C, n, k);
	getBucketEnds(C, B, k);

	for (i = 0; i < n; ++i)
		SA[i] = 0;

	b = &t;
	i = n - 1;
	j = n;
	m = 0;
	c0 = T[n - 1];


	for (;;) {
		do c1 = c0;
		while (0 <= --i && (c0 = T[i]) >= c1);

		if (0 > i)
			break;

		do c1 = c0;
		while (0 <= --i && (c0 = T[i]) <= c1);

		if (0 > i)
			break;

		*b = j;
		b = SA + --B[c1];
		j = i;
		++m;
	}

	if (1 < m) {
		if (flags & 48) {
			if (flags & 16) {
				if (D = malloc(k * 2 * sizeof *D), D == NULL)
					return (-2);
			} else
				D = B - k * 2;

			assert(j + 1 < n);
			++B[T[j + 1]];

			for(i = 0, j = 0; i < k; ++i) {
				j += C[i];
				if (B[i] != j) {
					assert(SA[B[i]] != 0);
					SA[B[i]] += n;
				}

				D[i] = D[i + k] = 0;
			}

			TSIZE(LMSsort2)(T, SA, C, B, D, n, k);
			*name = LMSpostproc2(SA, n, m);

			if (flags & 16)
				free(D);
		} else {
			TSIZE(LMSsort1)(T, SA, C, B, n, k);
			*name = TSIZE(LMSpostproc1)(T, SA, n, m);
		}
	} else if (m == 1) {
		*b = j + 1;
		*name = 1;
	} else
		*name = 0;

	return (m);
}

static int
TSIZE(sais_stage2)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type n, sais_index_type k, sais_index_type m,
    sais_index_type name, sais_index_type fs, int flags)
{
	sais_index_type *RA;
	sais_index_type i, j, c0, c1, newfs;

	newfs = (n + fs) - (m * 2);
	if ((flags & 13) == 0) {
		if (k + name <= newfs)
			newfs -= k;
		else
			flags |= 8;
	}

	assert(n >> 1 <= newfs + m);
	RA = SA + m + newfs;

	for (i = m + (n >> 1) - 1, j = m - 1; m <= i; --i)
		if (SA[i] != 0)
			RA[j--] = SA[i] - 1;

	if (sais_main_idx(RA, SA, newfs, m, name) != 0)
		return (-2);

	i = n - 1;
	j = m - 1;
	c0 = T[n - 1];


	for (;;) {
		do c1 = c0;
		while (0 <= --i && (c0 = T[i]) >= c1);

		if (0 > i)
			break;

		do c1 = c0;
		while (0 <= --i && (c0 = T[i]) <= c1);

		if (0 > i)
			break;

		RA[j--] = i + 1;
	}

	for (i = 0; i < m; ++i)
		SA[i] = RA[SA[i]];

	return (flags);
}

/* put all left-most S characters into their buckets */
static void
TSIZE(sais_stage3)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type *C, sais_index_type *B, sais_index_type n, sais_index_type k, sais_index_type m)
{
	sais_index_type i, j, p, q, c0, c1;

	getBucketEnds(C, B, k);

	i = m - 1;
	j = n;
	p = SA[m - 1];
	c1 = T[p];

	do {
		q = B[c0 = c1];

		while(q < j)
			SA[--j] = 0;

		do {
			SA[--j] = p;
			if (--i < 0)
				break;

			p = SA[i];
		} while ((c1 = T[p]) == c0);
	} while (0 <= i);

	while (0 < j)
		SA[--j] = 0;
}

/* find the suffix array SA of T[0..n-1] in {0..255}^n */
extern sais_index_type
TSIZE(sais_main)(const CHR *restrict T, sais_index_type *restrict SA,
    sais_index_type fs, sais_index_type n, sais_index_type k)
{
	sais_index_type *C, *B;
	sais_index_type m, name;
	int flags, newflags;
	/*
	 * flags:
	 *  1: C is malloc'ed
	 *  2: B is malloc'ed
	 *  4: C is malloc'ed and alias to B
	 *  8: B == C???
	 * 16: D is malloc'ed
	 * 32: D = B - k * 2
	 * none of 16, 32: LMSsort1 is used
	 */

	assert(T != NULL && SA != NULL);
	assert(0 <= fs && 0 < n && 1 <= k);

	/* step 0: decide on a memory allocation scheme, set flags */
	if (k <= MINBUCKETSIZE) {
		if (C = malloc(k * sizeof *C), C == NULL)
			return (-2);

		if (k <= fs) {
			B = SA + (n + fs - k);
			flags = 1;
		} else {
			if (B = malloc(k * sizeof *B), B == NULL) {
				free(C);
				return (-2);
			}

			flags = 3;
		}
	} else if (k <= fs) {
		C = SA + (n + fs - k);
		if (k <= fs - k) {
			B = C - k;
			flags = 0;
		} else if (k <= MINBUCKETSIZE * 4) {
			if (B = malloc(k * sizeof *B), B == NULL)
				return (-2);

			flags = 2;
		} else {
			B = C;
			flags = 8;
		}
	} else {
		if (C = B = malloc(k * sizeof *B), B == NULL)
			return (-2);

		flags = 4 | 8;
	}

	if (n <= SAIS_LMSSORT2_LIMIT && 2 <= n / k) {
		if (flags & 1)
			flags |= k * 2 <= fs - k ? 32 : 16;
		else if (flags == 0 && k * 2 <= fs - k * 2)
			flags |= 32;
	}

	m = TSIZE(sais_stage1)(T, SA, C, B, n, k, &name, flags);
	if (m < 0) {
		if (flags & 5)
			free(C);

		if (flags & 2)
			free(B);

		return (-2);
	}

	/* stage 2: solve the reduced problem, recurse if names are not yet unique */
	if (name < m) {
		if (flags & 4)
			free(C);

		if (flags & 2)
			free(B);

		newflags = TSIZE(sais_stage2)(T, SA, n, k, m, name, fs, flags);
		if (newflags < 0) {
			if (flags & 1)
				free(C);

			return (-2);
		} else
			flags = newflags;

		if (flags & 4)
			if (C = B = malloc(k * sizeof *B), B == NULL)
				return -2;

		if (flags & 2)
			if (B = malloc(k * sizeof *B), B == NULL) {
				if (flags & 1)
					free(C);

				return -2;
			}
	}

	/* stage 3: induce the result for the original problem */
	if (flags & 8)
		TSIZE(getCounts)(T, C, n, k);

	if (1 < m)
		TSIZE(sais_stage3)(T, SA, C, B, n, k, m);

	TSIZE(induceSA)(T, SA, C, B, n, k);

	if (flags & 5)
		free(C);

	if (flags & 2)
		free(B);

	return (0);
}
