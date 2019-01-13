/*
 * SAIS-x86 0.1 derived from sais-lite 2.4.1 by Yuta Mori.
 *
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
# error sais-chr.S must b included from sais.S, do not assemble alone!
#endif

	// getCounts(T, C, n, k)
	// count how often each character of T occurs.
	// Characters are in the range 0 .. k - 1
	.type TSIZE(getCounts), @function
	.align 16
TSIZE(getCounts):
	xor %eax, %eax
	mov %rsi, %r8	// r8 = C
	xchg %rdi, %rsi	// rsi = T, rdi = C
	rep stosl	// memset(rdi, 0, 4 * rcx); now i = ecx = 0

	test $1, %dl	// n odd?
	jz 2f		// if yes, make i odd, too

	MOVCHR (%rsi), %eax	// eax = T[0]
	incl (%r8, %rax, 4)	// C[eax]++

	inc %ecx		// i++
2:	cmp %edx, %ecx		// i < n?
	jge 0f

	.align 16
1:	MOVCHR (%rsi, %rcx, CHR), %eax	// eax = T[i]
	incl (%r8, %rax, 4)		// C[eax]++

	MOVCHR CHR(%rsi, %rcx, CHR), %eax	// eax = T[i + 1]
	incl (%r8, %rax, 4)			// C[eax]++

	add $2, %ecx		// i++
	cmp %edx, %ecx		// i < n?
	jl 1b			// if i < n, continue

0:	ret
	.size TSIZE(getCounts), .-TSIZE(getCounts)

	// LMSsort2(T, SA, C, B, D, n, k)
	// sort LMS substrings according to strategy 2
	// (for when n <= SAIS_LMSSORT2_LIMIT)
	.type TSIZE(LMSsort2), @function
	.align 16
TSIZE(LMSsort2):
	// T: r12 (orig. in rdi)
	// SA: r14 (orig. in rsi)
	// C: -32(%rbp) (orig. in rdx)
	// B: rbx (orig. in rcx)
	// D: r8
	// n: r9d
	// k: 16(%rbp)
	// i: esi
	// j: rdx
	// b: edi
	// d: r11d
	// c0, tmp: eax
	// c1: ecx
	// t: r10

	push %rbp		// establish ...
	mov %rsp, %rbp		// stack frame

	push %rbx		// - 8 save non-volatile registers
	push %r12		// -16
	push %r14		// -24

	// save T, SA, C, B, D, and n
	// in preparation of function call
	push %rdx		// -32 save C
	mov %rdi, %r12		// save T
	mov %rsi, %r14		// save SA
	mov %rcx, %rbx		// save B
	push %r8		// -40 save D
	push %r9		// -48 save n

	// call getBucketStarts(C, B, k)
	mov %rdx, %rdi
	mov %rcx, %rsi
	mov 16(%rbp), %edx
	call getBucketStarts

	// restore registers
	pop %r9			// -40 n
	pop %r8			// -32 D

	// initialise local variables
	MOVCHR -CHR(%r12, %r9, CHR), %ecx	// c1 = T[n - 1]
	mov (%rbx, %rcx, 4), %edi	// B[c1]
	lea -2(, %r9, 2), %edx		// j = 2 * n - 2

#if CHR == 1
	cmp %cl, -2*CHR(%r12, %r9, CHR)	// T[n - 2] < c1?
#else /* CHR == 4 */
	cmp %ecx, -2*CHR(%r12, %r9, CHR)// T[n - 2] < c1?
#endif
	sbb %eax, %eax			// tmp = T[n - 2] < c1 ? -1 : 0
	xor %edx, %eax			// tmp = T[n - 2] < c1 ? ~j : j
	mov %eax, (%r14, %rdi, 4)	// SA[b] = tmp
	inc %edi			// b++

	xor %esi, %esi			// i = 0
	xor %r11d, %r11d		// d = 0
	cmp %r9d, %esi			// i < n?
	jge 0f				// if i >= n, break

	// for (i = 0, d = 0, i < n; ++i)
	.align 16
1:	mov (%r14, %rsi, 4), %edx	// j = SA[i]

	test %edx, %edx			// j > 0?
	jle 2f				// if j > 0 ...
	cmp %edx, %r9d			// 	n <= j
	jg 3f				// 	if n <= j
	inc %r11d			//		d++
	sub %r9d, %edx			//		j -= n

3:	MOVCHR (%r12, %rdx, CHR), %eax	// 	c0 = T[j]
	cmp %eax, %ecx			//	c0 == c1?
	je 3f				//	if c0 != c1

	mov %edi, (%rbx, %rcx, 4)	//		B[c1] = b
	mov (%rbx, %rax, 4), %edi	//		b = B[c0]
	mov %eax, %ecx			//		c1 = c0

3:	dec %edx			//	j--
#if CHR == 1
	cmp %cl, (%r12, %rdx, CHR)	//	T[j] < c1
#else
	cmp %ecx, (%r12, %rdx, CHR)	//	T[j] < c1
#endif
	mov %eax, %r10d			//	t = c0
	adc %r10d, %r10d		//	t = (c0 << 1) + (T[j] < c1)

	cmp %r11d, (%r8, %r10, 4)	//	D[t] != d?
	je 3f				//	if D[t] != d
	add %r9d, %edx			//		j += n
	mov %r11d, (%r8, %r10, 4)	//		D[t] = d

3:	shr %r10d			//	t & 1?
	sbb %eax, %eax			//	tmp = t & 1 ? -1 : 0
	xor %edx, %eax			//	tmp = t & 1 ? ~j : j
	mov %eax, (%r14, %rdi, 4)	//	SA[b] = tmp
	inc %edi			//	b++

	movl $0, (%r14, %rsi, 4)	//	SA[i] = 0
	jmp 3f

2:	je 3f				// else if j < 0
	not %edx			//	j = ~j
	mov %edx, (%r14, %rsi, 4)	//	SA[i] = j

3:	inc %esi			// i++
	cmp %r9d, %esi			// i < n?
	jl 1b				// if i < n, continue

0:	dec %esi			// i-- resp. (i = n - 1)
	js 0f				// if 0 > i, break

	// for (i = n - 1; 0 <= i; --i)
	.align 16
1:	mov (%r14, %rsi, 4), %eax	// tmp = SA[i]
	test %eax, %eax			// 0 < SA[i]?
	jle 2f
	cmp %r9d, %eax			// SA[i] < n?
	jge 2f				// if (0 < SA[i] && SA[i] < n)

	add %r9d, %eax			//	tmp += n
	mov %eax, (%r14, %rsi, 4)	//	SA[i] = tmp

	// for (j = i - 1; SA[j] < n; --j) {}
	.align 16
3:	dec %esi			//	i-- (j = i - 1)
	cmp %r9d, (%r14, %rsi, 4)	//	SA[j] < n?
	jl 3b

	sub %r9d, (%r14, %rsi, 4)	//	SA[j] -= n

2:	dec %esi			// i--
	jns 1b				// if 0 > i, break

0:	pop %rdi			// -24 pop C off the stack
	push %r11			// -32 save d
	push %r8			// -40 save D
	push %r9			// -48 save n

	// call getBucketEnds(C, B, k)
	mov %rbx, %rsi
	mov 16(%rbp), %edx
	call getBucketEnds

	pop %r9				// -40 restore n
	pop %r8				// -32 restore D
	pop %r11			// -24 restore d

	// initialise variables for final loop
	mov %r9d, %esi			// i = n
	inc %r11d			// d++
	xor %ecx, %ecx			// c1 = 0
	mov (%rbx), %edi		// b = B[0]

	dec %esi			// i--
	js 0f				// if 0 > i break

	// for (...; 0 <= i; --i)
	.align 16

1:	mov (%r14, %rsi, 4), %edx	// j = SA[i]

	test %edx, %edx			// j > 0?
	jle 2f				// if j > 0 ...
	cmp %edx, %r9d			// 	n <= j
	jg 3f				// 	if n <= j
	inc %r11d			//		d++
	sub %r9d, %edx			//		j -= n

3:	MOVCHR (%r12, %rdx, CHR), %eax	// 	c0 = T[j]
	cmp %eax, %ecx			//	c0 == c1?
	je 3f				//	if c0 != c1
	mov %edi, (%rbx, %rcx, 4)	//		B[c1] = b
	mov (%rbx, %rax, 4), %edi	//		b = B[c0]
	mov %eax, %ecx			//		c1 = c0

3:	dec %edx			//	j--
#if CHR == 1
	cmp (%r12, %rdx, CHR), %cl	//	c1 < T[j]
#else
	cmp (%r12, %rdx, CHR), %ecx	//	c1 < T[j]
#endif
	mov %eax, %r10d			//	t = c0
	adc %r10d, %r10d		//	t = (c0 << 1) + (c1 < T[j])

	cmp %r11d, (%r8, %r10, 4)	//	D[t] != d?
	je 3f				//	if D[t] != d
	add %r9d, %edx			//		j += n
	mov %r11d, (%r8, %r10, 4)	//		D[t] = d

3:	shr %r10d			//	t & 1?
	sbb %eax, %eax			//	tmp = t & 1 ? -1 : 0
	sub %eax, %edx			//	j += t & 1
	xor %edx, %eax			//	tmp = t & 1 ? ~(j + 1) : j
	dec %edi			//	b--
	mov %eax, (%r14, %rdi, 4)	//	SA[b] = tmp
	movl $0, (%r14, %rsi, 4)	//	SA[i] = 0

2:	dec %esi			// i--
	jns 1b				// if i >= 0, continue

0:	pop %r14			// -16 restore non-volatile registers
	pop %r12			// - 8
	pop %rbx			// - 0
	pop %rbp			// tear down stack frame
	ret
	.size TSIZE(LMSsort2), .-TSIZE(LMSsort2)

	// LMSsort1(T, SA, C, B, n, k)
	// sort LMS substrings according to strategy 1
	// (for when n <= SAIS_LMSSORT2_LIMIT)
	.type TSIZE(LMSsort1), @function
	.align 16
TSIZE(LMSsort1):
	// T:  rdi
	// SA: rsi
	// C: -8(%rbp) (orig. in rdx)
	// B: rcx
	// n: r8d
	// k: -16(%rbp) (orig. in r9d)
	// b: r9
	// i: r10d
	// j: edx
	// c0, tmp: eax
	// c1: r11d

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

	push %rdx			// - 8 save C
	push %r9			// -16 save k

	push %rcx			// -24 save B
	push %rdi			// -32 save T
	push %rsi			// -40 save SA
	push %r8			// -48 save n

	cmp %rdx, %rcx			// B == C?
	jne 0f				// skip call to getCounts if B != C

	// getCounts(T, C, n, k)
					// setup argument 1: T (already set up)
	mov %rdx, %rsi			// setup argument 2: C
	mov %r8d, %edx			// setup argument 3: n
	mov %r9d, %ecx			// setup argument 4: k

	call TSIZE(getCounts)

0:	// getBucketStarts(C, B, k)
	mov -8(%rbp), %rdi		// setup argument 1: C
	mov -24(%rbp), %rsi		// setup argument 2: B
	mov -16(%rbp), %edx		// setup argument 3: k

	call getBucketStarts

	pop %r8				// -40 restore n
	pop %rsi			// -32 restore SA
	pop %rdi			// -24 restore T
	pop %rcx			// -16 restore B

	// initialise local variables
	lea -2(%r8), %edx		// j = n - 2
	MOVCHR -CHR(%rdi, %r8, CHR), %r11d	// c1 = T[n - 1]
	mov (%rcx, %r11, 4), %r9d	// b = B[c1]
#if CHR == 1
	cmp %r11b, (%rdi, %rdx, CHR)	// T[j] < c1
#else
	cmp %r11d, (%rdi, %rdx, CHR)	// T[j] < c1
#endif
	sbb %eax, %eax			// tmp = T[j] < c1 ? -1 : 0
	xor %edx, %eax			// tmp = T[j] < c1 ? ~j : j
	mov %eax, (%rsi, %r9, 4)	// SA[b] = tmp
	inc %r9d			// b++

	xor %r10d, %r10d		// i = 0
	cmp %r8d, %r10d			// i < n?
	jge 0f				// skip loop if i >= n

	// for (i = 0; i < n; i++)
	.align 16
1:	mov (%rsi, %r10, 4), %edx	// j = SA[i]
	test %edx, %edx			// edx > 0?
	jle 2f				// if 0 < j:

	MOVCHR (%rdi, %rdx, CHR), %eax	//	c0 = T[j]
	cmp %eax, %r11d			//	c0 == c1?
	je 3f				//	if c0 != c1:

	mov %r9d, (%rcx, %r11, 4)	//		B[c1] = b
	mov %eax, %r11d			//		c1 = c0
	mov (%rcx, %r11, 4), %r9d	//		b = B[c1]

3:	dec %edx			//	j--
#if CHR == 1
	cmp %r11b, (%rdi, %rdx, CHR)	// 	T[j] < c1
#else
	cmp %r11d, (%rdi, %rdx, CHR)	// 	T[j] < c1
#endif
	sbb %eax, %eax			// 	tmp = T[j] < c1 ? -1 : 0
	xor %edx, %eax			// 	tmp = T[j] < c1 ? ~j : j
	mov %eax, (%rsi, %r9, 4)	// 	SA[b] = tmp
	inc %r9d			// 	b++
	movl $0, (%rsi, %r10, 4)	//	SA[i] = 0
	jmp 3f

2:	je 3f				// else if 0 > j:

	not %edx			//	j = ~j
	mov %edx, (%rsi, %r10, 4)	//	SA[i] = j

3:	inc %r10d			// i++
	cmp %r8d, %r10d			// i < n?
	jl 1b				// continue if i < n

0:	push %rcx			// -24 save B
	push %rdi			// -32 save T
	push %rsi			// -40 save SA
	push %r8			// -48 save n

	cmp %rcx, -8(%rbp)		// C == B?
	jne 0f				// if yes, skip call to getCounts

	// getCounts(T, C, n, k)
					// setup argument 1: T (already set up)
	mov %rcx, %rsi			// setup argument 2: C (C == B)
	mov %r8d, %edx			// setup argument 3: n
	mov -16(%rbp), %ecx		// setup argument 4: k

	call TSIZE(getCounts)

0:	// getBucketEnds(C, B, k)
	mov -8(%rbp), %rdi		// setup argument 1: C
	mov -24(%rbp), %rsi		// setup argument 2: B
	mov -16(%rbp), %edx		// setup argument 3: k

	call getBucketEnds

	pop %r8				// -40 restore n
	pop %rsi			// -32 restore SA
	pop %rdi			// -24 restore T
	pop %rcx			// -16 restore B

	// initialise local variables for final loop
	mov %r8d, %r10d			// i = n
	xor %r11d, %r11d		// c1 = 0
	mov (%rcx), %r9d		// b = B[0]

	dec %r10d			// i = n - 1
	js 0f				// skip loop if i < 0

	// for (i = n - 1; 0 <= i; i--)
	.align 16
1:	mov (%rsi, %r10, 4), %edx	// j = SA[i]
	test %edx, %edx			// edx > 0?
	jle 2f				// if 0 < j:

	MOVCHR (%rdi, %rdx, CHR), %eax	//	c0 = T[j]
	cmp %eax, %r11d			//	c0 == c1?
	je 3f				//	if c0 != c1:

	mov %r9d, (%rcx, %r11, 4)	//		B[c1] = b
	mov %eax, %r11d			//		c1 = c0
	mov (%rcx, %r11, 4), %r9d	//		tmp = B[c1]

3:	dec %edx			//	j--
#if CHR == 1
	cmp (%rdi, %rdx, CHR), %r11b	// 	c1 < T[j]
#else
	cmp (%rdi, %rdx, CHR), %r11d	// 	c1 < T[j]
#endif
	sbb %eax, %eax			// 	tmp = c1 < T[j] ? -1 : 0
	sub %eax, %edx			//	j += c1 < T[j]
	xor %edx, %eax			// 	tmp = c1 < T[j] ? ~(j + 1) : j
	dec %r9d			// 	b--
	mov %eax, (%rsi, %r9, 4)	// 	SA[b] = tmp
	movl $0, (%rsi, %r10, 4)	//	SA[i] = 0

2:	dec %r10d			// i--
	jns 1b				// continue if i >= 0

0:	leave				// tear down the stack frame
	ret
	.size TSIZE(LMSsort1), .-TSIZE(LMSsort1)

	// LMSpostproc1(T, SA, n, m)
	// postprocess the result of LMSsort1
	.type TSIZE(LMSpostproc1), @function
	.align 16
TSIZE(LMSpostproc1):
	// T: rdi
	// SA: rsi
	// n: edx
	// m: ecx
	// name, tmp: eax
	// ebx: i
	// r14d: j
	// r8d: p
	// r9d: q
	// r10d: plen, tmp
	// r11d: qlen
	// r12: c0, &T[p]
	// r13: c1, &T[q]

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

	push %rbx			// save non-volatile registers
	push %r12
	push %r13
	push %r14

	xor %ebx, %ebx			// i = 0
	mov (%rsi), %r8d		// p = SA[0]
	xor $-1, %r8d			// ~p?
	js 0f				// if p < 0, skip loop

	// for (i = 0; (p = SA[i]) < 0; ++i)
	.align 16
1:	mov %r8d, (%rsi, %rbx, 4)	// SA[i] = ~p
	inc %ebx			// i++
	mov (%rsi, %rbx, 4), %r8d	// p = SA[++i]
	xor $-1, %r8d			// ~p?
	jns 1b				// if p < 0 continue

0:	cmp %ecx, %ebx			// i < m?
	jge 0f				// if i < m:

	mov %ebx, %r14d			//	j = i
	inc %ebx			//	i++

	// for (j = i++; ; ++i)
	.align 16
1:	mov (%rsi, %rbx, 4), %r8d	//	p = SA[i]
	xor $-1, %r8d			//	~p
	js 2f				//	if p < 0:

	mov %r8d, (%rsi, %r14, 4)	//		SA[j++] = ~p
	inc %r14d			//		j++
	movl $0, (%rsi, %rbx, 4)	//		SA[i] = 0
	cmp %ecx, %r14d			//		j == m?
	je 0f				//		if j == m break

2:	inc %ebx			//	i++
	jmp 1b				//	continue

0:	lea -1(%rdx), %ebx		// i = n - 1
	mov %ebx, %r14d			// j = n - 1
	MOVCHR (%rdi, %rbx, CHR), %r12d	// c0 = T[n - 1]

	// for (;;)
	// do c1 = c0; while (0 <= --i && (c0 = T[i]) >= c1);
	.align 16
1:	mov %r12d, %r13d		// c1 = c0
	dec %ebx			// i--
	jl 0f				// if i < 0 break (and skip the next loop, too)
	MOVCHR (%rdi, %rbx, CHR), %r12d	// c0 = T[i]
	cmp %r13d, %r12d		// c0 >= c1?
	jae 1b				// if c0 >= c1 continue

	// do c1 = c0; while (0 <= --i && (c0 = T[i]) <= c1);
	.align 16
2:	mov %r12d, %r13d		// c1 = c0
	dec %ebx			// i--
	jl 0f				// if i < 0 break (out of both loops)
	MOVCHR (%rdi, %rbx, CHR), %r12d	// c0 = T[i]
	cmp %r13d, %r12d		// c0 <= c1?
	jbe 2b				// if c0 <= c1 continue

	// no 0 <= i check here since we would have
	// taken the js 0f jump anyway if that was the case
	sub %ebx, %r14d			// j - i
	lea 1(%rbx), %eax		// tmp = i + 1
	shr %eax			// tmp = i + 1 >> 1
	add %ecx, %eax			// tmp = m + (i + 1 >> 1)
	mov %r14d, (%rsi, %rax, 4)	// SA[tmp] = j - i
	lea 1(%rbx), %r14d		// j = i + 1

	jmp 1b				// continue the outer loop

0:	xor %ebx, %ebx			// i = 0
	xor %eax, %eax			// name = 0
	mov %edx, %r9d			// q = n
	xor %r11d, %r11d		// qlen = 0

	cmp %ecx, %ebx			// i < m?
	jge 0f				// skip loop if i >= m

	// for (...; i < m; i++)
	.align 16
1:	mov (%rsi, %rbx, 4), %r8d	// p = SA[i]
	mov %r8d, %r10d			// p
	shr %r10d			// p >> 1
	add %ecx, %r10d			// m + (p >> 1)
	mov (%rsi, %r10, 4), %r10d	// plen = SA[m + (p >> 1)]

	cmp %r11d, %r10d		// plen == qlen?
	mov %r10d, %r11d		// qlen = plen (to free r10)
	jne 2f				// if plen != qlen, skip

	lea (%r9, %r11, 1), %r14d	// tmp = q + qlen
	cmp %edx, %r14d			// q + qlen < n?
	jge 2f				// if q + qlen >= n skip

					// if plen == qlen && q + qlen >= n:
	xor %r14d, %r14d		// 	j = 0
	lea (%rdi, %r8, CHR), %r12	// 	r12 = &T[p]
	lea (%rdi, %r9, CHR), %r13	// 	r13 = &T[q]

	cmp %r11d, %r14d		// 	j < qlen?
	jge 4f				//	if j == qlen, we have a match

	// for (j = 0; j < qlen && T[p + j] == T[q + j]; ++j);
	.align 16
3:	MOVCHR (%r12, %r14, CHR), %r10d	//	tmp = T[p + j]
#if CHR == 1
	cmp (%r13, %r14, CHR), %r10b	//	T[p + j] == T[q + j]?
#else
	cmp (%r13, %r14, CHR), %r10d	//	T[p + j] == T[q + j]?
#endif
	jne 2f				//	if T[p + j] != T[q + j] break
	inc %r14d			//	j++
	cmp %r11d, %r14d		//	j < qlen?
	jge 4f				//	if j == qlen, we have a match
	jmp 3b

2:	inc %eax			// name++
	mov %r8d, %r9d			// q = p
					// qlen = plen was already done above

4:	shr %r8d			// p >> 1
	add %ecx, %r8d			// m + (p >> 1)
	mov %eax, (%rsi, %r8, 4)	// SA[m + (p >> 1)] = name

	inc %ebx			// i++
	cmp %ecx, %ebx			// i < m?
	jl 1b				// if i < m continue

0:	pop %r14			// restore non-volatile registers
	pop %r13
	pop %r12
	pop %rbx

	pop %rbp			// tear down stack frame
	ret
	.size TSIZE(LMSpostproc1), .-TSIZE(LMSpostproc1)

	// induceSA(T, SA, C, B, n, k)
	// induce the suffix array from buckets B and the prepared array SA
	.type TSIZE(induceSA), @function
	.align 16
TSIZE(induceSA):
	// T: rdi
	// SA: r11 (orig. in rsi)
	// C: -8(%rbp) (orig. in rdx)
	// B: rcx
	// n: r8
	// k: -16(%rbp) (orig. in r9)
	// b: r10
	// i: r9
	// j: esi
	// c0, tmp: eax
	// c1: edx

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

	push %rdx			// - 8 save C
	push %r9			// -16 save k
	push %rdi			// -24 save T
	push %rsi			// -32 save SA
	push %rcx			// -40 save B
	push %r8			// -48 save n

	cmp %rcx, %rdx			// B == C?
	jne 0f				// skip call if equal

	// getCounts(T, C, n, k)
					// set up argument 1: T (already set up)
	mov %rdx, %rsi			// set up argument 2: C
	mov %r8d, %edx			// set up argument 3: n
	mov %r9d, %ecx			// set up argument 4: k

	call TSIZE(getCounts)

	// getBucketStarts(C, B, k)
0:	mov -8(%rbp), %rdi		// set up argument 1: C
	mov -40(%rbp), %rsi		// set up argument 2: B
	mov -16(%rbp), %edx		// set up argument 3: k

	call getBucketStarts

	pop %r8				// -40 restore n
	pop %rcx			// -32 restore B
	pop %r11			// -24 restore SA
	pop %rdi			// -16 restore T

	// initialise local variables
	lea -1(%r8), %esi		// j = n - 1
	MOVCHR (%rdi, %rsi, CHR), %edx	// c1 = T[j]
	mov (%rcx, %rdx, 4), %r10d	// b = B[c1]
	test %esi, %esi			// 0 < j?
	jle 0f				// if 0 < j:

#if CHR == 1
	cmp %dl, -CHR(%rdi, %rsi, CHR)	//	T[j - 1] < c1?
#else
	cmp %edx, -CHR(%rdi, %rsi, CHR)	//	T[j - 1] < c1?
#endif

	sbb %eax, %eax			//	tmp = T[j - 1] < c1 ? -1 : 0
	xor %eax, %esi			//	j = T[j - 1] < c1 ? ~j : j
0:	mov %esi, (%r11, %r10, 4)	// SA[b] = 0 < j && T[n - 2] < c1 ? ~j : j
	inc %r10d			// b++

	xor %r9d, %r9d			// i = 0
	cmp %r8d, %r9d			// i < n?
	jge 0f				// if i >= n skip loop

	// for (i = 0; i < n; i++)
	.align 16
1:	mov (%r11, %r9, 4), %esi	// j = SA[i]
	not %esi
	mov %esi, (%r11, %r9, 4)	// SA[i] = ~j

	xor $-1, %esi			// j = ~j, 0 < j?
	jle 2f				// if 0 < j:

	dec %esi			//	j--
	MOVCHR (%rdi, %rsi, CHR), %eax	//	c0 = T[j]
	cmp %edx, %eax			//	c0 == c1?
	je 3f				//	if c0 != c1:

	mov %r10d, (%rcx, %rdx, 4)	//		B[c1] = b
	mov %eax, %edx			//		c1 = c0
	mov (%rcx, %rdx, 4), %r10d	//		b = B[c1]

3:	inc %r10d			//	b++
	test %esi, %esi			//	0 < j?
	jle 4f				//	if 0 < j:

#if CHR == 1
	cmp %dl, -CHR(%rdi, %rsi, CHR)	//		T[j - 1] < c1?
#else
	cmp %edx, -CHR(%rdi, %rsi, CHR)	//		T[j - 1] < c1?
#endif
	sbb %eax, %eax			//		tmp = T[j - 1] < c1 ? -1 : 0
	xor %eax, %esi			//		tmp = T[j - 1] < c1 ? ~j : j
4:	mov %esi, -4(%r11, %r10, 4)	//		SA[b - 1] = j

2:	inc %r9d			// i++
	cmp %r8d, %r9d			// i < n?
	jl 1b				// if i < n continue

0:	push %rdi			// -24 save T
	push %r11			// -32 save SA
	push %rcx			// -40 save B
	push %r8			// -48 save n

	cmp %rcx, -8(%rbp)		// B == C?
	jne 0f				// if B != C skip call to getCounts

	// getCounts(T, C, n, k)
					// set up argument 1: T (already set up)
	mov %rcx, %rsi			// set up argument 2: C (C == B)
	mov %r8d, %edx			// set up argument 3: n
	mov -16(%rbp), %ecx		// set up argument 4: k

	call TSIZE(getCounts)

	// getBucketEnds(C, B, k)
0:	mov -8(%rbp), %rdi		// set up argument 1: C
	mov -40(%rbp), %rsi		// set up argument 2: B
	mov -16(%rbp), %edx		// set up argument 3: k

	call getBucketEnds

	pop %r8				// -40 restore n
	pop %rcx			// -32 restore B
	pop %r11			// -24 restore SA
	pop %rdi			// -16 restore T

	// initialise local variables for loop
	mov (%rcx), %r10d		// b = B[0]
	lea -1(%r8), %r9d		// i = n - 1
	xor %edx, %edx			// c1 = 0
	test %r9d, %r9d			// 0 <= i?
	jl 0f

	// for (...; 0 <= i; --i)
	.align 16
1:	mov (%r11, %r9, 4), %esi	// j = SA[i]
	test %esi, %esi			// 0 < j?
	jle 2f				// if 0 < j

	dec %esi			//	j--
	MOVCHR (%rdi, %rsi, CHR), %eax	//	c0 = T[j]
	cmp %edx, %eax			//	c0 == c1?
	je 3f				//	if c0 != c1:

	mov %r10d, (%rcx, %rdx, 4)	//		B[c1] = b
	mov %eax, %edx			//		c1 = c0
	mov (%rcx, %rdx, 4), %r10d	//		b = B[c1]

3:	dec %r10d			//	b--
	test %esi, %esi			//	j == 0?
	je 4f				//	if j == 0:

#if CHR == 1
	cmp -CHR(%rdi, %rsi, CHR), %dl	//		c1 < T[j - 1]?
#else
	cmp -CHR(%rdi, %rsi, CHR), %edx	//		c1 < T[j - 1]?
#endif
	cmc
	sbb %eax, %eax			//		tmp = c1 < T[j - 1] ? 0 : -1
	xor %eax, %esi			//		tmp = c1 < T[j - 1] ? j : ~j
4:	not %esi			//		tmp = c1 < T[j - 1] ? ~j : j
	mov %esi, (%r11, %r10, 4)	//		SA[b] = c1 < T[j - 1] ? ~j : j
	jmp 3f

2:					// else
	not %esi			//	~j
	mov %esi, (%r11, %r9, 4)	//	SA[i] = ~j

3:	dec %r9d			// i--
	jge 1b				// if 0 <= i continue

0:	leave				// tear down stack frame
	ret
	.size TSIZE(induceSA), .-TSIZE(induceSA)

	// sais_stage1(T, SA, C, B, n, k, name, flags)
	// reduce the problem by at least 1/2 and sort all LMS substrings
	.type TSIZE(sais_stage1), @function
	.align 16
TSIZE(sais_stage1):
	// T: r12 (orig. in rdi)
	// SA: rbx (orig. in rsi)
	// C: -24(%rbp), edi (orig. in rdx)
	// B: r13 (orig. in rcx)
	// n: -8(%rbp), r8d (orig. in r8)
	// k: -16(%rbp), r9d (orig. in r9)
	// name: 16(%rbp)
	// *name: eax
	// flags: 24(%rbp)
	// D: -56(%rbp), esi
	// b: r8
	// i: ecx
	// j: edx
	// m: -4(%rbp), eax
	// c0: esi
	// c1: edi
	// t: -12(%rbp)

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

	push %r8			// - 8 save n
	push %r9			// -16 save k
	push %rdx			// -24 save C

	push %rbx			// -32 preserve rbx
	push %r12			// -40 preserve r12
	push %r13			// -48 preserve r13

	mov %rdi, %r12			// set up T register
	mov %rsi, %rbx			// set up SA register
	mov %rcx, %r13			// set up B register

	// getCounts(T, C, n, k)
					// set up argument 1: T (already set up)
	mov %rdx, %rsi			// set up argument 2: C
	mov %r8d, %edx			// set up argument 3: n
	mov %r9d, %ecx			// set up argument 4: k

	call TSIZE(getCounts)

	// getBucketEnds(C, B, k)
	mov -24(%rbp), %rdi		// set up argument 1: C
	mov %r13, %rsi			// set up argument 2: B
	mov -16(%rbp), %edx		// set up argument 3: k

	call getBucketEnds

	// clear SA
	mov %rbx, %rdi			// SA
	mov -8(%rbp), %ecx		// n
	xor %eax, %eax			// tmp = 0
	rep stosl			// for (i = 0; i < n; i++) SA[i] = tmp;

	// initialise local variables
	mov -8(%rbp), %edx		// j = n
	lea -1(%rdx), %ecx		// i = n - 1
	lea -12(%rbp), %r8		// b = &t
	MOVCHR (%r12, %rcx, CHR), %esi	// c0 = T[i]
	xor %eax, %eax			// m = 0

	// do c1 = c0; while (0 <= --i && (c0 = T[i]) >= c1);
	.align 16
1:	mov %esi, %edi			// c1 = c0
	dec %ecx			// i--
	jl 0f				// if i < 0, skip all sorts of loops

	MOVCHR (%r12, %rcx, CHR), %esi	// c0 = T[i]
	cmp %edi, %esi			// c0 >= c1?
	jae 1b				// ... while (c0 >= c1)

	.align 16
	// do c1 = c0; while (0 <= --i && (c0 = T[i]) >= c1);
2:	mov %esi, %edi			// c1 = c0
	dec %ecx			// i--
	jl 0f				// if i < 0, break out of both loops

	MOVCHR (%r12, %rcx, CHR), %esi	// c0 = T[i]
	cmp %edi, %esi			// c0 <= c1?
	jbe 2b				// ... while (c0 <= c1)

	// 0 <= i has just been checked, so the
	// if statement vanishes
	mov %edx, (%r8)			// *b = j
	mov (%r13, %rdi, 4), %r8d	// B[c1]
	dec %r8d			// tmp--
	mov %r8d, (%r13, %rdi, 4)	// --B[c1]
	lea (%rbx, %r8, 4), %r8		// b = SA + --B[c1]
	mov %ecx, %edx			// j = i
	inc %eax			// m++

	jmp 1b				// continue

0:	mov %eax, -4(%rbp)		// save m to the stack
	cmp $1, %eax			// m > 1?
	jle 1f				// if m > 1:

	testb $48, 24(%rbp)		//	flags & 48?
	jz 2f				//	if flags & 48:


	MOVCHR CHR(%r12, %rdx, CHR), %eax	//	tmp = T[j + 1]
	incl (%r13, %rax, 4)		//		++B[t[j + 1]]

	testb $16, 24(%rbp)		//		flags & 16?
	jz 3f				//		if flags & 16:

	// malloc(k * 2 * sizeof *D)
	mov -16(%rbp), %edi
	shl $3, %rdi			//			set up arg 1: 8 * k

	call malloc@PLT

	mov %rax, %rsi			//			D = result
	test %rsi, %rsi			//			were we succesful?
	jnz 4f				//			if malloc failed:

	mov $-2, %eax			//				result = -2
	jmp 9f				//				go to return sequence

					//		else:
3:	imul $-2, -16(%rbp), %esi	//			-2 * k
	movslq %esi, %rsi		//			sign extend
	lea (%r13, %rsi, 4), %rsi	//			D = B - 2 * k

4:	mov -24(%rbp), %rdi		//		load C
	mov -8(%rbp), %r8d		//		load n
	mov -16(%rbp), %r9d		//		load k
	xor %ecx, %ecx			//		i = 0
	xor %edx, %edx			//		j = 0

	// for (i = 0, j = 0; i < k; ++i)
	.align 16
3:	add (%rdi, %rcx, 4), %edx	//		j += C[i]
	mov (%r13, %rcx, 4), %eax	//		tmp = B[i]
	cmp %edx, %eax			//		B[i] != j?
	je 4f				//		if B[i] != j:

	add %r8d, (%rbx, %rax, 4)	//			SA[B[i]] += n

4:	xor %eax, %eax			//		0
	lea (%rcx, %r9, 1), %r10	//		i + k
	mov %eax, (%rsi, %rcx, 4)	//		D[i] = 0
	mov %eax, (%rsi, %r10, 4)	//		D[i + k] = 0

	inc %ecx			//		i++
	cmp %r9d, %ecx			//		i < k?
	jl 3b				//		if i < k continue

	push %rsi			//		-56 save D

	// LMSsort2(T, SA, C, B, D, n, k)
	push %r9			//		-64 set up argument 7: k
	mov %rdi, %rdx			//		set up argument 3: C
	mov %r12, %rdi			//		set up argument 1: T
	mov %r8d, %r9d			//		set up argument 6: n
	mov %rsi, %r8			//		set up argument 5: D
	mov %rbx, %rsi			//		set up argument 2: SA
	mov %r13, %rcx			//		set up argument 4: B

	call TSIZE(LMSsort2)

	// LMSpostproc2(SA, n, m)
	mov %rbx, %rdi			//		set up argument 1: SA
	mov -8(%rbp), %esi		//		set up argument 2: n
	mov -4(%rbp), %edx		//		set up argument 3: m

	call LMSpostproc2		//		tmpname = result
	mov %eax, -64(%rbp)		//		stash tmpname away

	testb $16, 24(%rbp)		//		flags & 16?
	jz 3f				//		if flags & 16:

	// free(D)
	mov -56(%rbp), %rdi		//			set up argument 1: D

	call free@PLT

3:	pop %rax			//		-56 restore tmpname
	pop %rdi			//		-48 restore D (points to nowhere)
	jmp 0f				//	else

2:	// LMSsort1(T, SA, C, B, n, k)
	mov %r12, %rdi			//		set up argument 1: T
	mov %rbx, %rsi			//		set up argument 2: SA
	mov -24(%rbp), %rdx		//		set up argument 3: C
	mov %r13, %rcx			//		set up argument 4: B
	mov -8(%rbp), %r8d		//		set up argument 5: n
	mov -16(%rbp), %r9d		//		set up argument 6: k

	call TSIZE(LMSsort1)

	// tmpname = LMSpostproc1(T, SA, n, m)
	mov %r12, %rdi			//		set up argument 1: T
	mov %rbx, %rsi			//		set up argument 2: SA
	mov -8(%rbp), %rdx		//		set up argument 3: n
	mov -4(%rbp), %rcx		//		set up argument 4: m

	call TSIZE(LMSpostproc1)	//		tmpname = result

	jmp 0f				// else

1:	jl 1f				//	if m == 1:

	inc %edx			//		j++
	mov %edx, (%r8)			//		b = j + 1
	mov $1, %eax			//		tmpname = 1
	jmp 0f				//	else

1:	xor %eax, %eax			//		tmpname = 0

0:	mov 16(%rbp), %rdi		// name
	mov %eax, (%rdi)		// *name = tmpname

	mov -4(%rbp), %eax		// result = tmpname

9:	pop %r13			// -40 restore r13
	pop %r12			// -32 restore r12
	pop %rbx			// -24 restore rbx

	leave				// release stack frame
	ret
	.size TSIZE(sais_stage1), .-TSIZE(sais_stage1)

	// sais_stage2(T, SA, n, k, m, name, fs, flags)
	// solve the reduced problem by recursion
	.type TSIZE(sais_stage2), @function
	.align 16
TSIZE(sais_stage2):
	// T: rdi
	// SA: rsi
	// n: edx
	// k: ecx
	// m: r8d
	// name: r9d
	// fs: eax, 16(%rbp)
	// flags: 24(%rbp)
	// RA: r10
	// i: ecx, edx
	// j: r11d
	// c0: eax
	// c1: ecx

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

	mov 16(%rbp), %eax		// fs
	add %edx, %eax			// fs + n
	sub %r8d, %eax			// fs + n - m
	sub %r8d, %eax			// fs = fs + n - 2*m

	testb $13, 24(%rbp)		// flags & 13?
	jnz 0f				// if flags & 13:

	lea (%r9, %rcx, 1), %r10d	//	k + name
	cmp %eax, %r10d			//	k + name <= fs?
	jg 1f				//	if k + name <= fs:

	sub %ecx, %eax			//		fs -= k
	jmp 0f				//	else

1:	orb $8, 24(%rbp)		//		flags |= 8

0:	lea (%r8, %rax, 1), %r10d	// m + fs
	lea (%rsi, %r10, 4), %r10	// RA = &SA[m + fs]

	push %rdx			// - 8 save n

	mov %edx, %ecx			// n
	lea -1(%r8), %r11d		// j = m - 1
	shr %ecx			// n >> 1
	add %r11d, %ecx			// i = n + (n >> 1) - 1

	cmp %ecx, %r8d			// m <= i
	jg 0f				// if m > i, skip loop

	// for (i = m + (n >> 1) -- 1, j = m - 1; m <= i; --i)
	.align 16
1:	mov (%rsi, %rcx, 4), %edx	// SA[i]
	test %edx, %edx			// SA[i] != 0?
	jz 2f				// if SA[i] != 0:

	dec %edx			// SA[i] - 1
	mov %edx, (%r10, %r11, 4)	// RA[j] = SA[i] - 1
	dec %r11d			// j--

2:	dec %ecx
	cmp %ecx, %r8d			// m <= i?
	jle 1b				// if m <= i continue

0:	push %rdi			// -16 save T
	push %rsi			// -24 save SA
	push %r8			// -32 save m
	push %r10			// -40 save RA
	push %rcx			// -48 align stack to 16 bytes

	// sais_main_idx(RA, SA, fs, m, name)
	mov %r10, %rdi			// set up argument 1: RA
					// set up argument 2: SA (already set up)
	mov %eax, %edx			// set up argument 3: fs
	mov %r8d, %ecx			// set up argument 4: m
	mov %r9d, %r8d			// set up argument 5: name

	call sais_main_idx

	test %eax, %eax			// success?
	jz 0f				// on failure:

	leave				// tear down stack frame
	ret				// return whatever sais_main_idx returned

0:	pop %rcx			// -40 unalign stack
	pop %r10			// -32 restore RA
	pop %r8				// -24 restore m
	pop %rsi			// -16 restore SA
	pop %rdi			// - 8 restore T
	pop %rdx			//   0 restore n

	// set up variables for final loop
	lea -1(%r8), %r11d		// j = m - 1
	dec %edx			// i = n - 1
	MOVCHR (%rdi, %rdx, CHR), %eax	// c0 = T[n - 1]

	// for (;;)
	// do c1 = c0; while (0 <= --i && (c0 = T[i]) >= c1);
	.align 16
1:	mov %eax, %ecx			// c1 = c0
	dec %edx			// i--
	jl 0f				// if 0 > i, break both loops

	MOVCHR (%rdi, %rdx, CHR), %eax	// c0 = T[i]
	cmp %ecx, %eax			// c0 >= c1?
	jae 1b				// if c0 >= c1, continue

	// do c1 = c0; while (0 <= --i && (c0 = T[i]) <= c1);
	.align 16
2:	mov %eax, %ecx			// c1 = c0
	dec %edx			// i--
	jl 0f				// if 0 > i, break both loops

	MOVCHR (%rdi, %rdx, CHR), %eax	// c0 = T[i]
	cmp %ecx, %eax			// c0 <= c1?
	jbe 2b				// if c0 <= c1, continue

	lea 1(%rdx), %r9d		// i + 1
	mov %r9d, (%r10, %r11, 4)	// RA[j] = i + 1
	dec %r11d			// j--
	jmp 1b				// continue

0:	xor %edx, %edx			// i = 0

	test $1, %r8b			// is m odd?
	jz 0f				// if yes, make i odd, too

	mov (%rsi, %rdx, 4), %eax	// SA[i]
	mov (%r10, %rax, 4), %eax	// RA[SA[i]]
	mov %eax, (%rsi, %rdx, 4)	// SA[i] = RA[SA[i]]
	inc %edx

0:	cmp %r8d, %edx			// i < m?
	jge 0f

	// for (i = 0; i < m; i++)
	.align 16
1:	mov (%rsi, %rdx, 4), %eax	// SA[i]
	mov (%r10, %rax, 4), %eax	// RA[SA[i]]
	mov %eax, (%rsi, %rdx, 4)	// SA[i] = RA[SA[i]]

	mov 4(%rsi, %rdx, 4), %eax	// SA[i + 1]
	mov (%r10, %rax, 4), %eax	// RA[SA[i + 1]]
	mov %eax, 4(%rsi, %rdx, 4)	// SA[i + 1] = RA[SA[i + 1]]

	add $2, %edx			// i += 2
	cmp %r8d, %edx			// i < m?
	jl 1b

0:	mov 24(%rbp), %eax		// return flags

	leave				// tear down stack frame
	ret
	.size TSIZE(sais_stage2), .-TSIZE(sais_stage2)

	// sais_stage3(T, SA, B, n, m)
	// put all left-most S characters into their buckets
	.type TSIZE(sais_stage3), @function
	.align 16
TSIZE(sais_stage3):
	// T: rdi
	// SA: rsi
	// B: rdx
	// n: ecx
	// m: r8d
	// i: r8d
	// j: ecx
	// p: eax
	// q: r9d
	// c0: r10d
	// c1: r11d

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

					// j = n
	dec %r8d			// i = m - 1
	mov (%rsi, %r8, 4), %eax	// p = SA[i]
	MOVCHR (%rdi, %rax, CHR), %r11d	// c1 = T[p]

	// do ... while (0 <= i)
	.align 16
0:	mov %r11d, %r10d		// c0 = c1
	mov (%rdx, %r10, 4), %r9d	// q = B[c0]

	cmp %ecx, %r9d			// q < j?
	jge 1f

	// while (q < j)
	.align 16
2:	dec %ecx			// j--
	movl $0, (%rsi, %rcx, 4)	// SA[j] = 0
	cmp %ecx, %r9d			// q < j?
	jl 2b

	// do ... while ((c1 = T[p]) == c0);
	.align 16
1:	dec %ecx			// j--
	mov %eax, (%rsi, %rcx, 4)	// SA[j] = p
	dec %r8d			// i--
	jl 0f				// if (i < 0) break out of both loops

	mov (%rsi, %r8, 4), %eax	// p = SA[i]
	MOVCHR (%rdi, %rax, CHR), %r11d	// c1 = T[p]
	cmp %r10d, %r11d		// c1 == c0?
	je 1b				// if c1 == c0 continue

	jmp 0b				// continue

	test %ecx, %ecx			// 0 < j?
	jl 0b

	// while (0 < j) SA[--j] = 0
0:	mov %rsi, %rdi
	xor %eax, %eax
	rep stosl			// memset(SA, 0, j)

	pop %rbp			// release stack frame
	ret
	.size TSIZE(sais_stage3), .-TSIZE(sais_stage3)

	// sais_main(T, SA, fs, n, k)
	// find the suffix array SA of T[0..n-1] in {0..k}^n
	.type TSIZE(sais_main), @function
	.align 16
TSIZE(sais_main):
	// T: -8(%rbp) (orig in rdi)
	// SA: -16(%rbp) (orig in rsi)
	// fs: -24(%rbp) (orig in edx)
	// n: -32(%rbp) (orig in ecx)
	// k: -40(%rbp) (orig in r8d)
	// C: -56(%rbp)
	// B: -64(%rbp)
	// m: -28(%rbp)
	// name: -20(%rbp)
	// flags: bl
	// newflags: eax

	push %rbp			// establish ...
	mov %rsp, %rbp			// stack frame

	push %rdi			// - 8 save T
	push %rsi			// -16 save SA
	push %rdx			// -24 save fs
	push %rcx			// -32 save n
	push %r8			// -40 save k
	push %rbx			// -48 save volatile register rbx
	xor %ebx, %ebx			// break dependency on rbx

	cmp $MINBUCKETSIZE, %r8d	// k <= MINBUCKETSIZE?
	jl 1f				// if k <= MINBUCKETSIZE:

	// malloc(k * sizeof *C)
	mov %r8d, %edi			//	set up argument 1: k * sizeof *C
	shl $2, %rdi
	call malloc@PLT
	test %rax, %rax			//	did malloc fail?
	jz 9f

	push %rax			//	-56 save C

	mov -24(%rbp), %edx		//	fs
	mov -40(%rbp), %ecx		//	k
	cmp %edx, %ecx			//	k <= fs?
	jg 2f				//	if k <= fs:

	sub %ecx, %edx			//		fs - k
	add -32(%rbp), %edx		//		fs - k + n
	shl $2, %edx			//		fs - k + n << 2
	add -16(%rbp), %rdx		//		B = &SA[fs - k + n]
	push %rdx			//		-64 save B

	mov $1, %bl			//		flags = 1

	jmp 0f				//	else

2:	movb $3, %bl			//		flags = 3

	// malloc(k * sizeof *B)
	mov %ecx, %edi			//		set up argument 1: k * sizeof *B
	shl $2, %rdi
	call malloc@PLT
	test %rax, %rax			//		did malloc fail?
	jnz 2f				//		if malloc failed:

	// free(C)
	pop %rdi			//			set up argument 1: C
	call free@PLT
	jmp 9f

2:	push %rax			//		-64 save B

	jmp 0f				// else

1:	cmp %edx, %r8d			// k <= fs?
	jg 1f				// else if k <= fs:

	sub %r8d, %edx			//	fs -= k
	add %edx, %ecx			//	n + fs - k
	lea (%rsi, %rcx, 4), %rax	//	C = &SA[n + fs - k]
	push %rax			//	-56 save C
	cmp %edx, %r8d			//	k <= fs - k
	jg 2f				//	if k <= fs - k:

	neg %r8				//		-k
	lea (%rax, %r8, 4), %rax	//		B = C - k
	push %rax			//		-64 save B
	xor %ebx, %ebx			//		flags = 0
	jmp 0f				//	else

2:	cmp $MINBUCKETSIZE*4, %r8d	// k <= MINBUCKETSIZE * 4?
	jg 2f				//	else if k <= MINBUCKETSIZE * 4:

	// malloc(k * sizeof *B)
	mov %r8d, %edi			//		set up argument 1: k * sizeof *B
	shl $2, %rdi
	call malloc@PLT
	test %rax, %rax			//		did malloc fail?
	jz 9f

	push %rax			//		-64 save B
	mov $2, %bl			//		flags = 2
	jmp 0f				//	else

2:	push %rax			//		-64 save B = C
	movb $8, %bl			//		flags = 8
	jmp 0f				//	else

	// malloc(k * sizeof *B)
1:	mov %r8d, %edi			//		set up argument 1: k * sizeof *B
	shl $2, %rdi
	call malloc@PLT
	test %rax, %rax			//		did malloc fail?
	jz 9f

	push %rax			//		-56 save C
	push %rax			//		-64 save B = C

0:	mov -32(%rbp), %r8d		// n
	mov -40(%rbp), %r9d		// k
	cmp $SAIS_LMSSORT2_LIMIT, %r8d	// n <= SAIS_LMSSORT2_LIMIT?
	jg 0f

	lea (%r9, %r9, 1), %ecx		// 2*k
	cmp %r8d, %ecx			// 2*k <= n
	jg 0f				// if n <= SAIS_LMSSORT2_LIMIT && 2 <= n / k

	mov -24(%rbp), %edx		// 	fs
	sub %r9d, %edx			// 	fs - k

	test $1, %bl			//	flags & 1
	jz 2f				//	if flags & 1

	cmp %edx, %ecx			//		2*k <= fs - k
	jle 1f				//		if (2*k <= fs -k) goto 1f

	or $16, %bl			//		flags |= 16
	jmp 0f				//	else

2:	test %bl, %bl			//	flags == 0?
	jnz 0f				//	else if flags == 0:

	sub %r9d, %edx			//		fs - 2*k
	cmp %edx, %ecx			//		2*k <= fs - k
	jg 0f				//		if (2*k <= fs - k):

1:	or $32, %bl			//			flags |= 32

0:	// m = sais_stage1(T, SA, C, B, n, k, &name, flags)
	mov -8(%rbp), %rdi		//     set up argument 1: T
	mov -16(%rbp), %rsi		//     set up argument 2: SA
	mov -56(%rbp), %rdx		//     set up argument 3: C
	mov -64(%rbp), %rcx		//     set up argument 4: B
					//     set up argument 5: n (already set up)
					//     set up argument 6: k (already set up)
	push %rbx			// -72 set up argument 8: flags
	lea -20(%rbp), %rax
	push %rax			// -80 set up argument 7: name

	call TSIZE(sais_stage1)
	test %eax, %eax			// result < 0?
	jge 0f

	test $5, %bl			// flags & 5?
	jz 1f				// if flags & 5:

	// free(C)
	mov -56(%rbp), %rdi
	call free@PLT

1:	test $1, %bl			// flags & 1?
	jz 1f

	// free(B)
	mov -64(%rbp), %rdi
	call free@PLT

1:	jmp 9f

0:	mov %eax, -28(%rbp)		// m = result

	cmp %eax, -20(%rbp)		// name < m?
	jge 0f				// if name < m:

	test $4, %bl			//	flags & 4?
	jz 1f				//	if flags & 4:

	// free(C)
	mov -56(%rbp), %rdi		//		set up argument 1: C
	call free@PLT

1:	test $2, %bl			//	if flags & 2:
	jz 1f

	// free(B)
	mov -64(%rbp), %rdi		//		set up argument 1: B
	call free@PLT

1:	// sais_stage2(T, SA, n, k, m, name, fs, flags)
	mov -8(%rbp), %rdi		//	set up argument 1: T
	mov -16(%rbp), %rsi		//	set up argument 2: SA
	mov -32(%rbp), %edx		//	set up argument 3: n
	mov -40(%rbp), %ecx		//	set up argument 4: k
	mov -28(%rbp), %r8d		//	set up argument 5: m
	mov -20(%rbp), %r9d		//	set up argument 6: name
	push %rbx			//	set up argument 8: flags
	push -24(%rbp)			//	set up argument 7: fs

	call TSIZE(sais_stage2)		//	newflags = sais_stage2(...)
	test %eax, %eax			//	newflags < 0?
	jge 1f				//	if newflags < 0:

3:	test $1, %bl			//		flags & 1?
	jz 2f				//		if flags & 1:

	// free(C)
	mov -56(%rbp), %rdi		//			set up argument 1: C
	call free@PLT

2:	jmp 9f				//		goto failure

1:	mov %eax, %ebx			//	flags = newflags
	test $4, %bl			//	flags & 4?
	jz 1f				//	if flags & 4:

	// malloc(k * sizeof *B)
	mov -40(%rbp), %edi		//		k
	shl $2, %rdi			//		k * sizeof *B
	call malloc@PLT			//		tmp = malloc(k * sizeof *B)
	test %rax, %rax			//		was succesfull?
	jz 9f				//		if unsucessful fail

	mov %rax, -56(%rbp)		//		C = tmp
	mov %rax, -64(%rbp)		//		B = tmp

1:	test $2, %bl			//	if flags & 2:
	jz 0f

	// malloc(k * sizeof *B)
	mov -40(%rbp), %edi		//		k
	shl $2, %rdi			//		k * sizeof *B
	call malloc@PLT			//		tmp = malloc(k * sizeof *B)
	test %rax, %rax			//		was succesfull?
	jz 3b				//		if malloc failed, proceed as if
					//		sais_stage2 failed

	mov %rax, -64(%rbp)		//		B = tmp

0:	test $8, %bl			// flags & 8?
	jz 0f				// if flags & 8:

	// getCounts(T, C, n, k)
	mov -8(%rbp), %rdi		//	set up argument 1: T
	mov -56(%rbp), %rsi		//	set up argument 2: C
	mov -32(%rbp), %edx		//	set up argument 3: n
	mov -40(%rbp), %ecx		//	set up argument 4: k

	call TSIZE(getCounts)

0:	cmpl $1, -28(%rbp)		// m > 1?
	jle 0f				// if m > 1:

	// getBucketEnds(C, B, k)
	mov -56(%rbp), %rdi		//	set up argument 1: C
	mov -64(%rbp), %rsi		//	set up argument 2: B
	mov -40(%rbp), %edx		//	set up argument 3: k

	call getBucketEnds

	// sais_stage3(T, SA, B, n, m)
	mov -8(%rbp), %rdi		//	set up argument 1: T
	mov -16(%rbp), %rsi		//	set up argument 2: SA
	mov -64(%rbp), %rdx		//	set up argument 3: B
	mov -32(%rbp), %ecx		//	set up argument 4: n
	mov -28(%rbp), %r8d		//	set up argument 5: m

	call TSIZE(sais_stage3)

0:	// induceSA(T, SA, C, B, n, k)
	mov -8(%rbp), %rdi		// set up argument 1: T
	mov -16(%rbp), %rsi		// set up argument 2: SA
	mov -56(%rbp), %rdx		// set up argument 3: C
	mov -64(%rbp), %rcx		// set up argument 4: B
	mov -32(%rbp), %r8d		// set up argument 5: n
	mov -40(%rbp), %r9d		// set up argument 6: k

	call TSIZE(induceSA)

	test $5, %bl			// flags & 5?
	jz 0f				// if flags & 5:

	// free(C)
	mov -56(%rbp), %rdi		//	set up argument 1: C
	call free@PLT

0:	test $2, %bl			// flags & 2?
	jz 0f				// if flags & 2:

	mov -64(%rbp), %rdi		//	set up argument 1: B
	call free@PLT

	// success path
0:	mov -48(%rbp), %rbx		// restore volatile register rbx
	xor %eax, %eax			// return value: 0
	leave				// tear down stack frame
	ret

	// failure path
9:	mov -48(%rbp), %rbx		// restore volatile register rbx
	mov $-2, %eax			// return value: -2
	leave				// release stack frame
	ret
	.size TSIZE(sais_main), .-TSIZE(sais_main)
