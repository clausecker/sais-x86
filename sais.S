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

	// this is the one function we export
	.globl sais

	// Copyright notice
	.ident "Copyright (c) 2008--2010 Yuta Mori"
	.ident "Copyright (c) 2018--2019 Robert Clausecker"

	// useful constants
#define UCHAR_SIZE 256
#define MINBUCKETSIZE 256
#define SAIS_LMSSORT2_LIMIT 0x3fffffff

	.text

	// sais(T, SA, n)
	// compute suffix array SA of length n for byte array T
	.type sais, @function
	.align 16
sais:	push %rbp		// establish ...
	mov %rsp, %rbp		// stack frame

	test %rdi, %rdi		// sanity check: T == NULL?
	je 0f

	test %rsi, %rsi		// sanity check: SA == NULL?
	je 0f

	test %edx, %edx		// sanity check: n >= 0?
	js 0f

	cmp $1, %edx		// n <= 1?
	ja 1f

	je 2f			// if n == 1, SA[0] == 0 trivially
	movl $0, (%rsi)

2:	xor %eax, %eax		// return 0
	pop %rbp		// release stack frame
	ret

	// sanity check failed
0:	mov $-1, %eax		// return -1
	pop %rbp		// release stack frame
	ret

	// normal operation
1:	mov %edx, %ecx
	xor %edx, %edx
	mov $256, %r8d
	pop %rbp		// release stack frame
	jmp sais_main_chr	// return sais_main_chr(T, SA, 0, n, 256);
	.size sais, .-sais

	// getBucketStarts(C, B, k)
	// store in B where the buckets in C start
	// TODO: vectorize
	.type getBucketStarts, @function
	.align 16
getBucketStarts:
	push %rbp		// establish ...
	mov %rsp, %rbp		// stack frame

	xor %eax, %eax		// sum = 0
	xor %ecx, %ecx		// i = 0

	test $1, %dl		// is k odd?
	jz 2f			// if yes, make i odd, too

	mov (%rdi), %r8d	// tmp = C[0]
	mov %eax, (%rsi)	// B[0] = sum
	add %r8d, %eax		// sum += tmp

	inc %ecx			// i++
2:	cmp %edx, %ecx			// i < k?
	jge 0f

	.align 16
1:	mov (%rdi, %rcx, 4), %r8d	// tmp = C[i]
	mov %eax, (%rsi, %rcx, 4)	// B[i] = sum
	add %r8d, %eax			// sum += tmp

	mov 4(%rdi, %rcx, 4), %r8d	// tmp = C[i + 1]
	mov %eax, 4(%rsi, %rcx, 4)	// B[i + 1] = sum
	add %r8d, %eax			// sum += tmp

	add $2, %ecx			// i += 2
	cmp %edx, %ecx			// i < k?
	jb 1b

0:	pop %rbp		// release stack frame
	ret
	.size getBucketStarts, .-getBucketStarts


	// getBucketEnds(C, B, k)
	// store in B where the buckets in C end
	// TODO: vectorize
	.type getBucketEnds, @function
	.align 16
getBucketEnds:
	push %rbp		// establish ...
	mov %rsp, %rbp		// stack frame

	xor %eax, %eax		// sum = 0
	xor %ecx, %ecx		// i = 0

	test $1, %dl		// is k odd?
	jz 2f			// if yes, make i odd, too

	add (%rdi), %eax	// sum += C[0]
	mov %eax, (%rsi)	// B[0] = sum

	inc %ecx		// i++
2:	cmp %edx, %ecx		// i < k?
	jge 0f

	.align 16
1:	add (%rdi, %rcx, 4), %eax	// sum += C[i]
	mov %eax, (%rsi, %rcx, 4)	// B[i] = sum

	add 4(%rdi, %rcx, 4), %eax	// sum += C[i + 1]
	mov %eax, 4(%rsi, %rcx, 4)	// B[i + 1] = sum

	add $2, %ecx
	cmp %edx, %ecx			// i < k?
	jb 1b

0:	pop %rbp		// release stack frame
	ret
	.size getBucketEnds, .-getBucketEnds

	// LMSpostproc2(SA, n, m)
	// postprocess the result of LMSsort2
	.type LMSpostproc2, @function
	.align 16
LMSpostproc2:
	// rdi: SA
	// esi: n
	// edx: m
	// eax: name
	// ecx: i
	// r8: j
	// r9: d

	push %rbp		// establish ...
	mov %rsp, %rbp		// stack frame

	xor %ecx, %ecx		// i = 0
	xor %eax, %eax		// name = 0
	xor %r10d, %r10d	// scratch = 0
	mov (%rdi), %r8d	// j = SA[0]
	test %r8d, %r8d		// j < 0?
	jns 0f			// if (j >= 0) break;

	.align 16
1:	not %r8d		// j = ~j
	cmp %esi, %r8d		// j >= n?
	sbb $-1, %eax		// name -= (j >= n) - 1
	mov %r8d, (%rdi, %rcx, 4)	// SA[i] = j

	inc %ecx		// i++
	mov (%rdi, %rcx, 4), %r8d	// j = SA[i]
	test %r8d, %r8d
	js 1b			// if (j < 0) continue;

0:	cmp %edx, %ecx
	jge 2f			// skip if i >= m

	mov %ecx, %r9d		// d = i

	.align 16
1:	inc %ecx		// i++
	mov (%rdi, %rcx, 4), %r8d	// j = SA[i]
	test %r8d, %r8d		// j < 0?
	jns 1b			// if (j >= 0) continue

	not %r8d		// j = ~j
	cmp %esi, %r8d		// j >= n?
	sbb $-1, %eax		// name -= (j >= n) - 1
	mov %r8d, (%rdi, %r9, 4)	// SA[d] = j
	inc %r9d		// d++
	movl $0, (%rdi, %rcx, 4)	// SA[i] = 0
	cmp %r9d, %edx		// d == m?
	jne 1b			// if d == m break else continue

2:	cmp %edx, %eax		// name < m?
	jge 8f			// if name < m then ... else ...

	// if (name < m)
	lea -1(%rdx), %ecx	// i = m - 1
	lea 1(%rax), %r9d	// d = name + 1
	test %ecx, %ecx		// 0 <= i?
	js 9f

	.align 16
1:	mov (%rdi, %rcx, 4), %r8d	// j = SA[i]

	cmp %r8d, %esi		// n <= j?
	jg 2f			// if n > j skip

	sub %esi, %r8d		// j -= n
	dec %r9d		// d--

2:	shr %r8d		// j >>= 1
	add %edx, %r8d		// j += m
	mov %r9d, (%rdi, %r8, 4)	// SA[j] = d (i.e. SA[m + (j >> 1)] = d)
	dec %ecx		// i--
	jns 1b			// if (0 <= i) continue

	jmp 9f			// skip else part

	// else
8:	xor %ecx, %ecx		// i = 0
	test %edx, %edx		// m > 0?
	jle 9f			// skip if m <= 0

	.align 16
1:	mov (%rdi, %rcx, 4), %r8d	// j = SA[i]
	cmp %r8d, %esi		// n <= j?
	jg 2f			// if n > j skip

	sub %esi, %r8d		// j -= n
	mov %r8d, (%rdi, %rcx, 4)	// SA[i] = j

2:	inc %ecx
	cmp %edx, %edx		// i < m
	jl 1b			// if i < m continue

	// end of if-then-else
9:	pop %rbp		// release stack frame
	ret			// return (name)
	.size LMSpostproc2, .-LMSpostproc2

	// the code in sais-chr.S is included twice:
	// once with T being an array of bytes and once with
	// T being an array of ints.  The CHR macro lists the
	// type size and the TSIZE macro is used to decorate
	// identifiers.  The macro MOVCHR expands to movzbl
	// or mov depending on which of operand size is used.
#define TSIZE(ident) ident ## _idx
#define CHR 4
#define MOVCHR mov
#include "sais-chr.S"
#undef MOVCHR
#undef CHR
#undef TSIZE

#define TSIZE(ident) ident ## _chr
#define CHR 1
#define MOVCHR movzbl
#include "sais-chr.S"
#undef MOVCHR
#undef CHR
#undef TSIZE
