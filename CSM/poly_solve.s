	.file	"poly_solve.c"
	.section	.debug_abbrev,"",@progbits
.Ldebug_abbrev0:
	.section	.debug_info,"",@progbits
.Ldebug_info0:
	.section	.debug_line,"",@progbits
.Ldebug_line0:
	.text
.Ltext0:
	.section	.rodata
.LC1:
	.string	"bad args in zrhqr"
	.text
.globl zrhqr
	.type	zrhqr, @function
zrhqr:
.LFB5:
	.file 1 "poly_solve.c"
	.loc 1 17 0
	pushl	%ebp
.LCFI0:
	movl	%esp, %ebp
.LCFI1:
	subl	$56, %esp
.LCFI2:
	.loc 1 24 0
	movl	$50, 12(%esp)
	movl	$1, 8(%esp)
	movl	$50, 4(%esp)
	movl	$1, (%esp)
	call	dmatrix
	movl	%eax, -4(%ebp)
	.loc 1 25 0
	cmpl	$50, 12(%ebp)
	jg	.L2
	movl	12(%ebp), %eax
	sall	$3, %eax
	addl	8(%ebp), %eax
	fldl	(%eax)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jp	.L22
	je	.L2
.L22:
	jmp	.L4
.L2:
	movl	$.LC1, (%esp)
	call	nrerror
.L4:
	.loc 1 26 0
	movl	$1, -8(%ebp)
	jmp	.L6
.L7:
	.loc 1 27 0
	movl	-4(%ebp), %eax
	addl	$4, %eax
	movl	(%eax), %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-8(%ebp), %edx
	movl	12(%ebp), %eax
	subl	%edx, %eax
	sall	$3, %eax
	addl	8(%ebp), %eax
	fldl	(%eax)
	fchs
	movl	12(%ebp), %eax
	sall	$3, %eax
	addl	8(%ebp), %eax
	fldl	(%eax)
	fdivrp	%st, %st(1)
	fstpl	(%ecx)
	.loc 1 28 0
	movl	$2, -12(%ebp)
	jmp	.L8
.L9:
	movl	-12(%ebp), %eax
	sall	$2, %eax
	addl	-4(%ebp), %eax
	movl	(%eax), %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldz
	fstpl	(%eax)
	addl	$1, -12(%ebp)
.L8:
	movl	-12(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L9
	.loc 1 29 0
	movl	-8(%ebp), %eax
	cmpl	12(%ebp), %eax
	je	.L11
	movl	-4(%ebp), %edx
	addl	$4, %edx
	movl	-8(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fld1
	fstpl	(%eax)
.L11:
	.loc 1 26 0
	addl	$1, -8(%ebp)
.L6:
	movl	-8(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L7
	.loc 1 31 0
	movl	12(%ebp), %eax
	movl	%eax, 4(%esp)
	movl	-4(%ebp), %eax
	movl	%eax, (%esp)
	call	balanc
	.loc 1 32 0
	movl	20(%ebp), %eax
	movl	%eax, 12(%esp)
	movl	16(%ebp), %eax
	movl	%eax, 8(%esp)
	movl	12(%ebp), %eax
	movl	%eax, 4(%esp)
	movl	-4(%ebp), %eax
	movl	%eax, (%esp)
	call	hqr
	fstp	%st(0)
	.loc 1 33 0
	movl	$2, -12(%ebp)
	jmp	.L14
.L15:
	.loc 1 34 0
	movl	-12(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	(%eax)
	fstpl	-32(%ebp)
	.loc 1 35 0
	movl	-12(%ebp), %eax
	sall	$3, %eax
	addl	20(%ebp), %eax
	fldl	(%eax)
	fstpl	-24(%ebp)
	.loc 1 36 0
	movl	-12(%ebp), %eax
	subl	$1, %eax
	movl	%eax, -8(%ebp)
	jmp	.L16
.L17:
	.loc 1 37 0
	movl	-8(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	(%eax)
	fldl	-32(%ebp)
	fucompp
	fnstsw	%ax
	sahf
	jae	.L18
	.loc 1 38 0
	movl	16(%ebp), %edx
	addl	$8, %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	addl	%eax, %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	(%eax)
	fstpl	(%edx)
	.loc 1 39 0
	movl	20(%ebp), %edx
	addl	$8, %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	addl	%eax, %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	addl	20(%ebp), %eax
	fldl	(%eax)
	fstpl	(%edx)
	.loc 1 36 0
	subl	$1, -8(%ebp)
.L16:
	cmpl	$0, -8(%ebp)
	jg	.L17
.L18:
	.loc 1 41 0
	movl	16(%ebp), %edx
	addl	$8, %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	-32(%ebp)
	fstpl	(%eax)
	.loc 1 42 0
	movl	20(%ebp), %edx
	addl	$8, %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	-24(%ebp)
	fstpl	(%eax)
	.loc 1 33 0
	addl	$1, -12(%ebp)
.L14:
	movl	-12(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L15
	.loc 1 44 0
	movl	$50, 16(%esp)
	movl	$1, 12(%esp)
	movl	$50, 8(%esp)
	movl	$1, 4(%esp)
	movl	-4(%ebp), %eax
	movl	%eax, (%esp)
	call	free_dmatrix
	.loc 1 45 0
	leave
	ret
.LFE5:
	.size	zrhqr, .-zrhqr
	.section	.rodata
	.align 8
.LC4:
	.long	0
	.long	1074790400
	.align 8
.LC6:
	.long	0
	.long	1073741824
	.align 8
.LC8:
	.long	1717986918
	.long	1072588390
	.text
.globl balanc
	.type	balanc, @function
balanc:
.LFB6:
	.loc 1 50 0
	pushl	%ebp
.LCFI3:
	movl	%esp, %ebp
.LCFI4:
	subl	$64, %esp
.LCFI5:
	.loc 1 54 0
	fldl	.LC4
	fstpl	-24(%ebp)
	.loc 1 55 0
	movl	$0, -12(%ebp)
	.loc 1 56 0
	jmp	.L24
.L25:
	.loc 1 57 0
	movl	$1, -12(%ebp)
	.loc 1 58 0
	movl	$1, -4(%ebp)
	jmp	.L26
.L27:
	.loc 1 59 0
	fldz
	fstpl	-32(%ebp)
	fldl	-32(%ebp)
	fstpl	-56(%ebp)
	.loc 1 60 0
	movl	$1, -8(%ebp)
	jmp	.L28
.L29:
	.loc 1 61 0
	movl	-8(%ebp), %eax
	cmpl	-4(%ebp), %eax
	je	.L30
	.loc 1 62 0
	movl	-8(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-4(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	fldl	-32(%ebp)
	faddp	%st, %st(1)
	fstpl	-32(%ebp)
	.loc 1 63 0
	movl	-4(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	fldl	-56(%ebp)
	faddp	%st, %st(1)
	fstpl	-56(%ebp)
.L30:
	.loc 1 60 0
	addl	$1, -8(%ebp)
.L28:
	movl	-8(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L29
	.loc 1 65 0
	fldl	-32(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jne	.L35
	jp	.L35
	jmp	.L33
.L35:
	fldl	-56(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jne	.L37
	jp	.L37
	jmp	.L33
.L37:
	.loc 1 66 0
	fldl	-56(%ebp)
	fldl	.LC6
	fdivrp	%st, %st(1)
	fstpl	-48(%ebp)
	.loc 1 67 0
	fld1
	fstpl	-40(%ebp)
	.loc 1 68 0
	fldl	-32(%ebp)
	faddl	-56(%ebp)
	fstpl	-64(%ebp)
	.loc 1 69 0
	jmp	.L38
.L39:
	.loc 1 70 0
	fldl	-40(%ebp)
	fadd	%st(0), %st
	fstpl	-40(%ebp)
	.loc 1 71 0
	fldl	-32(%ebp)
	fmull	-24(%ebp)
	fstpl	-32(%ebp)
.L38:
	.loc 1 69 0
	fldl	-32(%ebp)
	fldl	-48(%ebp)
	fucompp
	fnstsw	%ax
	sahf
	ja	.L39
	.loc 1 73 0
	fldl	-56(%ebp)
	fadd	%st(0), %st
	fstpl	-48(%ebp)
	.loc 1 74 0
	jmp	.L41
.L42:
	.loc 1 75 0
	fldl	-40(%ebp)
	fldl	.LC6
	fdivrp	%st, %st(1)
	fstpl	-40(%ebp)
	.loc 1 76 0
	fldl	-32(%ebp)
	fdivl	-24(%ebp)
	fstpl	-32(%ebp)
.L41:
	.loc 1 74 0
	fldl	-32(%ebp)
	fldl	-48(%ebp)
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	ja	.L42
	.loc 1 78 0
	fldl	-32(%ebp)
	faddl	-56(%ebp)
	fdivl	-40(%ebp)
	fldl	-64(%ebp)
	fldl	.LC8
	fmulp	%st, %st(1)
	fucompp
	fnstsw	%ax
	sahf
	ja	.L45
	jmp	.L33
.L45:
	.loc 1 79 0
	movl	$0, -12(%ebp)
	.loc 1 80 0
	fld1
	fdivl	-40(%ebp)
	fstpl	-48(%ebp)
	.loc 1 81 0
	movl	$1, -8(%ebp)
	jmp	.L46
.L47:
	movl	-4(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-4(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-8(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-48(%ebp)
	fstpl	(%ecx)
	addl	$1, -8(%ebp)
.L46:
	movl	-8(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L47
	.loc 1 82 0
	movl	$1, -8(%ebp)
	jmp	.L49
.L50:
	movl	-8(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-4(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-8(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-4(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-40(%ebp)
	fstpl	(%ecx)
	addl	$1, -8(%ebp)
.L49:
	movl	-8(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L50
.L33:
	.loc 1 58 0
	addl	$1, -4(%ebp)
.L26:
	movl	-4(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L27
.L24:
	.loc 1 56 0
	cmpl	$0, -12(%ebp)
	je	.L25
	.loc 1 87 0
	leave
	ret
.LFE6:
	.size	balanc, .-balanc
	.section	.rodata
.LC12:
	.string	"Too many iterations in hqr"
	.align 8
.LC11:
	.long	0
	.long	1071644672
	.align 8
.LC13:
	.long	0
	.long	1072168960
	.align 8
.LC14:
	.long	0
	.long	-1076101120
	.text
.globl hqr
	.type	hqr, @function
hqr:
.LFB7:
	.loc 1 92 0
	pushl	%ebp
.LCFI6:
	movl	%esp, %ebp
.LCFI7:
	pushl	%ebx
.LCFI8:
	subl	$212, %esp
.LCFI9:
	.loc 1 94 0
	fldz
	fstpl	-80(%ebp)
	fldz
	fstpl	-72(%ebp)
	fldz
	fstpl	-64(%ebp)
	.loc 1 97 0
	fldz
	fstpl	-48(%ebp)
	.loc 1 99 0
	movl	8(%ebp), %eax
	addl	$4, %eax
	movl	(%eax), %eax
	addl	$8, %eax
	fldl	(%eax)
	fabs
	fstpl	-56(%ebp)
	.loc 1 100 0
	movl	$2, -16(%ebp)
	jmp	.L54
.L55:
	.loc 1 101 0
	movl	-16(%ebp), %eax
	subl	$1, %eax
	movl	%eax, -24(%ebp)
	jmp	.L56
.L57:
	.loc 1 102 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	fldl	-56(%ebp)
	faddp	%st, %st(1)
	fstpl	-56(%ebp)
	.loc 1 101 0
	addl	$1, -24(%ebp)
.L56:
	movl	-24(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L57
	.loc 1 100 0
	addl	$1, -16(%ebp)
.L54:
	movl	-16(%ebp), %eax
	cmpl	12(%ebp), %eax
	jle	.L55
	.loc 1 103 0
	movl	12(%ebp), %eax
	movl	%eax, -40(%ebp)
	.loc 1 104 0
	fldz
	fstpl	-96(%ebp)
	.loc 1 105 0
	jmp	.L60
.L61:
	.loc 1 106 0
	movl	$0, -20(%ebp)
.L62:
	.loc 1 108 0
	movl	-40(%ebp), %eax
	movl	%eax, -32(%ebp)
	jmp	.L63
.L64:
	.loc 1 109 0
	movl	8(%ebp), %edx
	subl	$4, %edx
	movl	-32(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-32(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	movl	-32(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-32(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	faddp	%st, %st(1)
	fstpl	-88(%ebp)
	.loc 1 110 0
	fldl	-88(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jp	.L140
	je	.L67
.L140:
	jmp	.L65
.L67:
	fldl	-56(%ebp)
	fstpl	-88(%ebp)
.L65:
	.loc 1 111 0
	movl	-32(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-32(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	faddl	-88(%ebp)
	fstpl	-48(%ebp)
	.loc 1 112 0
	fldl	-48(%ebp)
	fldl	-88(%ebp)
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jp	.L139
	je	.L68
.L139:
	.loc 1 108 0
	subl	$1, -32(%ebp)
.L63:
	cmpl	$1, -32(%ebp)
	jg	.L64
.L68:
	.loc 1 114 0
	movl	-40(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-128(%ebp)
	.loc 1 115 0
	movl	-32(%ebp), %eax
	cmpl	-40(%ebp), %eax
	jne	.L70
	.loc 1 116 0
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	-128(%ebp)
	faddl	-96(%ebp)
	fstpl	(%eax)
	.loc 1 117 0
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	20(%ebp), %eax
	fldz
	fstpl	(%eax)
	subl	$1, -40(%ebp)
	jmp	.L72
.L70:
	.loc 1 119 0
	movl	8(%ebp), %edx
	subl	$4, %edx
	movl	-40(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-136(%ebp)
	.loc 1 120 0
	movl	-40(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	movl	8(%ebp), %edx
	subl	$4, %edx
	movl	-40(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmulp	%st, %st(1)
	fstpl	-120(%ebp)
	.loc 1 121 0
	movl	-40(%ebp), %eax
	subl	$1, %eax
	cmpl	-32(%ebp), %eax
	jne	.L73
	.loc 1 122 0
	fldl	-136(%ebp)
	fsubl	-128(%ebp)
	fldl	.LC11
	fmulp	%st, %st(1)
	fstpl	-64(%ebp)
	.loc 1 123 0
	fldl	-64(%ebp)
	fmull	-64(%ebp)
	faddl	-120(%ebp)
	fstpl	-72(%ebp)
	.loc 1 124 0
	fldl	-72(%ebp)
	fstpl	-200(%ebp)
	movl	-200(%ebp), %edx
	movl	-196(%ebp), %ecx
	andl	$2147483647, %ecx
	movl	%edx, -176(%ebp)
	movl	%ecx, -172(%ebp)
	fldl	-176(%ebp)
	fsqrt
	fstpl	-184(%ebp)
	fldl	-184(%ebp)
	fucomp	%st(0)
	fnstsw	%ax
	sahf
	jp	.L138
	je	.L75
.L138:
	fldl	-176(%ebp)
	fstpl	(%esp)
	call	sqrt
	fstpl	-184(%ebp)
.L75:
	fldl	-184(%ebp)
	fstpl	-144(%ebp)
	.loc 1 125 0
	fldl	-128(%ebp)
	faddl	-96(%ebp)
	fstpl	-128(%ebp)
	.loc 1 126 0
	fldl	-72(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jae	.L78
	jmp	.L76
.L78:
	.loc 1 127 0
	fldl	-64(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jae	.L81
	jmp	.L79
.L81:
	fldl	-144(%ebp)
	fstpl	-200(%ebp)
	movl	-200(%ebp), %edx
	movl	-196(%ebp), %ecx
	andl	$2147483647, %ecx
	movl	%edx, -168(%ebp)
	movl	%ecx, -164(%ebp)
	jmp	.L82
.L79:
	fldl	-144(%ebp)
	fstpl	-200(%ebp)
	movl	-200(%ebp), %eax
	movl	-196(%ebp), %edx
	andl	$2147483647, %edx
	movl	%eax, %ecx
	movl	%edx, %ebx
	xorl	$-2147483648, %ebx
	movl	%ecx, -168(%ebp)
	movl	%ebx, -164(%ebp)
.L82:
	fldl	-168(%ebp)
	faddl	-64(%ebp)
	fstpl	-144(%ebp)
	.loc 1 128 0
	movl	16(%ebp), %edx
	subl	$8, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	%eax, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	-128(%ebp)
	faddl	-144(%ebp)
	fstpl	(%eax)
	fldl	(%eax)
	fstpl	(%edx)
	.loc 1 129 0
	fldl	-144(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jne	.L85
	jp	.L85
	jmp	.L83
.L85:
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	-120(%ebp)
	fdivl	-144(%ebp)
	fldl	-128(%ebp)
	fsubp	%st, %st(1)
	fstpl	(%eax)
.L83:
	.loc 1 130 0
	movl	20(%ebp), %edx
	subl	$8, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	%eax, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	20(%ebp), %eax
	fldz
	fstpl	(%eax)
	fldl	(%eax)
	fstpl	(%edx)
	jmp	.L86
.L76:
	.loc 1 132 0
	movl	16(%ebp), %edx
	subl	$8, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	%eax, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	16(%ebp), %eax
	fldl	-128(%ebp)
	faddl	-64(%ebp)
	fstpl	(%eax)
	fldl	(%eax)
	fstpl	(%edx)
	.loc 1 133 0
	movl	20(%ebp), %edx
	subl	$8, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	%eax, %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	addl	20(%ebp), %eax
	fldl	-144(%ebp)
	fstpl	(%eax)
	fldl	(%eax)
	fchs
	fstpl	(%edx)
.L86:
	.loc 1 135 0
	subl	$2, -40(%ebp)
	jmp	.L72
.L73:
	.loc 1 137 0
	cmpl	$30, -20(%ebp)
	jne	.L87
	movl	$.LC12, (%esp)
	call	nrerror
.L87:
	.loc 1 138 0
	cmpl	$10, -20(%ebp)
	je	.L89
	cmpl	$20, -20(%ebp)
	jne	.L91
.L89:
	.loc 1 139 0
	fldl	-96(%ebp)
	faddl	-128(%ebp)
	fstpl	-96(%ebp)
	.loc 1 140 0
	movl	$1, -16(%ebp)
	jmp	.L92
.L93:
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-16(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-16(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fsubl	-128(%ebp)
	fstpl	(%ecx)
	addl	$1, -16(%ebp)
.L92:
	movl	-16(%ebp), %eax
	cmpl	-40(%ebp), %eax
	jle	.L93
	.loc 1 141 0
	movl	-40(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	movl	8(%ebp), %edx
	subl	$4, %edx
	movl	-40(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	-16(%eax), %edx
	movl	-40(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	faddp	%st, %st(1)
	fstpl	-88(%ebp)
	.loc 1 142 0
	fldl	-88(%ebp)
	fldl	.LC13
	fmulp	%st, %st(1)
	fstpl	-128(%ebp)
	fldl	-128(%ebp)
	fstpl	-136(%ebp)
	.loc 1 143 0
	fldl	-88(%ebp)
	fldl	.LC14
	fmulp	%st, %st(1)
	fmull	-88(%ebp)
	fstpl	-120(%ebp)
.L91:
	.loc 1 145 0
	addl	$1, -20(%ebp)
	.loc 1 146 0
	movl	-40(%ebp), %eax
	subl	$2, %eax
	movl	%eax, -36(%ebp)
	jmp	.L95
.L96:
	.loc 1 147 0
	movl	-36(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-144(%ebp)
	.loc 1 148 0
	fldl	-128(%ebp)
	fsubl	-144(%ebp)
	fstpl	-80(%ebp)
	.loc 1 149 0
	fldl	-136(%ebp)
	fsubl	-144(%ebp)
	fstpl	-88(%ebp)
	.loc 1 150 0
	fldl	-80(%ebp)
	fmull	-88(%ebp)
	fsubl	-120(%ebp)
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-36(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fdivrp	%st, %st(1)
	movl	-36(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	faddp	%st, %st(1)
	fstpl	-64(%ebp)
	.loc 1 151 0
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-36(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fsubl	-144(%ebp)
	fsubl	-80(%ebp)
	fsubl	-88(%ebp)
	fstpl	-72(%ebp)
	.loc 1 152 0
	movl	8(%ebp), %edx
	addl	$8, %edx
	movl	-36(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-80(%ebp)
	.loc 1 153 0
	fldl	-64(%ebp)
	fabs
	fldl	-72(%ebp)
	fabs
	faddp	%st, %st(1)
	fldl	-80(%ebp)
	fabs
	faddp	%st, %st(1)
	fstpl	-88(%ebp)
	.loc 1 154 0
	fldl	-64(%ebp)
	fdivl	-88(%ebp)
	fstpl	-64(%ebp)
	.loc 1 155 0
	fldl	-72(%ebp)
	fdivl	-88(%ebp)
	fstpl	-72(%ebp)
	.loc 1 156 0
	fldl	-80(%ebp)
	fdivl	-88(%ebp)
	fstpl	-80(%ebp)
	.loc 1 157 0
	movl	-36(%ebp), %eax
	cmpl	-32(%ebp), %eax
	je	.L97
	.loc 1 158 0
	movl	-36(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	fldl	-72(%ebp)
	fabs
	fldl	-80(%ebp)
	fabs
	faddp	%st, %st(1)
	fmulp	%st, %st(1)
	fstpl	-104(%ebp)
	.loc 1 159 0
	fldl	-64(%ebp)
	fabs
	movl	8(%ebp), %edx
	subl	$4, %edx
	movl	-36(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	fldl	-144(%ebp)
	fabs
	faddp	%st, %st(1)
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-36(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-36(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fabs
	faddp	%st, %st(1)
	fmulp	%st, %st(1)
	fstpl	-112(%ebp)
	.loc 1 160 0
	fldl	-104(%ebp)
	faddl	-112(%ebp)
	fstpl	-48(%ebp)
	.loc 1 161 0
	fldl	-48(%ebp)
	fldl	-112(%ebp)
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jp	.L137
	je	.L97
.L137:
	.loc 1 146 0
	subl	$1, -36(%ebp)
.L95:
	movl	-36(%ebp), %eax
	cmpl	-32(%ebp), %eax
	jge	.L96
.L97:
	.loc 1 163 0
	movl	-36(%ebp), %eax
	addl	$2, %eax
	movl	%eax, -16(%ebp)
	jmp	.L100
.L101:
	.loc 1 164 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-16(%eax), %edx
	movl	-16(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldz
	fstpl	(%eax)
	.loc 1 165 0
	movl	-36(%ebp), %eax
	addl	$2, %eax
	cmpl	-16(%ebp), %eax
	je	.L102
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-24(%eax), %edx
	movl	-16(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldz
	fstpl	(%eax)
.L102:
	.loc 1 163 0
	addl	$1, -16(%ebp)
.L100:
	movl	-16(%ebp), %eax
	cmpl	-40(%ebp), %eax
	jle	.L101
	.loc 1 167 0
	movl	-36(%ebp), %eax
	movl	%eax, -28(%ebp)
	jmp	.L105
.L106:
	.loc 1 168 0
	movl	-28(%ebp), %eax
	cmpl	-36(%ebp), %eax
	je	.L107
	.loc 1 169 0
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-64(%ebp)
	.loc 1 170 0
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-72(%ebp)
	.loc 1 171 0
	fldz
	fstpl	-80(%ebp)
	.loc 1 172 0
	movl	-40(%ebp), %eax
	subl	$1, %eax
	cmpl	-28(%ebp), %eax
	je	.L109
	movl	8(%ebp), %edx
	addl	$8, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fstpl	-80(%ebp)
.L109:
	.loc 1 173 0
	fldl	-64(%ebp)
	fabs
	fldl	-72(%ebp)
	fabs
	faddp	%st, %st(1)
	fldl	-80(%ebp)
	fabs
	faddp	%st, %st(1)
	fstpl	-128(%ebp)
	fldl	-128(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jne	.L112
	jp	.L112
	jmp	.L107
.L112:
	.loc 1 174 0
	fldl	-64(%ebp)
	fdivl	-128(%ebp)
	fstpl	-64(%ebp)
	.loc 1 175 0
	fldl	-72(%ebp)
	fdivl	-128(%ebp)
	fstpl	-72(%ebp)
	.loc 1 176 0
	fldl	-80(%ebp)
	fdivl	-128(%ebp)
	fstpl	-80(%ebp)
.L107:
	.loc 1 179 0
	fldl	-64(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jae	.L115
	jmp	.L113
.L115:
	fldl	-64(%ebp)
	fmull	-64(%ebp)
	fldl	-72(%ebp)
	fmull	-72(%ebp)
	faddp	%st, %st(1)
	fldl	-80(%ebp)
	fmull	-80(%ebp)
	faddp	%st, %st(1)
	fstpl	(%esp)
	call	sqrt
	fstpl	-160(%ebp)
	jmp	.L116
.L113:
	fldl	-64(%ebp)
	fmull	-64(%ebp)
	fldl	-72(%ebp)
	fmull	-72(%ebp)
	faddp	%st, %st(1)
	fldl	-80(%ebp)
	fmull	-80(%ebp)
	faddp	%st, %st(1)
	fstpl	(%esp)
	call	sqrt
	fstpl	-200(%ebp)
	movl	-200(%ebp), %edx
	movl	-196(%ebp), %ecx
	xorl	$-2147483648, %ecx
	movl	%edx, -160(%ebp)
	movl	%ecx, -156(%ebp)
.L116:
	fldl	-160(%ebp)
	fstpl	-88(%ebp)
	fldl	-88(%ebp)
	fldz
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	jne	.L119
	jp	.L119
	jmp	.L117
.L119:
	.loc 1 180 0
	movl	-28(%ebp), %eax
	cmpl	-36(%ebp), %eax
	jne	.L120
	.loc 1 181 0
	movl	-32(%ebp), %eax
	cmpl	-36(%ebp), %eax
	je	.L124
	.loc 1 182 0
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fchs
	fstpl	(%ecx)
	jmp	.L124
.L120:
	.loc 1 184 0
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	-8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	-88(%ebp)
	fchs
	fmull	-128(%ebp)
	fstpl	(%eax)
.L124:
	.loc 1 185 0
	fldl	-64(%ebp)
	faddl	-88(%ebp)
	fstpl	-64(%ebp)
	.loc 1 186 0
	fldl	-64(%ebp)
	fdivl	-88(%ebp)
	fstpl	-128(%ebp)
	.loc 1 187 0
	fldl	-72(%ebp)
	fdivl	-88(%ebp)
	fstpl	-136(%ebp)
	.loc 1 188 0
	fldl	-80(%ebp)
	fdivl	-88(%ebp)
	fstpl	-144(%ebp)
	.loc 1 189 0
	fldl	-72(%ebp)
	fdivl	-64(%ebp)
	fstpl	-72(%ebp)
	.loc 1 190 0
	fldl	-80(%ebp)
	fdivl	-64(%ebp)
	fstpl	-80(%ebp)
	.loc 1 191 0
	movl	-28(%ebp), %eax
	movl	%eax, -24(%ebp)
	jmp	.L125
.L126:
	.loc 1 192 0
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-72(%ebp)
	faddp	%st, %st(1)
	fstpl	-64(%ebp)
	.loc 1 193 0
	movl	-40(%ebp), %eax
	subl	$1, %eax
	cmpl	-28(%ebp), %eax
	je	.L127
	.loc 1 194 0
	movl	8(%ebp), %edx
	addl	$8, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-80(%ebp)
	fldl	-64(%ebp)
	faddp	%st, %st(1)
	fstpl	-64(%ebp)
	.loc 1 195 0
	movl	8(%ebp), %edx
	addl	$8, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	8(%ebp), %edx
	addl	$8, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fldl	-64(%ebp)
	fmull	-144(%ebp)
	fsubrp	%st, %st(1)
	fstpl	(%ecx)
.L127:
	.loc 1 197 0
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	8(%ebp), %edx
	addl	$4, %edx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	leal	(%edx,%eax), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fldl	-64(%ebp)
	fmull	-136(%ebp)
	fsubrp	%st, %st(1)
	fstpl	(%ecx)
	.loc 1 198 0
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-28(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-24(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fldl	-64(%ebp)
	fmull	-128(%ebp)
	fsubrp	%st, %st(1)
	fstpl	(%ecx)
	.loc 1 191 0
	addl	$1, -24(%ebp)
.L125:
	movl	-24(%ebp), %eax
	cmpl	-40(%ebp), %eax
	jle	.L126
	.loc 1 200 0
	movl	-28(%ebp), %eax
	addl	$3, %eax
	movl	-40(%ebp), %edx
	movl	%edx, -192(%ebp)
	movl	%eax, -188(%ebp)
	movl	-192(%ebp), %ecx
	cmpl	%ecx, -188(%ebp)
	jle	.L130
	movl	-192(%ebp), %ebx
	movl	%ebx, -188(%ebp)
.L130:
	movl	-188(%ebp), %eax
	movl	%eax, -12(%ebp)
	.loc 1 201 0
	movl	-32(%ebp), %eax
	movl	%eax, -16(%ebp)
	jmp	.L131
.L132:
	.loc 1 202 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-128(%ebp)
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-136(%ebp)
	faddp	%st, %st(1)
	fstpl	-64(%ebp)
	.loc 1 203 0
	movl	-40(%ebp), %eax
	subl	$1, %eax
	cmpl	-28(%ebp), %eax
	je	.L133
	.loc 1 204 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	16(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fmull	-144(%ebp)
	fldl	-64(%ebp)
	faddp	%st, %st(1)
	fstpl	-64(%ebp)
	.loc 1 205 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	16(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	16(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fldl	-64(%ebp)
	fmull	-80(%ebp)
	fsubrp	%st, %st(1)
	fstpl	(%ecx)
.L133:
	.loc 1 207 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %eax
	leal	8(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fldl	-64(%ebp)
	fmull	-72(%ebp)
	fsubrp	%st, %st(1)
	fstpl	(%ecx)
	.loc 1 208 0
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %ecx
	movl	-16(%ebp), %eax
	sall	$2, %eax
	addl	8(%ebp), %eax
	movl	(%eax), %edx
	movl	-28(%ebp), %eax
	sall	$3, %eax
	leal	(%edx,%eax), %eax
	fldl	(%eax)
	fsubl	-64(%ebp)
	fstpl	(%ecx)
	.loc 1 201 0
	addl	$1, -16(%ebp)
.L131:
	movl	-16(%ebp), %eax
	cmpl	-12(%ebp), %eax
	jle	.L132
.L117:
	.loc 1 167 0
	addl	$1, -28(%ebp)
.L105:
	movl	-40(%ebp), %eax
	subl	$1, %eax
	cmpl	-28(%ebp), %eax
	jge	.L106
.L72:
	.loc 1 214 0
	movl	-40(%ebp), %eax
	subl	$1, %eax
	cmpl	-32(%ebp), %eax
	jg	.L62
.L60:
	.loc 1 105 0
	cmpl	$0, -40(%ebp)
	jg	.L61
	.loc 1 216 0
	fldl	-48(%ebp)
	.loc 1 217 0
	addl	$212, %esp
	popl	%ebx
	popl	%ebp
	ret
.LFE7:
	.size	hqr, .-hqr
	.local	sqrarg
	.comm	sqrarg,4,4
	.local	dsqrarg
	.comm	dsqrarg,8,8
	.local	dmaxarg1
	.comm	dmaxarg1,8,8
	.local	dmaxarg2
	.comm	dmaxarg2,8,8
	.local	dminarg1
	.comm	dminarg1,8,8
	.local	dminarg2
	.comm	dminarg2,8,8
	.local	maxarg1
	.comm	maxarg1,4,4
	.local	maxarg2
	.comm	maxarg2,4,4
	.local	minarg1
	.comm	minarg1,4,4
	.local	minarg2
	.comm	minarg2,4,4
	.local	lmaxarg1
	.comm	lmaxarg1,4,4
	.local	lmaxarg2
	.comm	lmaxarg2,4,4
	.local	lminarg1
	.comm	lminarg1,4,4
	.local	lminarg2
	.comm	lminarg2,4,4
	.local	imaxarg1
	.comm	imaxarg1,4,4
	.local	imaxarg2
	.comm	imaxarg2,4,4
	.local	iminarg1
	.comm	iminarg1,4,4
	.local	iminarg2
	.comm	iminarg2,4,4
	.section	.debug_frame,"",@progbits
.Lframe0:
	.long	.LECIE0-.LSCIE0
.LSCIE0:
	.long	0xffffffff
	.byte	0x1
	.string	""
	.uleb128 0x1
	.sleb128 -4
	.byte	0x8
	.byte	0xc
	.uleb128 0x4
	.uleb128 0x4
	.byte	0x88
	.uleb128 0x1
	.align 4
.LECIE0:
.LSFDE0:
	.long	.LEFDE0-.LASFDE0
.LASFDE0:
	.long	.Lframe0
	.long	.LFB5
	.long	.LFE5-.LFB5
	.byte	0x4
	.long	.LCFI0-.LFB5
	.byte	0xe
	.uleb128 0x8
	.byte	0x85
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI1-.LCFI0
	.byte	0xd
	.uleb128 0x5
	.align 4
.LEFDE0:
.LSFDE2:
	.long	.LEFDE2-.LASFDE2
.LASFDE2:
	.long	.Lframe0
	.long	.LFB6
	.long	.LFE6-.LFB6
	.byte	0x4
	.long	.LCFI3-.LFB6
	.byte	0xe
	.uleb128 0x8
	.byte	0x85
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI4-.LCFI3
	.byte	0xd
	.uleb128 0x5
	.align 4
.LEFDE2:
.LSFDE4:
	.long	.LEFDE4-.LASFDE4
.LASFDE4:
	.long	.Lframe0
	.long	.LFB7
	.long	.LFE7-.LFB7
	.byte	0x4
	.long	.LCFI6-.LFB7
	.byte	0xe
	.uleb128 0x8
	.byte	0x85
	.uleb128 0x2
	.byte	0x4
	.long	.LCFI7-.LCFI6
	.byte	0xd
	.uleb128 0x5
	.byte	0x4
	.long	.LCFI9-.LCFI7
	.byte	0x83
	.uleb128 0x3
	.align 4
.LEFDE4:
	.file 2 "nrutil.h"
	.text
.Letext0:
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
.LLST0:
	.long	.LFB5-.Ltext0
	.long	.LCFI0-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 4
	.long	.LCFI0-.Ltext0
	.long	.LCFI1-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 8
	.long	.LCFI1-.Ltext0
	.long	.LFE5-.Ltext0
	.value	0x2
	.byte	0x75
	.sleb128 8
	.long	0x0
	.long	0x0
.LLST1:
	.long	.LFB6-.Ltext0
	.long	.LCFI3-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 4
	.long	.LCFI3-.Ltext0
	.long	.LCFI4-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 8
	.long	.LCFI4-.Ltext0
	.long	.LFE6-.Ltext0
	.value	0x2
	.byte	0x75
	.sleb128 8
	.long	0x0
	.long	0x0
.LLST2:
	.long	.LFB7-.Ltext0
	.long	.LCFI6-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 4
	.long	.LCFI6-.Ltext0
	.long	.LCFI7-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 8
	.long	.LCFI7-.Ltext0
	.long	.LFE7-.Ltext0
	.value	0x2
	.byte	0x75
	.sleb128 8
	.long	0x0
	.long	0x0
	.section	.debug_info
	.long	0x55b
	.value	0x2
	.long	.Ldebug_abbrev0
	.byte	0x4
	.uleb128 0x1
	.long	.Ldebug_line0
	.long	.Letext0
	.long	.Ltext0
	.string	"GNU C 4.1.3 20070929 (prerelease) (Ubuntu 4.1.2-16ubuntu2)"
	.byte	0x1
	.string	"poly_solve.c"
	.string	"/home/user/src"
	.uleb128 0x2
	.long	.LASF0
	.byte	0x4
	.byte	0x7
	.uleb128 0x3
	.string	"unsigned char"
	.byte	0x1
	.byte	0x8
	.uleb128 0x3
	.string	"short unsigned int"
	.byte	0x2
	.byte	0x7
	.uleb128 0x3
	.string	"long unsigned int"
	.byte	0x4
	.byte	0x7
	.uleb128 0x3
	.string	"signed char"
	.byte	0x1
	.byte	0x6
	.uleb128 0x3
	.string	"short int"
	.byte	0x2
	.byte	0x5
	.uleb128 0x3
	.string	"int"
	.byte	0x4
	.byte	0x5
	.uleb128 0x3
	.string	"long long int"
	.byte	0x8
	.byte	0x5
	.uleb128 0x3
	.string	"long long unsigned int"
	.byte	0x8
	.byte	0x7
	.uleb128 0x3
	.string	"long int"
	.byte	0x4
	.byte	0x5
	.uleb128 0x2
	.long	.LASF0
	.byte	0x4
	.byte	0x7
	.uleb128 0x3
	.string	"char"
	.byte	0x1
	.byte	0x6
	.uleb128 0x3
	.string	"double"
	.byte	0x8
	.byte	0x4
	.uleb128 0x4
	.long	0x1b7
	.byte	0x1
	.string	"zrhqr"
	.byte	0x1
	.byte	0x11
	.byte	0x1
	.long	.LFB5
	.long	.LFE5
	.long	.LLST0
	.uleb128 0x5
	.string	"a"
	.byte	0x1
	.byte	0x10
	.long	0x1b7
	.byte	0x2
	.byte	0x91
	.sleb128 0
	.uleb128 0x5
	.string	"m"
	.byte	0x1
	.byte	0x10
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 4
	.uleb128 0x5
	.string	"rtr"
	.byte	0x1
	.byte	0x10
	.long	0x1b7
	.byte	0x2
	.byte	0x91
	.sleb128 8
	.uleb128 0x5
	.string	"rti"
	.byte	0x1
	.byte	0x10
	.long	0x1b7
	.byte	0x2
	.byte	0x91
	.sleb128 12
	.uleb128 0x6
	.string	"j"
	.byte	0x1
	.byte	0x15
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x6
	.string	"k"
	.byte	0x1
	.byte	0x15
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -16
	.uleb128 0x6
	.string	"hess"
	.byte	0x1
	.byte	0x16
	.long	0x1bd
	.byte	0x2
	.byte	0x91
	.sleb128 -12
	.uleb128 0x6
	.string	"xr"
	.byte	0x1
	.byte	0x16
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x6
	.string	"xi"
	.byte	0x1
	.byte	0x16
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.byte	0x0
	.uleb128 0x7
	.byte	0x4
	.long	0x11c
	.uleb128 0x7
	.byte	0x4
	.long	0x1b7
	.uleb128 0x4
	.long	0x26c
	.byte	0x1
	.string	"balanc"
	.byte	0x1
	.byte	0x32
	.byte	0x1
	.long	.LFB6
	.long	.LFE6
	.long	.LLST1
	.uleb128 0x5
	.string	"a"
	.byte	0x1
	.byte	0x31
	.long	0x1bd
	.byte	0x2
	.byte	0x91
	.sleb128 0
	.uleb128 0x5
	.string	"n"
	.byte	0x1
	.byte	0x31
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 4
	.uleb128 0x6
	.string	"last"
	.byte	0x1
	.byte	0x33
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x6
	.string	"j"
	.byte	0x1
	.byte	0x33
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -16
	.uleb128 0x6
	.string	"i"
	.byte	0x1
	.byte	0x33
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -12
	.uleb128 0x6
	.string	"s"
	.byte	0x1
	.byte	0x34
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x6
	.string	"r"
	.byte	0x1
	.byte	0x34
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x6
	.string	"g"
	.byte	0x1
	.byte	0x34
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.uleb128 0x6
	.string	"f"
	.byte	0x1
	.byte	0x34
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x6
	.string	"c"
	.byte	0x1
	.byte	0x34
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x6
	.string	"sqrdx"
	.byte	0x1
	.byte	0x34
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.byte	0x0
	.uleb128 0x8
	.long	0x3d0
	.byte	0x1
	.string	"hqr"
	.byte	0x1
	.byte	0x5c
	.byte	0x1
	.long	0x11c
	.long	.LFB7
	.long	.LFE7
	.long	.LLST2
	.uleb128 0x5
	.string	"a"
	.byte	0x1
	.byte	0x5b
	.long	0x1bd
	.byte	0x2
	.byte	0x91
	.sleb128 0
	.uleb128 0x5
	.string	"n"
	.byte	0x1
	.byte	0x5b
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 4
	.uleb128 0x5
	.string	"wr"
	.byte	0x1
	.byte	0x5b
	.long	0x1b7
	.byte	0x2
	.byte	0x91
	.sleb128 8
	.uleb128 0x5
	.string	"wi"
	.byte	0x1
	.byte	0x5b
	.long	0x1b7
	.byte	0x2
	.byte	0x91
	.sleb128 12
	.uleb128 0x6
	.string	"nn"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -48
	.uleb128 0x6
	.string	"m"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -44
	.uleb128 0x6
	.string	"l"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -40
	.uleb128 0x6
	.string	"k"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -36
	.uleb128 0x6
	.string	"j"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -32
	.uleb128 0x6
	.string	"its"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -28
	.uleb128 0x6
	.string	"i"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -24
	.uleb128 0x6
	.string	"mmin"
	.byte	0x1
	.byte	0x5d
	.long	0xcf
	.byte	0x2
	.byte	0x91
	.sleb128 -20
	.uleb128 0x6
	.string	"z"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -152
	.uleb128 0x6
	.string	"y"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -144
	.uleb128 0x6
	.string	"x"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -136
	.uleb128 0x6
	.string	"w"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -128
	.uleb128 0x6
	.string	"v"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -120
	.uleb128 0x6
	.string	"u"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -112
	.uleb128 0x6
	.string	"t"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -104
	.uleb128 0x6
	.string	"s"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -96
	.uleb128 0x6
	.string	"r"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -88
	.uleb128 0x6
	.string	"q"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x6
	.string	"p"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x3
	.byte	0x91
	.sleb128 -72
	.uleb128 0x6
	.string	"anorm"
	.byte	0x1
	.byte	0x5e
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -64
	.uleb128 0x6
	.string	"temp"
	.byte	0x1
	.byte	0x61
	.long	0x11c
	.byte	0x2
	.byte	0x91
	.sleb128 -56
	.byte	0x0
	.uleb128 0x6
	.string	"sqrarg"
	.byte	0x2
	.byte	0x4
	.long	0x3e4
	.byte	0x5
	.byte	0x3
	.long	sqrarg
	.uleb128 0x3
	.string	"float"
	.byte	0x4
	.byte	0x4
	.uleb128 0x6
	.string	"dsqrarg"
	.byte	0x2
	.byte	0x7
	.long	0x11c
	.byte	0x5
	.byte	0x3
	.long	dsqrarg
	.uleb128 0x6
	.string	"dmaxarg1"
	.byte	0x2
	.byte	0xa
	.long	0x11c
	.byte	0x5
	.byte	0x3
	.long	dmaxarg1
	.uleb128 0x6
	.string	"dmaxarg2"
	.byte	0x2
	.byte	0xa
	.long	0x11c
	.byte	0x5
	.byte	0x3
	.long	dmaxarg2
	.uleb128 0x6
	.string	"dminarg1"
	.byte	0x2
	.byte	0xe
	.long	0x11c
	.byte	0x5
	.byte	0x3
	.long	dminarg1
	.uleb128 0x6
	.string	"dminarg2"
	.byte	0x2
	.byte	0xe
	.long	0x11c
	.byte	0x5
	.byte	0x3
	.long	dminarg2
	.uleb128 0x6
	.string	"maxarg1"
	.byte	0x2
	.byte	0x12
	.long	0x3e4
	.byte	0x5
	.byte	0x3
	.long	maxarg1
	.uleb128 0x6
	.string	"maxarg2"
	.byte	0x2
	.byte	0x12
	.long	0x3e4
	.byte	0x5
	.byte	0x3
	.long	maxarg2
	.uleb128 0x6
	.string	"minarg1"
	.byte	0x2
	.byte	0x16
	.long	0x3e4
	.byte	0x5
	.byte	0x3
	.long	minarg1
	.uleb128 0x6
	.string	"minarg2"
	.byte	0x2
	.byte	0x16
	.long	0x3e4
	.byte	0x5
	.byte	0x3
	.long	minarg2
	.uleb128 0x6
	.string	"lmaxarg1"
	.byte	0x2
	.byte	0x1a
	.long	0x101
	.byte	0x5
	.byte	0x3
	.long	lmaxarg1
	.uleb128 0x6
	.string	"lmaxarg2"
	.byte	0x2
	.byte	0x1a
	.long	0x101
	.byte	0x5
	.byte	0x3
	.long	lmaxarg2
	.uleb128 0x6
	.string	"lminarg1"
	.byte	0x2
	.byte	0x1e
	.long	0x101
	.byte	0x5
	.byte	0x3
	.long	lminarg1
	.uleb128 0x6
	.string	"lminarg2"
	.byte	0x2
	.byte	0x1e
	.long	0x101
	.byte	0x5
	.byte	0x3
	.long	lminarg2
	.uleb128 0x6
	.string	"imaxarg1"
	.byte	0x2
	.byte	0x22
	.long	0xcf
	.byte	0x5
	.byte	0x3
	.long	imaxarg1
	.uleb128 0x6
	.string	"imaxarg2"
	.byte	0x2
	.byte	0x22
	.long	0xcf
	.byte	0x5
	.byte	0x3
	.long	imaxarg2
	.uleb128 0x6
	.string	"iminarg1"
	.byte	0x2
	.byte	0x26
	.long	0xcf
	.byte	0x5
	.byte	0x3
	.long	iminarg1
	.uleb128 0x6
	.string	"iminarg2"
	.byte	0x2
	.byte	0x26
	.long	0xcf
	.byte	0x5
	.byte	0x3
	.long	iminarg2
	.byte	0x0
	.section	.debug_abbrev
	.uleb128 0x1
	.uleb128 0x11
	.byte	0x1
	.uleb128 0x10
	.uleb128 0x6
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x25
	.uleb128 0x8
	.uleb128 0x13
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x1b
	.uleb128 0x8
	.byte	0x0
	.byte	0x0
	.uleb128 0x2
	.uleb128 0x24
	.byte	0x0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.byte	0x0
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x24
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.byte	0x0
	.byte	0x0
	.uleb128 0x4
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x1
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.uleb128 0x5
	.uleb128 0x5
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x6
	.uleb128 0x34
	.byte	0x0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0xa
	.byte	0x0
	.byte	0x0
	.uleb128 0x7
	.uleb128 0xf
	.byte	0x0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0x0
	.byte	0x0
	.uleb128 0x8
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x1
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0xc
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0xc
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x1
	.uleb128 0x40
	.uleb128 0x6
	.byte	0x0
	.byte	0x0
	.byte	0x0
	.section	.debug_pubnames,"",@progbits
	.long	0x2b
	.value	0x2
	.long	.Ldebug_info0
	.long	0x55f
	.long	0x126
	.string	"zrhqr"
	.long	0x1c3
	.string	"balanc"
	.long	0x26c
	.string	"hqr"
	.long	0x0
	.section	.debug_aranges,"",@progbits
	.long	0x1c
	.value	0x2
	.long	.Ldebug_info0
	.byte	0x4
	.byte	0x0
	.value	0x0
	.value	0x0
	.long	.Ltext0
	.long	.Letext0-.Ltext0
	.long	0x0
	.long	0x0
	.section	.debug_str,"",@progbits
.LASF0:
	.string	"unsigned int"
	.ident	"GCC: (GNU) 4.1.3 20070929 (prerelease) (Ubuntu 4.1.2-16ubuntu2)"
	.section	.note.GNU-stack,"",@progbits
