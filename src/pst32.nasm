; ---------------------------------------------------------
; Predizione struttura terziaria con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	e		resd		1

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global prova

input		equ		8

msg	db	'e:',32,0
nl	db	10,0



prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione
		
		; esempio: stampa input->e
		mov EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
        ; [EAX]      input->seq; 			// sequenza
		; [EAX + 4]  input->N; 			    // numero elementi sequenza
		; [EAX + 8]  input->sd;			    // seed
		; [EAX + 12] input->to;			    // temperatura
		; [EAX + 16] input->alpha;		    // tasso raffredamento
		; [EAX + 20] input->k; 		        // costante
		; [EAX + 24] input->hydrophobicity;	// hydrophobicity
		; [EAX + 28] input->volume;		    // volume
		; [EAX + 32] input->charge;		    // charge
		; [EAX + 36] input->phi;		    // vettore angoli phi
		; [EAX + 40] input->psi;		    // vettore angoli psi
		; [EAX + 44] input->e;			    // energy
		; [EAX + 48] input->dispaly;
		; [EAX + 52] input->silent;
		MOVSS XMM0, [EAX+44]
		MOVSS [e], XMM0 
		prints msg            
		printss e   
		prints nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

global matrix_product
matrix_product:
		;
		; sequenza di ingresso nella funzione
		;
		PUSH	EBP
		MOV		EBP, ESP
		PUSH 	EBX
		PUSH	ESI
		PUSH 	EDI
		;
		; lettura dei parametri dal record di attivazione
		;
		MOV		EDX, [EBP + n]
		;
		; corpo della funzione
		;
		XOR		EAX, EAX			; i=EAX=0
fori:	XOR		EBX, EBX			; j=EBX=0		
forj:  	MOV		EDI, EBX			; EDI=dim*j
		IMUL	EDI, EDX			; EDI=dim*j*n
		MOV		EDX, [EBP + C]		; XMM3=C
		ADD		EDX, EAX			; EDX = C + EAX
		MOVAPS	XMM0, [EDX + EDI]	; XMM0=C[i...i+p-1][j]
		MOV		EDX, [EBP + n]		; rimetto su EDX n
		XOR		ECX, ECX			; x=ECX=0
forx:	MOV		ESI, ECX			; ESI=dim*x*n
		IMUL	ESI, EDX			; ESI=dim*x*n
		MOV		EDX, [EBP + A]		; XMM4=A, carico su EDX A
		ADD		EDX, EAX			; EDX = A + EAX
		MOVAPS	XMM1, [EDX + ESI]	; XMM1=A[i...i+p-1][x]
		MOV		EDX, [EBP + B]		; XMM5 = B, EDX = B
		ADD		EDX, ECX			; EDX = B + ECX
		MOVSS	XMM2, [EDX + EDI]	; XMM2=B[x][j]
		MOV		EDX, [EBP + n]		; EDX = n
		SHUFPS	XMM2, XMM2, 0		
		MULPS	XMM1, XMM2			; XMM1=A[i...i+p-1][k]*B[k][j]
		ADDPS	XMM0, XMM1			; XMM0=XMM0+A[i...i+p-1][k]*B[k][j]
		ADD		ECX, dim			; x++, x=x+dim
		IMUL	EDX, EDX, dim		; EDX = n *dim
		CMP		ECX, EDX			; x<n?
		MOV		EDX, [EBP + n]		; EDX = n
		JL		forx
		MOV		EDX, [EBP + C]		; EDX = C
		ADD		EDX, EAX			; EDX = C + EAX
		MOVAPS	[EDX + EDI], XMM0	; C[i...i+p-1][j]=XMM0
		ADD		EBX, dim			; j++, j=j+dim
		MOV		EDX, [EBP + nn2]		; EDX = nn
		IMUL	EDX, EDX, dim		; EDX = nn * dim
		CMP		EBX, EDX			; j<nn?
		MOV		EDX, [EBP + n]		; EDX = n
		JL		forj
		ADD		EAX, dim*p			; i=i+p, i=i+dim*p
		IMUL	EDX, EDX, dim		; EDX = n*dim
		CMP		EAX, EDX			; i<n
		MOV		EDX, [EBP + n]		; EDX = n
		JL		fori
		;MOVD	EDX, XMM3			; EDX = C
 		;printps	EDX, n*nn/4
		;
		; sequenza di uscita della funzione [OPPURE stop]
		;
		pop		edi
		pop		esi
		pop		ebx
		pop		EAX
		pop		ECX
		POP		EDX
		mov		esp, ebp
		pop 	ebp
		ret