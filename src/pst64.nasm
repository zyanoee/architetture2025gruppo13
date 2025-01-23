; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli, F. Fassetti, S. Nisticò
; 12/11/2024
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mr9ante il packaging tool del sistema
; operativo; per esempio, su Ubuntu, mr9ante i comandi:
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
;     nasm -f elf64 pst64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
	dist0 dq 1.46
	dist1 dq 1.52
	dist2 dq 1.33
	val dq 1.0
	angle_cnca dq 2.124

section .bss			; Sezione contenente dati non inizializzati

alignb 32
e		resq		1

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in    rax
; l'indirizzo del primo bytes del blocco allocato
; (funziona mr9ante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mr9ante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block
extern alloc_vector
extern alloc_matrix
extern rotation
extern matrix_product


%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------
global prova

msg	db 'e:',32,0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		; rdi = indirizzo della struct input
		
		; esempio: stampa input->e
       	; [RDI] input->seq; 			    // sequenza
		; [RDI + 8]  input->N;			    // lunghezza della sequenza
		; [RDI + 12] input->sd; 		    // tasso raffredamento
		; [RDI + 16] input->to;			    // temperatura
		; [RDI + 24] input->alpha;		    // tasso raffredamento
		; [RDI + 32] input->k; 		        // numero di features da estrarre
		; [RDI + 40] input->hydrophobicity;	// hydrophobicity
		; [RDI + 48] input->volume;		    // volume
		; [RDI + 56] input->charge;		    // charge
		; [RDI + 64] input->phi;		    // vettore angoli phi
		; [RDI + 72] input-> psi;		    // vettore angoli psi
		; [RDI + 80] input->e:			    // energy
		; [RDI + 88] input->display;
		; [RDI + 92] input->silent;

		VMOVSD		XMM0, [RDI+80]
		VMOVSD		[e], XMM0
		prints 		msg
		printsd		e
		prints 		nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		    ; ripristina il Base Pointer
		ret				    ; torna alla funzione C chiamante



global distance_asm
align 32
distance_asm:
    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    push rbx						; salva i registri generali

    mov r8, rdi       ; r8 = i
    mov r9, rsi       ; r9 = j

    shl r8, 2                ; r8 = i * 4
    shl r9, 2                ; r9 = j * 4

    mov r10, rdx      ; r10 = c_alpha_coords



    vmovapd ymm0, [r10+r8*8] ; coords[i]
    vmovapd ymm1, [r10+r9*8] ; coord[j]
    vsubpd ymm1, ymm0         ; xmm1 = coords[j] - coords[i]


    vmulpd ymm1, ymm1


    vhaddpd ymm1, ymm1, ymm1         ; somma parziale
    vextractf128 xmm2, ymm1, 1
    vextractf128 xmm3, ymm1,0
    addpd xmm2, xmm3
    


    sqrtpd xmm2, xmm2
    xorps xmm0, xmm0
    movq  xmm0 , xmm2

    pop     rbx			; ripristina i registri generali
    mov		rsp, rbp	; ripristina lo Stack Pointer
    pop		rbp		    ; ripristina il Base Pointer
    ret				    ; torna alla funzione C chiamante

    ret


