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


global backbone_asm
align 32
backbone_asm: ;int N=rdi, VECTOR phi=rsi, VECTOR psi=rdx

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    vxorpd ymm10, ymm10
    vxorpd ymm11, ymm11
    vxorpd ymm12, ymm12
    vxorpd ymm13, ymm13

    vmovq   xmm10, [dist0]         ; dist0 = 1.46
    vmovq   xmm11, [dist1]         ; dist1 = 1.52
    vmovq   xmm12, [dist2]         ; dist2 = 1.33
    vmovq   xmm13, [angle_cnca]    ; angle_cnca = 2.124

    vinsertf128 ymm10, ymm10, xmm10, 0
    vinsertf128 ymm11, ymm11, xmm11, 0
    vinsertf128 ymm12, ymm12, xmm12, 0
    vinsertf128 ymm13, ymm13, xmm13, 0

    mov rbx, rdi
    mov r11, rsi
    mov r12, rdx


    ; Alloca memoria per coords

    sub rsp, 8
    mov [rsp], rbx
    imul    rbx, 3


    mov rdi, rbx
    mov rsi,   4
    call    alloc_matrix
    mov     r10, rax

    mov rdi, 4
    call alloc_matrix
    mov rcx, rax

    mov rbx, [rsp]
    add rsp, 8


    ; Imposta i primi valori di coords
    vxorpd   ymm4, ymm4
    vmovapd  [r10], ymm4         ; coords[0] = 0, coords[1] = 0, coords[2] = 0, coords[3] = 0 (Padding)
	vaddpd  ymm4, ymm10          ; xmm4 = dist0
    vmovapd  [r10 + 32], ymm4    ; coords[4] = dist0, coords[5] = 0, coords[6] = 0, coords[7] = 0 (Padding)
	

	
    ; Prepara il ciclo per i calcoli
    xor     r9, r9            ; i = 0

.loop:
    ; Calcola idx
    mov     r8, r9
    imul    r8, 12      ; idx = i * 12

    ; Salta il calcolo di N e CA per la prima iterazione
    cmp     r9, 0
    je      .skip_c

    ; Calcola N (tmp[0..3] = coords[4*(idx-1)] - coords[4*(idx-2)])

    vmovapd  ymm5, [r10 + r8*8 - 32]
    vmovapd  ymm6, [r10 + r8*8 - 64]
    vsubpd    ymm5, ymm6
    vmovapd  [rcx], ymm5         ; tmp = coords[4*(idx-1)] - coords[4*(idx-2)]



    sub rsp, 8
    mov [rsp], rcx

    push r8
    push r9
    push r10
    push r11
    push r12

    ; Ruota tmp con angle_cnca

    vextractf128 xmm0, ymm13, 0x0
    mov rdi, rcx
    call rotation

    pop r12
    pop r11
    pop r10
    pop r9
    pop r8

    mov [rsp], rcx
    add rsp, 8

    sub rsp, 8
    mov [rsp],    rax

    ; Crea nuovo vettore (tmp e dist2)
    vxorpd   ymm5, ymm5

    vextractf128 xmm7, ymm12, 0x0
    vinsertf128 ymm5, ymm5, xmm7, 0x10

    vmovapd  [rcx], xmm5

    mov rax, [rsp]
    mov [rsp], rcx

    push r8
    push r9
    push r10
    push r11
    push r12

    mov rsi,    rax
    mov rdi, rcx
    call matrix_product

    pop r12
    pop r11
    pop r10
    pop r9
    pop r8

    mov rcx, [rsp]
    mov [rsp], rax


    mov     r8, r9
    imul    r8, 12


    mov rax, [rsp]
    add rsp, 8

    vmovapd ymm5, [rax]
    vmovapd ymm4, [r10 + r8*8 - 32]
    vaddpd ymm4, ymm5
    vmovapd [r10 + r8*8], ymm4



    ; Calcola CA (tmp[0..3] = coords[4*(idx+1)] - coords[4*(idx)])
    vmovapd  ymm5, [r10 + r8*8]
    vmovapd  ymm6, [r10 + r8*8 - 32]
    vsubpd    ymm5, ymm6
    vmovapd  [rcx], ymm5


    ; Ruota tmp con phi[i]

    mov    rax, [r11]
    mov    rax, [rax+r9*8]

    sub rsp, 8
    mov [rsp], rcx

    push r8
    push r9
    push r10
    push r11
    push r12


    mov rsi,    rax
    mov rdi,    rcx     ; phi[i]
    call    rotation

    pop r12
    pop r11
    pop r10
    pop r9
    pop r8



    mov rcx, [rsp]
    mov [rsp], rax

    ; Crea nuovo vettore (tmp e dist0)
    vxorpd ymm5, ymm5


    vextractf128 xmm7, ymm10, 0x0
    vinsertf128 ymm5, ymm5, xmm7, 0x10


    vmovapd  [rcx], ymm5


    mov rax, [rsp]
    mov [rsp], rcx

    push r8
    push r9
    push r10
    push r11
    push r12


    mov rsi, rax
    mov rdi, rcx
    call matrix_product

    pop r12
    pop r11
    pop r10
    pop r8

    mov rcx, [rsp]
    mov [rsp], rax

    mov     r8, r9
    imul    r8, 12

    mov rax, [rsp]
    add rsp, 8

    vmovapd ymm5, [   rax]
    vmovapd ymm4, [r10 + r8*8]
    vaddpd ymm5, ymm4

    vmovapd  [r10 + r8*8 + 32], ymm5 ; Assegna risultato in coords[idx+1]

.skip_c:

    ; Calcola C (tmp[0..3] = coords[4*(idx+2)] - coords[4*(idx+1)])
    vmovapd  ymm5, [r10 + r8*8 + 32]
    vmovapd  ymm6, [r10 + r8*8]
    vsubpd    ymm5, ymm6
    vmovapd  [rcx], ymm5

    ; Ruota tmp con psi[i]

    mov    rax, [r12]
    mov    rax, [rax+r9*8]


    push rcx
    push r8
    push r9
    push r10
    push r11
    push r12

    push       rax
    push    rcx       ; psi[i]
    call    rotation
    add     esp, 8

    pop r12
    pop r11
    pop r10
    pop r9
    pop r8
    pop rcx


    ; Crea nuovo vettore (tmp e dist0)
    vxorpd ymm5, ymm5

    vextractf128 xmm7, ymm11, 0x0
    vinsertf128 ymm5, ymm5, xmm7, 0x10

    vmovapd  [rcx], ymm5



    push rcx
    push r8
    push r9
    push r10
    push r11
    push r12

    push    rax
    push rcx
    call matrix_product
    add esp,8

    pop r12
    pop r11
    pop r10
    pop r9
    pop r8
    pop rcx

    mov     r8, r9
    imul    r8, 12

    vmovapd ymm5, [rax]
    vmovapd ymm4, [r10 + r8*8 + 32]
    vaddpd ymm5, ymm4

    vmovapd  [r10 + r8*8 + 64], ymm5 ; Assegna risultato in coords[idx+1]

    ; Incrementa i registri e continua il ciclo
    add     r9, 1
    cmp     r9, rbx
    jl      .loop

    ; Fine ciclo
    ; Dealloca memoria e ripristina i registri
    fremem rcx
    mov    rax, r10
    popaq				; ripristina i registri generali
    mov		rsp, rbp	; ripristina lo Stack Pointer
    pop		rbp		    ; ripristina il Base Pointer
    ret				    ; torna alla funzione C chiamante


global distance_asm
align 32
distance_asm:
    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

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
    vhaddpd ymm1, ymm1, ymm1         ; somma totale


    vsqrtpd ymm1, ymm1         ; sqrt(xmm1)
    vextractf128 xmm1, ymm1, 0   ; Estrae la parte inferiore di ymm1 (128 bit) in xmm1
    movq    rax, xmm1
    xorps xmm0, xmm0
    movq xmm0, rax

    popaq				; ripristina i registri generali
    mov		rsp, rbp	; ripristina lo Stack Pointer
    pop		rbp		    ; ripristina il Base Pointer
    ret				    ; torna alla funzione C chiamante

    ret


