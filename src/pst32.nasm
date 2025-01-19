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
	dist0 dd 1.46
	dist1 dd 1.52
	dist2 dd 1.33
	angle_cnca dd 2.124

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
global backbone

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



backbone:
    ; Prologo della funzione
    push    ebp
    mov     ebp, esp
    sub     esp, 32  ; Alloca spazio per variabili locali

    ; Salva i registri da preservare
    push    ebx
    push    esi
    push    edi

    ; Carica i parametri
    mov     eax, [ebp + 8]   ; eax = s (char*)
    mov     ebx, [ebp + 12]  ; ebx = phi (VECTOR)
    mov     ecx, [ebp + 16]  ; ecx = psi (VECTOR)

    ; Calcola la lunghezza della stringa s
    mov     edi, eax
    xor     eax, eax
    cld
    repne scasb
    not     ecx
    dec     ecx
    mov     [esp], ecx  ; N = strlen(s)

    ; Alloca la matrice coords
    mov     eax, ecx
    imul    eax, 3
    push    eax
    push    3
    call    alloc_matrix
    add     esp, 8
    mov     esi, eax  ; esi = coords

    ; Inizializza le costanti
    movss   xmm0, dword [dist0]
    movss   xmm1, dword [dist1]
    movss   xmm2, dword [dist2]
    movss   xmm3, dword [angle_cnca]

    ; Inizializza il primo N
    mov     dword [esi], 0
    mov     dword [esi + 4], 0
    mov     dword [esi + 8], 0

    ; Alloca il vettore tmp
    push    3
    call    alloc_vector
    add     esp, 4
    mov     edi, eax  ; edi = tmp

    ; Loop principale
    xor     eax, eax  ; i = 0
.loop_start:
    cmp     eax, [esp]
    jge     .loop_end

    imul    edx, eax, 3
    mov     ebp, edx  ; idx = i * 3

    ; Se i > 0
    test    eax, eax
    jz      .next_iteration

    ; Calcola tmp per N
    mov     ecx, ebp
    sub     ecx, 3
    movss   xmm4, dword [esi + ecx * 4]
    movss   xmm5, dword [esi + ecx * 4 + 4]
    movss   xmm6, dword [esi + ecx * 4 + 8]
    sub     ecx, 3
    subss   xmm4, dword [esi + ecx * 4]
    subss   xmm5, dword [esi + ecx * 4 + 4]
    subss   xmm6, dword [esi + ecx * 4 + 8]
    movss   dword [edi], xmm4
    movss   dword [edi + 4], xmm5
    movss   dword [edi + 8], xmm6

    ; Chiama la funzione rotation
    push    dword [angle_cnca]
    push    edi
    call    rotation
    add     esp, 8
    mov     ebx, eax  ; ebx = R

    ; Calcola tmp per CA
    movss   xmm4, dword [dist2]
    movss   dword [edi + 4], xmm4
    movss   dword [edi], xmm0
    movss   dword [edi + 8], xmm0

    ; Chiama la funzione matrix_product
    push    ebx
    push    edi
    call    matrix_product
    add     esp, 8
    mov     ebx, eax  ; ebx = neww

    ; Aggiorna coords per N
    movss   xmm4, dword [esi + ecx * 4]
    movss   xmm5, dword [esi + ecx * 4 + 4]
    movss   xmm6, dword [esi + ecx * 4 + 8]
    addss   xmm4, dword [ebx]
    addss   xmm5, dword [ebx + 4]
    addss   xmm6, dword [ebx + 8]
    movss   dword [esi + ebp * 4], xmm4
    movss   dword [esi + ebp * 4 + 4], xmm5
    movss   dword [esi + ebp * 4 + 8], xmm6

    ; Calcola tmp per CA
    movss   xmm4, dword [esi + ebp * 4]
    movss   xmm5, dword [esi + ebp * 4 + 4]
    movss   xmm6, dword [esi + ebp * 4 + 8]
    subss   xmm4, dword [esi + ecx * 4]
    subss   xmm5, dword [esi + ecx * 4 + 4]
    subss   xmm6, dword [esi + ecx * 4 + 8]
    movss   dword [edi], xmm4
    movss   dword [edi + 4], xmm5
    movss   dword [edi + 8], xmm6

    ; Chiama la funzione rotation
    push    dword [phi + eax * 4]
    push    edi
    call    rotation
    add     esp, 8
    mov     ebx, eax  ; ebx = R

    ; Calcola tmp per CA
    movss   xmm4, dword [dist0]
    movss   dword [edi + 4], xmm4
    movss   dword [edi], xmm0
    movss   dword [edi + 8], xmm0

    ; Chiama la funzione matrix_product
    push    ebx
    push    edi
    call    matrix_product
    add     esp, 8
    mov     ebx, eax  ; ebx = neww

    ; Aggiorna coords per CA
    movss   xmm4, dword [esi + ebp * 4]
    movss   xmm5, dword [esi + ebp * 4 + 4]
    movss   xmm6, dword [esi + ebp * 4 + 8]
    addss   xmm4, dword [ebx]
    addss   xmm5, dword [ebx + 4]
    addss   xmm6, dword [ebx + 8]
    movss   dword [esi + (ebp + 3) * 4], xmm4
    movss   dword [esi + (ebp + 3) * 4 + 4], xmm5
    movss   dword [esi + (ebp + 3) * 4 + 8], xmm6

.next_iteration:
    inc     eax
    jmp     .loop_start

.loop_end:
    ; Ripristina i registri da preservare
    pop     edi
    pop     esi
    pop     ebx

    ; Epilogo della funzione
    mov     esp, ebp
    pop     ebp
    ret