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

    ; Carica i parametris
    mov     eax, [ebp + 8]   ; eax = N (int)
    mov     ebx, [ebp + 12]  ; ebx = phi (VECTOR)
    mov     ecx, [ebp + 16]  ; ecx = psi (VECTOR)
    mov     edx, [ebp + 20]  ; edx = coords (MATRIX)

    xorps xmm0, xmm0
    xorps xmm1, xmm1
    xorps xmm2, xmm2
    xorps xmm3, xmm3

    ; Inizializza le costanti
    movss   xmm0, [dist0]
    movss   xmm1, [dist1]
    movss   xmm2, [dist2]
    movss   xmm3, [angle_cnca]

    ; Inizializziamo degli indici
    xor edi, edi ;i=0
  
    ; Inizializza il primo N e il primo CA
    movaps xmm4, [edx+edi*4]
    xorps xmm4, xmm4
    movaps [edx+edi*4], xmm4
    add edi, 4
    movaps xmm4, [edx+edi*4]
    addps xmm4, xmm0
    movaps [edx+edi*4], xmm4

    push 4
    call alloc_vector
    add esp, 4
    mov esi, eax

    movaps xmm5, [esi]
    xorps xmm5, xmm5 ;XMM5 = TMP
    mov eax, [ebp + 8] ; EAX = N
    
    ; Loop principale
    xor     edi, edi  ; i = 0
.loop_start:
    test     edi, eax
    jge     .loop_end

    ; Se i > 0
    cmp    edi, edi
    jz      .next_iteration

    ; CALCOLA N
    movaps xmm5, [edx+edi*4 - 16]
    subps xmm5, [edx+edi*4 - 32]
    movaps [esi], xmm5
    
    push esi
    push dword [angle_cnca]
    call rotation
    add esp, 8
    
    xorps xmm5, xmm5
    insertps xmm5, xmm2, 0x10
    push esi
    push eax
    call matrix_product
    add esp, 8

    xorps xmm5, xmm5
    movaps xmm5, [edx+edi*4 - 16]
    xorps xmm6, xmm6
    movaps xmm6, [eax]
    addps xmm5, xmm6
    movaps [edx+edi*4], xmm5
    xorps xmm5, xmm5

    ; CALCOLA CA
    movaps xmm5, [edx+edi*4]
    subps xmm5, [edx+edi*4 - 16]
    movaps [esi], xmm5
    mov eax, [ebx + edi * 4]
    
    push esi
    push eax
    call rotation
    add esp, 8
    
    xorps xmm5, xmm5
    insertps xmm5, xmm0, 0x10
    push esi
    push eax
    call matrix_product
    add esp, 8

    xorps xmm5, xmm5
    movaps xmm5, [edx+edi*4]
    xorps xmm6, xmm6
    movaps xmm6, [eax]
    addps xmm5, xmm6
    movaps [edx+edi*4 + 16], xmm5

.next_iteration
    ; CALCOLA C 
    movaps xmm5, [edx+edi*4 + 16]
    subps xmm5, [edx+edi*4 ]
    movaps [esi], xmm5
    mov eax, [ecx + edi * 4]
    
    push esi
    push eax
    call rotation
    add esp, 8
    
    xorps xmm5, xmm5
    insertps xmm5, xmm1, 0x10
    push esi
    push eax
    call matrix_product
    add esp, 8

    xorps xmm5, xmm5
    movaps xmm5, [edx+edi*4 + 16]
    xorps xmm6, xmm6
    movaps xmm6, [eax]
    addps xmm5, xmm6
    movaps [edx+edi*4 + 32], xmm5

    inc edi
    mov eax, [ebp + 8]
    jmp .loop_start
.loop_end
    fremem esi

    pop	edi		; ripristina i registri da preservare
    pop	esi
    pop	ebx
    mov	esp, ebp	; ripristina lo Stack Pointer
    pop	ebp		; ripristina il Base Pointer
    ret			; torna alla funzione C chiamante