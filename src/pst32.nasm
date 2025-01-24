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
	ten dd 10.0
    temp dd 0.0


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
extern alloc_vector
extern alloc_matrix
extern rotation
extern matrix_product
extern distance

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

global backbone_asm
align 16
backbone_asm:
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

    ; Inizializza le costanti
    movss   xmm0, [dist0]         ; dist0 = 1.46
    movss   xmm1, [dist1]         ; dist1 = 1.52
    movss   xmm2, [dist2]         ; dist2 = 1.33
    movss   xmm3, [angle_cnca]    ; angle_cnca = 2.124

    mov ebx, [ebp + 8]

    ; Alloca memoria per coords
    push    ebx
    imul    ebx, 3


    push    ebx
    push    4                   ; Dimensione di ogni elemento (4 float)
    call    alloc_matrix
    add     esp, 8
    mov     esi, eax            ; coords = alloc_matrix(N*3, 4)

    push 4
    call alloc_matrix
    add esp, 4
    mov ecx, eax
    pop ebx

    ; Imposta i primi valori di coords
    xorps   xmm4, xmm4
    movaps  [esi], xmm4         ; coords[0] = 0, coords[1] = 0, coords[2] = 0, coords[3] = 0 (Padding)
    addps  xmm4, xmm0          ; xmm4 = dist0
    movaps  [esi + 16], xmm4    ; coords[4] = dist0, coords[5] = 0, coords[6] = 0, coords[7] = 0 (Padding)

    ; Prepara il ciclo per i calcoli
    xor     edi, edi            ; i = 0

.loop:
    ; Calcola idx
    mov     edx, edi
    imul    edx, 12      ; idx = i * 12

    ; Salta il calcolo di N e CA per la prima iterazione
    cmp     edi, 0
    je      .skip_c

    ; Calcola N (tmp[0..3] = coords[4*(idx-1)] - coords[4*(idx-2)])

    movaps  xmm5, [esi + edx*4 - 16]
    movaps  xmm6, [esi + edx*4 - 32]
    subps    xmm5, xmm6
    movaps  [ecx], xmm5         ; tmp = coords[4*(idx-1)] - coords[4*(idx-2)]



    push ecx

    movd eax, xmm3 ; Ruota tmp con angle_cnca

    push    eax
    push    ecx
    call    rotation
    add     esp, 8

    pop ecx


    push eax

    ; Crea nuovo vettore (tmp e dist2)
    xorps   xmm5, xmm5
    insertps xmm5, xmm2, 0x10         ; dist2
    movaps  [ecx], xmm5

    pop eax
    push ecx



    push eax
    push ecx
    call matrix_product
    add esp, 8

    pop ecx
    push eax


    mov     edx, edi
    imul    edx, 12
    pop eax
    movaps xmm5, [eax]
    movaps xmm4, [esi + edx*4 - 16]
    addps xmm4, xmm5
    movaps [esi + edx*4], xmm4



    ; Calcola CA (tmp[0..3] = coords[4*(idx+1)] - coords[4*(idx)])
    movaps  xmm5, [esi + edx*4]
    movaps  xmm6, [esi + edx*4 - 16]
    subps    xmm5, xmm6
    movaps  [ecx], xmm5


    ; Ruota tmp con phi[i]

    mov eax, [ebp+12]
    mov eax, [eax+edi*4]

    push ecx

    push    eax
    push    ecx     ; phi[i]
    call    rotation
    add     esp, 8


    pop ecx
    push eax

    ; Crea nuovo vettore (tmp e dist0)
    xorps xmm5, xmm5
    insertps  xmm5, xmm0, 0x10          ; dist0
    movaps  [ecx], xmm5


    pop eax
    push ecx

    push eax
    push ecx
    call matrix_product
    add esp,8

    pop ecx

    push eax

    mov     edx, edi
    imul    edx, 12

    pop eax

    movaps xmm5, [eax]
    movaps xmm4, [esi + edx*4]
    addps xmm5, xmm4

    movaps  [esi + edx*4 + 16], xmm5 ; Assegna risultato in coords[idx+1]

.skip_c:

    ; Calcola C (tmp[0..3] = coords[4*(idx+2)] - coords[4*(idx+1)])
    movaps  xmm5, [esi + edx*4 + 16]
    movaps  xmm6, [esi + edx*4]
    subps    xmm5, xmm6
    movaps  [ecx], xmm5

    ; Ruota tmp con psi[i]

    mov eax, [ebp+16]
    mov eax, [eax+edi*4]


    push ecx

    push    eax
    push    ecx       ; psi[i]
    call    rotation
    add     esp, 8


    pop ecx


    ; Crea nuovo vettore (tmp e dist0)
    xorps xmm5, xmm5
    insertps  xmm5, xmm1, 0x10          ; dist1
    movaps  [ecx], xmm5



    push ecx

    push eax
    push ecx
    call matrix_product
    add esp,8

    pop ecx

    mov     edx, edi
    imul    edx, 12

    movaps xmm5, [eax]
    movaps xmm4, [esi + edx*4 + 16]
    addps xmm5, xmm4

    movaps  [esi + edx*4 + 32], xmm5 ; Assegna risultato in coords[idx+1]

    ; Incrementa i registri e continua il ciclo
    add     edi, 1
    cmp     edi, ebx
    jl      .loop

    ; Fine ciclo
    ; Dealloca memoria e ripristina i registri
    fremem ecx
    mov eax, esi
    pop     edi
    pop     esi
    pop     ebx
    mov     esp, ebp
    pop     ebp
    ret



global distance_asm
align 16
distance_asm:
    push		ebp		; salva il Base Pointer
    mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
    push		ebx		; salva i registri da preservare
    push		esi
    push		edi

    mov edx, [ebp + 8]        ; edx = i
    mov edi, [ebp + 12]       ; edi = j

    shl edx, 2                ; edx = i * 4
    shl edi, 2                ; edi = j * 4 

    mov esi, [ebp + 16]       ; esi = c_alpha_coords



    movaps xmm0, [esi+edx*4] ; coords[i]
    movaps xmm1, [esi+edi*4] ; coord[j]
    subps xmm1, xmm0         ; xmm1 = coords[j] - coords[i]


    mulps xmm1, xmm1


    haddps xmm1, xmm1         ; somma parziale
    haddps xmm1, xmm1         ; somma totale


    sqrtss xmm1, xmm1         ; sqrt(xmm1)
    movd eax, xmm1            ; eax = sqrt(xmm1)

    push eax
    fld dword [esp]
    add esp, 4


    pop edi                   ; ripristina i registri da preservare
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp

    ret













