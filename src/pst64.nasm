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
;     nasm -f elf64 pst64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

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
backbone_asm:

    push		ebp		; salva il Base Pointer
    mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
    push		ebx		; salva i registri da preservare
    push		esi
    push		edi

    movss   ymm0, [dist0]         ; dist0 = 1.46
    movss   ymm1, [dist1]         ; dist1 = 1.52
    movss   ymm2, [dist2]         ; dist2 = 1.33
    movss   ymm3, [angle_cnca]    ; angle_cnca = 2.124

    mov ebx, [ebp + 8]

    ; Alloca memoria per coords
    push    ebx
    imul    ebx, 3


    push    ebx
    push    4                   
    call    alloc_matrix
    add     esp, 8
    mov     esi, eax            

    push 4
    call alloc_matrix
    add esp, 4
    mov ecx, eax
    pop ebx

    ; Imposta i primi valori di coords
    xorps   ymm4, ymm4
    vmovapd  [esi], ymm4         ; coords[0] = 0, coords[1] = 0, coords[2] = 0, coords[3] = 0 (Padding)
	addps  ymm4, ymm0          ; xmm4 = dist0
    vmovapd  [esi + 32], ymm4    ; coords[4] = dist0, coords[5] = 0, coords[6] = 0, coords[7] = 0 (Padding)
	

	
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

    vmovapd  ymm5, [esi + edx*8 - 32]
    vmovapd  ymm6, [esi + edx*8 - 64]
    vsubpd    ymm5, ymm6
    vmovapd  [ecx], ymm5         ; tmp = coords[4*(idx-1)] - coords[4*(idx-2)]



    push ecx

    ; Ruota tmp con angle_cnca

    vextractf128 xmm1, ymm3, 0x0
    movd eax, xmm1
    push eax
    push ecx
    call rotation
    add esp, 8
    pop ecx


    push eax

    ; Crea nuovo vettore (tmp e dist2)
    vxorpd   ymm5, ymm5

    vextractf128 xmm7, ymm2, 0x0        
    movsd xmm7, xmm7                   
    vshufpd xmm7, xmm7, xmm7, 0x1      
    vblendpd ymm5, ymm5, xmm7, 0b0010  

    vmovapd  [ecx], xmm5

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
    vmovapd ymm5, [eax]
    vmovapd ymm4, [esi + edx*8 - 32]
    vaddpd ymm4, ymm5
    vmovapd [esi + edx*8], ymm4



    ; Calcola CA (tmp[0..3] = coords[4*(idx+1)] - coords[4*(idx)])
    vmovapd  ymm5, [esi + edx*8]
    vmovapd  ymm6, [esi + edx*8 - 32]
    vsubpd    ymm5, ymm6
    vmovapd  [ecx], ymm5


    ; Ruota tmp con phi[i]

    mov eax, [ebp+12]
    mov eax, [eax+edi*8]

    push ecx

    push    eax
    push    ecx     ; phi[i]
    call    rotation
    add     esp, 8


    pop ecx
    push eax

    ; Crea nuovo vettore (tmp e dist0)
    vxorpd ymm5, ymm5

    vextractf128 xmm7, ymm0, 0x0        
    movsd xmm7, xmm7                   
    vshufpd xmm7, xmm7, xmm7, 0x1      
    vblendpd ymm5, ymm5, xmm7, 0b0010  

    vmovapd  [ecx], ymm5


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

    vmovapd ymm5, [eax]
    vmovapd ymm4, [esi + edx*8]
    vaddpd ymm5, ymm4

    vmovapd  [esi + edx*8 + 32], ymm5 ; Assegna risultato in coords[idx+1]

.skip_c:

    ; Calcola C (tmp[0..3] = coords[4*(idx+2)] - coords[4*(idx+1)])
    vmovapd  ymm5, [esi + edx*8 + 32]
    vmovapd  ymm6, [esi + edx*8]
    vsubpd    ymm5, ymm6
    vmovapd  [ecx], ymm5

    ; Ruota tmp con psi[i]

    mov eax, [ebp+16]
    mov eax, [eax+edi*8]


    push ecx

    push    eax
    push    ecx       ; psi[i]
    call    rotation
    add     esp, 8


    pop ecx


    ; Crea nuovo vettore (tmp e dist0)
    vxorpd ymm5, ymm5

    vextractf128 xmm7, ymm1, 0x0        
    movsd xmm7, xmm7                   
    vshufpd xmm7, xmm7, xmm7, 0x1      
    vblendpd ymm5, ymm5, xmm7, 0b0010  

    vmovapd  [ecx], ymm5



    push ecx

    push eax
    push ecx
    call matrix_product
    add esp,8

    pop ecx

    mov     edx, edi
    imul    edx, 12

    vmovapd ymm5, [eax]
    vmovapd ymm4, [esi + edx*8 + 32]
    vaddpd ymm5, ymm4

    vmovapd  [esi + edx*8 + 64], ymm5 ; Assegna risultato in coords[idx+1]

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
align 32
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



    vmovapd ymm0, [esi+edx*8] ; coords[i]
    vmovapd ymm1, [esi+edi*8] ; coord[j]
    vsubpd ymm1, ymm0         ; xmm1 = coords[j] - coords[i]


    vmulpd ymm1, ymm1


    haddpd ymm1, ymm1         ; somma parziale
    haddpd ymm1, ymm1         ; somma totale


    sqrtsd ymm1, ymm1         ; sqrt(xmm1)
    vextractf128 xmm1, ymm1, 0   ; Estrae la parte inferiore di ymm1 (128 bit) in xmm1
    movd eax, xmm1               ; Sposta il valore di xmm1 (double) in eax

    push eax
    fld dword [esp]
    add esp, 4


    pop edi                   ; ripristina i registri da preservare
    pop esi
    pop ebx
    pop ebp

    ret


