/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-64 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 pst64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o pst64.o pst64c.c -o pst64c -lm && ./pst64c $pars
* 
* oppure
* 
* ./runpst64
* 
*/

#include <stdlib.h>

#include <stdint.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

typedef struct {
	char* seq;		// sequenza di amminoacidi
	int N;			// lunghezza sequenza
	unsigned int sd; 	// seed per la generazione casuale
	type to;		// temperatura INIZIALE
	type alpha;		// tasso di raffredamento
	type k;		// costante
	VECTOR hydrophobicity; // hydrophobicity
	VECTOR volume;		// volume
	VECTOR charge;		// charge
	VECTOR phi;		// vettore angoli phi
	VECTOR psi;		// vettore angoli psi
	type e;		// energy
	int display;
	int silent;

} params;


/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

VECTOR alloc_vector(int n) {
	return (VECTOR) get_block(sizeof(type),n);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

void dealloc_vector(void* mat) {
	free_block(mat);
}



/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
* 
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
* 
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY
extern void prova(params* input);
extern MATRIX backbone_asm(int N, VECTOR phi, VECTOR psi);
extern type distance_asm(int i, int j, MATRIX coords);

// PROCEDURE IN C
	//TRIGONOMETRICHE
type cosine(type theta){
	type t0 = 1;
	type t1 = -theta*theta / 2;
	type t2 = theta*theta*theta*theta / 24; //fact(4)
	type t3 = -theta*theta*theta*theta*theta*theta / 720;  //fact(6)

	return t0 + t1 + t2 + t3;
}
type sine(type theta){
	type t0 = theta;
	type t1 = -theta*theta*theta / 6; //fact(3)
	type t2 = theta*theta*theta*theta*theta / 120; //fact(5)
	type t3 = -theta*theta*theta*theta*theta*theta*theta / 5040; //fact(7)

	return t0 + t1 + t2 + t3;
}

//funzione per la normalizzazione dell'axis
void normalize_vector(VECTOR axis){
	type magn = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
	if(magn != 0){
		axis[0] /= magn;
		axis[1] /= magn;
		axis[2] /= magn;
	}
}
//ROTAZIONE
MATRIX rotation(VECTOR axis, type theta){
	MATRIX R = alloc_matrix(3,3);
	normalize_vector(axis);
	type scalar_prod = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
	axis[0] = axis[0]/scalar_prod;
	axis[1] = axis[1]/scalar_prod;
	axis[2] = axis[2]/scalar_prod;


	type a = cosine(theta/2);
	type stheta = sine(theta/2);

	type b = -1 * axis[0] * stheta;
	type c = -1 * axis[1] * stheta;
	type d = -1 * axis[2] * stheta;

	R[0] = a*a + b*b - c*c - d*d;
	R[1] = 2*b*c + 2*a*d;
	R[2] = 2*b*d - 2*a*c;
	
	R[3] = 2*b*c - 2*a*d;
	R[4] = a*a + c*c - b*b - d*d;
	R[5] = 2*c*d + 2*a*b;

	R[6] = 2*b*d + 2*a*c;
	R[7] = 2*c*d - 2*a*b;
	R[8] = a*a + d*d - b*b - c*c;
	return R;
}

VECTOR matrix_product(VECTOR v, MATRIX m){
	VECTOR res = alloc_vector(4);
	res[0] = v[0]*m[0] + v[1]*m[3] + v[2]*m[6];
	res[1] = v[0]*m[1] + v[1]*m[4] + v[2]*m[7];
	res[2] = v[0]*m[2] + v[1]*m[5] + v[2]*m[8];
	res[3] = 0;

	return res;
}

type distance_angles(type phi, type psi, type a_phi, type a_psi){
	return sqrt(pow(phi-a_phi, 2) + pow(psi-a_psi, 2));
}

type min(type a, type b){
	if(a<b)
		return a;
	return b;
}

MATRIX backbone(int N, VECTOR phi, VECTOR psi){

	MATRIX coords = alloc_matrix(N*3,4);


	type dist0=1.46; //rcan
	type dist1=1.52; //rcac
	type dist2=1.33; //rcn
	type angle_cnca = 2.124;

	//Primo N
	coords[0] = 0;
	coords[1] = 0;
	coords[2] = 0;
	coords[3] = 0; //PADDING PER SSE

	//Primo CA
	coords[4] = dist0;
	coords[5] = 0;
	coords[6] = 0;
	coords[7] = 0; //PADDING PER SSE

	VECTOR tmp = alloc_vector(4);
	for(int i = 0; i<N; i++){
		int idx = i*3;
		MATRIX R;
		VECTOR neww;
		if(i>0){

			//N
			tmp[0] = coords[4*(idx-1)] - coords[4*(idx-2)];
			tmp[1] = coords[4*(idx-1)+1] - coords[4*(idx-2)+1];
			tmp[2] = coords[4*(idx-1)+2] - coords[4*(idx-2)+2];
			tmp[3] = 0;
			R = rotation(tmp, angle_cnca);
			tmp[0] = 0;
			tmp[1]= dist2;
			tmp[2]= 0;
			tmp[3]= 0;
			neww = matrix_product(tmp, R);
			coords[idx*4]=coords[(idx-1)*4]+neww[0];
			coords[idx*4+1]=coords[(idx-1)*4+1]+neww[1];
			coords[idx*4+2]=coords[(idx-1)*4+2]+neww[2];
			coords[idx*4+3]=0; //PADDING PER SSE

			//CA
			tmp[0] = coords[idx*4] - coords[4*(idx-1)];
			tmp[1] = coords[idx*4+1] - coords[4*(idx-1)+1];
			tmp[2] = coords[idx*4+2] - coords[4*(idx-1)+2];
			tmp[3] = 0;
			R = rotation(tmp, phi[i]);
			tmp[0] = 0;
			tmp[1]= dist0;
			tmp[2]= 0;
			tmp[3]= 0;
			neww = matrix_product(tmp, R);
			coords[4*(idx+1)]=coords[4*idx]+neww[0];
			coords[4*(idx+1)+1]=coords[4*idx+1]+neww[1];
			coords[4*(idx+1)+2]=coords[4*idx+2]+neww[2];
			coords[4*(idx+1)+3]=0; //PADDING PER SSE
		}	

		//C
		tmp[0] = coords[4*(idx+1)] - coords[4*(idx)];
		tmp[1] = coords[4*(idx+1)+1] - coords[4*(idx)+1];
		tmp[2] = coords[4*(idx+1)+2] - coords[4*(idx)+2];
		tmp[3] = 0;
		R = rotation(tmp, psi[i]);
		tmp[0] = 0;
		tmp[1]= dist1;
		tmp[2]= 0;
		tmp[3]= 0;
		neww = matrix_product(tmp, R);
		coords[4*(idx+2)]=coords[4*(idx+1)]+neww[0];
		coords[4*(idx+2)+1]=coords[4*(idx+1)+1]+neww[1];
		coords[4*(idx+2)+2]=coords[4*(idx+1)+2]+neww[2];
		coords[4*(idx+2)+3]=0; 

	}
	dealloc_vector(tmp);
	return coords;
}

MATRIX c_alpha_coords(MATRIX coords, int N){
	MATRIX c_alpha_coords = alloc_matrix(N,4);
	for(int i=0; i<N; i++){
		int inx = i*12;
		int idx = i*4;
		c_alpha_coords[idx] = coords[inx+4];
		c_alpha_coords[idx+1] = coords[inx+5];
		c_alpha_coords[idx+2] = coords[inx+6];
		c_alpha_coords[idx+3] = coords[inx+7];
	}
	
	return c_alpha_coords;
}

type distance(int i, int j, MATRIX coords){
	int real_i = i*4;
	int real_j = j*4;
	return sqrt(pow(coords[real_j]-coords[real_i], 2) +
				pow(coords[real_j+1]-coords[real_i+1], 2) +
				pow(coords[real_j+2]-coords[real_i+2], 2) +
				pow(coords[real_j+3] -coords[real_i+3], 2)
	);
}

//ENERGIE
type rama_energy(VECTOR phi, VECTOR psi, int N){
	type alpha_phi = -57.8;
	type alpha_psi = -47.0;
	type beta_phi = -119.0;
	type beta_psi = 113.0;
	type energy = 0;
	for (int i = 0; i<N; i++){
		type a_dist = distance_angles(phi[i], psi[i], alpha_phi, alpha_psi);
		type b_dist = distance_angles(phi[i], psi[i], beta_phi, beta_psi);
		energy = energy + 0.5*min(a_dist, b_dist);
	}
	return energy;
}

type hydrophobic_energy(char* s, MATRIX c_alpha_coords, int N){
	type energy = 0;
	for(int i = 0; i<N; i++){
		for(int j = i+1; j<N; j++ ){
			type dist = (type)distance_asm(i,j, c_alpha_coords);
			if( dist < 10.0){
				energy = energy + (hydrophobicity[s[i]-'A']*hydrophobicity[s[j]-'A'])/dist;
			}
		}
	}
	return energy;
}

type electrostatic_energy(char* s, MATRIX c_alpha_coords, int N){
	type energy = 0;
	for(int i = 0; i<N; i++){
		for(int j = i+1; j<N; j++ ){
			type dist = (type)distance_asm(i,j, c_alpha_coords);
			if( dist < 10.0 && charge[s[i]-'A'] != 0 && charge[s[j]-'A'] != 0){
				energy = energy + (charge[s[i]-'A']*charge[s[j]-'A'])/(dist*4);
			}
		}
	}
	return energy;
}

type packing_energy(char* s, MATRIX c_alpha_coords, int N){
	type energy = 0;
	for(int i = 0; i<N; i++){
		type density = 0;
		for(int j = 0; j<N; j++ ){
			type dist = (type)distance_asm(i,j, c_alpha_coords);
			if(i!=j && dist < 10.0){
				type d = (volume[s[j]-'A']/(dist*dist*dist));
				density = density + d;
			}
		}
		type voldens = (volume[s[i]-'A']-density)*(volume[s[i]-'A']-density);
		energy = energy + voldens;
	}
	return energy;
}

type energy(char* s, int N, VECTOR phi, VECTOR psi){

	//MATRIX coords = backbone(N, phi, psi);
	MATRIX coords = (MATRIX) backbone_asm(N, phi, psi);
	MATRIX c_alpha = c_alpha_coords(coords, N);

	type rama = rama_energy(phi, psi, N);
	type hydrophobic = hydrophobic_energy(s, c_alpha, N) ;
	type electrostatic = electrostatic_energy(s, c_alpha, N);
	type packing = packing_energy(s, c_alpha, N);

	type w_rama = 1.0;
	type w_hydrophobic = 0.5;
	type w_electrostatic = 0.2;
	type w_packing = 0.3;

	dealloc_matrix(coords);
	dealloc_matrix(c_alpha);

	return w_rama*rama + w_hydrophobic*hydrophobic + w_electrostatic*electrostatic + w_packing*packing;
}

void simulated_annealing(char* s, int N, VECTOR phi, VECTOR psi, type T0, type alpha, type k){
	printf("Simulated Annealing start...");
	type E = energy(s,N,phi,psi);
	type T = T0;
	int t = 0;
	printf("SA: t > %i - ENERGIA> %f - DELTA-E > NA\n", t, E);

	while (T>0){
		int i = (int) (random()*N);
		type dphi = (random()*2 * M_PI) - M_PI;
		type dpsi = (random()*2 * M_PI) - M_PI;
		phi[i]+= dphi;
		psi[i]+= dpsi;


		type new_energy = energy(s,N,phi,psi);
		type deltaE = new_energy - E;
		printf("SA: t > %i - ENERGIA> %f - DELTA-E > %f \n", t, E, deltaE);
		if(deltaE<= 0){
			printf("SA: t>%i - Accettazione nuova configurazione per indice %i [dE <= 0] \n", t,i);
			E = new_energy;
		}else{

			type sp = -deltaE/(k*T);
			type P = exp(sp);
			type r = random();
			if(r<= P){
				printf("SA: t>%i - Accettazione nuova configurazione per indice %i [r <= P] \n", t,i);
				E = energy(s,N,phi,psi);
			}else{
				printf("SA: t>%i -Rifiuto nuova configurazione per indice %i \n",t,i);
				phi[i]-= dphi;
				psi[i]-= dpsi;
				
			}
			
		}
		printf("\n ------------------ \n");
		t++;
		T = T0 - sqrt(alpha*t);
	}
	
	printf("SA: Energia = %f\n", E);

	return;
}

void pst(params* input){
	simulated_annealing(input->seq, input->N, input->phi, input->psi, input->to,input->alpha, input->k);
}

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
	int d;
	
	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;	
	input->N = -1;			
	input->to = -1;
	input->alpha = -1;
	input->k = -1;		
	input->sd = -1;		
	input->phi = NULL;		
	input->psi = NULL;
	input->silent = 0;
	input->display = 0;
	input->e = -1;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);

	
	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);
	input->psi = alloc_matrix(input->N, 1);
	// Impostazione seed 
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	t = clock();
	pst(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out32_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}
