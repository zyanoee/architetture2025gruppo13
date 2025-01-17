#include <stdlib.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)


//----------------------------------------------------------------------

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

//----------------------------------------------------------------------

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

int main(int argc, char** argv) {

    if (argc != 2) {
        printf("Usage: %s <filename.ds2>\n", argv[0]);
        return 1;
    }

    char* filename = argv[1];
    int n, k;
    char* data = load_seq(filename, &n, &k);

    // Create output filename
    char output_filename[256];
    strcpy(output_filename, filename);
    char* dot = strrchr(output_filename, '.');
    if (dot != NULL) {
        strcpy(dot, ".txt");
    } else {
        strcat(output_filename, ".txt");
    }

    // Write data to output file
    FILE* output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
        printf("Error opening output file '%s'\n", output_filename);
        free(data);
        return 1;
    }

    for (int i = 0; i < n; i++) {
       fprintf(output_file, "%f", data[i]);
       fprintf(output_file, "\n");
    }

    fclose(output_file);
    free(data);

    return 0;

}