#include <stdio.h>  
#include <stdlib.h>
#include <unistd.h>  
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>

#define MAX(a, b) ( (a > b) ? a : b )
#define MIN(a, b) ( (a < b) ? a : b )
#define MAX_SEQ_LEN 8192

int seq[MAX_SEQ_LEN];
float v[MAX_SEQ_LEN+1][MAX_SEQ_LEN+1];
int gap_p;
int mismatch_p;
int match_score;

float getMax3(float a, float b, float c)
{
    return ( MAX( MAX(a,b),c));
}

float getMax4(float a, float b, float c, float d)
{
    return ( MAX( MAX(a,b),MAX(c,d)));
}

char decode( int a) 
{
    switch(a)
    {
        case 0: return 'A'; 
        case 1: return 'C'; 
        case 2: return 'G'; 
        case 3: return 'T'; 
        case -1: return '-';
    }
}

int encode( char c) 
{
    switch(c)
    {
        case 'A': return 0; 
        case 'C': return 1; 
        case 'G': return 2; 
        case 'T': return 3; 
        case '-': return -1; 
    }
}

int getScore(int x, int y)
{
    if(x == y)
    {
        return match_score;
    }
    else
    {
        if( x == -1 || y == -1) return gap_p;
        else return mismatch_p;
    }
}

// score for aligning x with column j
float getS(int** aln_matrix , float** profile, int x, int j)
{
    float sum = 0;

    for( int i = 0; i < 5; i++)
    {
        sum += profile[i][j-1] * getScore(x, i);
        
    }
    return sum;
}

void generateOutput(int seq1len, int* seq1algn, char* output_fn, char* seq_fn)
{
    FILE* output_file;
    FILE* seq_file;
    char c;
    seq_file = fopen(seq_fn, "r");
    output_file = fopen(output_fn, "w");
    
    if (!seq_file || !output_file) 
    {
        printf("ERROR: File could not open.");
        return 1;
    }

    while( (c = fgetc(seq_file)) != EOF ) fputc(c, output_file);
    fprintf(output_file, "\nsequence: ");
    for(int i = 0; i < seq1len; i++)
    {
        if( seq1algn[i] == -1) fprintf(output_file, "-");
        else fprintf(output_file, "%c",  decode(seq1algn[i]));
    }
    fclose(output_file);
    fclose(seq_file);
}

void naiveAlignmentMethod( int** aln_matrix, float** profile, int t, int n, int seq_len, char* output_fn, char* seq_fn)
{
    float del_case, mmatch_case, ins_case;
    int traceback_table[seq_len][n]; // 1: diagonal, 2: vertical, 3: horizontal 
    memset(traceback_table, 0, sizeof traceback_table);
    memset(v, 0, sizeof v);

    // init matrix's first row and col with gap penalties
    // base cases
    for( int i = 0; i <= seq_len; i++)
    {
        v[i][0] = i * gap_p;
    }
    for( int i = 0; i <= n; i++)
    {
        v[0][i] = i * gap_p;
    }
    // calculate maximum score using DP
    for( int i = 1; i <= seq_len; i++)
    {
        for( int j = 1; j <= n; j++)
        {
            ins_case = v[i][j-1] + getS(aln_matrix, profile, -1, j);
            mmatch_case = v[i-1][j-1] + getS(aln_matrix, profile, seq[i-1], j);
            del_case = v[i-1][j] + getScore(seq[i-1], -1);
            v[i][j] = getMax3( mmatch_case, // mismatch
                               del_case, // deletion
                               ins_case); // insertion
                               
            if( v[i][j] == mmatch_case) traceback_table[i-1][j-1] = 1;
            else if( v[i][j] == del_case) traceback_table[i-1][j-1] = 2;
            else if( v[i][j] == ins_case) traceback_table[i-1][j-1] = 3;
        }
    }

    // TRACEBACK (does not use a traceback table)

    int seq1pos = seq_len-1;
    int seq2pos = n-1;
    int seq1algn[n];
    int x = n-1;

    
    /*for(int i = 0; i < seq_len; i++)
    {
        for(int j = 0; j<n; j++)
        {
            printf("%f ",v[i][j]);
        }
        printf("\n");
        
    }
    printf("---\n");
    for(int i = 0; i < seq_len; i++)
    {
        for(int j = 0; j<n; j++)
        {
            printf("%d ",traceback_table[i][j]);
        }
        printf("\n");   
    }
    printf("---\n");
    */

    while( !(seq1pos == -1 || seq2pos == -1))
    {
        // traceback diagonal
        if( traceback_table[seq1pos][seq2pos] == 1)
        {
            seq1algn[x] = seq[seq1pos];
            x -= 1;
            seq1pos -= 1;
            seq2pos -= 1;
        }
        // traceback vertical
        else if( traceback_table[seq1pos][seq2pos] == 2)
        {
            seq1algn[x] = seq[seq1pos];
            x -= 1;
            seq1pos -= 1;
        }
        // traceback horizontal
        else if( traceback_table[seq1pos][seq2pos] == 3)
        {
            seq1algn[x] = -1;
            x -= 1;
            seq2pos -= 1;
        }
    }
    // output alignment result
    generateOutput(n, seq1algn, output_fn, seq_fn);
}

int main(int argc, char **argv) {
    
    // read inputs
    int option;
	FILE* fp;
    char* output_fn;
	char* seq_fn;
	char* aln_fn;

    static struct option long_options[] =
    {
        {"fasta", required_argument, NULL, 'fasta'},
        {"aln", required_argument, NULL, 'aln'},
        {"out", required_argument, NULL, 'out'},
        {"gap", required_argument, NULL, 'gap'},
        {"mismatch", required_argument, NULL, 'mismatch'},
        {"match", required_argument, NULL, 't'},
        {NULL, 0, NULL, 0}
    };

	while(( option = getopt_long( argc, argv, "fasta:aln:out:gap:mismatch:t:", long_options, NULL)) != -1)
	{
		switch(option)
		{   
			case 'fasta':
                seq_fn = optarg;
            	break;
			case 'aln':
                aln_fn = optarg;
            	break;
            case 'out':
                output_fn = optarg;
                break;
            case 'gap':
                gap_p = atoi(optarg);
                break;
            case 'mismatch':
                mismatch_p = atoi(optarg);
                break;
            case 't':
                match_score = atoi(optarg);
                break;
		}
	}

	fp = fopen(seq_fn,"r");
    if( fp == NULL)
	{
		printf("Input file is NULL\n");
		return 0;
	}
    
    // read sequence
    char c;
    int seq_len = 0;

	int i = 0;
	while(c != EOF)
    {   
        while((c = getc(fp)) != '\n' && c != EOF)
        {
            if( c == '>')
            {
                while((c = getc(fp)) != '\n' && c != EOF);
            }
            else
            {
                if (c != '\n') 
                {
                    seq[seq_len] = encode(c);
                    seq_len++;
                }  
            }
        }
    }
    fclose(fp);
	fp = fopen(aln_fn,"r");
    if( fp == NULL)
	{
		printf("Input file is NULL\n");
		return 0;
	}

 	int t = 0; // number of DNA sequences
 	int n = 0; // length of each DNA sequence
    char nucleotide;
 	int j;

    // get N value
	while((nucleotide = getc(fp)) != '\n')
    {
        n++;
        if( nucleotide == ' ') 
        {
            n = 0;
        }
    }
    //printf("n: %d\n", n);

    // get T value
    while(nucleotide != EOF)
    {   
        while((nucleotide = getc(fp)) != '\n' && nucleotide != EOF);
        t++;
    }
    t++;
    //printf("t: %d\n", t); 

    fclose(fp);
	fp = fopen(aln_fn,"r");

    // read alignment file to a matrix, then form the profile
    int *aln_matrix[t];
    for (i=0; i<t; i++) 
        aln_matrix[i] = (int *)malloc(n * sizeof(int));

    for( i = 0; i < t; i++)
    {
        while((nucleotide = getc(fp)) != ' ' && nucleotide != EOF);
        for( j = 0; j < n; j++)
        {
            nucleotide = getc(fp);
            if(nucleotide != "\n")
            {
                aln_matrix[i][j] =  encode(nucleotide);
            }
        }
    }
/*
    for( i = 0; i < t; i++)
    {
        for( j = 0; j < n; j++)
        {
            printf("%i ",aln_matrix[i][j]);
        }
        printf("\n");
    }
    printf("-----------------\n");
*/
    // calculate frequencies
    float* profile[5];
    for (i=0; i<5; i++) 
        profile[i] = (float *)malloc(n * sizeof(float)); 

    int* counts = (int*)calloc(5, sizeof(int));

    for( j = 0; j < n; j++)
    {
        counts = (int*)calloc(5, sizeof(int));
        for( i = 0; i < t; i++)
        {
            if( aln_matrix[i][j] == 0) counts[0] += 1;
            else if( aln_matrix[i][j] == 1) counts[1] += 1;
            else if( aln_matrix[i][j] == 2) counts[2] += 1;
            else if( aln_matrix[i][j] == 3) counts[3] += 1;
            else if( aln_matrix[i][j] == -1) counts[4] += 1;
        }
        profile[0][j] = counts[0]/(float)t; // A
        profile[1][j] = counts[1]/(float)t; // C
        profile[2][j] = counts[2]/(float)t; // G
        profile[3][j] = counts[3]/(float)t; // T
        profile[4][j] = counts[4]/(float)t; // -
    }
    free(counts);

    naiveAlignmentMethod( aln_matrix, profile, t, n, seq_len, output_fn, aln_fn);

    // free aln_matrix and profile
    for(i=0;i<t;i++)
    {
        free(aln_matrix[i]);
    }
    for(i=0;i<5;i++)
    {
        free(profile[i]);
    }

    printf("\nfile generated.\n");
    return 0;
}
