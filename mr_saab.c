/*

   1) The mismatchPenalty, indelPenalty and matchScore should be tested such that a best one is picked
   2) Check the greedy method to see how good it does
   3) Check the maximum spaning tree vcersion and see which one is working better

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h> 
#include <sys/param.h>

#define MISMATCH_PENALTY 2
#define INDEL_PENALTY 4
#define MATCH_SCORE 1

#define MAX_READ_LENGTH 250


int clusterId_Min, clusterId_Max;
const int scoreVariance=3;

int OEASeqMaxLen=10000;
int read_number;
int histrogram[1000];

struct oea_read
{
	char* read_name;
	char* sequence;
	int chromosome_number;
	int mapping_position;
	char orientation;
	float score;
	int cluster_id;
	int mark; //0 if it is not yet picked. 1 if it is picked (it should not be picked again);
};

float **graph_weights_matrix;

//CKOCKAN: I changed the cluster file formatting. Make the necessary changes to the other parts that are effected
int main( int argc, char** argv)
{
	/* Variables */
	FILE* oea_cluster_file;
	FILE* histogram;
	int cluster_id;
	int max_num_reads;
	int is_forward;
	char orientation;
	char* result_contig;
	struct oea_read* oea_reads;
	int i;
	    
	/* Check for the correct number of arguments */
    if( argc != 5)
    {
        fprintf( stderr, "Incorrect number of arguments entered.\n");
        fprintf( stderr, "Correct usage: ./mr_saab <histogram> <cluster_file> <max_num_reads> <orientation>\n");
        exit( 1);
    }

	/* Open the histogram file for reading in text mode. */
	/* The histogram file contains the insert sizes of the paired end reads. */
	histogram = safe_fopen( argv[1], "r");

	/* Open the cluster file which contains the OEA read sequences, mapping positions of the mapped ends, the strand, and the cluster id. */
	oea_cluster_file = safe_fopen( argv[2], "r");

	/* Parse the remaining command line arguments */
	max_num_reads = atoi( argv[3]);
	is_forward = atoi( argv[4]);
 
	/* Allocate memory for the OEA reads */
	oea_reads = ( struct oea_read*) malloc( max_num_reads * sizeof( struct oea_read));

	/* Process the histogram */
	get_len_hist( histogram, MAX_READ_LENGTH);
	
	/* Parse the OEA read information from the cluster file */
	read_oea_reads( oea_cluster_file, &oea_reads);

	/* Create the Overlap Graph */
	create_the_graph();

	/* Use the greedy Overlap Graph Traversal algorithm to obtain the result contig */
	result_contig = greedy_graph_traversal( graph_weights_matrix, read_number);

	// CKOCKAN: REMOVE AFTER TESTING
	fprintf( stderr, "Result OEA Contig: %s\n", result_contig);

	/* Free memory that was previously allocated for the graph */
	for( i = 0; i < read_number; i++)
	{
		free( graph_weights_matrix[i]);
	}

	free( graph_weights_matrix);
}
