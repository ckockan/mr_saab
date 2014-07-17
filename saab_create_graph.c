int create_the_graph()
{
	int i;
	int j;
	graph_weights_matrix = ( float**) malloc( read_number * sizeof( float*));
	for( i = 0; i < read_number; i++)
	{
		graph_weights_matrix[i] = ( float*) malloc( read_number * sizeof( float));
	}

	/* i and j represent Nodes A and B */
	for( i = 0; i < read_number; i++)
	{
		for( j = 0; j < read_number; j++)
		{
			graph_weights_matrix[i][j] = prob( read[i], read[j]);
		}
	}
}

float prob( struct oea_read x, struct oea_read y)
{
	int dist1;
	int overlap;
	float result;
	float* overlap_score_based_alignment = ( float*) malloc( sizeof( float) * ( strlen( y.seq) + 1));

	dist1 = y.pos - x.pos;
	overlap_score_based_alignment = suffix_prefix_alignments( x.seq, y.seq);

	result = 0;
	/* overlap is the length of the overlap between the end of suffix-prefix */
	for( overlap = 0; overlap < min( strlen( y.seq) + 1, strlen( x.seq) + 1); overlap++) 
	{
		if( strlen( x.seq) - overlap >= 0)
		{ 
			if( result < overlap_score_based_alignment[overlap])
			{
				result = overlap_score_based_alignment[overlap];
			}
		}
	}

	free( overlap_score_based_alignment);
	return result;
}

/* aligns suffix of str1 with prefix of str2 */
float* suffix_prefix_alignments( char* str1, char* str2) 
{
	int row_size; /* row_size is strlen( str2) */
	int col_size; /* col_size is strlen( str1). The last column should be returned. */
	float* return_col;
	float** align_matrix; // alignMatrix[x][y] is aligning str2[x] and str1[y]
	int i;
	int j;

	rowSize = strlen( str2);
	colSize = strlen( str1);
	return_col = ( float*) malloc( ( row_size + 1) * sizeof( float));
	align_matrix=( float**) malloc( ( row_size + 1) * sizeof( float*));

	for( i = 0; i < row_size + 1; i++)
	{
		align_matrix[i] = ( float*) malloc( ( col_size + 1) * sizeof( float));
	}

	for( i = 0; i < row_size + 1; i++)
	{
		for( j = 0; j < col_size + 1; j++)
		{
			align_matrix[i][j] = 0;
		}
	}

	for( i = 0; i < row_size + 1; i++)
	{
		align_matrix[i][0] = ( -1) * INDEL_PENALTY * i;
	}

	for( i = 1; i < row_size + 1; i++)
	{
		for( j = 1; j < col_size + 1; j++)
		{
			if( str1[j-1] != str2[i-1])
			{
				align_matrix[i][j] = max3way( align_matrix[i-1][j-1] - MISMATCH_PENALTY,
									align_matrix[i-1][j] - INDEL_PENALTY, align_matrix[i][j-1] - INDEL_PENALTY);
			}
			else
			{
				align_matrix[i][j] = max3way( align_matrix[i-1][j-1] + MATCH_SCORE,
									align_matrix[i-1][j] - INDEL_PENALTY, align_matrix[i][j-1] - INDEL_PENALTY);
			}
		}	
	}

	for( i = 0; i < row_size + 1; i++)
	{
		return_col[i] = align_matrix[i][col_size];
	}

	for( i = 0; i < row_size + 1; i++)
	{
		free( align_matrix[i]);
	}
	free( align_matrix);

	return return_col;
}
