struct edge
{
	int head_node;
	int tail_node;
	float score;
	struct edge* next;
};

struct edge* egde_list;

char* greedy_graph_traversal( float** graph_weights_matrix, int num_nodes)
{
	int** edges_picked;
	float* align_score;
	int overlap_count;
	float overlap_score;
	char* temp_result_seq;
	char* result_seq;
	float min_score_for_adding = 15;
	int node_chosen;
	float max_value = -1000;
	float min_value_in_list = -1000;
	int start_node;
	int stop_node;
	int prev_start_node = stop_node;	
	int prev_stop_node = start_node;
	int first_added;
	int second_added;
	int i;
	int j;

	temp_result_seq = ( char*) malloc( OEASeqMaxLen * sizeof( char));
	egde_list = NULL;

	for( i = 0; i < num_nodes; i++)
	{
		for( j = 0; j < num_nodes; j++)
		{
			if( graph_weights_matrix[i][j] > max_value && i != j)
			{
				add_max_edge( i, j);
				max_value = graph_weights_matrix[i][j];
			}
			else if( graph_weights_matrix[i][j] > max_value - score_variance && i != j)
			{
				add_edge( i, j);
			}
		}
	}

	best_edge_to_pick( &start_node, &stop_node);
	read[start_node].mark = 1;
	read[stop_node].mark = 1;	
	align_score = suffix_prefix_alignments( read[start_node].seq, read[stop_node].seq);
	overlap_count = -1;
	overlap_score = -1;

	for( i = 0; i < strlen( read[stop_node].seq); i++)
	{
		if( align_score[i] > overlap_score)
		{
			overlap_score = align_score[i];
			overlap_count = i;	
		}
	}

	free( align_score);

	prev_start_node = start_node;
	prev_stop_node = stop_node;
	result_seq = add_seq( read[start_node].seq, read[stop_node].seq, overlap_count);

	free_list();
	edge_list = NULL;

	do
	{
		max_value = -1000;
		node_chosen = 0;

		for( i = 0; i < num_nodes; i++)
		{
			first_added = 0;
			second_added = 0;

			if( read[i].mark == 0)
			{	      
				if( max_value < graph_weights_matrix[i][prev_start_node])
				{  
					max_value = graph_weights_matrix[i][prev_start_node];
					add_max_edge( i, prev_start_node);
					first_added = 1;
				}

				if( max_value < graph_weights_matrix[prev_stop_node][i])
				{	  
					max_value = graph_weights_matrix[prev_stop_node][i];
					add_max_edge( prev_stop_node, i);
					second_added = 1; 
				}

				if( first_added == 0 && max_value - score_variance < graph_weights_matrix[i][prev_start_node])
				{
					add_edge( i, prev_start_node);
				}

				if( second_added == 0 &&  max_value - score_variance < graph_weights_matrix[prev_stop_node][i])
				{
					add_edge( prev_stop_node, i);
				}
			}
		}

		best_edge_to_pick( &start_node, &stop_node);
		align_score = suffix_prefix_alignments( read[start_node].seq, read[stop_node].seq);

		if( max_value > min_score_for_adding)
		{
			if( start_node == prev_stop_node)
			{
				free( align_score);
				align_score = suffix_prefix_alignments( read[prev_stop_node].seq, read[stop_node].seq);
				overlap_count = -1;
				overlap_score = -100;
				for( i = 0; i < strlen( read[stop_node].seq) + 1; i++)
				{
					if( align_score[i] > overlap_score)
					{
						overlap_score = align_score[i];
						overlap_count = i;
					}
				}

				strcpy( temp_result_seq, result_seq);
				result_seq = add_seq( result_seq, read[stop_node].seq, overlap_count);

				if( strlen( result_seq) > OEASeqMaxLen)
				{
					return temp_result_seq;
				}

				read[stop_node].mark = 1;
				prev_stop_node = stop_node;
				node_chosen = 1;
				free( align_score);

			}
			else if( stop_node == prev_start_node)
			{
				free( align_score);
				align_score = suffix_prefix_alignments( read[start_node].seq, read[prev_start_node].seq);
				overlap_count = -1;
				overlap_score = -1;

				for( i = 0; i < strlen( read[prev_start_node].seq) + 1; i++)
				{
					if( align_score[i] > overlap_score)
					{
						overlap_score = align_score[i];
						overlap_count = i;
					}
				}

				strcpy( temp_result_seq, result_seq);
				result_seq = add_seq( read[start_node].seq, result_seq, overlap_count); 
                                                      
				if( strlen( result_seq) > OEASeqMaxLen)
				{
					return temp_result_seq;
				}

				read[start_node].mark = 1;
				prev_start_node = start_node;
				node_chosen = 1;
				free( align_score);
			}
		}
		else 
		{
			free( align_score);
		}

		free_list();
		edge_list = NULL;
	} while( max_value > min_score_for_adding && strlen( result_seq) < 400);

	return result_seq;
}

char* add_seq( char* seq1, char* seq2, int overlap_count)
{
	char* result_seq = ( char*) malloc( ( strlen( seq1) + strlen( seq2) + 1 - overlap_count) * sizeof( char));

	int i;
	for( i = 0; i < strlen( seq1); i++)
	{
		resultSeq[i] = seq1[i];
	}

	for( i = overlap_count; i < strlen( seq2); i++)
	{
		result_seq[i + strlen( seq1) - overlap_count] = seq2[i];
	}

	result_seq[( strlen( seq1) + strlen( seq2) - overlap_count)] = '\0';
	return result_seq;
}		

int add_edge( int head, int tail)
{
	struct edge* new_edge;
	new_edge = ( struct edge*) malloc( sizeof( struct edge));
	new_edge->head_hode = head;
	new_edge->tail_node = tail;
	new_edge->score = graph_weights_matrix[head][tail];

	new_edge->next = edge_list;
	edge_list = new_edge;
}

int add_max_edge( int start, int stop)
{
	struct edge* ptList;
	struct edge* headList;
	float max_score;

	max_score = graph_weights_matrix[start][stop];
  
	if( edge_list != NULL)
	{
		while( edge_list != NULL && edge_list->score < max_score - score_variance)
		{
			ptList = edge_list->next;
			edge_list = ptList;
		}

		while( edge_list != NULL && edge_list->next != NULL)
		{
			if( edge_list->next->score < max_score - score_variance)
			{
				ptList = edge_list->next;
				edge_list->next = ptList->next;
			}
			else
			{
				edge_list = edge_list->next;
			}
		}
		add_edge( start, stop);
	}
	else
	{
		add_edge( start, stop);
	}
}


int best_edge_to_pick( int* source, int* dest)
{
	struct edge* ptr_edge;
	float min_value = 1000;
	int min_value_source = -1;
	int min_value_dest = -1;

	ptr_edge = edge_list;
	while( ptr_edge != NULL)
	{
		if( abs( read[ptr_edge->head_node].pos - read[ptr_edge->tail_node].pos) < min_value)
		{
			min_value = abs( read[ptr_edge->head_node].pos - read[ptr_edge->tail_node].pos);
			min_value_source = ptr_edge->head_node;
			min_value_dest = ptr_edge->tail_node;
		}
		ptr_edge = ptr_edge->next;
	}

	*source = min_value_source;
	*dest = min_value_dest;
}


void free_list()
{
	struct edge* temp;
	while( edge_list != NULL)
	{
		temp = edge_list->next;
		free( edge_list);
		edge_list = temp;
	}
	edge_list = NULL;
}
