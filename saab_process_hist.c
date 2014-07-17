void get_len_hist( FILE* histogram, int max_read_length)
{
	int delta_min;
	int delta_max;
	float* len_hist;
	float** prob_overlap_based_pos;
	int len_base;
	float prob;
	int i;
	int j;

	fscanf( histogram, "%d %d\n", &delta_min, &delta_max);
	len_hist = ( float*) malloc( ( delta_max + 1) * sizeof( float));
	prob_overlap_based_pos = ( float**) malloc( ( 2 * deltaMax + 1) * sizeof( float*));

	for( i = 0; i < delta_max; i++)
	{
		len_hist[i] = 0;
	}

	for( i = 0; i < 2 * deltaMax + 1; i++)
	{
		prob_overlap_based_pos[i] = ( float*) malloc( max_read_length * sizeof( float));
	}

	for( i = 0; i < 2 * delta_max + 1; i++)
	{
		for( j = 0; j < max_read_length; j++)
		{
			prob_overlap_based_pos[i][j] = 0;
		}
	}

	while( fscanf( histogram, "%d %f\n", &len_base, &prob) == 2)
	{
		len_hist[len_base] = prob;
	}

	calculate_prob_hist( &prob_overlap_based_pos);
}

void calculate_prob_hist( float*** prob_overlap_based_pos)
{
	int dist1;
	int dist2;
	int paired_end_len;

	/* dist1 is the distance between the ends of the mapped paired-ends if Forward strand */
	for( dist1 = ( -1) * delta_max; dist1 < delta_max; dist1++) 
	{
		/* dist2 is the distance between the ends of the unmapped paired-end if Reverse strand */
		for( dist2 = 0; dist2 < max_read_length; dist2++) 
		{
			/* the distance between the paired-ends of the first paired-end */
			for( paired_end_len = delta_min; paired_end_len < delta_max; paired_end_len++) 
			{
				if( paired_end_len + dist2 - dist1 > 0 && paired_end_len - dist1 + dist2 < deltaMax)
				{
					( *prob_overlap_based_pos)[dist1 + delta_max][dist2] += len_hist[paired_end_len] * len_hist[paired_end_len - dist1 + dist2]; 
				}
			}
		}
	}
}
