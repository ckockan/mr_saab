void read_oea_reads( FILE* oea_cluster_file, struct oea_read** oea_reads)
{
	char read_name_buffer[256];
	char sequence_buffer[256];
	int chromosome_number;
	int mapping_position;
	char orientation;
	float score;
	int cluster_id;
	int i;
	
	i = 0;
	while( fscanf( oea_cluster_file, "%s\t%s\t%d\t%d\t%c\t%.3f\t%d\n", 
				   read_name_buffer, sequence_buffer, &chromosome_number, &mapping_position, &orientation, &score, &cluster_id) == 7)
	{
		set_str( ( *oea_reads)[i]->read_name, read_name_buffer);
		set_str( ( *oea_reads)[i]->sequence, sequence_buffer);
		( *oea_reads)[i]->chromosome_number = chromosome_number;
		( *oea_reads)[i]->mapping_position = mapping_position;
		( *oea_reads)[i]->orientation = orientation;
		( *oea_reads)[i]->score = score;
		( *oea_reads)[i]->cluster_id = cluster_id;
		i++;	
	}
}

