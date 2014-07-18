float max( float x, float y)
{
	if( x < y)
	{
		return y;
	}
	else
	{
		return x;
	}
}

float min( float x, float y)
{
	if( x < y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

float max3way( float a, float b, float c)
{
	if( ( a >= b) && ( a >= c))
	{
		return a;
	}
	if( ( b >= a) && ( b >= c))
	{
		return b;
	}
	if( ( c >= b) && ( c >= a))
	{
		return c;
	}
}
