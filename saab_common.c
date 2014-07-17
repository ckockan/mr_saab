float max3way(float a, float b, float c)
{
   if ((a>=b)&&(a>=c))
      return a;
   if ((b>=a)&&(b>=c))
      return b;
   if ((c>=b)&&(c>=a))
      return c;
}

float max (float a, float b)
{
   if (a<b)
      return b;
   else return a;
}
float min(float a, float b)
{

   if (a<b)
      return a;
   else return b;
}
