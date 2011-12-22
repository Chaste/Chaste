#include <stdio.h>

#define NUM 40

int main(int argc, char *argv[])
{
  int i, j, k, n;
  n=0;
  printf("%i\t3\t0\t0\n", (NUM+1)*(NUM+1)*(NUM+1));
  for (i=0;i<=NUM;i++)
    {
      for (j=0;j<=NUM;j++)
	{
	  for (k=0;k<=NUM;k++)
	    {
	      printf("%i\t%f\t%f\t%f\n", n++, i*0.1/NUM, j*0.1/NUM, k*0.1/NUM);
	    }
	}
    }
}
