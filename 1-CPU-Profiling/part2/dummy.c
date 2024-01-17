void dummy_(double *a, double *b, double *c, double *d);
void dummy(double *a, double *b, double *c, double *d);
void dummy_(double *a, double *b, double *c, double *d)
{
  dummy(a, b, c, d);
}

void dummy(double *a, double *b, double *c, double *d)
{
  a[0] = b[0];
}

