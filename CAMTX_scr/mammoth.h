/*
typedef struct{
 int id[2];
 double dist;
} DLIST;
*/
int dlist_comp(const void *_a, const void *_b)
{
  const DLIST *a = (const DLIST *)_a;
  const DLIST *b = (const DLIST *)_b;
  
  if (a->dist > b->dist)
   return 1;
  else if (a->dist < b->dist)
   return -1;
  else
   return 0;
}

