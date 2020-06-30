/*
#define int8
#ifdef int8
#define ANUMKIND 8
#define NATMKIND 8
#else
#define ANUMKIND 1
#define NATMKIND 2
#endif
*/
#ifndef ANUMKIND
#define ANUMKIND 1
#endif

#ifndef NATMKIND
#define NATMKIND 2
#endif

#define MAXATOM 99 
#define MAXBAS 150
#define NUMELS 2
#define DEBUG
