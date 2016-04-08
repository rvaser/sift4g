/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* Sept. 7, 2000 queue for amino acids (type char; NOT type Residue!!) */

#ifndef _Queue_h
#define _Queue_h

#include <stdio.h>

extern FILE* errorfp;

struct QueueRecord;
typedef struct QueueRecord *Queue;

int IsEmpty( Queue Q );
int IsFull( Queue Q );
Queue CreateQueue( int MaxElements );
void DisposeQueue( Queue Q );
void MakeEmpty( Queue Q );
void Enqueue( char X, Queue Q );
char Front( Queue Q );
void Dequeue( Queue Q );
char FrontAndDequeue( Queue Q );

#endif  /* _Queue_h */
