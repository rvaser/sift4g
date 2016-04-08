/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */  

#ifndef _QUEUE_C_
#define _QUEUE_C_
#include "queue.h"
#include <stdlib.h>

        struct QueueRecord
        {
            int Capacity;
            int Front;   	/* index to front of queue */
            int Rear; 		/* index to end of queue */
            int Size; 		/* how many elements */
            char *Array;
        };

        int
        IsEmpty( Queue Q )
        {
            return Q->Size == 0;
        }

        int
        IsFull( Queue Q )
        {
            return Q->Size == Q->Capacity;
        }

Queue
CreateQueue( int MaxElements )
{
      Queue Q;


      Q = malloc( sizeof( struct QueueRecord ) );
      if( Q == NULL ) {
          fprintf ( errorfp, "Out of space for creating queue!!!" );
	  exit (-1);
	}

      Q->Array = malloc( sizeof( char ) * MaxElements );
      if( Q->Array == NULL) {
	fprintf (errorfp, "Out of space for creating queue!!");
	exit (-1);
	}
      Q->Capacity = MaxElements;
      MakeEmpty( Q );

      return Q;
}

        void
        MakeEmpty( Queue Q )
        {
            Q->Size = 0;
            Q->Front = 1;
            Q->Rear = 0;
        }

        void
        DisposeQueue( Queue Q )
        {
            if( Q != NULL )
            {
                free( Q->Array );
                free( Q );
            }
        }

        static int
        Succ( int Value, Queue Q )
        {
            if( ++Value == Q->Capacity )
                Value = 0;
            return Value;
        }

        void
        Enqueue( char X, Queue Q )
        {
            if( IsFull( Q ) ) {
                fprintf (errorfp, "Full queue" );
            	exit (-1);
	    }
        	else
            {
                Q->Size++;
                Q->Rear = Succ( Q->Rear, Q );
                Q->Array[ Q->Rear ] = X;
            }
        }



        char 
        Front( Queue Q )
        {
            if( !IsEmpty( Q ) ) 
                return Q->Array[ Q->Front ];
            return '-';  /* Return value used to avoid warning */
        }

        void
        Dequeue( Queue Q )
        {
            if( IsEmpty( Q ) )
                printf ("Error, trying to dequeue empty queue\n");
            else
            {
                Q->Size--;
                Q->Front = Succ( Q->Front, Q );
            }
        }

        char
        FrontAndDequeue( Queue Q )
        {
            char X;

	    X = '-';

            if( IsEmpty( Q ) ) {
                printf ("queye was empty and tried to front & deque it\n");
		return  X;
            } else
            {
                Q->Size--;
                X = Q->Array[ Q->Front ];
                Q->Front = Succ( Q->Front, Q );
            }
            return X;
        }

#endif
