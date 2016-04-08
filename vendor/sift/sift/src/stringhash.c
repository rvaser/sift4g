        /* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

#ifndef _STRINGHASH_C_
#define _STRINGHASH_C_

#include "stringhash.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MinTableSize (10)
#define TRUE 1
#define FALSE 0

        /* Return next prime; assume N >= 10 */

        static int
        NextPrime( int N )
        {
            int i;

            if( N % 2 == 0 )
                N++;
            for( ; ; N += 2 )
            {
                for( i = 3; i * i <= N; i += 2 )
                    if( N % i == 0 )
                        goto ContOuter;  /* Sorry about this! */
                return N;
              ContOuter: ;
            }
        }

        Index
	Hash ( char* Key, int TableSize )
        {
            int i =  0;
   	    int j = 0;
	    char* temp;
	    char ch;

	    while (Key[j] != '\0') {
		ch = Key[j];
		i += (int) ch;
		j++;
	    }

	    return i % TableSize;
	}

        HashTable
        InitializeTable( int TableSize )
        {
            HashTable H;
            int i;
            if( TableSize < MinTableSize )
            {
/* 2*/          printf( "Table size too small" );
/* 3*/          return NULL;
            }

            /* Allocate table */
/* 4*/      H = (struct HashTbl* ) malloc( sizeof( struct HashTbl ) );
/* 5*/      if( H == NULL ) {
/* 6*/          fprintf(errorfp, "Out of space to make hash table of %d !!!", TableSize );
		exit (-1);
	    }
/* 7*/      H->TableSize = NextPrime( TableSize );
	    H->no_of_elements = 0;

            /* Allocate array of Cells */
/* 8*/      H->TheLists = (List *) malloc( sizeof( List ) * H->TableSize );
/* 9*/      if( H->TheLists == NULL ) {
/*10*/          fprintf (errorfp, "Out of space to make hash table of %d!!!", TableSize);
		exit (-1);
	    }
/*11*/      for( i = 0; i < H->TableSize; i++ ) {
                H->TheLists[ i ] = (struct ListNode*) malloc
						(sizeof (struct ListNode));
		H->TheLists[i]->key[0] = '\0';
		if (H->TheLists[i] == NULL) {
			fprintf (errorfp, "out of space in hash table making nodes!\n");
			exit (-1);
		} else {
			H->TheLists[i]->Next = NULL;
		}
	    }

/*13*/      return H;
        }

        Position
        Find( char* Key, HashTable H )
        {
            Position P;
	    List L;

	    L = H->TheLists[ Hash (Key, H->TableSize )];
	    P = L->Next;
	    while (P != NULL && (strcmp (P->key, Key) != 0) ) {
	        P = P->Next;
	    }
	return P;
	}

        void
        Insert( char* Key, HashTable H )
        {
            Position Pos, NewCell;
	    List L;

	    Pos = Find( Key, H );
	    if (Pos == NULL) { /* key not found */
                NewCell = (struct ListNode*) malloc (sizeof (struct ListNode ));
		if (NewCell == NULL) {
			fprintf (errorfp, "out of space in insert hash table \n");
			exit (-1);
		} else {
			L = H->TheLists [Hash (Key, H->TableSize) ];
			NewCell->Next = L->Next; /* insert it at the beginning*/
						/* of the list */
			strcpy (NewCell->key, Key);
			L->Next = NewCell;
			H->no_of_elements++;
		}
	    } /* end of if key not found */
	} /*end of insert */


	int
	Exists (char* Key, HashTable H)
	{
		Position P;

      		P = Find (Key, H);
		if (P == NULL) {
			return FALSE;
		} else {
			return TRUE;
		}

	}


        char*
        Retrieve( Position P, HashTable H )
        {
	    return P->key;
        }

        void
        DestroyTable( HashTable H )
        {
            Position P, Temp;
	    int i;
		for (i = 0; i < H->TableSize; i++) {
			P = H->TheLists[i];
			while (P != NULL) {
				Temp = P->Next;
				free (P);
				P = Temp;
			}
		} /* end of for loop */
		free (H->TheLists);
		free(H);

	 }

#endif
