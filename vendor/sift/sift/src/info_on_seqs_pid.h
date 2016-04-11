#ifndef INFO_ON_SEQS_PID_H_
#define INFO_ON_SEQS_PID_H_

#include "blimps/blocksprogs.h"

void generate_predictions(Sequence** sequences, int sequences_length, FILE* polymorphism_fp,
	int sequence_identity, FILE* out_fp);

#endif
