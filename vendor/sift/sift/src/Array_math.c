/* COPYRIGHT 2000 , Fred Hutchinson Cancer Research Center */

#include <assert.h>

#include "Array_math.h"

void
sort_double_list (double array[], int array_length)
{
        int j, P;
        double tmp;

        for (P = 1; P < array_length; P++) {
                tmp = array[P];
                for (j = P; j > 0 && array[j-1] < tmp; j--) {
                        array[j] = array[j-1];
                }
                array[j] = tmp;
        }
} /* end of sort_list */

double
median_of_array (double array[], int array_length)
{
        double median;
        int upper_index, lower_index;

	if (array_length == 1) {
		return array[0];
	}

        sort_double_list (array, array_length);
        if ((array_length % 2) == 0) {
                median = array[(int) array_length/2];
                return median;
        } else {
                lower_index = floor ((double) array_length / 2.0);
                upper_index = ceil ((double) array_length / 2.0);
                median = array[lower_index] + array[upper_index];
                median /= 2;
                return median;
        }

}

double
mean_of_array (double array[], int array_length)
{
	int i;
	double sum;

	assert (array_length != 0);
	sum = 0.0;
	for (i = 0; i < array_length; i++) {
		sum += array[i];
	}
	return (sum/array_length);
}
