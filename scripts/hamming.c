#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int hamming_distance(const char *str1, const char *str2) {
    int length = strlen(str1);
    int distance = 0;
    for (int i = 0; i < length; i++) {
        if (str1[i] != str2[i]) {
            distance++;
        }
    }

    return distance;
}

int* hamming_comparison(const char *str1, const char *str2, int length) {
    int *result = (int*) malloc(length * sizeof(int));
    
    if (result == NULL) {
        printf("Memory allocation fail\n");
        return NULL;
    }

    for (int i = 0; i < length; i++) {
        if (str1[i] == str2[i]) {
            result[i] = 0;
        } else {
            result[i] = 1;
        }
    }
    
    return result;
}

// free memory function for python
void free_result(int *result) {
    free(result);
}

