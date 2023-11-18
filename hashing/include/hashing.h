#include <linked_list.h>

#ifndef HASHING_H
#define HASHING_H

typedef struct Bucket
{
    int capacity;
    Node **head;
} Bucket;

typedef struct HashTable
{
    int size;
    Bucket **buckets;
} HashTable;

Bucket *create_bucket();

/**
 * @param size The table will have 2^size buckets
 * @brief Creates a new database
 */
HashTable *create_hash_table(int size);

#endif