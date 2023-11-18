#include <hashing.h>

static const int MAX_BUCKET_SIZE = 10;

Bucket *create_bucket()
{
    Bucket *bucket = (Bucket *)malloc(sizeof(Bucket));
    bucket->capacity = MAX_BUCKET_SIZE;

    Node *head = (Node *)malloc(sizeof(Node));

    bucket->head = &head;

    return bucket;
}

HashTable *create_hash_table(int size)
{
    int num_addresses = 1 << size;
    HashTable *hash_table = (HashTable *)malloc(sizeof(HashTable));

    hash_table->size = size;
    hash_table->buckets = (Bucket **)malloc(sizeof(Bucket *));
    hash_table->buckets[0] = create_bucket();

    return hash_table;
}

HashTable *double_address_table(HashTable *table)
{
}

/**
 * Gets the first depth bits in the binary representation of a hash
 */
char *to_binary(unsigned long hash, int depth)
{
    char *binary = malloc(65);
    binary[64] = '\0';

    for (int i = depth - 1; i >= 0; --i)
    {
        binary[i] = (hash & 1) ? '1' : '0';
        hash >>= 1;
    }

    return binary;
}

/**
 * Hash function from https://stackoverflow.com/a/7666577/13219239
 */
unsigned long hash_employee(Employee *employee, int table_size)
{
    unsigned long hash = 5381;
    int c;

    while ((c = *employee->name++))
    {
        hash = ((hash << 5) + hash) + c;
    }

    return hash % table_size;
}