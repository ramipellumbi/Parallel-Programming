#include <stdio.h>
#include <stdlib.h>

struct list_node_s
{
    int data;
    struct list_node_s *next;
};

int member(int value, struct list_node_s **head_p)
{
    struct list_node_s *curr_p = *head_p;

    while (curr_p != NULL && curr_p->data < value)
    {
        curr_p = curr_p->next;
    }

    if (curr_p == NULL || curr_p->data > value)
    {
        return 0;
    }

    return 1;
}

int insert(int value, struct list_node_s **head_pp)
{
    struct list_node_s *curr_p = *head_pp;
    struct list_node_s **pp = head_pp;

    while (curr_p != NULL && curr_p->data < value)
    {
        pp = &curr_p->next;
        curr_p = curr_p->next;
    }

    if (curr_p == NULL || curr_p->data > value)
    {
        struct list_node_s *new_p = malloc(sizeof(struct list_node_s));
        new_p->data = value;
        new_p->next = curr_p;

        *pp = new_p;

        return 1;
    }

    return 0;
}

int delete(int value, struct list_node_s **head_pp)
{
    struct list_node_s *curr_p = *head_pp;
    struct list_node_s **pp = head_pp;

    while (curr_p)
    {
        if (curr_p->data == value)
        {
            *pp = curr_p->next;
            free(curr_p);
            return 1;
        }

        pp = &curr_p->next;
        curr_p = curr_p->next;
    }

    return 0;
}

void print_list(struct list_node_s **head_pp)
{
    struct list_node_s *curr_p = *head_pp;
    while (curr_p)
    {
        printf("%d\t", curr_p->data);
        curr_p = curr_p->next;
    }
    printf("\n");
}

int main(int argc, char *argv[])
{
    struct list_node_s *list = malloc(sizeof(struct list_node_s));

    for (int i = 0; i < 10; ++i)
    {
        insert(i, &list);
    }

    print_list(&list);
    delete (1, &list);
    delete (8, &list);
    print_list(&list);
}