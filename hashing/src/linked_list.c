#include <linked_list.h>
#include <stdio.h>
#include <stdlib.h>

int is_member(Employee employee_to_check, Node **head_p)
{
    Node *curr_p = *head_p;

    while (curr_p != NULL && curr_p->employee.id < employee_to_check.id)
    {
        curr_p = curr_p->next;
    }

    return curr_p != NULL && curr_p->employee.id == employee_to_check.id;
}

int insert(Employee new_employee, Node **head_pp)
{
    Node *curr_p = *head_pp;
    Node **pp = head_pp;

    while (curr_p != NULL && curr_p->employee.id < new_employee.id)
    {
        pp = &curr_p->next;
        curr_p = curr_p->next;
    }

    if (curr_p == NULL || curr_p->employee.id > new_employee.id)
    {
        Node *new_p = malloc(sizeof(Node));
        new_p->employee = new_employee;
        new_p->next = curr_p;

        *pp = new_p;

        return 1;
    }

    return 0;
}

int remove(Employee employee, struct Node **head_pp)
{
    struct Node *curr_p = *head_pp;
    struct Node **pp = head_pp;

    while (curr_p)
    {
        if (curr_p->employee.id == employee.id)
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

void print_list(struct Node **head_pp)
{
    struct Node *curr_p = *head_pp;
    while (curr_p)
    {
        printf("%d\t", curr_p->employee.id);
        curr_p = curr_p->next;
    }
    printf("\n");
}