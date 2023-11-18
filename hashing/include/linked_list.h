#ifndef LINKED_LIST_H
#define LINKED_LIST_H

typedef struct Employee
{
    unsigned int id;
    char *name;
    int age;
    double salary;
    char sex;
    int ssn;
} Employee;

typedef struct Node
{
    Employee employee;
    struct Node *next;
} Node;

int is_member(Employee employee_to_check, Node *head_p);
int insert(Employee new_employee, Node **head_pp);
int remove(Employee employee, Node **head_pp);

#endif