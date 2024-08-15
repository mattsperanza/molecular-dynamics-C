#ifndef MOLECULAR_DYNAMICS_C_LIST_H
#define MOLECULAR_DYNAMICS_C_LIST_H

#include <stdbool.h>

/**
 * Generic double-linked list using void pointers.
 * <p>
 * C++ creator said these are slow af, so only use if absolutely necessary. Since every
 * data pointer points to data away from where nodes are defined.
 */
typedef struct ListNode{
    void* data;
    struct ListNode* next;
    struct ListNode* previous;
} ListNode;

typedef void(*CallbackFree)(void *);
typedef int(*CallbackCompare)(void *a, void *b);

typedef struct LinkedList{
    int size;
    int bytesPerData;
    ListNode* head;
    ListNode* tail;
    CallbackFree callbackFree;
    CallbackCompare callbackCompare;
} LinkedList;

// Construct & Destruct
LinkedList* listCreate(int dataSize, CallbackCompare cbCompare, CallbackFree cbFree);
void listDestroy(LinkedList *list);

// Modifiers
void listAppend(LinkedList *list, void *data);
bool listRemove(LinkedList *list, int index);

// Search
bool listContains(LinkedList *list, void *data);
ListNode* listGetNode(LinkedList *list, int index);

// Extraction
void* listGetNodeData(LinkedList *list, int index);
void listToArray(LinkedList *list, void** ret, int size);

bool hasNext(ListNode* node);
bool hasPrevious(ListNode* node);

#endif //MOLECULAR_DYNAMICS_C_LIST_H
